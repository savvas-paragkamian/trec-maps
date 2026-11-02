"""forge3d path-traced terrain plate, rendered in tiles, plus the camera maths for overlays.

Engine conventions (verified on forge3d 1.44.0, Metal, in ~/projects/crete3d):
  * grid column -> +X (east), grid row -> +Z (south), height -> +Y;
  * the plate is centred on the origin, in units of ``spacing`` (metres here);
  * ``sun_azimuth_deg`` runs counter-clockwise relative to compass bearings;
  * every native call has a hard 512 MiB budget (about 50 bytes per cell plus the tile);
  * tile sizes must be the frame size divided by a power of two (float32 sensor_rect check).
"""
from __future__ import annotations

import math
import time

import numpy as np
from scipy import ndimage

from .config import hex_to_rgb


def srgb_to_linear(rgb: np.ndarray) -> np.ndarray:
    rgb = np.clip(np.asarray(rgb, dtype=np.float32), 0.0, 1.0)
    return np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4).astype(np.float32)


def compass_to_engine_azimuth(compass_deg: float) -> float:
    return (360.0 - float(compass_deg)) % 360.0


def compass_unit(compass_deg: float) -> tuple[float, float]:
    """Compass bearing -> engine (x, z) unit vector (east = +x, north = -z)."""
    a = math.radians(compass_deg)
    return math.sin(a), -math.cos(a)


def make_camera(dem_shape: tuple[int, int], cell_m: float, size_px: tuple[int, int], cfg: dict, exaggeration: float) -> dict:
    """Pinhole camera on a compass bearing, placed so the plate's near edge spans `fill` of the frame width."""
    rows, cols = dem_shape
    plate_w, plate_d = cols * cell_m, rows * cell_m
    el = math.radians(cfg["elevation_deg"])
    fov_y = math.radians(cfg["fov_y_deg"])
    aspect = size_px[0] / size_px[1]
    tan_half_x = math.tan(fov_y / 2) * aspect
    distance = (plate_w / 2) / (cfg["fill"] * tan_half_x) + (plate_d / 2) * math.cos(el)
    ux, uz = compass_unit(cfg["azimuth_deg"])
    target = (cfg.get("target_shift_east_frac", 0.0) * plate_w, cfg.get("target_lift_m", 0.0) * exaggeration,
              -cfg.get("target_shift_north_frac", 0.0) * plate_d)
    horiz = distance * math.cos(el)
    origin = (target[0] + ux * horiz, target[1] + distance * math.sin(el), target[2] + uz * horiz)
    return {"model": "pinhole", "origin": origin, "look_at": target, "up": (0.0, 1.0, 0.0), "fov_y": cfg["fov_y_deg"]}


def tile_grid(width: int, height: int, max_tile: int) -> tuple[int, int]:
    """Tile sizes W/2^j, H/2^k so every sensor_rect boundary is an exact binary fraction."""
    sizes = []
    for extent in (width, height):
        j = 0
        while extent / 2 ** j > max_tile:
            j += 1
        if extent % 2 ** j:
            raise ValueError(f"frame size {extent} is not divisible by 2^{j}; pick a size like 3840, 2560 or 1600")
        sizes.append(extent // 2 ** j)
    return sizes[0], sizes[1]


def render_plate(dem: np.ndarray, albedo_rgba: np.ndarray, *, cell_m: float, size_px: tuple[int, int], frames: int,
                 camera_cfg: dict, sun_cfg: dict, exaggeration: float, env_intensity: float, tile: int = 960,
                 albedo_sampling: str = "bilinear", turbidity: float = 1.0, seed: int = 11,
                 earth_model: str | None = None, refraction_model: str | None = None, log=print) -> dict:
    """Render the full frame as tiles; returns rgba (uint8) plus the depth/albedo/normal AOVs and the camera."""
    from forge3d.path_tracing import hybrid_render_terrain_reference  # needs the GPU wheel: import late

    dem = np.ascontiguousarray(dem, dtype=np.float32)
    albedo_rgba = np.ascontiguousarray(albedo_rgba, dtype=np.float32)
    assert albedo_rgba.shape[:2] == dem.shape, (albedo_rgba.shape, dem.shape)
    width, height = size_px
    camera = make_camera(dem.shape, cell_m, size_px, camera_cfg, exaggeration)
    sun_az = compass_to_engine_azimuth(sun_cfg["azimuth_deg"])
    extra = {}
    if earth_model:
        extra["earth_model"] = earth_model
    if refraction_model:
        extra["refraction_model"] = refraction_model

    rgba = np.zeros((height, width, 4), np.uint8)
    depth = np.full((height, width), np.nan, np.float32)
    albedo = np.zeros((height, width, 3), np.float32)
    normal = np.zeros((height, width, 3), np.float32)
    tw, th = tile_grid(width, height, tile)
    n_x, n_y = width // tw, height // th
    t0 = time.time()
    for iy in range(n_y):
        for ix in range(n_x):
            x0, y0 = ix * tw, iy * th
            rect = (x0 / width, y0 / height, (x0 + tw) / width, (y0 + th) / height)
            res = hybrid_render_terrain_reference(
                dem, tw, th, camera, spacing=(cell_m, cell_m), exaggeration=exaggeration,
                albedo_map=albedo_rgba, albedo_sampling=albedo_sampling,
                sensor_rect=rect, full_width=width, full_height=height, pixel_offset=(x0, y0),
                sun_azimuth_deg=sun_az, sun_elevation_deg=sun_cfg["elevation_deg"],
                sun_intensity=sun_cfg["intensity"], env_intensity=env_intensity,
                turbidity=turbidity, spp=1, min_frames=max(2, frames), max_frames=max(2, frames),
                variance_threshold=1e9, seed=seed, **extra,
            )
            rgba[y0:y0 + th, x0:x0 + tw] = res["rgba"]
            depth[y0:y0 + th, x0:x0 + tw] = res["depth"]
            albedo[y0:y0 + th, x0:x0 + tw] = res["albedo"]
            normal[y0:y0 + th, x0:x0 + tw] = res["normal"]
            log(f"tile {iy * n_x + ix + 1}/{n_x * n_y} ({tw}x{th}) {time.time() - t0:5.1f}s "
                f"host {res.get('peak_host_visible_bytes', 0) / 2**20:.0f} MiB")
    return {"rgba": rgba, "depth": depth, "albedo": albedo, "normal": normal, "camera": camera,
            "seconds": time.time() - t0}


# ----------------------------------------------------------------------------- albedo

def ramp_rgb(values: np.ndarray, stops: list) -> np.ndarray:
    """Per-channel linear interpolation of sRGB hex stops [[height_m, '#rrggbb'], ...] -> floats (..., 3)."""
    stops = sorted(stops, key=lambda s: float(s[0]))
    xs = np.array([float(s[0]) for s in stops])
    cols = np.array([hex_to_rgb(s[1]) for s in stops])      # (n, 3)
    v = np.asarray(values, dtype=np.float32)
    return np.stack([np.interp(v, xs, cols[:, c]) for c in range(3)], axis=-1).astype(np.float32)


def soften_hex(color: str, soften) -> np.ndarray:
    """Scale the HSV saturation and value of a colour (trec-coastal-microbiome's `landcover_palette`)."""
    import colorsys
    h, s_, v = colorsys.rgb_to_hsv(*hex_to_rgb(color))
    return np.array(colorsys.hsv_to_rgb(h, s_ * float(soften[0]), v * float(soften[1])), np.float32)


def landcover_lut(classes: dict, soften, water: str, fallback: str) -> np.ndarray:
    """(256, 3) sRGB look-up table: class code -> softened legend colour; 80 (water) and unknown codes aside."""
    lut = np.tile(hex_to_rgb(fallback), (256, 1)).astype(np.float32)
    for code, (_, colour) in classes.items():
        lut[int(code)] = soften_hex(colour, soften)
    lut[80] = hex_to_rgb(water)
    return lut


def make_albedo(z: np.ndarray, land: np.ndarray, cfg: dict, landcover: np.ndarray | None = None,
                lc_cfg: dict | None = None) -> np.ndarray:
    """Land cover (or hypsometric tint) on land, depth tint at sea, a light shoreline ring; linear RGB, alpha 1."""
    use_classes = (cfg.get("land_mode", "hypsometric") == "landcover" and landcover is not None
                   and lc_cfg is not None and bool((landcover > 0).any()))
    if use_classes:
        classes = {int(k): (v[0], v[1]) for k, v in lc_cfg["classes"].items()}
        lut = landcover_lut(classes, lc_cfg.get("soften_3d", [0.55, 0.92]), lc_cfg.get("water_color", "#a9c6d8"),
                            lc_cfg.get("fallback_color", "#cfc6b4"))
        land_rgb = lut[landcover]
    else:
        land_rgb = ramp_rgb(z, cfg["land"])
    rgb = np.where(land[..., None], land_rgb, ramp_rgb(z, cfg["sea"]))
    shore = cfg.get("shoreline") or {}
    if int(shore.get("width_cells", 0)) > 0:
        ring = ndimage.binary_dilation(land, iterations=int(shore["width_cells"])) & ~land
        rgb[ring] = hex_to_rgb(shore["color"])
    out = np.ones(z.shape + (4,), np.float32)
    out[..., :3] = srgb_to_linear(rgb)
    return out


def legend_gradient(stops: list, lo: float, hi: float, n: int = 256) -> np.ndarray:
    """(1, n, 3) sRGB strip of a ramp between two heights, for the figure legend."""
    return ramp_rgb(np.linspace(lo, hi, n), stops)[None, :, :]


# ----------------------------------------------------------------------------- projection for overlays

def grid_to_world(x, y, z, meta: dict, exaggeration: float) -> np.ndarray:
    """Map coordinates (metres in the grid CRS) and height -> engine world coordinates (N, 3) = (X east, Y up, Z south)."""
    cell = meta["cell_m"]
    col = (np.asarray(x, float) - meta["minx"]) / cell
    row = (meta["maxy"] - np.asarray(y, float)) / cell
    X = (col - meta["width"] / 2.0) * cell
    Z = (row - meta["height"] / 2.0) * cell
    Y = np.asarray(z, float) * exaggeration
    return np.column_stack([X, Y, Z])


def camera_basis(camera: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    origin = np.asarray(camera["origin"], float)
    forward = np.asarray(camera["look_at"], float) - origin
    forward /= np.linalg.norm(forward)
    right = np.cross(forward, np.asarray(camera["up"], float))
    right /= np.linalg.norm(right)
    up = np.cross(right, forward)
    return origin, forward, right, up


def world_to_pixel(points: np.ndarray, camera: dict, width: int, height: int) -> dict:
    """Inverse of the engine's pinhole ray generation: continuous pixel coordinates (centre of pixel i at i),
    the ray distance and the camera-axis depth of each point, and whether it lies in front of the camera."""
    origin, forward, right, up = camera_basis(camera)
    d = np.asarray(points, float) - origin
    zc, xc, yc = d @ forward, d @ right, d @ up
    half_h = math.tan(math.radians(float(camera["fov_y"])) / 2.0)
    half_w = half_h * width / height
    with np.errstate(divide="ignore", invalid="ignore"):
        ndc_x = xc / (zc * half_w)
        ndc_y = yc / (zc * half_h)
    px = (ndc_x + 1.0) / 2.0 * width - 0.5
    py = (1.0 - ndc_y) / 2.0 * height - 0.5
    return {"px": px, "py": py, "dist": np.linalg.norm(d, axis=1), "zc": zc, "in_front": zc > 0}


def depth_kind(proj: dict, depth: np.ndarray, mask: np.ndarray | None = None) -> tuple[str, float]:
    """Does the depth AOV hold the ray distance or the camera-axis depth? Returns the better match and its median ratio."""
    h, w = depth.shape
    ix, iy = np.rint(proj["px"]).astype(int), np.rint(proj["py"]).astype(int)
    ok = proj["in_front"] & (ix >= 0) & (ix < w) & (iy >= 0) & (iy < h)
    if mask is not None:
        ok &= mask
    dd = depth[iy[ok], ix[ok]]
    ok2 = np.isfinite(dd)
    if not ok2.any():
        return "dist", float("nan")
    r_dist = np.median(dd[ok2] / proj["dist"][ok][ok2])
    r_zc = np.median(dd[ok2] / proj["zc"][ok][ok2])
    return ("dist", float(r_dist)) if abs(r_dist - 1) <= abs(r_zc - 1) else ("zc", float(r_zc))


def visible(proj: dict, depth: np.ndarray, kind: str = "dist", rel_tol: float = 0.02) -> np.ndarray:
    """True where the projected point is inside the frame and not behind the rendered surface."""
    h, w = depth.shape
    ix, iy = np.rint(proj["px"]).astype(int), np.rint(proj["py"]).astype(int)
    ok = proj["in_front"] & (ix >= 0) & (ix < w) & (iy >= 0) & (iy < h)
    out = np.zeros(len(ix), bool)
    dd = depth[iy[ok], ix[ok]]
    ref = proj[kind][ok]
    out[ok] = np.isfinite(dd) & (ref <= dd * (1.0 + rel_tol))
    return out
