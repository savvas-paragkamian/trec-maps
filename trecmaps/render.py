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

def world_y(z, exaggeration: float, sea_factor: float = 1.0) -> np.ndarray:
    """Height -> engine Y: `exaggeration` on land, `exaggeration * sea_factor` below sea level (continuous at 0)."""
    z = np.asarray(z, float)
    return np.where(z < 0, z * exaggeration * sea_factor, z * exaggeration)


def exaggerate_sea(z: np.ndarray, sea_factor: float) -> np.ndarray:
    """The heightmap handed to the tracer: depths stretched by `sea_factor` (the tracer applies `exaggeration` on top)."""
    out = np.array(z, dtype=np.float32, copy=True)
    out[out < 0] *= float(sea_factor)
    return np.ascontiguousarray(out)


def pad_block(z: np.ndarray, albedo: np.ndarray, pad: int, base_z: float, cut_rgb: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Surround the heightmap with `pad` cells at `base_z`: the heightfield then drops to a base along its
    rim, which renders as cut walls (the sea floor and the land as a solid block). The rim takes the cut colour."""
    zp = np.pad(z, pad, mode="constant", constant_values=float(base_z)).astype(np.float32)
    ap = np.pad(albedo, ((pad, pad), (pad, pad), (0, 0)), mode="edge")
    ring = np.ones(zp.shape, bool)
    ring[pad:-pad, pad:-pad] = False
    ap[ring, :3] = srgb_to_linear(cut_rgb)
    ap[ring, 3] = 1.0
    return np.ascontiguousarray(zp), np.ascontiguousarray(ap)


def grid_to_world(x, y, z, meta: dict, exaggeration: float, sea_factor: float = 1.0) -> np.ndarray:
    """Map coordinates (metres in the grid CRS) and height -> engine world coordinates (N, 3) = (X east, Y up, Z south).
    The plate is centred on the origin, so symmetric padding does not move the real cells."""
    cell = meta["cell_m"]
    col = (np.asarray(x, float) - meta["minx"]) / cell
    row = (meta["maxy"] - np.asarray(y, float)) / cell
    X = (col - meta["width"] / 2.0) * cell
    Z = (row - meta["height"] / 2.0) * cell
    return np.column_stack([X, world_y(z, exaggeration, sea_factor), Z])


def block_outline(meta: dict, pad: int) -> dict:
    """World X/Z of the block's outer faces (the rim's outer edge) and of the real data's edge."""
    cell = meta["cell_m"]
    W, H = meta["width"], meta["height"]
    return {"x_west": -(W / 2.0 + pad) * cell, "x_east": (W / 2.0 + pad) * cell,
            "z_north": -(H / 2.0 + pad) * cell, "z_south": (H / 2.0 + pad) * cell,
            "xs": (np.arange(W) + 0.5 - W / 2.0) * cell, "zs": (np.arange(H) + 0.5 - H / 2.0) * cell}


def waterline(meta: dict, pad: int) -> np.ndarray:
    """Sea level (Y = 0) along the four walls, as one closed loop of world points (N, 3)."""
    o = block_outline(meta, pad)
    south = np.column_stack([o["xs"], np.zeros_like(o["xs"]), np.full_like(o["xs"], o["z_south"])])
    east = np.column_stack([np.full_like(o["zs"], o["x_east"]), np.zeros_like(o["zs"]), o["zs"][::-1]])
    north = np.column_stack([o["xs"][::-1], np.zeros_like(o["xs"]), np.full_like(o["xs"], o["z_north"])])
    west = np.column_stack([np.full_like(o["zs"], o["x_west"]), np.zeros_like(o["zs"]), o["zs"]])
    loop = np.vstack([south, east, north, west, south[:1]])
    return loop


def depth_axis(meta: dict, pad: int, base_z: float, top_z: float, exaggeration: float, sea_factor: float,
               step_land_m: float, step_sea_m: float, offset_cells: float = 2.0) -> dict:
    """A vertical axis just outside the south-west corner: its line (2, 3) and ticks [(world point, label)]."""
    o = block_outline(meta, pad)
    cell = meta["cell_m"]
    x, zc = o["x_west"] - offset_cells * cell, o["z_south"] + offset_cells * cell
    values = [0.0]
    v = -step_sea_m
    while v >= base_z:
        values.append(v)
        v -= step_sea_m
    v = step_land_m
    while v <= top_z:
        values.append(v)
        v += step_land_m
    ticks = [(np.array([x, float(world_y(val, exaggeration, sea_factor)), zc]),
              "0 (sea level)" if val == 0 else f"{val:+.0f}") for val in sorted(values)]
    line = np.array([[x, float(world_y(base_z, exaggeration, sea_factor)), zc],
                     [x, float(world_y(max(top_z, 0.0), exaggeration, sea_factor)), zc]])
    return {"line": line, "ticks": ticks, "tick_dx": 0.8 * cell}


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


def visible(proj: dict, depth: np.ndarray, kind: str = "dist", rel_tol: float = 0.02, miss_ok: bool = False) -> np.ndarray:
    """True where the projected point is inside the frame and not behind the rendered surface.
    `miss_ok` counts pixels where the ray hits nothing (depth NaN) as visible: points in the air."""
    h, w = depth.shape
    ix, iy = np.rint(proj["px"]).astype(int), np.rint(proj["py"]).astype(int)
    ok = proj["in_front"] & (ix >= 0) & (ix < w) & (iy >= 0) & (iy < h)
    out = np.zeros(len(ix), bool)
    dd = depth[iy[ok], ix[ok]]
    ref = proj[kind][ok]
    hit = np.isfinite(dd)
    out[ok] = np.where(hit, ref <= np.where(hit, dd, np.inf) * (1.0 + rel_tol), miss_ok)
    return out


def pixel_rays(camera: dict, width: int, height: int) -> np.ndarray:
    """Unit ray direction of every pixel centre, (H, W, 3): the engine's pinhole ray generation."""
    origin, forward, right, up = camera_basis(camera)
    half_h = math.tan(math.radians(float(camera["fov_y"])) / 2.0)
    half_w = half_h * width / height
    ndc_x = (np.arange(width) + 0.5) / width * 2.0 - 1.0
    ndc_y = 1.0 - (np.arange(height) + 0.5) / height * 2.0
    d = (ndc_x[None, :, None] * half_w * right[None, None, :] + ndc_y[:, None, None] * half_h * up[None, None, :]
         + forward[None, None, :])
    return d / np.linalg.norm(d, axis=-1, keepdims=True)


def water_alpha(depth_m, cfg: dict) -> np.ndarray:
    """Opacity of the water tint for a true depth below sea level (metres)."""
    d = np.clip(np.asarray(depth_m, float), 0.0, None)
    return float(cfg.get("max_alpha", 0.75)) * (1.0 - np.exp(-d / float(cfg.get("depth_scale_m", 40.0))))


def water_tint(rgb: np.ndarray, depth: np.ndarray, camera: dict, exaggeration: float, sea_factor: float, cfg: dict) -> np.ndarray:
    """Tint every rendered point below sea level (sea floor and cut walls alike) by its true depth, so the
    water column reads as one blue body and its cross-section shows on the walls. Needs the depth AOV
    (ray distance) and the camera to recover each pixel's hit point."""
    h, w = depth.shape
    origin = np.asarray(camera["origin"], float)
    dirs = pixel_rays(camera, w, h)
    hit = np.isfinite(depth)
    y_hit = origin[1] + dirs[..., 1] * np.where(hit, depth, 0.0)
    depth_m = np.where(hit & (y_hit < 0.0), -y_hit / (exaggeration * sea_factor), 0.0)
    alpha = water_alpha(depth_m, cfg)
    colour = hex_to_rgb(cfg.get("color", "#1b4f7a")) * 255.0
    out = rgb.astype(np.float32) * (1.0 - alpha[..., None]) + colour[None, None, :] * alpha[..., None]
    return np.clip(np.rint(out), 0, 255).astype(np.uint8)


def sea_legend_gradient(stops: list, lo: float, hi: float, water_cfg: dict | None, n: int = 256) -> np.ndarray:
    """(1, n, 3) strip of the sea floor albedo with the water tint applied, as it appears in the render."""
    depths = np.linspace(lo, hi, n)
    strip = ramp_rgb(depths, stops)
    if water_cfg and water_cfg.get("enabled", True):
        a = water_alpha(-depths, water_cfg)[:, None]
        strip = strip * (1.0 - a) + hex_to_rgb(water_cfg.get("color", "#1b4f7a"))[None, :] * a
        S = water_cfg.get("surface") or {}
        veil = float(S.get("base_alpha", 0.16)) + float(S.get("fresnel_f0", 0.02))
        strip = strip * (1.0 - veil) + hex_to_rgb(S.get("sky_color", "#e4ebf0"))[None, :] * veil
    return strip[None, :, :]


# ----------------------------------------------------------------------------- water surface (post-processing)

def sun_vector(sun_cfg: dict) -> np.ndarray:
    """Unit vector towards the sun in engine axes (+X east, +Y up, +Z south) from a compass bearing and elevation."""
    b, e = math.radians(float(sun_cfg["azimuth_deg"])), math.radians(float(sun_cfg["elevation_deg"]))
    return np.array([math.sin(b) * math.cos(e), math.sin(e), -math.cos(b) * math.cos(e)], np.float32)


def wave_slopes(X: np.ndarray, Z: np.ndarray, outline: dict, cfg: dict) -> tuple[np.ndarray, np.ndarray]:
    """Surface slopes (dh/dx, dh/dz) of a few directional swells plus a smoothed random ripple field,
    evaluated at world (X, Z). Wavelengths in metres, scaled by `wave_scale`; amplitudes are slopes."""
    scale = float(cfg.get("wave_scale", 1.0))
    sx = np.zeros(X.shape, np.float32)
    sz = np.zeros(X.shape, np.float32)
    for i, (wavelength, slope, direction) in enumerate(cfg.get("waves", [])):
        k = 2.0 * math.pi / (float(wavelength) * scale)
        d = math.radians(float(direction))
        dx, dz = math.sin(d), -math.cos(d)
        arg = k * (X * dx + Z * dz) + 1.7 * i
        c = (float(slope) * np.cos(arg)).astype(np.float32)
        sx += c * dx
        sz += c * dz
    ripple = cfg.get("ripple") or {}
    if ripple.get("slope", 0):
        spacing = float(ripple.get("scale_m", 60.0)) * scale / 4.0
        rng = np.random.default_rng(int(cfg.get("seed", 3)))
        nx = int((outline["x_east"] - outline["x_west"]) / spacing) + 4
        nz = int((outline["z_south"] - outline["z_north"]) / spacing) + 4
        field = ndimage.gaussian_filter(rng.standard_normal((nz, nx)).astype(np.float32), 3.0)
        field /= max(float(field.std()), 1e-6)
        height = field * float(ripple["slope"]) * float(ripple.get("scale_m", 60.0)) * scale
        gz, gx = np.gradient(height, spacing)
        col = (X - outline["x_west"]) / spacing + 1.0
        row = (Z - outline["z_north"]) / spacing + 1.0
        sx += ndimage.map_coordinates(gx, [row, col], order=3, mode="nearest").astype(np.float32)
        sz += ndimage.map_coordinates(gz, [row, col], order=3, mode="nearest").astype(np.float32)
    return sx, sz


def water_surface(rgb: np.ndarray, depth: np.ndarray, camera: dict, meta: dict, pad: int, exaggeration: float,
                  sea_factor: float, sun_cfg: dict, cfg: dict) -> np.ndarray:
    """Composite a semi-transparent, wavy sea surface at Y = 0 over everything seen through it.

    Per pixel: the ray's crossing of sea level inside the block (the surface point), the hit point
    below it (the floor, attenuated towards the water colour by its true depth and displaced a little
    by the wave slope: refraction), a wave normal, Fresnel reflection of the sky, a sun glint. Cut
    walls below sea level keep the depth tint only (the water column in cross-section)."""
    h, w = depth.shape
    S = cfg.get("surface") or {}
    origin = np.asarray(camera["origin"], np.float32)
    dirs = pixel_rays(camera, w, h).astype(np.float32)
    hit = np.isfinite(depth)
    dd = np.where(hit, depth, 0.0).astype(np.float32)
    y_hit = origin[1] + dirs[..., 1] * dd
    with np.errstate(divide="ignore", invalid="ignore"):
        t0 = np.where(dirs[..., 1] < 0, -origin[1] / dirs[..., 1], np.inf).astype(np.float32)
    o = block_outline(meta, pad)
    Xs = origin[0] + dirs[..., 0] * t0
    Zs = origin[2] + dirs[..., 2] * t0
    inside = (Xs >= o["x_west"]) & (Xs <= o["x_east"]) & (Zs >= o["z_north"]) & (Zs <= o["z_south"])
    under = hit & (y_hit < 0.0)
    water = under & inside & (t0 < dd)
    depth_m = np.where(under, -y_hit / (exaggeration * sea_factor), 0.0)

    # the floor as seen through the water: refracted a little, then attenuated towards the water colour
    body = hex_to_rgb(cfg.get("color", "#1b4f7a")) * 255.0
    sx, sz = wave_slopes(np.where(water, Xs, 0.0), np.where(water, Zs, 0.0), o, S)
    shift = float(S.get("refraction_gain", 150.0)) * w / 3840.0      # pixels per unit slope, at 3840 wide
    rows, cols = np.mgrid[0:h, 0:w].astype(np.float32)
    r2 = np.clip(rows + np.where(water, sz, 0.0) * shift, 0, h - 1)
    c2 = np.clip(cols + np.where(water, sx, 0.0) * shift, 0, w - 1)
    base = rgb.astype(np.float32)
    floor = np.stack([ndimage.map_coordinates(base[..., c], [r2, c2], order=1, mode="nearest") for c in range(3)], axis=-1)
    floor = np.where(water[..., None], floor, base)
    a_body = water_alpha(depth_m, cfg)[..., None]
    seen = floor * (1.0 - a_body) + body[None, None, :] * a_body     # applies to the walls too (cross-section)

    # the surface: wave normal, Fresnel, sky reflection shaded by the sun, glint
    n = np.stack([-sx, np.ones_like(sx), -sz], axis=-1)
    n /= np.linalg.norm(n, axis=-1, keepdims=True)
    v = -dirs
    cos_v = np.clip(np.sum(n * v, axis=-1), 0.0, 1.0)
    f0 = float(S.get("fresnel_f0", 0.02))
    fresnel = f0 + (1.0 - f0) * (1.0 - cos_v) ** 5
    # A mirror-like surface reflects a sky gradient: a small tilt of the normal swings the reflected
    # direction a lot, so the tilt towards the sun is amplified (`wave_shading`) rather than used as
    # a Lambert term. Crests brighter than `crest` get a touch of white.
    L = sun_vector(sun_cfg)
    tilt = np.sum(n * L[None, None, :], axis=-1) - L[1]
    shade = np.clip(1.0 + float(S.get("wave_shading", 8.0)) * tilt, 0.45, 1.7)
    crest = float(S.get("crest", 1.3))
    foam = np.clip(shade - crest, 0.0, None) * float(S.get("crest_white", 0.5))
    sky = hex_to_rgb(S.get("sky_color", "#e4ebf0")) * 255.0
    A = np.clip(float(S.get("base_alpha", 0.16)) + fresnel, 0.0, 0.95)[..., None]
    surface = seen * (1.0 - A) + (sky[None, None, :] * shade[..., None]) * A + 255.0 * foam[..., None]
    out = np.where(water[..., None], surface, seen)
    return np.clip(np.rint(out), 0, 255).astype(np.uint8)
