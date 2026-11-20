"""One metric grid for a transect: Copernicus GLO-30 on land, EMODnet bathymetry at sea.

Nothing here goes online. The DEM tiles and the EMODnet extract are fetched by
scripts/00_get_dem.py and scripts/01_get_bathymetry.py; this module only decides which
files they need (`dem_tiles_for`, `wcs_request`) and merges them (`build_grid`).
The metric CRS is the path's UTM zone when `grid.crs` is `auto`.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject
from scipy import ndimage

from .config import ROOT, read_json

RESAMPLING = {"nearest": Resampling.nearest, "bilinear": Resampling.bilinear,
              "cubic": Resampling.cubic, "average": Resampling.average}
SIDES = ("north", "south", "east", "west")


# ----------------------------------------------------------------------------- grid definition

def utm_crs(points: list[dict]) -> str:
    """The WGS 84 / UTM zone of the path's centre (EPSG:326xx north, 327xx south)."""
    lon = sum(p["lon"] for p in points) / len(points)
    lat = sum(p["lat"] for p in points) / len(points)
    zone = int((lon + 180) // 6) + 1
    return f"EPSG:{(32600 if lat >= 0 else 32700) + zone}"


def metric_bbox(points: list[dict], crs: str, margin, cell_m: float, max_cells: int | None = None) -> dict:
    """Bounding box of the path in the metric CRS plus margins, snapped outward to whole cells.
    `crs: auto` = the path's UTM zone; above `max_cells` the cell is coarsened to a multiple of 5 m."""
    if str(crs).lower() == "auto":
        crs = utm_crs(points)
    tr = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    xs, ys = tr.transform([p["lon"] for p in points], [p["lat"] for p in points])
    m = {k: float(margin.get(k, 0)) for k in SIDES} if isinstance(margin, dict) else {k: float(margin) for k in SIDES}
    minx = math.floor((min(xs) - m["west"]) / cell_m) * cell_m
    maxx = math.ceil((max(xs) + m["east"]) / cell_m) * cell_m
    miny = math.floor((min(ys) - m["south"]) / cell_m) * cell_m
    maxy = math.ceil((max(ys) + m["north"]) / cell_m) * cell_m
    width, height = int(round((maxx - minx) / cell_m)), int(round((maxy - miny) / cell_m))
    if max_cells and width * height > max_cells:
        coarser = math.ceil(cell_m * math.sqrt(width * height / max_cells) / 5) * 5
        return metric_bbox(points, crs, margin, coarser)
    return {"crs": crs, "cell_m": float(cell_m), "minx": minx, "miny": miny, "maxx": maxx, "maxy": maxy,
            "width": width, "height": height}


def grid_box(G: dict, points: list[dict]) -> tuple[dict, tuple]:
    """(metric bbox, lon/lat box) of a transect's grid from the `grid` section of params.yml."""
    bbox = metric_bbox(points, G["crs"], G["margin_m"], G["cell_m"], G.get("max_cells"))
    return bbox, bbox_lonlat(bbox)


def grid_transform(bbox: dict):
    return from_origin(bbox["minx"], bbox["maxy"], bbox["cell_m"], bbox["cell_m"])


def bbox_lonlat(bbox: dict, pad_deg: float = 0.003) -> tuple[float, float, float, float]:
    """(minlon, minlat, maxlon, maxlat) enclosing the metric box: its edges sampled, then padded."""
    tr = Transformer.from_crs(bbox["crs"], "EPSG:4326", always_xy=True)
    xs = np.linspace(bbox["minx"], bbox["maxx"], 25)
    ys = np.linspace(bbox["miny"], bbox["maxy"], 25)
    ex = np.concatenate([xs, xs, np.full(25, bbox["minx"]), np.full(25, bbox["maxx"])])
    ey = np.concatenate([np.full(25, bbox["miny"]), np.full(25, bbox["maxy"]), ys, ys])
    lon, lat = tr.transform(ex, ey)
    return (float(lon.min()) - pad_deg, float(lat.min()) - pad_deg, float(lon.max()) + pad_deg, float(lat.max()) + pad_deg)


# ----------------------------------------------------------------------------- sources

def tiles_for(lonlat: tuple, size_deg: int = 1) -> list[tuple[str, str]]:
    """Tile names by south-west corner ('N37', 'E023') for `size_deg`-degree tiles covering a lon/lat box."""
    minlon, minlat, maxlon, maxlat = lonlat
    corner = lambda v: math.floor(v / size_deg) * size_deg
    tiles = []
    for lat in range(corner(minlat), corner(maxlat) + 1, size_deg):
        for lon in range(corner(minlon), corner(maxlon) + 1, size_deg):
            tiles.append((f"{'N' if lat >= 0 else 'S'}{abs(lat):02d}", f"{'E' if lon >= 0 else 'W'}{abs(lon):03d}"))
    return tiles


def dem_tiles_for(lonlat: tuple) -> list[tuple[str, str]]:
    """GLO-30 tile names (1-degree tiles) covering a lon/lat box."""
    return tiles_for(lonlat, 1)


def landcover_tiles_for(params: dict, lonlat: tuple) -> list[tuple[str, str]]:
    return tiles_for(lonlat, int(params["sources"]["landcover"].get("tile_deg", 3)))


def landcover_url(params: dict, lat_s: str, lon_s: str) -> str:
    return params["sources"]["landcover"]["url"].format(lat=lat_s, lon=lon_s)


def landcover_path(params: dict, name: str, lat_s: str, lon_s: str) -> Path:
    return params["paths"]["raw_dir"] / "landcover" / f"worldcover_2021_{name}_{lat_s}{lon_s}.tif"


def landcover_paths(params: dict, name: str) -> list[Path]:
    return sorted((params["paths"]["raw_dir"] / "landcover").glob(f"worldcover_2021_{name}_*.tif"))


def dem_tile_url(params: dict, lat_s: str, lon_s: str) -> str:
    return params["sources"]["dem"]["url"].format(lat=lat_s, lon=lon_s)


def dem_path(params: dict, name: str, lat_s: str, lon_s: str) -> Path:
    """The window of one GLO-30 tile over a transect's grid (scripts/00_get_dem.py reads it from the COG)."""
    return params["paths"]["raw_dir"] / "dem" / f"glo30_{name}_{lat_s}{lon_s}.tif"


def dem_paths(params: dict, name: str) -> list[Path]:
    return sorted((params["paths"]["raw_dir"] / "dem").glob(f"glo30_{name}_*.tif"))


def bathy_path(params: dict, name: str) -> Path:
    return params["paths"]["raw_dir"] / "bathy" / f"emodnet_mean_{name}.tif"


def coastal_zones_path(params: dict, name: str) -> Path:
    """The Coastal Zones 2018 polygons cut to a transect's grid box (scripts/03_get_coastal_zones.py)."""
    return params["paths"]["raw_dir"] / "coastal_zones" / f"cz_2018_{name}.geojson.gz"


def wcs_request(lonlat: tuple, cfg: dict) -> dict:
    """EMODnet WCS GetCoverage snapped to the native grid (1/960 degree), so nothing is resampled server side."""
    per_deg = int(round(1.0 / float(cfg["native_deg"])))
    pad = int(cfg.get("margin_cells", 0))
    i0, i1 = math.floor(lonlat[0] * per_deg) - pad, math.ceil(lonlat[2] * per_deg) + pad
    j0, j1 = math.floor(lonlat[1] * per_deg) - pad, math.ceil(lonlat[3] * per_deg) + pad
    bbox = (i0 / per_deg, j0 / per_deg, i1 / per_deg, j1 / per_deg)
    query = {"service": "WCS", "version": "1.0.0", "request": "GetCoverage", "coverage": cfg["coverage"],
             "crs": "EPSG:4326", "bbox": ",".join(f"{v:.10f}" for v in bbox), "width": i1 - i0, "height": j1 - j0,
             "format": "GeoTIFF", "interpolation": "nearest"}
    return {"url": cfg["wcs"] + "?" + urlencode(query), "bbox": bbox, "width": i1 - i0, "height": j1 - j0}


def tile_extent(lat_s: str, lon_s: str, size_deg: int) -> tuple[float, float, float, float]:
    """(minlon, minlat, maxlon, maxlat) of a tile named by its south-west corner ('N37', 'E023')."""
    lat = int(lat_s[1:]) * (1 if lat_s[0] == "N" else -1)
    lon = int(lon_s[1:]) * (1 if lon_s[0] == "E" else -1)
    return float(lon), float(lat), float(lon + size_deg), float(lat + size_deg)


def window_is_short(path: Path, lonlat: tuple, extent: tuple) -> bool:
    """Does an existing window file fall short of the box (clipped to its tile's extent)? -> fetch it again.
    One cell is allowed: GLO-30 tiles start half a cell off the whole degree, so a window at a tile edge
    is always that much short of the nominal extent."""
    clipped = (max(lonlat[0], extent[0]), max(lonlat[1], extent[1]), min(lonlat[2], extent[2]), min(lonlat[3], extent[3]))
    return (path.exists() and clipped[2] > clipped[0] and clipped[3] > clipped[1]
            and not covers(path, clipped, tol_cells=1.0))


def covers(raster_path: Path, lonlat: tuple, tol_deg: float = 1e-6, tol_cells: float = 0.0) -> bool:
    """Does the EPSG:4326 raster's extent contain the (minlon, minlat, maxlon, maxlat) box?"""
    with rasterio.open(raster_path) as src:
        b, tol = src.bounds, tol_deg + tol_cells * max(src.res)
    return (b.left <= lonlat[0] + tol and b.bottom <= lonlat[1] + tol
            and b.right >= lonlat[2] - tol and b.top >= lonlat[3] - tol)


def polygons_cover(path: Path, lonlat: tuple, tol_deg: float = 1e-6) -> bool:
    """Was the polygon file (its `bbox` member) cut to a box that contains this one?"""
    b = read_json(path).get("bbox") or (0, 0, 0, 0)
    return (b[0] <= lonlat[0] + tol_deg and b[1] <= lonlat[1] + tol_deg
            and b[2] >= lonlat[2] - tol_deg and b[3] >= lonlat[3] - tol_deg)


def clip_polygons(features: list[dict], lonlat: tuple, field: str) -> list[dict]:
    """GeoJSON polygon features (EPSG:4326) cut to a lon/lat box, as features with `code` (the class in
    `field`) and `area` (square degrees), the largest first. The sea polygons of Coastal Zones span
    whole basins; only the piece inside the box is kept."""
    import shapely
    from shapely.geometry import mapping, shape

    rect = shapely.box(*lonlat)
    out = []
    for f in features:
        if not f.get("geometry"):
            continue
        g = shape(f["geometry"])
        g = (g if g.is_valid else shapely.make_valid(g)).intersection(rect)
        # a polygon that only touches the box leaves lines and points: keep the polygons
        parts = [q for p in shapely.get_parts(g) for q in shapely.get_parts(p) if q.geom_type == "Polygon" and q.area > 0]
        if parts:
            g = shapely.MultiPolygon(parts)
            out.append({"type": "Feature", "properties": {"code": int(f["properties"][field]), "area": float(g.area)},
                        "geometry": mapping(g)})
    return sorted(out, key=lambda f: -f["properties"]["area"])


# ----------------------------------------------------------------------------- reprojection

def _reproject(source, bbox: dict, resampling: str, *, src_transform=None, src_crs=None, src_nodata=None) -> np.ndarray:
    dst = np.full((bbox["height"], bbox["width"]), np.nan, np.float32)
    reproject(source, dst, src_transform=src_transform, src_crs=src_crs, src_nodata=src_nodata,
              dst_transform=grid_transform(bbox), dst_crs=bbox["crs"], dst_nodata=np.nan,
              resampling=RESAMPLING[resampling])
    return dst


def read_dem(paths: list[Path], bbox: dict, resampling: str = "bilinear") -> np.ndarray:
    """DEM tiles onto the grid; later tiles only fill cells the earlier ones left NaN."""
    dem = np.full((bbox["height"], bbox["width"]), np.nan, np.float32)
    for path in paths:
        with rasterio.open(path) as src:
            tmp = _reproject(rasterio.band(src, 1), bbox, resampling, src_nodata=src.nodata)
        fill = np.isnan(dem) & np.isfinite(tmp)
        dem[fill] = tmp[fill]
    return dem


def read_classes(paths: list[Path], bbox: dict) -> np.ndarray:
    """Class rasters (uint8) onto the grid by nearest neighbour; 0 where nothing maps."""
    out = np.zeros((bbox["height"], bbox["width"]), np.uint8)
    for path in paths:
        with rasterio.open(path) as src:
            tmp = _reproject(rasterio.band(src, 1), bbox, "nearest", src_nodata=src.nodata)
        fill = (out == 0) & np.isfinite(tmp)
        out[fill] = tmp[fill].astype(np.uint8)
    return out


def read_polygons(path: Path, bbox: dict) -> np.ndarray:
    """Class polygons (`clip_polygons` features) burnt onto the grid as int32 codes, in the file's order
    (the largest first, so a small polygon is never hidden); 0 where nothing maps. Cells within two
    cells of a polygon take the nearest class: the generalised outlines leave slivers one cell wide."""
    from rasterio.features import rasterize
    from rasterio.warp import transform_geom

    features = read_json(path)["features"]
    shape = (bbox["height"], bbox["width"])
    if not features:
        return np.zeros(shape, np.int32)
    geoms = transform_geom("EPSG:4326", bbox["crs"], [f["geometry"] for f in features])
    codes = rasterize(zip(geoms, (int(f["properties"]["code"]) for f in features)), out_shape=shape,
                      transform=grid_transform(bbox), fill=0, dtype="int32")
    if (codes > 0).any():
        dist, idx = ndimage.distance_transform_edt(codes == 0, return_indices=True)
        codes = np.where(dist <= 2.0, codes[tuple(idx)], codes)
    return codes


def read_bathy(path: Path, bbox: dict, resampling: str = "cubic") -> tuple[np.ndarray, dict]:
    """EMODnet extract onto the grid. Fill values, land and zero cells become nodata before resampling,
    and the sign is checked: the grid must hold negative metres below the datum."""
    with rasterio.open(path) as src:
        raw = src.read(1).astype(np.float32)
        src_transform, src_crs, nodata = src.transform, src.crs, src.nodata
    valid = np.isfinite(raw) & (np.abs(raw) < 12000)
    if nodata is not None and np.isfinite(nodata):
        valid &= raw != nodata
    info = {"raw_shape": list(raw.shape), "raw_nodata": None if nodata is None else float(nodata),
            "raw_valid_fraction": float(valid.mean()),
            "raw_min": float(raw[valid].min()) if valid.any() else None,
            "raw_max": float(raw[valid].max()) if valid.any() else None}
    # EMODnet stores elevations: negative below the datum at sea, positive on land. A grid with no
    # negative value at all holds depths as positive numbers and is flipped.
    flipped = bool(valid.any() and raw[valid].min() >= 0)
    if flipped:
        raw = -raw
    info["sign_flipped"] = flipped
    clean = np.where(valid & (raw < 0), raw, np.nan).astype(np.float32)
    bathy = _reproject(clean, bbox, resampling, src_transform=src_transform, src_crs=src_crs, src_nodata=np.nan)
    return bathy, info


# ----------------------------------------------------------------------------- land / sea merge

def land_mask(dem: np.ndarray, threshold_m: float, min_cells: int) -> np.ndarray:
    land = np.isfinite(dem) & (dem > threshold_m)
    if min_cells > 1 and land.any():
        labels, n = ndimage.label(land)
        sizes = ndimage.sum(land, labels, index=np.arange(1, n + 1))
        small = np.zeros(n + 1, bool)
        small[1:] = sizes < min_cells
        land &= ~small[labels]
    return land


def nearest_fill(arr: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Every cell takes the value of the nearest valid cell."""
    if valid.all():
        return arr
    if not valid.any():
        raise ValueError("nothing to fill from")
    idx = ndimage.distance_transform_edt(~valid, return_indices=True)[1]
    return arr[tuple(idx)]


def merge_land_sea(dem: np.ndarray, bathy: np.ndarray, land: np.ndarray, cfg: dict) -> tuple[np.ndarray, dict]:
    """Land cells take the DEM. Sea cells take the EMODnet floor where it exists; where it does not
    (the 115 m grid stops short of the 30 m coastline) the floor slopes from the cap at the shoreline
    down to the nearest known depth, weighted by the two distances, and the whole sea floor gets a
    light Gaussian so the coarse source cells do not read as terraces."""
    cap = float(cfg["sea_floor_cap_m"])
    cell = float(cfg.get("cell_m", 10.0))
    sea = ~land
    dem_filled = nearest_fill(dem, np.isfinite(dem))
    sigma = float(cfg.get("land_smooth_sigma_cells", 0) or 0)
    if sigma > 0:
        dem_filled = ndimage.gaussian_filter(dem_filled, sigma)
    valid_b = np.isfinite(bathy) & sea
    if valid_b.any():
        d_valid, idx = ndimage.distance_transform_edt(~valid_b, return_indices=True)
        nearest = np.minimum(bathy[tuple(idx)], cap)
        d_shore = ndimage.distance_transform_edt(sea)       # sea cells: distance in cells to the nearest land cell
        w = d_shore / np.maximum(d_shore + d_valid, 1e-6)   # 0 at the shoreline, 1 at a known cell
        floor = np.where(valid_b, nearest, cap + (nearest - cap) * w)
    else:
        floor = np.full_like(dem, cap)
    seam = float(cfg.get("seam_cells", 0) or 0)
    if seam > 0:                                            # soften the first cells off the beach as well
        d_shore = ndimage.distance_transform_edt(sea)
        floor = cap + (floor - cap) * np.clip(d_shore / seam, 0.0, 1.0)
    sigma_sea = float(cfg.get("bathy_smooth_sigma_m", 0) or 0) / cell
    if sigma_sea > 0:
        smooth = ndimage.gaussian_filter(floor, sigma_sea, mode="nearest")
        floor = np.where(sea, np.minimum(smooth, cap), floor)
    z = np.ascontiguousarray(np.where(land, dem_filled, floor).astype(np.float32))
    assert np.isfinite(z).all()
    stats = {"z_min": float(z.min()), "z_max": float(z.max()), "land_fraction": float(land.mean()),
             "sea_cells_filled_from_neighbours": int((sea & ~valid_b).sum()),
             "dem_nan_cells": int((~np.isfinite(dem)).sum())}
    return z, stats


# ----------------------------------------------------------------------------- build / load

def build_grid(params: dict, tcfg: dict, points: list[dict], *, force: bool = False, log=print) -> dict:
    G = params["grid"]
    name = tcfg["name"]
    derived = params["paths"]["derived_dir"] / name
    derived.mkdir(parents=True, exist_ok=True)
    npz, tif, meta_path = derived / "grid.npz", derived / "z.tif", derived / "grid.json"
    bbox, lonlat = grid_box(G, points)
    cz_path = coastal_zones_path(params, name)
    cz_file = cz_path.name if cz_path.exists() else None
    if npz.exists() and meta_path.exists() and not force:
        old = json.loads(meta_path.read_text())
        same_box = all(old.get(k) == bbox[k] for k in ("crs", "cell_m", "minx", "miny", "maxx", "maxy"))
        if same_box and old.get("coastal_zones_file") == cz_file:
            log(f"have {npz.relative_to(ROOT)} (use --force to rebuild)")
            return load_grid(derived)
        log(f"the grid box changed ({old.get('width')} x {old.get('height')} -> {bbox['width']} x {bbox['height']} cells); rebuilding"
            if not same_box else f"the grid has no Coastal Zones layer and {cz_file} is there; rebuilding")
    log(f"grid {bbox['width']} x {bbox['height']} cells of {bbox['cell_m']:g} m in {bbox['crs']} "
        f"({bbox['width'] * bbox['height'] / 1e6:.2f} M cells)")
    present = dem_paths(params, name)
    if not present:
        raise SystemExit(f"no DEM windows for {name} in data/raw/dem; run scripts/00_get_dem.py first")
    dem = read_dem(present, bbox, G.get("dem_resampling", "bilinear"))
    bpath = bathy_path(params, name)
    if not bpath.exists():
        raise SystemExit(f"missing {bpath.relative_to(ROOT)}; run scripts/01_get_bathymetry.py first")
    if not covers(bpath, lonlat):
        raise SystemExit(f"{bpath.relative_to(ROOT)} does not cover the grid; run scripts/01_get_bathymetry.py again")
    bathy, binfo = read_bathy(bpath, bbox, G.get("bathy_resampling", "cubic"))
    log(f"bathymetry raw {binfo['raw_min']:.1f}..{binfo['raw_max']:.1f} m, valid {binfo['raw_valid_fraction']:.0%}"
        + (", sign flipped" if binfo["sign_flipped"] else ""))
    land = land_mask(dem, G["land_threshold_m"], G["min_land_component_cells"])
    lc_paths = landcover_paths(params, name)
    landcover = read_classes(lc_paths, bbox) if lc_paths else np.zeros_like(land, dtype=np.uint8)
    if not lc_paths:
        log("no land cover extract (scripts/02_get_landcover.py); the land keeps the hypsometric tint")
    elif G.get("landcover_as_land", False):
        # WorldCover is a 10 m land/sea mask too: it keeps polders, marshes and flats below the
        # DEM threshold on land (the Dutch, German and Danish coasts); 0 = no data, 80 = water
        extra = np.isin(landcover, [c for c in np.unique(landcover) if c not in (0, 80)]) & ~land
        land |= extra
        log(f"land cover adds {int(extra.sum())} land cells under the {G['land_threshold_m']:g} m threshold")
    coastal = read_polygons(cz_path, bbox) if cz_file else np.zeros(land.shape, np.int32)
    if cz_file and not polygons_cover(cz_path, lonlat):
        raise SystemExit(f"{cz_path.relative_to(ROOT)} does not cover the grid; run scripts/03_get_coastal_zones.py again")
    if not cz_file:
        log("no Coastal Zones extract (scripts/03_get_coastal_zones.py); the figures go without that layer")
    z, stats = merge_land_sea(dem, bathy, land, G)
    log(f"z {stats['z_min']:.1f}..{stats['z_max']:.1f} m, land {stats['land_fraction']:.0%}, "
        f"{stats['sea_cells_filled_from_neighbours']} sea cells filled from neighbours")
    lc_fractions = {}
    if lc_paths:
        codes, counts = np.unique(landcover[land], return_counts=True)
        lc_fractions = {int(c): float(n / max(1, land.sum())) for c, n in zip(codes, counts)}
        log("land cover on land: " + ", ".join(f"{c}: {f:.0%}" for c, f in sorted(lc_fractions.items(), key=lambda kv: -kv[1])))
    codes, counts = np.unique(coastal[land], return_counts=True)
    cz_fractions = {int(c): float(n / max(1, land.sum())) for c, n in zip(codes, counts) if c > 0}
    if cz_file:
        log(f"Coastal Zones: {len(np.unique(coastal[coastal > 0]))} classes, {(coastal > 0).mean():.0%} of the grid, "
            f"{sum(cz_fractions.values()):.0%} of the land; largest on land: "
            + ", ".join(f"{c}: {f:.0%}" for c, f in sorted(cz_fractions.items(), key=lambda kv: -kv[1])[:6]))

    np.savez_compressed(npz, z=z, land=land, dem=dem, bathy=bathy, landcover=landcover, coastal_zones=coastal)
    with rasterio.open(tif, "w", driver="GTiff", width=bbox["width"], height=bbox["height"], count=2, dtype="float32",
                       crs=bbox["crs"], transform=grid_transform(bbox), compress="deflate", tiled=True) as dst:
        dst.write(z, 1)
        dst.write(land.astype(np.float32), 2)
        dst.set_band_description(1, "z_m (DEM on land, EMODnet sea floor at sea)")
        dst.set_band_description(2, "land mask")
    meta = {**bbox, "transform": list(grid_transform(bbox))[:6], "lonlat_bbox": list(lonlat),
            "dem_files": [p.name for p in present],
            "bathymetry_file": bpath.name, "bathymetry": binfo, **stats,
            "landcover_files": [p.name for p in lc_paths], "landcover_fractions_on_land": lc_fractions,
            "coastal_zones_file": cz_file, "coastal_zones_fractions_on_land": cz_fractions,
            "built": time.strftime("%Y-%m-%d %H:%M")}
    meta_path.write_text(json.dumps(meta, indent=2))
    return {"z": z, "land": land, "dem": dem, "bathy": bathy, "landcover": landcover, "coastal_zones": coastal, "meta": meta}


def load_grid(derived: Path) -> dict:
    data = np.load(derived / "grid.npz")
    meta = json.loads((derived / "grid.json").read_text())
    landcover = data["landcover"] if "landcover" in data.files else np.zeros(data["z"].shape, np.uint8)
    coastal = data["coastal_zones"] if "coastal_zones" in data.files else np.zeros(data["z"].shape, np.int32)
    return {"z": data["z"], "land": data["land"].astype(bool), "dem": data["dem"], "bathy": data["bathy"],
            "landcover": landcover, "coastal_zones": coastal, "meta": meta}


def xy_to_rowcol(x, y, meta: dict) -> tuple[np.ndarray, np.ndarray]:
    """Continuous (row, col) array indices; the centre of cell (i, j) is (i, j)."""
    col = (np.asarray(x, float) - meta["minx"]) / meta["cell_m"] - 0.5
    row = (meta["maxy"] - np.asarray(y, float)) / meta["cell_m"] - 0.5
    return row, col
