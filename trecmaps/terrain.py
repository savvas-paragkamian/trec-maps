"""One metric grid for a transect: Copernicus GLO-30 on land, EMODnet bathymetry at sea.

Nothing here goes online. The DEM tiles and the EMODnet extract are fetched by
scripts/00_get_dem.py and scripts/01_get_bathymetry.py; this module only decides which
files they need (`dem_tiles_for`, `wcs_request`) and merges them (`build_grid`).
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

from .config import ROOT

RESAMPLING = {"nearest": Resampling.nearest, "bilinear": Resampling.bilinear,
              "cubic": Resampling.cubic, "average": Resampling.average}
SIDES = ("north", "south", "east", "west")


# ----------------------------------------------------------------------------- grid definition

def metric_bbox(points: list[dict], crs: str, margin, cell_m: float) -> dict:
    """Bounding box of the path in the metric CRS plus margins, snapped outward to whole cells."""
    tr = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    xs, ys = tr.transform([p["lon"] for p in points], [p["lat"] for p in points])
    m = {k: float(margin.get(k, 0)) for k in SIDES} if isinstance(margin, dict) else {k: float(margin) for k in SIDES}
    minx = math.floor((min(xs) - m["west"]) / cell_m) * cell_m
    maxx = math.ceil((max(xs) + m["east"]) / cell_m) * cell_m
    miny = math.floor((min(ys) - m["south"]) / cell_m) * cell_m
    maxy = math.ceil((max(ys) + m["north"]) / cell_m) * cell_m
    return {"crs": crs, "cell_m": float(cell_m), "minx": minx, "miny": miny, "maxx": maxx, "maxy": maxy,
            "width": int(round((maxx - minx) / cell_m)), "height": int(round((maxy - miny) / cell_m))}


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

def dem_tiles_for(lonlat: tuple) -> list[tuple[str, str]]:
    """GLO-30 tile names (south-west corner, e.g. ('N37', 'E023')) covering a lon/lat box."""
    minlon, minlat, maxlon, maxlat = lonlat
    tiles = []
    for lat in range(math.floor(minlat), math.floor(maxlat) + 1):
        for lon in range(math.floor(minlon), math.floor(maxlon) + 1):
            tiles.append((f"{'N' if lat >= 0 else 'S'}{abs(lat):02d}", f"{'E' if lon >= 0 else 'W'}{abs(lon):03d}"))
    return tiles


def dem_tile_url(params: dict, lat_s: str, lon_s: str) -> str:
    return params["sources"]["dem"]["url"].format(lat=lat_s, lon=lon_s)


def dem_tile_path(params: dict, lat_s: str, lon_s: str) -> Path:
    return params["paths"]["raw_dir"] / "dem" / Path(dem_tile_url(params, lat_s, lon_s)).name


def bathy_path(params: dict, name: str) -> Path:
    return params["paths"]["raw_dir"] / "bathy" / f"emodnet_mean_{name}.tif"


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


def covers(raster_path: Path, lonlat: tuple, tol_deg: float = 1e-6) -> bool:
    """Does the EPSG:4326 raster's extent contain the (minlon, minlat, maxlon, maxlat) box?"""
    with rasterio.open(raster_path) as src:
        b = src.bounds
    return (b.left <= lonlat[0] + tol_deg and b.bottom <= lonlat[1] + tol_deg
            and b.right >= lonlat[2] - tol_deg and b.top >= lonlat[3] - tol_deg)


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
    """Land cells take the DEM, sea cells the sea floor, capped below the surface and ramped to the
    cap over `seam_cells` off the shoreline so the 115 m bathymetry meets the 30 m coastline."""
    cap = float(cfg["sea_floor_cap_m"])
    seam = float(cfg.get("seam_cells", 0) or 0)
    sea = ~land
    dem_filled = nearest_fill(dem, np.isfinite(dem))
    sigma = float(cfg.get("land_smooth_sigma_cells", 0) or 0)
    if sigma > 0:
        dem_filled = ndimage.gaussian_filter(dem_filled, sigma)
    valid_b = np.isfinite(bathy) & sea
    floor = nearest_fill(bathy, valid_b) if valid_b.any() else np.full_like(dem, cap)
    floor = np.minimum(floor, cap)
    if seam > 0:
        dist = ndimage.distance_transform_edt(sea)      # sea cells: distance in cells to the nearest land cell
        floor = cap + (floor - cap) * np.clip(dist / seam, 0.0, 1.0)
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
    if npz.exists() and meta_path.exists() and not force:
        log(f"have {npz.relative_to(ROOT)} (use --force to rebuild)")
        return load_grid(derived)

    bbox = metric_bbox(points, G["crs"], G["margin_m"], G["cell_m"])
    log(f"grid {bbox['width']} x {bbox['height']} cells of {bbox['cell_m']:g} m in {bbox['crs']} "
        f"({bbox['width'] * bbox['height'] / 1e6:.2f} M cells)")
    lonlat = bbox_lonlat(bbox)
    tile_paths = [dem_tile_path(params, la, lo) for la, lo in dem_tiles_for(lonlat)]
    present = [p for p in tile_paths if p.exists()]
    if not present:
        raise SystemExit("no DEM tiles in data/raw/dem; run scripts/00_get_dem.py first")
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
    z, stats = merge_land_sea(dem, bathy, land, G)
    log(f"z {stats['z_min']:.1f}..{stats['z_max']:.1f} m, land {stats['land_fraction']:.0%}, "
        f"{stats['sea_cells_filled_from_neighbours']} sea cells filled from neighbours")

    np.savez_compressed(npz, z=z, land=land, dem=dem, bathy=bathy)
    with rasterio.open(tif, "w", driver="GTiff", width=bbox["width"], height=bbox["height"], count=2, dtype="float32",
                       crs=bbox["crs"], transform=grid_transform(bbox), compress="deflate", tiled=True) as dst:
        dst.write(z, 1)
        dst.write(land.astype(np.float32), 2)
        dst.set_band_description(1, "z_m (DEM on land, EMODnet sea floor at sea)")
        dst.set_band_description(2, "land mask")
    meta = {**bbox, "transform": list(grid_transform(bbox))[:6], "lonlat_bbox": list(lonlat),
            "dem_tiles": [p.name for p in present], "dem_tiles_missing": [p.name for p in tile_paths if not p.exists()],
            "bathymetry_file": bpath.name, "bathymetry": binfo, **stats, "built": time.strftime("%Y-%m-%d %H:%M")}
    meta_path.write_text(json.dumps(meta, indent=2))
    return {"z": z, "land": land, "dem": dem, "bathy": bathy, "meta": meta}


def load_grid(derived: Path) -> dict:
    data = np.load(derived / "grid.npz")
    meta = json.loads((derived / "grid.json").read_text())
    return {"z": data["z"], "land": data["land"].astype(bool), "dem": data["dem"], "bathy": data["bathy"], "meta": meta}


def xy_to_rowcol(x, y, meta: dict) -> tuple[np.ndarray, np.ndarray]:
    """Continuous (row, col) array indices; the centre of cell (i, j) is (i, j)."""
    col = (np.asarray(x, float) - meta["minx"]) / meta["cell_m"] - 0.5
    row = (meta["maxy"] - np.asarray(y, float)) / meta["cell_m"] - 0.5
    return row, col
