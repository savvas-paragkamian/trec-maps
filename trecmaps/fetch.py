"""Downloads: the only code in trecmaps that touches the network (used by scripts/0*_get_*.py)."""
from __future__ import annotations

import shutil
import urllib.request
from pathlib import Path

USER_AGENT = "trec-maps/0.1 (https://github.com/savvas-paragkamian; TREC transect maps)"


def head_size(url: str) -> int | None:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            size = resp.headers.get("Content-Length")
            return int(size) if size else None
    except Exception:
        return None


def download(url: str, dest: Path, *, dry_run: bool = False, log=print) -> bool:
    """Stream `url` to `dest` via `dest.part`; skip when `dest` exists. Returns True when the file is present."""
    dest = Path(dest)
    if dest.exists():
        log(f"have {dest.name} ({dest.stat().st_size / 2**20:.1f} MiB), skipped")
        return True
    if dry_run:
        size = head_size(url)
        log(f"would download {url}" + (f" ({size / 2**20:.1f} MiB)" if size else ""))
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=300) as resp, open(part, "wb") as f:
            shutil.copyfileobj(resp, f, length=1 << 20)
    except urllib.error.HTTPError as exc:
        part.unlink(missing_ok=True)
        if exc.code == 404:
            log(f"404 for {url} (an all-ocean DEM tile is not published), skipped")
            return False
        raise
    part.rename(dest)
    log(f"saved {dest.name} ({dest.stat().st_size / 2**20:.1f} MiB)")
    return True


def read_cog_window(url: str, lonlat: tuple, dest: Path, *, dry_run: bool = False, log=print) -> bool:
    """Save the lon/lat window of a cloud-optimised GeoTIFF: GDAL fetches only the blocks it needs."""
    import rasterio
    from rasterio.windows import Window, from_bounds, intersection

    dest = Path(dest)
    if dest.exists():
        log(f"have {dest.name} ({dest.stat().st_size / 2**20:.1f} MiB), skipped")
        return True
    if dry_run:
        size = head_size(url)
        log(f"would read the window {', '.join(f'{v:.4f}' for v in lonlat)} of {url}"
            + (f" ({size / 2**20:.0f} MiB tile; only the blocks under the window travel)" if size else ""))
        return False
    env = rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                       GDAL_HTTP_USERAGENT=USER_AGENT)
    with env, rasterio.open("/vsicurl/" + url) as src:
        window = from_bounds(*lonlat, transform=src.transform).round_offsets().round_lengths()
        window = intersection(window, Window(0, 0, src.width, src.height))
        if window.width <= 0 or window.height <= 0:
            log(f"the window lies outside {Path(url).name}, skipped")
            return False
        data = src.read(1, window=window)
        profile = {"driver": "GTiff", "width": data.shape[1], "height": data.shape[0], "count": 1, "dtype": data.dtype,
                   "crs": src.crs, "transform": src.window_transform(window), "nodata": src.nodata, "compress": "deflate"}
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    with rasterio.open(part, "w", **profile) as dst:
        dst.write(data, 1)
    part.rename(dest)
    log(f"saved {dest.name} ({data.shape[1]} x {data.shape[0]} cells)")
    return True
