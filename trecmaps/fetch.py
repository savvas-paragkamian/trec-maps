"""Downloads: the only code in trecmaps that touches the network (used by scripts/0*_get_*.py)."""
from __future__ import annotations

import gzip
import json
import math
import shutil
import time
import urllib.parse
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


def get_json(url: str, params: dict, *, tries: int = 3, timeout: float = 600) -> dict:
    """GET `url?params` and parse the JSON answer (gzip accepted). An ArcGIS map service reports a failed
    query as an `error` member of a 200 answer; that and network errors are tried again after a pause."""
    req = urllib.request.Request(url + "?" + urllib.parse.urlencode(params),
                                 headers={"User-Agent": USER_AGENT, "Accept-Encoding": "gzip"})
    for attempt in range(1, tries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                body = resp.read()
                if resp.headers.get("Content-Encoding") == "gzip":
                    body = gzip.decompress(body)
            answer = json.loads(body)
            if not (isinstance(answer, dict) and "error" in answer):
                return answer
            problem = f"the service answered {answer['error']}"
        except (OSError, ValueError) as exc:               # URLError and timeouts are OSErrors
            problem = repr(exc)
        if attempt < tries:
            time.sleep(10 * attempt)
    raise SystemExit(f"{url}: {problem} ({tries} tries)")


def read_cog_window(url: str, lonlat: tuple, dest: Path, *, refresh: bool = False, dry_run: bool = False, log=print) -> bool:
    """Save the lon/lat window of a cloud-optimised GeoTIFF: GDAL fetches only the blocks it needs.
    `refresh` reads a window that exists again; the old file stays until the new one is saved."""
    import rasterio
    from rasterio.windows import Window, from_bounds, intersection

    dest = Path(dest)
    if dest.exists() and not refresh:
        log(f"have {dest.name} ({dest.stat().st_size / 2**20:.1f} MiB), skipped")
        return True
    if dry_run:
        size = head_size(url)
        log(f"would read the window {', '.join(f'{v:.4f}' for v in lonlat)} of {url}"
            + (f" ({size / 2**20:.0f} MiB tile; only the blocks under the window travel)" if size else ""))
        return False
    from rasterio.errors import RasterioIOError, WindowError

    env = rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                       GDAL_HTTP_USERAGENT=USER_AGENT)
    try:
        with env, rasterio.open("/vsicurl/" + url) as src:
            # whole pixels, rounded outward: rounded to the nearest pixel the saved window fell up to half a
            # pixel short of the box, so `terrain.window_is_short` had it fetched again on every run
            w = from_bounds(*lonlat, transform=src.transform)
            col, row = math.floor(w.col_off), math.floor(w.row_off)
            window = Window(col, row, math.ceil(w.col_off + w.width) - col, math.ceil(w.row_off + w.height) - row)
            window = intersection(window, Window(0, 0, src.width, src.height))
            if window.width <= 0 or window.height <= 0:
                log(f"the window lies outside {Path(url).name}, skipped")
                return False
            data = src.read(1, window=window)
            profile = {"driver": "GTiff", "width": data.shape[1], "height": data.shape[0], "count": 1, "dtype": data.dtype,
                       "crs": src.crs, "transform": src.window_transform(window), "nodata": src.nodata, "compress": "deflate"}
    except RasterioIOError as exc:
        if "404" not in str(exc):                          # a network problem is not a missing tile
            raise SystemExit(f"{Path(url).name}: {exc}")
        log(f"no tile {Path(url).name} (an all-ocean tile is not published), skipped")
        return False
    except WindowError:
        log(f"the window lies outside {Path(url).name}, skipped")
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    with rasterio.open(part, "w", **profile) as dst:
        dst.write(data, 1)
    part.rename(dest)
    log(f"saved {dest.name} ({data.shape[1]} x {data.shape[0]} cells)")
    return True
