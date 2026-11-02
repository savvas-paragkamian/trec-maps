#!/usr/bin/env python
"""Fetch the EMODnet Bathymetry DTM over a transect's grid through the WCS service (network step)."""
import json
import sys
import time

from _common import C, parser, setup
from trecmaps import fetch, sites, terrain


def main() -> None:
    ap = parser(__doc__)
    ap.add_argument("--dry-run", action="store_true", help="print the request without downloading")
    args = ap.parse_args()
    params, tcfg = setup(args)
    log = lambda m: C.log("bathy", m)
    points, _, _ = sites.load_path(params, tcfg)
    _, lonlat = terrain.grid_box(params["grid"], points)
    req = terrain.wcs_request(lonlat, params["sources"]["bathymetry"])
    dest = terrain.bathy_path(params, tcfg["name"])
    log(f"{req['width']} x {req['height']} native cells, bbox {', '.join(f'{v:.5f}' for v in req['bbox'])}")
    log(req["url"])
    if dest.exists():
        if terrain.covers(dest, req["bbox"]):
            log(f"have {dest.name} and it covers the grid, skipped")
            return
        log(f"{dest.name} does not cover the current grid box; fetching it again")
        if not args.dry_run:
            dest.unlink()
    if not fetch.download(req["url"], dest, dry_run=args.dry_run, log=log):
        return
    magic = dest.read_bytes()[:4]
    if magic not in (b"II*\x00", b"MM\x00*"):
        bad = dest.with_suffix(".bad")
        dest.rename(bad)
        log(f"the response is not a TIFF (moved to {bad.name}); first bytes: {bad.read_bytes()[:400]!r}")
        sys.exit(1)
    import rasterio
    with rasterio.open(dest) as src:
        band = src.read(1)
        log(f"saved {src.width} x {src.height}, {src.crs}, nodata {src.nodata}, values {band.min():.1f}..{band.max():.1f}")
    dest.with_suffix(".json").write_text(json.dumps({**req, "fetched": time.strftime("%Y-%m-%d %H:%M"),
                                                      "citation": params["sources"]["bathymetry"]["citation"]}, indent=2))


if __name__ == "__main__":
    main()
