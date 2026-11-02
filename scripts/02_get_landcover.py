#!/usr/bin/env python
"""Save the ESA WorldCover 2021 window over a transect's grid (network step; GDAL range requests on the AWS tiles)."""
from _common import C, parser, setup
from trecmaps import fetch, sites, terrain


def main() -> None:
    ap = parser(__doc__)
    ap.add_argument("--dry-run", action="store_true", help="list the tiles without reading them")
    args = ap.parse_args()
    params, tcfg = setup(args)
    log = lambda m: C.log("landcover", m)
    points, _, _ = sites.load_path(params, tcfg)
    G = params["grid"]
    lonlat = terrain.bbox_lonlat(terrain.metric_bbox(points, G["crs"], G["margin_m"], G["cell_m"]))
    for lat_s, lon_s in terrain.landcover_tiles_for(params, lonlat):
        fetch.read_cog_window(terrain.landcover_url(params, lat_s, lon_s), lonlat,
                              terrain.landcover_path(params, tcfg["name"], lat_s, lon_s), dry_run=args.dry_run, log=log)


if __name__ == "__main__":
    main()
