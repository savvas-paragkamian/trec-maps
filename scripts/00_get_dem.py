#!/usr/bin/env python
"""Download the Copernicus GLO-30 tiles that cover a transect's grid (the only network step for the DEM)."""
from _common import C, parser, setup
from trecmaps import fetch, sites, terrain


def main() -> None:
    ap = parser(__doc__)
    ap.add_argument("--dry-run", action="store_true", help="list the tiles and their sizes without downloading")
    args = ap.parse_args()
    params, tcfg = setup(args)
    log = lambda m: C.log("dem", m)
    points, _, _ = sites.load_path(params, tcfg)
    G = params["grid"]
    lonlat = terrain.bbox_lonlat(terrain.metric_bbox(points, G["crs"], G["margin_m"], G["cell_m"]))
    log(f"grid box lon {lonlat[0]:.4f}..{lonlat[2]:.4f}, lat {lonlat[1]:.4f}..{lonlat[3]:.4f}")
    for lat_s, lon_s in terrain.dem_tiles_for(lonlat):
        fetch.download(terrain.dem_tile_url(params, lat_s, lon_s), terrain.dem_tile_path(params, lat_s, lon_s),
                       dry_run=args.dry_run, log=log)


if __name__ == "__main__":
    main()
