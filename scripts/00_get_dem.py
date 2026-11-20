#!/usr/bin/env python
"""Save the Copernicus GLO-30 windows over a transect's grid (network step; GDAL range requests on the AWS COG tiles)."""
from _common import C, parser, setup
from trecmaps import fetch, sites, terrain


def main() -> None:
    ap = parser(__doc__)
    ap.add_argument("--dry-run", action="store_true", help="list the tiles without reading them")
    args = ap.parse_args()
    params, tcfg = setup(args)
    log = lambda m: C.log("dem", m)
    points, _, _ = sites.load_path(params, tcfg)
    bbox, lonlat = terrain.grid_box(params["grid"], points)
    log(f"grid box lon {lonlat[0]:.4f}..{lonlat[2]:.4f}, lat {lonlat[1]:.4f}..{lonlat[3]:.4f} ({bbox['crs']})")
    for lat_s, lon_s in terrain.dem_tiles_for(lonlat):
        dest = terrain.dem_path(params, tcfg["name"], lat_s, lon_s)
        short = terrain.window_is_short(dest, lonlat, terrain.tile_extent(lat_s, lon_s, 1))
        if short:
            log(f"{dest.name} does not cover the current grid box")
        fetch.read_cog_window(terrain.dem_tile_url(params, lat_s, lon_s), lonlat, dest, refresh=short,
                              dry_run=args.dry_run, log=log)


if __name__ == "__main__":
    main()
