#!/usr/bin/env python
"""Merge the DEM and the bathymetry onto one metric grid (data/derived/<transect>/grid.npz, z.tif, grid.json)."""
import numpy as np

from _common import C, parser, setup
from trecmaps import sites, terrain


def main() -> None:
    ap = parser(__doc__)
    ap.add_argument("--force", action="store_true", help="rebuild even when the grid exists")
    args = ap.parse_args()
    params, tcfg = setup(args)
    log = lambda m: C.log("grid", m)
    points, land, stations = sites.load_path(params, tcfg)
    for p in points:
        log(f"path point {p['name']:<12} lon {p['lon']:.5f} lat {p['lat']:.5f}")
    grid = terrain.build_grid(params, tcfg, points, force=args.force, log=log)
    meta = grid["meta"]
    # the DEM at the land site against the elevation recorded in the site export, when there is one
    site = sites.load_sites(params).loc[tcfg["land_site"]]
    if np.isfinite(site["elevation_m"]):
        from pyproj import Transformer
        x, y = Transformer.from_crs("EPSG:4326", meta["crs"], always_xy=True).transform(site["lon"], site["lat"])
        row, col = terrain.xy_to_rowcol(x, y, meta)
        r, c = int(round(float(row))), int(round(float(col)))
        if 0 <= r < meta["height"] and 0 <= c < meta["width"]:
            log(f"DEM at site {tcfg['land_site']}: {grid['dem'][r, c]:.1f} m (site export says {site['elevation_m']:.1f} m)")
    for _, st in stations.iterrows():
        from pyproj import Transformer
        x, y = Transformer.from_crs("EPSG:4326", meta["crs"], always_xy=True).transform(st.lon, st.lat)
        row, col = terrain.xy_to_rowcol(x, y, meta)
        r, c = int(round(float(row))), int(round(float(col)))
        if 0 <= r < meta["height"] and 0 <= c < meta["width"]:
            log(f"sea floor at {st.alias}: grid {grid['z'][r, c]:.1f} m, EMODnet cell {grid['bathy'][r, c]:.1f} m, "
                f"recorded {-abs(st.bathymetry_m):.1f} m")


if __name__ == "__main__":
    main()
