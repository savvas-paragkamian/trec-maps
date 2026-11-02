#!/usr/bin/env python
"""Draw the TREC sampling-concept schematic with the real transect data (needs 20_profile_2d.py's tables)."""
import json

import pandas as pd
from pyproj import Transformer
from shapely.geometry import Point

from _common import C, parser, setup
from trecmaps import profile as P, schematic, sites, terrain


def main() -> None:
    args = parser(__doc__).parse_args()
    params, tcfg = setup(args)
    log = lambda m: C.log("schematic", m)
    derived, out = C.transect_dirs(params, tcfg["name"])
    prof = pd.read_csv(derived / "profile.csv")
    land_samples = pd.read_csv(derived / "land_samples.csv")
    stations = pd.read_csv(derived / "stations.csv")
    if "landcover" not in prof.columns or (prof["landcover"] == 0).all():
        log("profile.csv has no land cover; run 02_get_landcover.py, 10_build_grid.py --force and 20_profile_2d.py")
    water = sites.load_shoreline_water(params, tcfg["land_site"])
    if water:
        path_file = derived / "path.json"                       # the path as 20_profile_2d.py routed it
        points = json.loads(path_file.read_text()) if path_file.exists() else sites.load_path(params, tcfg)[0]
        crs = terrain.load_grid(derived)["meta"]["crs"]
        line, _, s0 = P.build_line(points, crs)
        x, y = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(water["lon"], water["lat"])
        water["s"] = float(line.project(Point(x, y))) - s0
        log(f"shoreline water: {water['n_samples']} samples, {water['depth_min_m']:g}-{water['depth_max_m']:g} m, "
            f"at s = {water['s']:.0f} m")
    schematic.plot_schematic(prof, land_samples, stations, water, params, tcfg, out / "schematic.png", out / "schematic.pdf")
    log(f"wrote {out / 'schematic.png'} and .pdf")


if __name__ == "__main__":
    main()
