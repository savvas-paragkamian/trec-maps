#!/usr/bin/env python
"""Sample the grid along the transect path, place the samples and stations, draw the 2D profile."""
from _common import C, parser, setup
from trecmaps import profile as P, sites, terrain


def main() -> None:
    args = parser(__doc__).parse_args()
    params, tcfg = setup(args)
    log = lambda m: C.log("profile", m)
    derived, out = C.transect_dirs(params, tcfg["name"])
    points, land, stations = sites.load_path(params, tcfg)
    grid = terrain.load_grid(derived)
    crs = grid["meta"]["crs"]
    line, _, s0 = P.build_line(points, crs)
    prof = P.sample_line(line, grid, s0, params["profile"]["step_m"])
    land_samples = P.place_land_samples(land, prof, line, s0, crs)
    placed = P.place_stations(stations, prof, line, s0, crs)
    log(f"path {line.length / 1000:.2f} km, {len(prof)} samples every {params['profile']['step_m']:g} m; "
        f"s from {prof.s.min():.0f} to {prof.s.max():.0f} m")
    for m in P.check_profile(prof, placed):
        log(m)
    prof.to_csv(derived / "profile.csv", index=False)
    land_samples.to_csv(derived / "land_samples.csv", index=False)
    placed.to_csv(derived / "stations.csv", index=False)
    P.plot_profile(prof, land_samples, placed, params, tcfg, out / "profile.png", out / "profile.pdf")
    log(f"wrote {out / 'profile.png'} and .pdf")


if __name__ == "__main__":
    main()
