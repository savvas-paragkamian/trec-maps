#!/usr/bin/env python
"""Path-trace the land-sea block with forge3d and draw the transect path and stations on it."""
import json
import math

import numpy as np
import pandas as pd

from _common import C, parser, setup
from trecmaps import compose, profile as P, render, terrain


def main() -> None:
    ap = parser(__doc__)
    ap.add_argument("--quality", default="preview", help="key of render.quality in params.yml")
    ap.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE",
                    help="override render settings, e.g. camera.azimuth_deg=200 exaggeration=4")
    ap.add_argument("--tag", default="", help="suffix for the output files, e.g. _try")
    ap.add_argument("--recompose", action="store_true", help="reuse the cached render, redo the page only")
    ap.add_argument("--no-overlay", action="store_true", help="skip the path and station overlay")
    ap.add_argument("--dpi", type=int, default=200, help="dpi of the page (changes the physical size only)")
    args = ap.parse_args()
    params, tcfg = setup(args)
    R = C.apply_overrides(params["render"], args.set)
    log = lambda m: C.log("render", m)
    derived, out = C.transect_dirs(params, tcfg["name"])
    grid = terrain.load_grid(derived)
    meta, z, land = grid["meta"], grid["z"], grid["land"]
    q = R["quality"][args.quality]
    size, frames = (int(q["size"][0]), int(q["size"][1])), int(q["frames"])
    exag = float(R["exaggeration"])
    cache = derived / f"render_{args.quality}{args.tag}.npz"
    cam_file = derived / f"camera_{args.quality}{args.tag}.json"

    if args.recompose and cache.exists():
        data = np.load(cache)
        rgba, depth = data["rgba"], data["depth"]
        camera = json.loads(cam_file.read_text())
        log(f"reusing {cache.name}")
    else:
        albedo = render.make_albedo(z, land, R["albedo"])
        log(f"{size[0]}x{size[1]} px, {frames} frames, {z.shape[1]}x{z.shape[0]} cells, exaggeration {exag:g}")
        res = render.render_plate(z, albedo, cell_m=meta["cell_m"], size_px=size, frames=frames,
                                  camera_cfg=R["camera"], sun_cfg=R["sun"], exaggeration=exag,
                                  env_intensity=R["env_intensity"], tile=int(R["tile_px"]),
                                  albedo_sampling=R["albedo_sampling"], turbidity=R["turbidity"], seed=int(R["seed"]),
                                  earth_model=R.get("earth_model"), refraction_model=R.get("refraction_model"), log=log)
        rgba, depth, camera = res["rgba"], res["depth"], res["camera"]
        np.savez_compressed(cache, rgba=rgba, depth=depth)
        cam_file.write_text(json.dumps({**camera, "size": size, "exaggeration": exag}, indent=2))
        log(f"rendered in {res['seconds']:.0f} s; misses {np.isnan(depth).mean():.1%} of the frame")

    page_style = R["page"]
    rgb = compose.vignette(compose.key_background(rgba, depth, page_style["background"]), float(page_style["vignette"]))
    page = compose.Page(rgb, dpi=args.dpi)
    W, H = size

    if not args.no_overlay:
        prof = pd.read_csv(derived / "profile.csv")
        stations = pd.read_csv(derived / "stations.csv")
        land_samples = pd.read_csv(derived / "land_samples.csv")
        lift = 2.0 * exag   # metres above the surface so the path does not z-fight with the terrain
        world = render.grid_to_world(prof.x, prof.y, prof.z, meta, exag)
        world[:, 1] += lift
        proj = render.world_to_pixel(world, camera, W, H)
        kind, ratio = render.depth_kind(proj, depth)
        log(f"depth AOV matches the {'ray distance' if kind == 'dist' else 'camera-axis depth'} "
            f"(median rendered/projected = {ratio:.4f}; 1.0000 means the overlay is calibrated)")
        vis = render.visible(proj, depth, kind)
        log(f"path: {vis.mean():.0%} of {len(vis)} points visible")
        page.polyline(proj["px"], proj["py"], vis, color=page_style["path"], width=5)

        for _, st in stations.iterrows():
            casts = P.casts_of(st.cast_depths_m)
            pts = render.grid_to_world([st.x] * (2 + len(casts)), [st.y] * (2 + len(casts)),
                                       [0.0, float(st.z_grid)] + [-d for d in casts], meta, exag)
            pr = render.world_to_pixel(pts, camera, W, H)
            v = render.visible(pr, depth, kind, rel_tol=0.03)
            sx, sy, fx, fy = pr["px"][0], pr["py"][0], pr["px"][1], pr["py"][1]
            page.segment(sx, sy, fx, fy, color=page_style["station"], width=3, ls=(0, (2, 2)), alpha=0.9 if v[1] else 0.4)
            page.marker(fx, fy, size=16, face="white", edge=page_style["station"], edge_width=2.5, shape="s")
            for i, d in enumerate(casts):
                page.marker(pr["px"][2 + i], pr["py"][2 + i], size=22, face=page_style["station"], edge="white")
            depth_txt = f"sea floor {abs(float(st.bathymetry_m)):.0f} m" if not math.isnan(float(st.bathymetry_m)) else ""
            cast_txt = ("sampled at " + ", ".join(f"{d:g}" for d in casts) + " m") if casts else ""
            label = "\n".join(t for t in (f"{st.alias}  {st.station_name}", depth_txt, cast_txt) if t)
            page.text(sx + page.px(26), sy - page.px(26), label, size=34, color=page_style["station"], va="bottom",
                      linespacing=1.15, weight="bold")

        pts = render.grid_to_world(land_samples.x, land_samples.y, land_samples.z + 1.0, meta, exag)
        pr = render.world_to_pixel(pts, camera, W, H)
        for i in range(len(pr["px"])):
            page.marker(pr["px"][i], pr["py"][i], size=16, face=page_style["sample"], edge="white", shape="v")
        if len(pr["px"]):
            i = int(np.argmax(-np.asarray(land_samples.shore_distance_m)))   # the shoreline position
            page.text(pr["px"][i] - page.px(26), pr["py"][i] - page.px(40),
                      f"{tcfg['land_site']} soil transect\n{len(land_samples)} positions, "
                      f"{land_samples.shore_distance_m.min():g}-{land_samples.shore_distance_m.max():g} m from the shoreline",
                      size=34, color=page_style["sample"], ha="right", va="bottom", linespacing=1.15, weight="bold")

    lo_sea = math.floor(float(z.min()) / 25) * 25
    hi_land = max(25.0, math.ceil(float(z[land].max()) / 25) * 25) if land.any() else 25.0
    legend = {"land": (render.legend_gradient(R["albedo"]["land"], 0.0, hi_land), 0.0, hi_land),
              "sea": (render.legend_gradient(R["albedo"]["sea"], 0.0, lo_sea), 0.0, lo_sea),
              "note": f"Vertical exaggeration x{exag:g}. Grid {meta['cell_m']:g} m, {meta['crs']}."}
    credits = [f"Land: {params['sources']['dem']['citation']}",
               f"Sea floor: {params['sources']['bathymetry']['citation']}",
               "Samples: TREC (EMBL) data portal. Rendered with forge3d."]
    compose.draw_furniture(page, title=tcfg["title"], subtitle=tcfg.get("subtitle", ""), credits=credits,
                           legend=legend, style=page_style)
    out_png = out / f"block_{args.quality}{args.tag}.png"
    page.save(out_png, out_png.with_suffix(".pdf"), title=f"{tcfg['title']} land-sea block")
    log(f"wrote {out_png} and .pdf")


if __name__ == "__main__":
    main()
