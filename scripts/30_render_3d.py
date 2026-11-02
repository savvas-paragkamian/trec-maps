#!/usr/bin/env python
"""Path-trace the land-sea block with forge3d and draw the transect path, the stations and a depth axis on it.

The block is a solid: the grid is rimmed with cells at a base depth, so the heightfield drops into cut
walls along its edges; sea level is drawn along the walls and a vertical axis stands at the south-west
corner. Depths are exaggerated `exaggeration * sea_exaggeration_factor`, land `exaggeration` (legend note).
"""
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
    ap.add_argument("--no-overlay", action="store_true", help="skip the path, station and axis overlay")
    ap.add_argument("--dpi", type=int, default=200, help="dpi of the page (changes the physical size only)")
    args = ap.parse_args()
    params, tcfg = setup(args)
    R = C.apply_overrides(params["render"], args.set)
    log = lambda m: C.log("render", m)
    derived, out = C.transect_dirs(params, tcfg["name"])
    grid = terrain.load_grid(derived)
    meta, z, land = grid["meta"], grid["z"], grid["land"]
    camera_cfg = R["camera"]
    if str(camera_cfg.get("azimuth_deg", "auto")).lower() == "auto":
        # the camera sits over the sea on the transect's seaward bearing (inland -> shoreline), snapped to a plate edge
        from pyproj import Geod
        ls = pd.read_csv(derived / "land_samples.csv").sort_values("shore_distance_m")
        shore, inland = ls.iloc[0], ls.iloc[-1]
        bearing, _, _ = Geod(ellps="WGS84").inv(float(inland.lon), float(inland.lat), float(shore.lon), float(shore.lat))
        snap = float(camera_cfg.get("azimuth_snap_deg", 90) or 0)
        az = (round(bearing / snap) * snap) % 360 if snap else bearing % 360
        camera_cfg["azimuth_deg"] = float(az)
        log(f"camera azimuth {az:.0f} deg (seaward bearing {bearing % 360:.0f} deg of the land transect)")
    if str(R["sun"].get("azimuth_deg", "auto")).lower() == "auto":
        R["sun"]["azimuth_deg"] = (float(camera_cfg["azimuth_deg"]) + float(R["sun"].get("offset_from_camera_deg", -55))) % 360
    q = R["quality"][args.quality]
    size, frames = (int(q["size"][0]), int(q["size"][1])), int(q["frames"])
    W, H = size
    exag = float(R["exaggeration"])
    sea_f = float(R.get("sea_exaggeration_factor", 1.0))
    B = R.get("block") or {}
    pad = int(B.get("pad_cells", 0) or 0)
    base_z = B.get("base_depth_m")
    if base_z is None:
        base_z = math.floor((float(z.min()) - 25.0) / 50.0) * 50.0
    base_z = float(base_z)
    cache = derived / f"render_{args.quality}{args.tag}.npz"
    cam_file = derived / f"camera_{args.quality}{args.tag}.json"

    if args.recompose and cache.exists():
        data = np.load(cache)
        rgba, depth = data["rgba"], data["depth"]
        cut = data["cut"] if "cut" in data.files else None
        saved = json.loads(cam_file.read_text())
        camera = {k: saved[k] for k in ("model", "origin", "look_at", "up", "fov_y")}
        exag, sea_f, pad, base_z = saved.get("exaggeration", exag), saved.get("sea_factor", sea_f), \
            saved.get("pad", pad), saved.get("base_z", base_z)
        if "azimuth_deg" in saved:
            camera_cfg["azimuth_deg"] = float(saved["azimuth_deg"])
        log(f"reusing {cache.name}")
    else:
        albedo = render.make_albedo(z, land, R["albedo"], grid["landcover"], params["sources"]["landcover"])
        z_render = render.exaggerate_sea(z, sea_f)
        if pad > 0:
            z_render, albedo = render.pad_block(z_render, albedo, pad, base_z * sea_f,
                                                C.hex_to_rgb(B.get("cut_face_color", "#c9bfae")))
        log(f"{W}x{H} px, {frames} frames, {z_render.shape[1]}x{z_render.shape[0]} cells, exaggeration x{exag:g} "
            f"(sea x{exag * sea_f:g}), block base {base_z:g} m")
        res = render.render_plate(z_render, albedo, cell_m=meta["cell_m"], size_px=size, frames=frames,
                                  camera_cfg=R["camera"], sun_cfg=R["sun"], exaggeration=exag,
                                  env_intensity=R["env_intensity"], tile=int(R["tile_px"]),
                                  albedo_sampling=R["albedo_sampling"], turbidity=R["turbidity"], seed=int(R["seed"]),
                                  earth_model=R.get("earth_model"), refraction_model=R.get("refraction_model"), log=log)
        rgba, depth, camera = res["rgba"], res["depth"], res["camera"]
        cut_lin = render.srgb_to_linear(C.hex_to_rgb(B.get("cut_face_color", "#c9bfae")))
        cut = np.all(np.abs(res["albedo"] - cut_lin[None, None, :]) < 0.02, axis=-1) if pad > 0 else np.zeros(depth.shape, bool)
        np.savez_compressed(cache, rgba=rgba, depth=depth, cut=cut)
        cam_file.write_text(json.dumps({**camera, "size": size, "exaggeration": exag, "sea_factor": sea_f,
                                        "pad": pad, "base_z": base_z, "azimuth_deg": float(camera_cfg["azimuth_deg"])}, indent=2))
        log(f"rendered in {res['seconds']:.0f} s; misses {np.isnan(depth).mean():.1%} of the frame")

    page_style = R["page"]
    rgb = rgba[..., :3]
    water = R.get("water") or {}
    if water.get("enabled", True):
        rgb = render.water_surface(rgb, depth, camera, meta, pad, exag, sea_f, R["sun"], water)
    rgb = compose.vignette(compose.key_background(np.dstack([rgb, rgba[..., 3]]), depth, page_style["background"]),
                           float(page_style["vignette"]))
    page = compose.Page(rgb, dpi=args.dpi)

    def project(points):
        return render.world_to_pixel(np.asarray(points, float), camera, W, H)

    if not args.no_overlay:
        prof = pd.read_csv(derived / "profile.csv")
        stations = pd.read_csv(derived / "stations.csv")
        land_samples = pd.read_csv(derived / "land_samples.csv")
        lift = 2.0 * exag   # metres above the surface so the path does not z-fight with the terrain
        world = render.grid_to_world(prof.x, prof.y, prof.z, meta, exag, sea_f)
        world[:, 1] += lift
        proj = project(world)
        kind, ratio = render.depth_kind(proj, depth)
        log(f"depth AOV matches the {'ray distance' if kind == 'dist' else 'camera-axis depth'} "
            f"(median rendered/projected = {ratio:.4f}; 1.0000 means the overlay is calibrated)")
        vis = render.visible(proj, depth, kind)
        log(f"path: {vis.mean():.0%} of {len(vis)} points visible")

        # sea level along the walls, and the vertical axis at the south-west corner
        if pad > 0:
            az = float(camera_cfg["azimuth_deg"])
            wl = project(render.front_edge(meta, pad, az))
            wl_vis = render.visible(wl, depth, kind, rel_tol=0.03, miss_ok=True)
            page.polyline(wl["px"], wl["py"], wl_vis, color=B.get("waterline_color", "#2f6fa3"), width=3.0, occluded_alpha=0.0)
            A = B.get("axis") or {}
            ax_geo = render.depth_axis(meta, pad, base_z, float(z[land].max()) if land.any() else 0.0, exag, sea_f,
                                       float(A.get("step_land_m", 100)), float(A.get("step_sea_m", 50)), azimuth_deg=az)
            ln = project(ax_geo["line"])
            axis_col = A.get("color", page_style["text"])
            page.segment(ln["px"][0], ln["py"][0], ln["px"][1], ln["py"][1], color=axis_col, width=3)
            for pt, label in ax_geo["ticks"]:
                t = project([pt])
                page.segment(t["px"][0] - page.px(12), t["py"][0], t["px"][0], t["py"][0], color=axis_col, width=2.5)
                page.text(t["px"][0] - page.px(18), t["py"][0], label, size=26, color=axis_col, ha="right", va="center")
            top = project([ax_geo["line"][1]])
            page.text(top["px"][0] - page.px(14), top["py"][0] - page.px(40),
                      f"metres above /\nbelow sea level\n(x{exag:g} land, x{exag * sea_f:g} sea)",
                      size=26, color=page_style["muted"], ha="right", va="bottom", linespacing=1.15)

        page.polyline(proj["px"], proj["py"], vis, color=page_style["path"], width=5)

        for _, st in stations.iterrows():
            casts = P.casts_of(st.cast_depths_m)
            pts = render.grid_to_world([st.x] * (2 + len(casts)), [st.y] * (2 + len(casts)),
                                       [0.0, float(st.z_grid)] + [-d for d in casts], meta, exag, sea_f)
            pr = project(pts)
            v = render.visible(pr, depth, kind, rel_tol=0.03, miss_ok=True)
            sx, sy, fx, fy = pr["px"][0], pr["py"][0], pr["px"][1], pr["py"][1]
            # the pin: sea level (hollow circle) -> sampling depth (filled) -> sea floor (square)
            page.segment(sx, sy, fx, fy, color=page_style["station"], width=3, ls=(0, (2, 2)), alpha=0.9 if v[1] else 0.4)
            page.marker(fx, fy, size=16, face="white", edge=page_style["station"], edge_width=2.5, shape="s")
            page.marker(sx, sy, size=18, face="white", edge=page_style["station"], edge_width=2.5)
            for i, d in enumerate(casts):
                page.marker(pr["px"][2 + i], pr["py"][2 + i], size=22, face=page_style["station"], edge="white")
            depth_txt = P.floor_label(st)
            cast_txt = ("sampled at " + ", ".join(f"{d:g}" for d in casts) + " m") if casts else ""
            label = "\n".join(t for t in (f"{st.alias}  {st.station_name}", cast_txt, depth_txt,
                                          f"{int(st.n_samples)} samples") if t)
            if sx > 0.72 * W:   # near the right edge the land is behind the pin: label below-left, over the water
                page.text(sx - page.px(30), sy + page.px(40), label, size=34, color=page_style["station"],
                          ha="right", va="top", linespacing=1.15, weight="bold")
            else:
                page.text(sx + page.px(30), sy - page.px(30), label, size=34, color=page_style["station"],
                          va="bottom", linespacing=1.15, weight="bold")

        pts = render.grid_to_world(land_samples.x, land_samples.y, land_samples.z + 1.0, meta, exag, sea_f)
        pr = project(pts)
        for i in range(len(pr["px"])):
            page.marker(pr["px"][i], pr["py"][i], size=16, face=page_style["sample"], edge="white", shape="v")
        if len(pr["px"]):
            i = int(np.argmax(np.asarray(land_samples.shore_distance_m)))   # the inland-most position: label above it
            page.text(pr["px"][i], pr["py"][i] - page.px(44),
                      f"{tcfg['land_site']} soil transect\n{len(land_samples)} positions, "
                      f"{land_samples.shore_distance_m.min():g}-{land_samples.shore_distance_m.max():g} m from the shore",
                      size=34, color=page_style["sample"], ha="center", va="bottom", linespacing=1.15, weight="bold")

    lo_sea = math.floor(float(z.min()) / 25) * 25
    hi_land = max(25.0, math.ceil(float(z[land].max()) / 25) * 25) if land.any() else 25.0
    note = f"Vertical exaggeration x{exag:g} on land" + (f", x{exag * sea_f:g} below sea level" if sea_f != 1 else "") \
        + f". Block base {base_z:g} m. Grid {meta['cell_m']:g} m, {meta['crs']}."
    legend = {"land": (render.legend_gradient(R["albedo"]["land"], 0.0, hi_land), 0.0, hi_land),
              "sea": (render.sea_legend_gradient(R["albedo"]["sea"], 0.0, lo_sea, water), 0.0, lo_sea),
              "note": note,
              "symbols": [("Tara station: sea level, sampling depth, sea floor", "pin", page_style["station"]),
                          ("soil / sediment sampling position", "v", page_style["sample"])]}
    LC = params["sources"]["landcover"]
    use_classes = R["albedo"].get("land_mode") == "landcover" and bool((grid["landcover"] > 0).any())
    if use_classes:
        classes = C.landcover_classes(params)
        lut = render.landcover_lut(classes, LC.get("soften_3d", [0.55, 0.92]), LC.get("water_color", "#a9c6d8"),
                                   LC.get("fallback_color", "#cfc6b4"))
        fractions = {int(k): v for k, v in meta.get("landcover_fractions_on_land", {}).items()}
        present = [c for c, f in sorted(fractions.items(), key=lambda kv: -kv[1]) if c in classes and f >= 0.01]
        legend["classes"] = [(classes[c][0], tuple(float(v) for v in lut[c])) for c in present]
        legend["classes_title"] = "Land cover (ESA WorldCover 2021)"
    credits = [f"Land: {params['sources']['dem']['citation']}"]
    if use_classes:
        credits.append(f"Land cover: {LC['citation']}")
    credits += [f"Sea floor: {params['sources']['bathymetry']['citation']}",
                "Samples: TREC (EMBL) data portal. Rendered with forge3d."]
    compose.draw_furniture(page, title=tcfg["title"], subtitle=tcfg.get("subtitle", ""), credits=credits,
                           legend=legend, style=page_style)
    out_png = out / f"block_{args.quality}{args.tag}.png"
    page.save(out_png, out_png.with_suffix(".pdf"), title=f"{tcfg['title']} land-sea block")
    log(f"wrote {out_png} and .pdf")


if __name__ == "__main__":
    main()
