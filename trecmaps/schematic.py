"""The TREC sampling-concept figure drawn with real data: one land-sea section on a symmetric-log
distance axis around the shoreline, land cover along the surface, soil and sediment tubes on land,
water samples in the water column at their sampling depths, and the inland / coastal transition /
offshore band. Vertical scales differ above and below sea level (two stacked axes)."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from .config import hex_to_rgb
from .profile import casts_of, floor_label, station_title
from .render import landcover_lut

MUTED = "#5d6670"

# Marker paths (unit-ish coordinates; matplotlib scales them to the marker size, origin kept):
# a sample tube standing on its base, and a sailing boat sitting on its waterline.
TUBE = [(-0.28, 0.0), (0.28, 0.0), (0.28, 1.45), (0.18, 1.6), (-0.18, 1.6), (-0.28, 1.45), (-0.28, 0.0)]
BOAT_HULL = [(-0.55, 0.0), (0.55, 0.0), (0.38, -0.22), (-0.38, -0.22), (-0.55, 0.0)]
BOAT_SAIL = [(0.03, 0.08), (0.03, 0.95), (0.5, 0.12), (0.03, 0.08)]


def _marker(points) -> "matplotlib.path.Path":
    from matplotlib.path import Path
    return Path(points, closed=True)


def _boat() -> "matplotlib.path.Path":
    from matplotlib.path import Path
    verts = BOAT_HULL + BOAT_SAIL
    codes = [Path.MOVETO] + [Path.LINETO] * (len(BOAT_HULL) - 2) + [Path.CLOSEPOLY] \
        + [Path.MOVETO] + [Path.LINETO] * (len(BOAT_SAIL) - 2) + [Path.CLOSEPOLY]
    return Path(verts, codes)


def symlog_ticks(lo: float, hi: float) -> list[float]:
    """0 and a 1-3-10 sequence on both sides of the shoreline, within [lo, hi]."""
    vals = {0.0}
    for sign, lim in ((-1, -lo), (1, hi)):
        for e in range(1, 6):
            for m in (1, 3):
                v = m * 10 ** e
                if 10 <= v <= lim:
                    vals.add(sign * v)
    return sorted(vals)


def fmt_distance(v: float) -> str:
    a = abs(v)
    if a == 0:
        return "0"
    txt = f"{a / 1000:g} km" if a >= 1000 else f"{a:g} m"
    return ("−" if v < 0 else "") + txt


def plot_schematic(profile: pd.DataFrame, land_samples: pd.DataFrame, stations: pd.DataFrame,
                   shoreline_water: dict | None, params: dict, tcfg: dict, out_png, out_pdf) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

    S = params["schematic"]
    C, F = S["colors"], S["figure"]
    LC = params["sources"]["landcover"]
    classes = {int(k): (v[0], v[1]) for k, v in LC["classes"].items()}
    lut = landcover_lut(classes, LC.get("soften", [0.72, 0.96]), LC.get("water_color", "#a9c6d8"), C["land"])
    plt.rcParams.update({"pdf.fonttype": 42, "font.family": "sans-serif",
                         "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
                         "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": C["text"]})

    s = profile["s"].to_numpy()
    z = profile["z"].to_numpy()
    land = profile["land"].to_numpy(bool)
    lc = profile["landcover"].to_numpy(int) if "landcover" in profile.columns else np.zeros(len(s), int)
    s_min, s_max = float(s.min()), float(s.max())
    ztop = max(float(z[land].max()) * 1.6, 12.0) if land.any() else 12.0
    zlo = min(float(z.min()) * 1.12, -10.0)

    fig = plt.figure(figsize=(F["width_in"], F["height_in"]))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.7, 0.14], hspace=0.0, left=0.06, right=0.985, top=0.86, bottom=0.12)
    ax_up = fig.add_subplot(gs[0])
    ax_dn = fig.add_subplot(gs[1], sharex=ax_up)
    ax_band = fig.add_subplot(gs[2], sharex=ax_up)
    ax_up.set_xscale("symlog", linthresh=float(S["linthresh_m"]), linscale=float(S["linscale"]))
    to_axis = ax_up.xaxis.get_transform()                  # room after the last station, measured on the symlog axis
    t0, t1 = to_axis.transform(np.array([s_min - 5.0, s_max]))
    x_right = float(to_axis.inverted().transform(np.array([t1 + float(S.get("end_padding_frac", 0.0)) * (t1 - t0)]))[0])
    ax_up.set_xlim(s_min - 5.0, x_right)

    # Above sea level: the land, its cover, the soil tubes, the boats.
    zu = np.clip(z, 0.0, None)
    land_step = land | np.r_[False, land[:-1]] | np.r_[land[1:], False]   # the steps where land turns to sea too: no white wedge
    ax_up.fill_between(s, zu, 0.0, where=land_step, color=C["land"], lw=0, zorder=2)
    pts = np.column_stack([s, zu])
    segs = np.stack([pts[:-1], pts[1:]], axis=1)
    keep = land[:-1] & land[1:]
    ax_up.add_collection(LineCollection(segs[keep], colors=lut[lc[:-1]][keep], linewidths=float(S["ribbon_pt"]),
                                        capstyle="butt", zorder=3))
    ax_up.axhline(0.0, color=C["sea_level"], lw=0.9, zorder=1)
    ax_up.set_ylim(0.0, ztop)
    for sp in ("top", "right", "bottom"):
        ax_up.spines[sp].set_visible(False)
    ax_up.tick_params(bottom=False, labelbottom=False, labelsize=7)
    ax_up.set_ylabel("Elevation (m)", fontsize=8)

    # Below sea level: the water column, the ground under the sea floor, the land body on top
    # (zs is 0 on land, so the water and the sea-floor fills need no mask and leave no white wedge).
    zs = np.where(land, 0.0, np.minimum(z, 0.0))
    ax_dn.fill_between(s, zs, 0.0, color=C["water"], lw=0, zorder=1)
    ax_dn.fill_between(s, zs, zlo, color=C["sea_floor"], lw=0, zorder=2)
    ax_dn.fill_between(s, 0.0, zlo, where=land, color=C["land"], lw=0, zorder=2.5)
    ax_dn.plot(s, np.where(land, np.nan, np.minimum(z, 0.0)), color=C["ground_line"], lw=0.8, zorder=3)
    ax_dn.set_ylim(zlo, 0.0)
    for sp in ("top", "right"):
        ax_dn.spines[sp].set_visible(False)
    ax_dn.tick_params(bottom=False, labelbottom=False, labelsize=7)
    ax_dn.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(v):g}"))
    ax_dn.set_ylabel("Depth (m)", fontsize=8)

    tube = _marker(TUBE)
    boat = _boat()
    tube_pt = float(S["tube_pt"])
    placements = [(0, 4, "center"), (4, 12, "left"), (-4, 12, "right")]
    for i, (_, r) in enumerate(land_samples.sort_values("s", ascending=False).iterrows()):
        colour = C["sediment"] if "sediment" in str(r.medium).lower() else C["soil"]
        y = max(float(r.z), 0.0)
        ax_up.plot(r.s, y, marker=tube, ms=tube_pt, mfc=colour, mec=C["tube_edge"], mew=0.6, ls="none", zorder=6)
        dx, dy, ha = placements[i % len(placements)]
        if i == 0:
            ha, dx = "left", 4
        ax_up.annotate(f"{r.zone} {r.shore_distance_m:g} m\nn = {int(r.n_samples)}", (r.s, y),
                       xytext=(dx, tube_pt + dy), textcoords="offset points", ha=ha, va="bottom", fontsize=5.6,
                       color=C["text"], linespacing=1.1, zorder=7)

    if shoreline_water:
        d0, d1 = float(shoreline_water["depth_min_m"]), float(shoreline_water["depth_max_m"])
        s_sw = float(shoreline_water.get("s", 0.0))
        floor_here = float(np.interp(s_sw, s, np.where(land, 0.0, np.minimum(z, 0.0))))
        y_sw = max(-(d0 + d1) / 2, floor_here * 0.5)   # the filled floor may be shallower than the real water here
        ax_dn.plot(s_sw, y_sw, marker="v", ms=7, mfc=C["water_sample"], mec=C["water_sample_edge"],
                   mew=0.8, ls="none", zorder=6)
        ax_dn.annotate(f"shoreline water {d0:g}-{d1:g} m\nn = {shoreline_water['n_samples']}", (s_sw, y_sw),
                       xytext=(6, -2), textcoords="offset points", ha="left", va="top", fontsize=5.6, color=C["text"],
                       linespacing=1.1, zorder=7)

    fig.canvas.draw()

    def align(ax, x):
        x0, x1 = ax.get_window_extent().x0, ax.get_window_extent().x1
        f = (ax.transData.transform((x, 0.0))[0] - x0) / (x1 - x0)
        return ("left", 4) if f < 0.08 else ("right", -4) if f > 0.9 else ("center", 0)

    for _, st in stations.iterrows():
        casts = casts_of(st.cast_depths_m)
        ax_dn.plot([st.s, st.s], [0.0, float(st.z_grid)], color=C["cast_line"], lw=0.8, ls=(0, (3, 3)), zorder=4)
        ax_dn.plot([st.s] * len(casts), [-d for d in casts], marker="v", ms=7, mfc=C["water_sample"],
                   mec=C["water_sample_edge"], mew=0.8, ls="none", zorder=6)
        ax_dn.plot(st.s, float(st.z_grid), marker="_", ms=10, mew=1.4, color=C["cast_line"], ls="none", zorder=6)
        ax_up.plot(st.s, 0.0, marker=boat, ms=22, mfc=C["boat"], mec=C["boat"], mew=0.4, ls="none", zorder=6)
        depth_txt = floor_label(st)
        cast_txt = "sampled at " + ", ".join(f"{d:g}" for d in casts) + " m" if casts else ""
        label = "\n".join(t for t in (station_title(st), cast_txt, depth_txt, f"n = {int(st.n_samples)}") if t)
        ha, dx = align(ax_up, st.s)
        ax_up.annotate(label, (st.s, 0.0), xytext=(dx, 24), textcoords="offset points", ha=ha, va="bottom",
                       fontsize=5.8, color=C["text"], linespacing=1.15, zorder=7)
        for d in casts:
            ha_c = "right" if ha == "right" else "left"
            ax_dn.annotate(f"{d:g} m", (st.s, -d), xytext=(-7 if ha_c == "right" else 7, 0), textcoords="offset points",
                           ha=ha_c, va="center", fontsize=5.6, color=C["cast_line"], zorder=7)

    # The band: inland / coastal transition zone / offshore, labelled at the visual midpoints.
    ticks = symlog_ticks(s_min, s_max)
    ax_band.set_ylim(0, 1)
    ax_band.set_yticks([])
    for sp in ("top", "right", "left"):
        ax_band.spines[sp].set_visible(False)
    fig.canvas.draw()
    for zone in S["zones"]:
        a = s_min - 5.0 if zone["from"] is None else float(zone["from"])
        b = x_right if zone["to"] is None else float(zone["to"])
        ax_band.axvspan(a, b, color=zone["color"], lw=0, zorder=0)
        xa, xb = ax_band.transData.transform([(a, 0.5), (b, 0.5)])[:, 0]
        xm = ax_band.transData.inverted().transform([((xa + xb) / 2, 0.0)])[0, 0]
        ax_band.text(xm, 0.5, zone["name"], ha="center", va="center", fontsize=7, color=C["text"], fontweight="bold")
    ax_band.xaxis.set_major_locator(FixedLocator(ticks))
    ax_band.xaxis.set_minor_locator(NullLocator())
    ax_band.xaxis.set_major_formatter(FuncFormatter(lambda v, _: fmt_distance(v)))
    ax_band.tick_params(labelsize=6.5, length=3)
    ax_band.set_xlabel("Distance from the shoreline along the transect path (symmetric-log scale; inland negative)",
                       fontsize=7.5, color=MUTED)

    # Legend in the sky over the sea: land cover classes along the land, then the sample symbols.
    present = [c for c, _ in sorted(((int(c), n) for c, n in zip(*np.unique(lc[land], return_counts=True))),
                                    key=lambda cn: -cn[1]) if c in classes] if land.any() else []
    handles = [Patch(facecolor=lut[c], edgecolor="none", label=classes[c][0]) for c in present]
    handles += [
        Line2D([], [], marker=tube, ms=9, mfc=C["soil"], mec=C["tube_edge"], mew=0.6, ls="none", label="topsoil samples"),
        Line2D([], [], marker=tube, ms=9, mfc=C["sediment"], mec=C["tube_edge"], mew=0.6, ls="none", label="superficial sediments"),
        Line2D([], [], marker="v", ms=7, mfc=C["water_sample"], mec=C["water_sample_edge"], mew=0.8, ls="none",
               label="water column samples (sampling depth)"),
        Line2D([], [], marker=boat, ms=14, mfc=C["boat"], mec=C["boat"], ls="none", label="Tara station"),
    ]
    ax_up.legend(handles=handles, loc="upper right", fontsize=6, frameon=False, ncol=2, handletextpad=0.6,
                 columnspacing=1.4, borderaxespad=0.2, title="ESA WorldCover 2021 along the land  |  samples",
                 title_fontsize=6.2, alignment="left")

    fig.text(0.01, 0.975, tcfg["title"], fontsize=11, fontweight="bold", va="top", color=C["text"])
    fig.text(0.01, 0.925, tcfg.get("subtitle", ""), fontsize=7.5, va="top", color=MUTED)
    fig.text(0.985, 0.008, "Vertical scales differ above and below sea level. Land: " + params["sources"]["dem"]["citation"]
             + ". Sea floor: " + params["sources"]["bathymetry"]["citation"] + ". Land cover: " + LC["citation"] + ".",
             fontsize=4.8, ha="right", va="bottom", color=MUTED, wrap=True)
    fig.savefig(out_png, dpi=F["dpi"])
    fig.savefig(out_pdf)
    plt.close(fig)
