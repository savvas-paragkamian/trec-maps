"""Elevation / depth profile along the transect path: sampling, sample placement, the 2D figure."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
from pyproj import Transformer
from scipy import ndimage
from shapely.geometry import LineString, Point

from .terrain import xy_to_rowcol


def build_line(points: list[dict], crs: str) -> tuple[LineString, list[tuple[float, float]], float]:
    """Metric polyline through the path points, their (x, y), and the along-path distance of the shoreline."""
    tr = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    xy = [tuple(map(float, tr.transform(p["lon"], p["lat"]))) for p in points]
    line = LineString(xy)
    s0 = 0.0
    for p, (x, y) in zip(points, xy):
        if p["kind"] == "shoreline":
            s0 = float(line.project(Point(x, y)))
    return line, xy, s0


def sample_line(line: LineString, grid: dict, s0: float, step_m: float) -> pd.DataFrame:
    """Terrain height and land flag every `step_m` along the line; `s` is 0 at the shoreline, negative inland."""
    meta = grid["meta"]
    s = np.arange(0.0, line.length, step_m)
    if line.length - s[-1] > 1e-6:
        s = np.append(s, line.length)
    pts = [line.interpolate(float(d)) for d in s]
    x = np.array([p.x for p in pts])
    y = np.array([p.y for p in pts])
    row, col = xy_to_rowcol(x, y, meta)
    z = ndimage.map_coordinates(grid["z"], [row, col], order=1, mode="nearest")
    land = ndimage.map_coordinates(grid["land"].astype(np.uint8), [row, col], order=0, mode="nearest").astype(bool)
    landcover = ndimage.map_coordinates(grid["landcover"], [row, col], order=0, mode="nearest").astype(int)
    lon, lat = Transformer.from_crs(meta["crs"], "EPSG:4326", always_xy=True).transform(x, y)
    return pd.DataFrame({"s": s - s0, "x": x, "y": y, "lon": lon, "lat": lat, "z": z, "land": land, "landcover": landcover})


def z_at(profile: pd.DataFrame, s) -> np.ndarray:
    return np.interp(np.asarray(s, float), profile["s"].to_numpy(), profile["z"].to_numpy())


def place_land_samples(land: pd.DataFrame, profile: pd.DataFrame, line: LineString, s0: float, crs: str) -> pd.DataFrame:
    """Land positions at s = -shore_distance (the TREC design distance), with the terrain height there."""
    out = land.copy()
    out["s"] = -out["shore_distance_m"].astype(float)
    out["z"] = z_at(profile, out["s"])
    tr = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    x, y = tr.transform(out["lon"].to_numpy(), out["lat"].to_numpy())
    out["x"], out["y"] = x, y
    out["s_projected"] = [float(line.project(Point(px, py))) - s0 for px, py in zip(x, y)]
    return out


def place_stations(stations: pd.DataFrame, profile: pd.DataFrame, line: LineString, s0: float, crs: str) -> pd.DataFrame:
    """Marine stations projected onto the path, with the grid's sea floor height at that point."""
    out = stations.copy()
    tr = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    x, y = tr.transform(out["lon"].to_numpy(), out["lat"].to_numpy())
    out["x"], out["y"] = x, y
    out["s"] = [float(line.project(Point(px, py))) - s0 for px, py in zip(x, y)]
    out["z_grid"] = z_at(profile, out["s"])
    out["cast_depths_m"] = out["cast_depths_m"].map(lambda v: ";".join(f"{d:g}" for d in v) if isinstance(v, (list, tuple)) else str(v))
    return out


def casts_of(value) -> list[float]:
    """Sampling depths from a list, a 'd1;d2' string or a single number (a one-cast station read back from CSV)."""
    if isinstance(value, (list, tuple)):
        return [float(v) for v in value]
    if isinstance(value, str) and value.strip():
        return [float(v) for v in value.split(";") if v.strip()]
    try:
        f = float(value)
        return [] if math.isnan(f) else [f]
    except (TypeError, ValueError):
        return []


def check_profile(profile: pd.DataFrame, stations: pd.DataFrame) -> list[str]:
    """Sanity messages: the height at the shoreline, the grid versus the recorded bathymetry, land on marine legs."""
    msgs = [f"z at the shoreline (s = 0): {float(z_at(profile, 0.0)):.1f} m"]
    for _, st in stations.iterrows():
        rec = -abs(float(st.bathymetry_m)) if not math.isnan(float(st.bathymetry_m)) else math.nan
        msgs.append(f"{st.alias} at s = {st.s / 1000:.2f} km: grid {st.z_grid:.1f} m, recorded bathymetry {rec:.1f} m")
    marine = profile[profile["s"] > 10.0]
    on_land = marine[marine["land"]]
    if len(on_land):
        runs = []
        start = prev = None
        for s in on_land["s"]:
            if start is None:
                start = prev = s
            elif s - prev > 10.0 * (profile["s"].diff().median() or 1):
                runs.append((start, prev))
                start = prev = s
            else:
                prev = s
        runs.append((start, prev))
        where = ", ".join(f"{a / 1000:.2f}-{b / 1000:.2f} km" for a, b in runs)
        msgs.append(f"WARNING: {len(on_land)} profile points on land along the marine legs ({where}); "
                    f"add a waypoint in config/transects.yml to keep the path at sea")
    else:
        msgs.append("marine legs stay at sea (no land cells)")
    return msgs


# ----------------------------------------------------------------------------- figure

def _draw_terrain(ax, s, z, land, ylo, C):
    ax.fill_between(s, np.minimum(z, 0.0), 0.0, where=z < 0, color=C["water"], lw=0, zorder=1)
    ax.fill_between(s, z, ylo, where=land, color=C["land"], lw=0, interpolate=True, zorder=2)
    ax.fill_between(s, z, ylo, where=~land, color=C["sea_floor"], lw=0, interpolate=True, zorder=2)
    ax.plot(s, z, color=C["ground_line"], lw=0.9, zorder=3)
    ax.axhline(0.0, color=C["sea_level"], lw=0.7, ls=(0, (4, 3)), zorder=3)


def _draw_zones(ax, zones, a, b, s_min, s_max, ytop, yspan, C, *, min_frac=0.08, fontsize=6, stagger=0.0):
    """Shade the zones; label those wide enough for the panel, alternating two rows when `stagger` > 0."""
    n = 0
    for zone in zones:
        z0 = s_min if zone["from"] is None else float(zone["from"])
        z1 = s_max if zone["to"] is None else float(zone["to"])
        lo, hi = max(z0, a), min(z1, b)
        if hi <= lo:
            continue
        ax.axvspan(lo, hi, color=C["zone"], alpha=0.07, lw=0, zorder=0)
        if (hi - lo) / (b - a) >= min_frac:
            y = ytop - (0.015 + stagger * (n % 2)) * yspan
            ax.text((lo + hi) / 2, y, zone["name"], ha="center", va="top", fontsize=fontsize,
                    color=C["zone"], zorder=5, linespacing=1.1)
            n += 1


def _align(s, a, b) -> str:
    frac = (s - a) / (b - a)
    return "left" if frac < 0.12 else "right" if frac > 0.88 else "center"


def plot_profile(profile: pd.DataFrame, land_samples: pd.DataFrame, stations: pd.DataFrame, params: dict,
                 tcfg: dict, out_png, out_pdf) -> None:
    """Two panels with a broken x axis (metres near the shore, kilometres offshore) sharing the y axis,
    plus a zoom inset of the land transect, where the sampling positions are a few tens of metres apart."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    P = params["profile"]
    C, F = P["colors"], P["figure"]
    muted = "#5d6670"
    plt.rcParams.update({"pdf.fonttype": 42, "font.family": "sans-serif",
                         "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
                         "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": muted,
                         "xtick.color": muted, "ytick.color": muted, "axes.labelcolor": "#20262e"})
    s, z, land = profile["s"].to_numpy(), profile["z"].to_numpy(), profile["land"].to_numpy(bool)
    step = float(P["step_m"])
    brk = float(P["break_at_m"])
    s_min, s_max = float(s.min()), float(s.max())
    ranges = [(s_min, brk), (brk, s_max)] if s_max > 2 * brk else [(s_min, s_max)]
    ylo = min(float(z.min()) * 1.12, -10.0)
    yhi = max(float(z.max()) * 1.3, 48.0)   # headroom for the station labels under the zone labels
    yspan = yhi - ylo

    fig, axes = plt.subplots(1, len(ranges), sharey=True, figsize=(F["width_in"], F["height_in"]),
                             gridspec_kw={"width_ratios": [1.0, 1.6][:len(ranges)], "wspace": 0.03})
    axes = np.atleast_1d(axes)
    for k, (ax, (a, b)) in enumerate(zip(axes, ranges)):
        m = (s >= a - 2 * step) & (s <= b + 2 * step)
        _draw_terrain(ax, s[m], z[m], land[m], ylo, C)
        _draw_zones(ax, P["zones"], a, b, s_min, s_max, yhi, yspan, C)
        for _, r in land_samples.iterrows():
            if a <= r.s <= b:
                ax.plot(r.s, r.z, marker="v", ms=4.5, color=C["sample"], mec="white", mew=0.5, ls="none", zorder=6)
        for _, st in stations.iterrows():
            if not (a <= st.s <= b):
                continue
            casts = casts_of(st.cast_depths_m)
            ax.plot([st.s, st.s], [0.0, st.z_grid], color=C["station"], lw=0.8, ls=":", zorder=4)
            ax.plot([st.s] * len(casts), [-d for d in casts], marker="o", ms=4.5, mfc=C["station"], mec="white",
                    mew=0.6, ls="none", zorder=6)
            ax.plot(st.s, st.z_grid, marker="_", ms=9, mew=1.4, color=C["station"], ls="none", zorder=6)
            depth_txt = f"sea floor {abs(float(st.bathymetry_m)):.0f} m" if not math.isnan(float(st.bathymetry_m)) else ""
            cast_txt = "sampled at " + ", ".join(f"{d:g}" for d in casts) + " m" if casts else ""
            label = "\n".join(t for t in (str(st.alias), str(st.station_name), depth_txt, cast_txt) if t)
            ha = _align(st.s, a, b)
            dx = {"left": 4, "right": -4, "center": 0}[ha]
            ax.annotate(label, (st.s, 0.0), xytext=(dx, 6), textcoords="offset points", ha=ha, va="bottom",
                        fontsize=5.8, color=C["station"], zorder=7, linespacing=1.15)
        ax.set_xlim(a, b)
        ax.set_ylim(ylo, yhi)
        if b <= brk * 1.001:
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g} m"))
        else:
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v / 1000:g} km"))
        ax.tick_params(labelsize=7)
        if k > 0:
            ax.spines["left"].set_visible(False)
            ax.tick_params(left=False)
    if len(axes) == 2:
        d = 0.012
        axes[0].plot((1 - d, 1 + d), (-d, d), transform=axes[0].transAxes, color=muted, clip_on=False, lw=0.8)
        axes[1].plot((-d, d), (-d, d), transform=axes[1].transAxes, color=muted, clip_on=False, lw=0.8)
    axes[0].set_ylabel("Elevation (m)", fontsize=8)
    sec = axes[-1].secondary_yaxis("right", functions=(lambda v: -v, lambda v: -v))
    sec.set_ylabel("Depth (m)", fontsize=8)
    sec.tick_params(labelsize=7, colors=muted)
    sec.spines["right"].set_color(muted)

    # Zoom inset of the land transect, in the empty water column of the last panel.
    zoom = P.get("land_zoom", {})
    if zoom.get("enabled", True) and len(land_samples):
        za, zb = float(zoom.get("from_m", s_min)), float(zoom.get("to_m", 120.0))
        m = (s >= za - step) & (s <= zb + step)
        ztop = max(float(z[m].max()) * 1.6, 12.0)
        zlo = -max(6.0, 0.35 * ztop)
        ins = axes[-1].inset_axes(zoom.get("box", [0.04, 0.40, 0.40, 0.46]))
        _draw_terrain(ins, s[m], z[m], land[m], zlo, C)
        _draw_zones(ins, P["zones"], za, zb, s_min, s_max, ztop, ztop - zlo, C, min_frac=0.0, fontsize=5.2, stagger=0.13)
        placements = [(0, 6, "center"), (5, 14, "left"), (-5, 14, "right")]
        for i, (_, r) in enumerate(land_samples.sort_values("s", ascending=False).iterrows()):
            if za <= r.s <= zb:
                dx, dy, ha = placements[i % len(placements)]
                ins.plot(r.s, r.z, marker="v", ms=5, color=C["sample"], mec="white", mew=0.5, ls="none", zorder=6)
                ins.annotate(f"{r.zone} {r.shore_distance_m:g} m", (r.s, r.z), xytext=(dx, dy),
                             textcoords="offset points", ha=ha, va="bottom", fontsize=5.2, color=C["sample"], zorder=7)
        ins.set_xlim(za, zb)
        ins.set_ylim(zlo, ztop)
        ins.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g} m"))
        ins.tick_params(labelsize=5.5, length=2, pad=1.5)
        ins.set_title("land transect, zoom", fontsize=6, color=muted, pad=2, loc="left")
        for sp in ("top", "right"):
            ins.spines[sp].set_visible(False)
        ins.patch.set_alpha(0.9)

    fig.supxlabel("Distance along the transect path from the shoreline (inland is negative)", fontsize=8, color=muted, y=0.055)
    fig.text(0.01, 0.975, tcfg["title"], fontsize=11, fontweight="bold", va="top", color="#20262e")
    fig.text(0.01, 0.925, tcfg.get("subtitle", ""), fontsize=7.5, va="top", color=muted)
    fig.text(0.99, 0.008, f"Land: {params['sources']['dem']['citation']}.  Sea floor: {params['sources']['bathymetry']['citation']}.",
             fontsize=5, ha="right", va="bottom", color=muted)
    fig.subplots_adjust(left=0.06, right=0.94, bottom=0.16, top=0.86, wspace=0.03)
    fig.canvas.draw()
    for ax, (a, b) in zip(axes, ranges):
        bb = ax.get_window_extent()
        ve = ((b - a) / bb.width) / (yspan / bb.height)
        ax.text(0.99, 0.02, f"vertical exaggeration x{ve:.0f}", transform=ax.transAxes, ha="right", va="bottom",
                fontsize=5.5, color=muted)
    fig.savefig(out_png, dpi=F["dpi"])
    fig.savefig(out_pdf)
    plt.close(fig)
