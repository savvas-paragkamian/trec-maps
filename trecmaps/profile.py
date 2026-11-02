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


def floor_label(st) -> str:
    """'sea floor 131 m', or both numbers when the grid's floor disagrees with the recorded bathymetry."""
    rec = abs(float(st.bathymetry_m)) if not math.isnan(float(st.bathymetry_m)) else math.nan
    grid = abs(float(st.z_grid)) if "z_grid" in st and not math.isnan(float(st.z_grid)) else math.nan
    if math.isnan(rec) and math.isnan(grid):
        return ""
    if math.isnan(rec):
        return f"sea floor {grid:.0f} m (EMODnet)"
    if not math.isnan(grid) and abs(grid - rec) > max(5.0, 0.25 * rec):
        return f"sea floor {rec:.0f} m recorded, {grid:.0f} m in EMODnet here"
    return f"sea floor {rec:.0f} m"


def check_profile(profile: pd.DataFrame, stations: pd.DataFrame) -> list[str]:
    """Sanity messages: the height at the shoreline, the grid versus the recorded bathymetry, land on marine legs."""
    msgs = [f"z at the shoreline (s = 0): {float(z_at(profile, 0.0)):.1f} m"]
    for _, st in stations.iterrows():
        rec = -abs(float(st.bathymetry_m)) if not math.isnan(float(st.bathymetry_m)) else math.nan
        msgs.append(f"{st.alias} at s = {st.s / 1000:.2f} km: grid {st.z_grid:.1f} m, recorded bathymetry {rec:.1f} m")
    marine = profile[profile["s"] > 50.0]          # the first 50 m off the shoreline event may still be beach in a 30 m DEM
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
    """Three panels on one distance line: the land transect in detail with its own vertical scale, then
    the nearshore (metres) and the offshore (kilometres) panels sharing one elevation/depth axis."""
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
    zoom = P.get("land_zoom", {})
    use_zoom = bool(zoom.get("enabled", True)) and len(land_samples) > 0
    za, zb = float(zoom.get("from_m", s_min)), float(zoom.get("to_m", 120.0))
    start = zb if use_zoom else s_min
    ranges = [(start, brk), (brk, s_max)] if s_max > 2 * brk else [(start, s_max)]
    ylo = min(float(z.min()) * 1.12, -10.0)
    in_main = s >= start
    z_main_max = float(z[in_main].max()) if in_main.any() else 0.0
    yhi = max(z_main_max * 1.3, 0.3 * -ylo, 5.0)   # headroom for the station and zone labels
    yspan = yhi - ylo

    widths = ([0.9, 0.12] if use_zoom else []) + [1.0, 1.5][:len(ranges)]     # a spacer column after the detail panel
    fig = plt.figure(figsize=(F["width_in"], F["height_in"]))
    gs = fig.add_gridspec(1, len(widths), width_ratios=widths, wspace=0.03,
                          left=0.06, right=0.94, bottom=0.16, top=0.86)
    first = 2 if use_zoom else 0
    axes = [fig.add_subplot(gs[first])]
    for k in range(1, len(ranges)):
        axes.append(fig.add_subplot(gs[first + k], sharey=axes[0]))

    def draw_stations(ax, a, b):
        for _, st in stations.iterrows():
            if not (a <= st.s <= b):
                continue
            casts = casts_of(st.cast_depths_m)
            ax.plot([st.s, st.s], [0.0, st.z_grid], color=C["station"], lw=0.8, ls=":", zorder=4)
            ax.plot([st.s] * len(casts), [-d for d in casts], marker="o", ms=4.5, mfc=C["station"], mec="white",
                    mew=0.6, ls="none", zorder=6)
            ax.plot(st.s, st.z_grid, marker="_", ms=9, mew=1.4, color=C["station"], ls="none", zorder=6)
            depth_txt = floor_label(st)
            cast_txt = "sampled at " + ", ".join(f"{d:g}" for d in casts) + " m" if casts else ""
            label = "\n".join(t for t in (str(st.alias), str(st.station_name), depth_txt, cast_txt) if t)
            ha = _align(st.s, a, b)
            dx = {"left": 4, "right": -4, "center": 0}[ha]
            ax.annotate(label, (st.s, 0.0), xytext=(dx, 6), textcoords="offset points", ha=ha, va="bottom",
                        fontsize=5.8, color=C["station"], zorder=7, linespacing=1.15)

    for k, (ax, (a, b)) in enumerate(zip(axes, ranges)):
        m = (s >= a - 2 * step) & (s <= b + 2 * step)
        _draw_terrain(ax, s[m], z[m], land[m], ylo, C)
        _draw_zones(ax, P["zones"], a, b, s_min, s_max, yhi, yspan, C)
        for _, r in land_samples.iterrows():
            if a <= r.s <= b:
                ax.plot(r.s, r.z, marker="v", ms=4.5, color=C["sample"], mec="white", mew=0.5, ls="none", zorder=6)
        draw_stations(ax, a, b)
        ax.set_xlim(a, b)
        ax.set_ylim(ylo, yhi)
        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(7, prune="upper" if k < len(ranges) - 1 else "lower"))
        if b <= brk * 1.001:
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g} m"))
        else:
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v / 1000:g} km"))
        ax.tick_params(labelsize=7)
        if k > 0:
            ax.spines["left"].set_visible(False)
            ax.tick_params(left=False, labelleft=False)
    axes[0].set_ylabel("Elevation (m)", fontsize=8)
    sec = axes[-1].secondary_yaxis("right", functions=(lambda v: -v, lambda v: -v))
    sec.set_ylabel("Depth (m)", fontsize=8)
    sec.tick_params(labelsize=7, colors=muted)
    sec.spines["right"].set_color(muted)

    panels = list(zip(axes, ranges))
    if use_zoom:
        ax0 = fig.add_subplot(gs[0])
        m = (s >= za - step) & (s <= zb + step)
        ztop = max(float(z[m].max()) * 1.6, 12.0)
        zlo = -max(6.0, 0.35 * ztop)
        _draw_terrain(ax0, s[m], z[m], land[m], zlo, C)
        _draw_zones(ax0, P["zones"], za, zb, s_min, s_max, ztop, ztop - zlo, C, min_frac=0.0, fontsize=5.6, stagger=0.11)
        placements = [(0, 6, "center"), (5, 14, "left"), (-5, 14, "right")]
        for i, (_, r) in enumerate(land_samples.sort_values("s", ascending=False).iterrows()):
            if za <= r.s <= zb:
                dx, dy, ha = placements[i % len(placements)]
                ax0.plot(r.s, r.z, marker="v", ms=5, color=C["sample"], mec="white", mew=0.5, ls="none", zorder=6)
                ax0.annotate(f"{r.zone} {r.shore_distance_m:g} m", (r.s, r.z), xytext=(dx, dy), textcoords="offset points",
                             ha=ha, va="bottom", fontsize=5.6, color=C["sample"], zorder=7)
        draw_stations(ax0, za, zb)
        ax0.set_xlim(za, zb)
        ax0.set_ylim(zlo, ztop)
        from matplotlib.ticker import MaxNLocator
        ax0.xaxis.set_major_locator(MaxNLocator(5))
        ax0.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g} m"))
        ax0.tick_params(labelsize=7)
        ax0.set_ylabel("Elevation (m), detail scale", fontsize=8)
        ax0.set_title("land transect (own vertical scale)", fontsize=7, color=muted, loc="left", pad=3)
        panels.insert(0, (ax0, (za, zb)))

    # axis-break marks between neighbouring panels
    d = 0.012
    for (left, _), (right, _) in zip(panels[:-1], panels[1:]):
        left.plot((1 - d, 1 + d), (-d, d), transform=left.transAxes, color=muted, clip_on=False, lw=0.8)
        right.plot((-d, d), (-d, d), transform=right.transAxes, color=muted, clip_on=False, lw=0.8)

    fig.supxlabel("Distance along the transect path from the shoreline (inland is negative)", fontsize=8, color=muted, y=0.055)
    fig.text(0.01, 0.975, tcfg["title"], fontsize=11, fontweight="bold", va="top", color="#20262e")
    fig.text(0.01, 0.925, tcfg.get("subtitle", ""), fontsize=7.5, va="top", color=muted)
    fig.text(0.99, 0.008, f"Land: {params['sources']['dem']['citation']}.  Sea floor: {params['sources']['bathymetry']['citation']}.",
             fontsize=5, ha="right", va="bottom", color=muted)
    fig.canvas.draw()
    for ax, (a, b) in panels:
        bb = ax.get_window_extent()
        y0, y1 = ax.get_ylim()
        ve = ((b - a) / bb.width) / ((y1 - y0) / bb.height)
        ax.text(0.99, 0.02, f"vertical exaggeration x{ve:.0f}", transform=ax.transAxes, ha="right", va="bottom",
                fontsize=5.5, color=muted)
    fig.savefig(out_png, dpi=F["dpi"])
    fig.savefig(out_pdf)
    plt.close(fig)
