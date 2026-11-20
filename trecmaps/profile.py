"""Elevation / depth profile along the transect path: sampling, sample placement, the 2D figure."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
from pyproj import Transformer
from scipy import ndimage
from shapely.geometry import LineString, Point

from .render import soften_hex
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


def route_around_land(points: list[dict], grid: dict, cfg: dict, log=print) -> list[dict]:
    """Marine legs that cross land are replaced by the shortest route through sea cells: Dijkstra on a
    coarse copy of the land mask (`cell_m`) dilated by `land_buffer_m`; endpoints on land (the beach, a
    station on a flat) snap to the nearest free cell within `snap_m`. A leg keeps its straight line when
    no route exists inside the grid or the route is longer than `max_detour` x the straight leg."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import dijkstra

    meta, land = grid["meta"], grid["land"]
    cell = float(meta["cell_m"])
    f = max(1, int(round(float(cfg.get("cell_m", 50)) / cell)))
    H, W = land.shape
    h, w = H // f, W // f
    coarse = land[:h * f, :w * f].reshape(h, f, w, f).any(axis=(1, 3))
    buf = int(math.ceil(float(cfg.get("land_buffer_m", 100)) / (f * cell)))
    tr = Transformer.from_crs("EPSG:4326", meta["crs"], always_xy=True)
    inv = Transformer.from_crs(meta["crs"], "EPSG:4326", always_xy=True)
    xy = [tuple(map(float, tr.transform(p["lon"], p["lat"]))) for p in points]
    idx = np.arange(h * w).reshape(h, w)
    snap_cells = float(cfg.get("snap_m", 1500)) / (f * cell)

    def sea_graph(buffer_cells):
        """(free mask, sparse 8-connected graph of the free cells) with the land dilated by `buffer_cells`."""
        free = ~(ndimage.binary_dilation(coarse, iterations=buffer_cells) if buffer_cells > 0 else coarse)
        src, dst, wts = [], [], []
        for dr, dc, wgt in ((0, 1, 1.0), (1, 0, 1.0), (1, 1, math.sqrt(2)), (1, -1, math.sqrt(2))):
            a = idx[0:h - dr, max(0, -dc):w - max(0, dc)]
            b = idx[dr:h, max(0, -dc) + dc:w - max(0, dc) + dc]
            ok = free.ravel()[a] & free.ravel()[b]
            src += [a[ok], b[ok]]
            dst += [b[ok], a[ok]]
            wts += [np.full(int(ok.sum()), wgt * f * cell)] * 2
        graph = coo_matrix((np.concatenate(wts), (np.concatenate(src), np.concatenate(dst))), shape=(h * w, h * w)).tocsr()
        return free, graph

    graphs = {}                                              # built on demand: the buffered mask first, the bare one as a retry

    def cell_of(x, y):
        r, c = xy_to_rowcol(np.array([x]), np.array([y]), meta)
        return int(min(max(r[0] // f, 0), h - 1)), int(min(max(c[0] // f, 0), w - 1))

    def snap(free, r0, c0):
        if free[r0, c0]:
            return r0, c0
        free_rc = np.argwhere(free)
        d2 = (free_rc[:, 0] - r0) ** 2 + (free_rc[:, 1] - c0) ** 2
        k = int(np.argmin(d2))
        return (int(free_rc[k, 0]), int(free_rc[k, 1])) if d2[k] <= snap_cells ** 2 else None

    def shortest(a, b):
        """Coarse-cell chain of the shortest sea route, or (None, reason)."""
        reason = "no sea cell near an endpoint"
        for buffer_cells in (buf, 0) if buf > 0 else (0,):
            if buffer_cells not in graphs:
                graphs[buffer_cells] = sea_graph(buffer_cells)
            free, graph = graphs[buffer_cells]
            s, t = snap(free, *cell_of(*a)), snap(free, *cell_of(*b))
            if s is None or t is None:
                continue
            dist, pred = dijkstra(graph, indices=idx[s], return_predecessors=True)
            if np.isfinite(dist[idx[t]]):
                node, chain = idx[t], []
                while node != idx[s] and node >= 0:
                    chain.append(node)
                    node = pred[node]
                return list(reversed(chain)), ("" if buffer_cells == buf else "hugging the shore")
            reason = "no sea route inside the grid"
        return None, reason

    def land_cells_on(a, b, skip_m):
        n = max(2, int(math.hypot(b[0] - a[0], b[1] - a[1]) / cell))
        xs, ys = np.linspace(a[0], b[0], n), np.linspace(a[1], b[1], n)
        r, c = xy_to_rowcol(xs, ys, meta)
        r = np.clip(np.rint(r), 0, H - 1).astype(int)
        c = np.clip(np.rint(c), 0, W - 1).astype(int)
        return int((land[r, c] & (np.hypot(xs - a[0], ys - a[1]) > skip_m)).sum())

    shore_i = next((i for i, p in enumerate(points) if p["kind"] == "shoreline"), 0)
    if points[-1]["kind"] == "sea" and len(points) == shore_i + 2:
        # a transect without a station: when the offshore leg hits land, swing its bearing until it is at sea
        a, b = xy[shore_i], xy[-1]
        length, bearing = math.hypot(b[0] - a[0], b[1] - a[1]), math.degrees(math.atan2(b[0] - a[0], b[1] - a[1]))
        moved = None
        for turn in (0, 15, -15, 30, -30, 45, -45):
            cand = (a[0] + length * math.sin(math.radians(bearing + turn)), a[1] + length * math.cos(math.radians(bearing + turn)))
            if land_cells_on(a, cand, 50.0) == 0:
                moved = (cand, f"offshore end turned {turn:+d} deg so the {length / 1000:g} km leg stays at sea") if turn else None
                break
        else:
            # land all round (an estuary, a ria): end on the nearest sea cell and let the leg be routed to it
            r0, c0 = cell_of(*b)
            if coarse[r0, c0]:
                if 0 not in graphs:
                    graphs[0] = sea_graph(0)
                cell_rc = snap(graphs[0][0], r0, c0)
                if cell_rc is not None:
                    cand = (meta["minx"] + (cell_rc[1] * f + f / 2) * cell, meta["maxy"] - (cell_rc[0] * f + f / 2) * cell)
                    moved = (cand, f"offshore end was on land; moved {math.hypot(cand[0] - b[0], cand[1] - b[1]):.0f} m to the nearest sea cell")
        if moved:
            cand, why = moved
            lon, lat = inv.transform(*cand)
            points = points[:-1] + [{**points[-1], "lon": float(lon), "lat": float(lat)}]
            xy = xy[:-1] + [cand]
            log(why)
    out = list(points[:shore_i + 1])
    max_detour = float(cfg.get("max_detour", 2.5))
    for i in range(shore_i, len(points) - 1):
        a, b = xy[i], xy[i + 1]
        leg = f"{points[i]['name']} -> {points[i + 1]['name']}"
        straight = math.hypot(b[0] - a[0], b[1] - a[1])
        n_land = land_cells_on(a, b, 50.0 if points[i]["kind"] == "shoreline" else 0.0)
        route = None
        if n_land:
            chain, note = shortest(a, b)
            if chain is None:
                log(f"leg {leg}: {n_land} land cells on the straight line, {note}; kept straight")
            else:
                cells = [(meta["minx"] + (n % w * f + f / 2) * cell, meta["maxy"] - (n // w * f + f / 2) * cell) for n in chain]
                line = LineString([a] + cells[:-1] + [b]).simplify(0.75 * f * cell)
                if line.length > max_detour * straight:
                    log(f"leg {leg}: the sea route is {line.length / 1000:.1f} km for a {straight / 1000:.1f} km leg "
                        f"(> {max_detour:g} x); kept straight")
                else:
                    route = list(line.coords)[1:-1]
                    rr, cc = np.array([n // w for n in chain]), np.array([n % w for n in chain])
                    if (rr.min() < 2) or (cc.min() < 2) or (rr.max() > h - 3) or (cc.max() > w - 3):
                        note = (note + ", " if note else "") + "along the grid edge: wider grid.margin_m would give a natural route"
                    log(f"leg {leg}: {n_land} land cells on the straight line; routed through {len(route)} waypoints, "
                        f"{line.length / 1000:.2f} km instead of {straight / 1000:.2f} km" + (f" ({note})" if note else ""))
        for x, y in route or []:
            lon, lat = inv.transform(x, y)
            out.append({"name": "sea route", "lon": float(lon), "lat": float(lat), "kind": "waypoint", "routed": True})
        out.append(points[i + 1])
    return out


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
    coastal = ndimage.map_coordinates(grid["coastal_zones"], [row, col], order=0, mode="nearest").astype(int)
    lon, lat = Transformer.from_crs(meta["crs"], "EPSG:4326", always_xy=True).transform(x, y)
    return pd.DataFrame({"s": s - s0, "x": x, "y": y, "lon": lon, "lat": lat, "z": z, "land": land, "landcover": landcover,
                         "coastal_zone": coastal})


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


def station_title(st) -> str:
    """'TE_060  Samnanger', or just the alias when the export has no station name (or repeats the alias)."""
    name = st.station_name if isinstance(st.station_name, str) else ""
    name = name.strip()
    return f"{st.alias}  {name}" if name and name.lower() not in ("nan", str(st.alias).lower()) else str(st.alias)


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
    """Water down to the floor, the ground under the whole line in the sea-floor colour, the land on top
    (a land mask alone leaves a white wedge in each step where land turns to sea). The water stops one
    sample into the land, so land below sea level (polders, saltpans: Texel 2, Split 4) stays dry."""
    sea = ~np.asarray(land, bool)
    wet = sea | np.r_[False, sea[:-1]] | np.r_[sea[1:], False]
    ax.fill_between(s, np.minimum(z, 0.0), 0.0, where=wet, color=C["water"], lw=0, zorder=1)
    ax.fill_between(s, z, ylo, color=C["sea_floor"], lw=0, zorder=2)
    ax.fill_between(s, z, ylo, where=land, color=C["land"], lw=0, zorder=2.5)
    ax.plot(s, z, color=C["ground_line"], lw=0.9, zorder=3)
    ax.axhline(0.0, color=C["sea_level"], lw=0.7, ls=(0, (4, 3)), zorder=3)


def draw_class_strip(ax, s, codes, colours, frac, edge) -> list[tuple[float, float, int]]:
    """A strip along the foot of the panel, `frac` of its height, coloured by the class at each distance;
    a white hairline where the class changes (several legend colours are close). Returns the runs
    drawn as (from, to, class)."""
    runs = []
    change = np.flatnonzero(np.diff(codes)) + 1
    for i, j in zip(np.r_[0, change], np.r_[change, len(codes)]):
        x0 = s[i] if i == 0 else (s[i - 1] + s[i]) / 2
        x1 = s[j - 1] if j == len(codes) else (s[j - 1] + s[j]) / 2
        if int(codes[i]) in colours:
            ax.axvspan(x0, x1, ymin=0.0, ymax=frac, color=colours[int(codes[i])], lw=0, zorder=4.4)
            runs.append((float(x0), float(x1), int(codes[i])))
        if i > 0:
            ax.axvline(x0, ymin=0.0, ymax=frac, color="white", lw=0.5, zorder=4.5)
    ax.plot([0, 1], [frac, frac], transform=ax.transAxes, color=edge, lw=0.5, zorder=4.6)
    return runs


def class_legend(fig, codes, colours, classes, title, subtitle, edge, right=0.94) -> None:
    """The classes of a strip as a legend at the top right of the figure, in the order the path meets them
    (land to sea), read row by row (a legend fills by column); fewer columns when it would reach the subtitle."""
    from matplotlib.patches import Patch
    order = [c for c in dict.fromkeys(np.asarray(codes).tolist()) if c in colours]
    renderer = fig.canvas.get_renderer()
    for ncol in (min(len(order), 3), 2, 1):
        rows = -(-len(order) // ncol)
        cells = [order[r * ncol + c] if r * ncol + c < len(order) else None for c in range(ncol) for r in range(rows)]
        handles = [Patch(facecolor=colours[c], edgecolor=edge, lw=0.3, label=classes[c][0]) if c is not None
                   else Patch(alpha=0.0, label="") for c in cells]
        leg = fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(right, 0.995), ncol=ncol, fontsize=6,
                         frameon=False, handlelength=1.5, handleheight=0.8, handletextpad=0.5, columnspacing=1.3,
                         labelspacing=0.3, borderaxespad=0.0, title=title, title_fontsize=6.2, alignment="left")
        if ncol == 1 or leg.get_window_extent(renderer).x0 > subtitle.get_window_extent(renderer).x1 + 10:
            break
        leg.remove()


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
    # Copernicus Coastal Zones along the path: a strip at the foot of every panel, when the transect has the layer
    CZ = params["sources"]["coastal_zones"]
    cz = profile["coastal_zone"].to_numpy(int) if "coastal_zone" in profile.columns else np.zeros(len(s), int)
    cz_colours = {int(c): soften_hex(v[1], CZ.get("soften", [1, 1])) for c, v in CZ["classes"].items()}
    strip = float(P.get("coastal_zones_strip_frac", 0.0)) if np.isin(cz, list(cz_colours)).any() else 0.0
    step = float(P["step_m"])
    brk = float(P["break_at_m"])
    s_min, s_max = float(s.min()), float(s.max())
    zoom = P.get("land_zoom", {})
    use_zoom = bool(zoom.get("enabled", True)) and len(land_samples) > 0
    za, zb = float(zoom.get("from_m", s_min)), float(zoom.get("to_m", 120.0))
    if use_zoom:                                   # every land sample in the detail panel (Helsingor 2 has one 3 km inland)
        za = min(za, -float(land_samples["shore_distance_m"].max()) - 20.0)
    start = zb if use_zoom else s_min
    ranges = [(start, brk), (brk, s_max)] if s_max > 2 * brk else [(start, s_max)]
    s_end = s_max + float(P.get("end_padding_frac", 0.0)) * (s_max - ranges[-1][0])   # room after the last station
    ranges[-1] = (ranges[-1][0], s_end)
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
            label = "\n".join(t for t in (station_title(st).replace("  ", "\n"), depth_txt, cast_txt) if t)
            ha = _align(st.s, a, b)
            dx = {"left": 4, "right": -4, "center": 0}[ha]
            ax.annotate(label, (st.s, 0.0), xytext=(dx, 6), textcoords="offset points", ha=ha, va="bottom",
                        fontsize=5.8, color=C["station"], zorder=7, linespacing=1.15)

    for k, (ax, (a, b)) in enumerate(zip(axes, ranges)):
        m = (s >= a - 2 * step) & (s <= b + 2 * step)
        _draw_terrain(ax, s[m], z[m], land[m], ylo, C)
        if strip:
            draw_class_strip(ax, s[m], cz[m], cz_colours, strip, C["ground_line"])
        _draw_zones(ax, P["zones"], a, b, s_min, s_end, yhi, yspan, C)
        for _, r in land_samples.iterrows():
            if a <= r.s <= b:
                ax.plot(r.s, r.z, marker="v", ms=4.5, color=C["sample"], mec="white", mew=0.5, ls="none", zorder=6)
        draw_stations(ax, a, b)
        ax.set_xlim(a, b)
        ax.set_ylim(ylo, yhi)
        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(7, prune="upper" if k < len(ranges) - 1 else "lower"))
        if min(b, s_max) <= brk * 1.001:           # the end padding does not turn a metre panel into kilometres
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
        if strip:
            draw_class_strip(ax0, s[m], cz[m], cz_colours, strip, C["ground_line"])
        wide = zb - za > 1.2 * (float(zoom.get("to_m", 120.0)) - float(zoom.get("from_m", -320.0)))   # widened for a far sample
        _draw_zones(ax0, P["zones"], za, zb, s_min, s_end, ztop, ztop - zlo, C, min_frac=0.05 if wide else 0.0,
                    fontsize=5.6, stagger=0.11)
        ax0.set_xlim(za, zb)
        ax0.set_ylim(zlo, ztop)
        placements = [(0, 6, "center"), (5, 14, "left"), (-5, 14, "right")]
        rows = [(i, r) for i, (_, r) in enumerate(land_samples.sort_values("s", ascending=False).iterrows())
                if za <= r.s <= zb]
        marks = [ax0.plot(r.s, r.z, marker="v", ms=5, color=C["sample"], mec="white", mew=0.5, ls="none", zorder=6)[0]
                 for _, r in rows]
        renderer = fig.canvas.get_renderer()       # the markers and the drawn labels: a label over one moves up
        placed = [m.get_window_extent(renderer) for m in marks]
        for i, r in rows:
            dx, dy, ha = placements[i % len(placements)]
            frac = (r.s - za) / (zb - za)
            if frac > 0.88 or frac < 0.12:     # near a panel edge: the label opens towards the panel
                dx, ha = (-5, "right") if frac > 0.88 else (5, "left")
            ann = ax0.annotate(f"{r.zone} {r.shore_distance_m:g} m", (r.s, r.z), xytext=(dx, dy),
                               textcoords="offset points", ha=ha, va="bottom", fontsize=5.6, color=C["sample"], zorder=7)
            for _ in range(8):
                if not any(ann.get_window_extent(renderer).overlaps(b) for b in placed):
                    break
                dy += 8
                ann.set_position((dx, dy))
            placed.append(ann.get_window_extent(renderer))
        draw_stations(ax0, za, zb)
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
    subtitle = fig.text(0.01, 0.925, tcfg.get("subtitle", ""), fontsize=7.5, va="top", color=muted)
    credit = f"Land: {params['sources']['dem']['citation']}.  Sea floor: {params['sources']['bathymetry']['citation']}."
    if strip:
        class_legend(fig, cz, cz_colours, CZ["classes"], "Copernicus Coastal Zones 2018 along the path, land to sea "
                     "(the strip at the foot of each panel)", subtitle, muted)
        credit = f"Land cover / land use: {CZ['citation']}.\n" + credit
    fig.text(0.99, 0.008, credit, fontsize=5, ha="right", va="bottom", color=muted, linespacing=1.3)
    fig.canvas.draw()
    for ax, (a, b) in panels:
        bb = ax.get_window_extent()
        y0, y1 = ax.get_ylim()
        ve = ((b - a) / bb.width) / ((y1 - y0) / bb.height)
        ax.text(0.99, 0.02 + strip, f"vertical exaggeration x{ve:.0f}", transform=ax.transAxes, ha="right", va="bottom",
                fontsize=5.5, color=muted)
    fig.savefig(out_png, dpi=F["dpi"])
    fig.savefig(out_pdf)
    plt.close(fig)
