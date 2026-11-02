"""Sites, land sample positions, marine stations and auto-generated transects from the TREC data exports.

The xlsx exports carry a title in row 1 and the header in row 2. The attribute columns
("Site specific attributes", "Event specific attributes") hold JSON lists of
{"attribute", "value", "unit"} records. Both exports are sample tables that repeat the site
columns ("Site ID", "Site alias", "Site latitude (DD)", ...), so the site table is derived
from them. Each export is cached as a pickle under data/derived/_cache (parsing the marine
export takes 11 s, reading the pickle 0.1 s).
"""
from __future__ import annotations

import json
import math
import re

import numpy as np
import pandas as pd
from pyproj import Geod

GEOD = Geod(ellps="WGS84")
_MEMO: dict[str, pd.DataFrame] = {}


def read_export(params: dict, key: str) -> pd.DataFrame:
    """`paths.<key>` as a DataFrame (header in row 2), cached as a pickle under the derived dir."""
    path = params["paths"][key]
    if str(path) in _MEMO:
        return _MEMO[str(path)]
    if not path.exists():
        raise SystemExit(f"missing {path} (see data/README.md)")
    cache = params["paths"]["derived_dir"] / "_cache" / f"{path.stem}.pkl"
    if cache.exists() and cache.stat().st_mtime >= path.stat().st_mtime:
        df = pd.read_pickle(cache)
    else:
        df = pd.read_excel(path, header=1)
        df.columns = [str(c).strip() for c in df.columns]
        if "Site ID" in df.columns:
            df["Site ID"] = df["Site ID"].astype(str).str.strip()
        cache.parent.mkdir(parents=True, exist_ok=True)
        df.to_pickle(cache)
    _MEMO[str(path)] = df
    return df


def parse_attributes(cell) -> dict:
    """JSON list of {attribute, value, unit} -> {attribute: value}; {} when empty or invalid."""
    if not isinstance(cell, str) or not cell.strip():
        return {}
    try:
        obj = json.loads(cell)
    except ValueError:
        return {}
    if isinstance(obj, dict):
        return obj
    out = {}
    if isinstance(obj, list):
        for rec in obj:
            if isinstance(rec, dict) and "attribute" in rec:
                out[str(rec["attribute"])] = rec.get("value")
    return out


def to_float(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _numeric(df: pd.DataFrame, column: str) -> pd.Series:
    if column not in df.columns:
        return pd.Series(math.nan, index=df.index, dtype=float)
    return pd.to_numeric(df[column], errors="coerce")


def _split_ids(value) -> list[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()] if isinstance(value, str) else []


def load_sites(params: dict) -> pd.DataFrame:
    """One row per TREC site (LSI and TARA), indexed by site ID, from the site columns of both exports:
    alias, locality, site type, country, lat, lon, elevation, bathymetry, station name, related site IDs."""
    df = pd.concat([read_export(params, "marine_xlsx"), read_export(params, "samples_xlsx")], ignore_index=True)
    first = df.groupby("Site ID").first()            # first non-null value per column
    attrs = first["Site specific attributes"].map(parse_attributes)
    related = df.groupby("Site ID")["Site related sites"].agg(
        lambda s: sorted({x for v in s.dropna() for x in _split_ids(v)}))
    sites = pd.DataFrame({
        "site_id": first.index,
        "alias": first["Site alias"].astype(str).str.strip(),
        "locality": first["Site locality"].astype(str).str.strip(),
        "site_type": first["Site type"].astype(str).str.strip(),
        "country": first["Event country"].astype(str).str.strip(),
        "lat": pd.to_numeric(first["Site latitude (DD)"], errors="coerce"),
        "lon": pd.to_numeric(first["Site longitude (DD)"], errors="coerce"),
        "elevation_m": attrs.map(lambda a: to_float(a.get("elevation"))),
        "bathymetry_m": attrs.map(lambda a: to_float(a.get("bathymetry"))),
        "station_name": attrs.map(lambda a: a.get("station_name")),
        "related": related,
    }, index=first.index)
    return sites


def load_land_positions(params: dict, site_id: str) -> pd.DataFrame:
    """Distinct sampling positions of a land transect, sorted by distance from the shoreline.

    Columns: context, zone (first word of the context), shore_distance_m, lat, lon,
    n_samples, medium. Event coordinates are used when present, else the site's. A missing
    shore distance is replaced by the zone's nominal one (`profile.nominal_shore_distance_m`).
    """
    df = read_export(params, "samples_xlsx")
    rows = df[df["Site ID"] == site_id].copy()
    if rows.empty:
        raise SystemExit(f"no soil/sediment samples for site {site_id} in {params['paths']['samples_xlsx'].name}")
    attrs = rows["Event specific attributes"].map(parse_attributes)
    rows["shore_distance_m"] = attrs.map(lambda a: to_float(a.get("shore_distance")))
    rows["context"] = rows["Local environmental context"].astype(str).str.strip()
    nominal = params["profile"].get("nominal_shore_distance_m") or {}
    by_zone = rows["context"].map(lambda c: to_float(nominal.get(c.split()[0] if c else "", math.nan))).astype(float)
    rows["shore_distance_m"] = rows["shore_distance_m"].fillna(by_zone)
    rows["lat"] = _numeric(rows, "Event latitude start (DD)").fillna(_numeric(rows, "Site latitude (DD)"))
    rows["lon"] = _numeric(rows, "Event longitude start (DD)").fillna(_numeric(rows, "Site longitude (DD)"))
    medium_col = "Environmental medium" if "Environmental medium" in rows.columns else "Sample type"
    grouped = (rows.groupby(["context", "shore_distance_m"], dropna=False)
               .agg(lat=("lat", "median"), lon=("lon", "median"), n_samples=("lat", "size"),
                    medium=(medium_col, lambda s: ", ".join(sorted(set(map(str, s))))))
               .reset_index())
    grouped["zone"] = grouped["context"].str.split().str[0]
    grouped = grouped.sort_values("shore_distance_m", na_position="last").reset_index(drop=True)
    return grouped


def load_marine_stations(params: dict, site_ids: list[str]) -> pd.DataFrame:
    """Tara stations in the given order: position, bathymetry, sampled depths (cable_out), sounder depths."""
    sites = load_sites(params)
    df = read_export(params, "marine_xlsx")
    records = []
    for sid in site_ids:
        if sid not in sites.index:
            raise SystemExit(f"site {sid} is in neither export (paths.marine_xlsx, paths.samples_xlsx)")
        site = sites.loc[sid]
        rows = df[df["Site ID"] == sid]
        events = rows["Event specific attributes"].map(parse_attributes) if "Event specific attributes" in rows.columns else []
        casts, sounder = set(), []
        for a in events:
            cable = to_float(a.get("cable_out"))
            if not math.isnan(cable) and cable >= 0:
                casts.add(round(cable, 1))
            for key in ("sounder_in", "sounder_out"):
                v = to_float(a.get(key))
                if not math.isnan(v) and v > 0:
                    sounder.append(v)
        if not casts:
            casts = {round(float(v), 1) for v in _numeric(rows, "Depth (m)").dropna() if v >= 0}
        bathy = to_float(site["bathymetry_m"])
        if math.isnan(bathy) and sounder:
            bathy = float(pd.Series(sounder).median())
        records.append({
            "site_id": sid, "alias": site["alias"],
            "station_name": site["station_name"] if isinstance(site["station_name"], str) else site["alias"],
            "lat": float(site["lat"]), "lon": float(site["lon"]), "bathymetry_m": bathy,
            "cast_depths_m": sorted(casts),
            "sounder_min_m": min(sounder) if sounder else math.nan,
            "sounder_max_m": max(sounder) if sounder else math.nan,
            "n_samples": int(len(rows)),
        })
    columns = ["site_id", "alias", "station_name", "lat", "lon", "bathymetry_m", "cast_depths_m",
               "sounder_min_m", "sounder_max_m", "n_samples"]
    return pd.DataFrame(records, columns=columns)


def auto_transects(params: dict, rules: dict, taken: set[str] = frozenset()) -> dict[str, dict]:
    """A transect per LSI site of the samples export that no explicit entry covers (`taken` land sites).

    Stations: the site's related Tara stations within `max_station_km` of it, nearest first, at most
    `max_stations`; when none is that close, the nearest Tara station within `max_unrelated_km`; else
    none (the profile continues `profile.sea_extension_m` offshore). Names are the site aliases in
    lower case (Roscoff_LSI_1 -> roscoff_lsi_1); `overrides` in the auto section are merged on top.
    """
    sites = load_sites(params)
    tara = sites[sites.site_type == "TARA"].dropna(subset=["lat", "lon"])
    max_km = float(rules.get("max_station_km", 20))
    near_km = float(rules.get("max_unrelated_km", max_km))
    max_n = int(rules.get("max_stations", 2))
    offshore_km = float(params["profile"].get("sea_extension_m", 0) or 0) / 1000
    out = {}
    for sid in sorted(set(read_export(params, "samples_xlsx")["Site ID"])):
        if sid in taken or sid not in sites.index or sites.loc[sid, "site_type"] != "LSI":
            continue
        site = sites.loc[sid]
        if not (np.isfinite(site.lat) and np.isfinite(site.lon)):
            continue
        n = len(tara)
        _, _, dist = GEOD.inv(np.full(n, site.lon), np.full(n, site.lat), tara.lon.to_numpy(), tara.lat.to_numpy())
        within = sorted((d / 1000, t) for d, t in zip(dist, tara.index) if d / 1000 <= max_km)
        related = [(d, t) for d, t in within if t in set(site.related)]
        chosen = (related or [(d, t) for d, t in within if d <= near_km][:1])[:max_n]
        name = re.sub(r"[^a-z0-9]+", "_", str(site.alias).lower()).strip("_")
        where = f"{str(site.locality).replace('_', ' ')}, {site.country}"
        if chosen:
            what = " and ".join(f"{tara.loc[t, 'alias']} ({d:.1f} km)" for d, t in chosen)
            subtitle = f"TREC land-sea transect {sid} ({where}) with Tara station{'s' if len(chosen) > 1 else ''} {what}"
        else:
            subtitle = (f"TREC land-sea transect {sid} ({where}); no Tara station nearby, "
                        f"the profile continues {offshore_km:g} km offshore")
        out[name] = {"title": str(site.alias).replace("_", " "), "subtitle": subtitle, "land_site": sid,
                     "marine_sites": [t for _, t in chosen], "waypoints": [], "land_path": "from_samples",
                     "auto": True}
    return out


def extend_point(from_lonlat, to_lonlat, distance_m: float) -> tuple[float, float]:
    """The point `distance_m` beyond `to` on the geodesic from `from` through `to`."""
    az, _, _ = GEOD.inv(from_lonlat[0], from_lonlat[1], to_lonlat[0], to_lonlat[1])
    lon, lat, _ = GEOD.fwd(to_lonlat[0], to_lonlat[1], az, distance_m)
    return float(lon), float(lat)


def path_points(tcfg: dict, land: pd.DataFrame, stations: pd.DataFrame, land_extension_m: float,
                sea_extension_m: float = 0.0) -> list[dict]:
    """Ordered path: inland end -> shoreline -> [waypoints] -> stations, each {name, lon, lat, kind}.
    Without stations the path continues `sea_extension_m` offshore along the transect's bearing."""
    land_path = tcfg.get("land_path", "from_samples")
    if land_path == "from_samples":
        shore, inland = land.iloc[0], land.iloc[-1]
        shore_ll, inland_ll = (float(shore.lon), float(shore.lat)), (float(inland.lon), float(inland.lat))
        if abs(shore_ll[0] - inland_ll[0]) < 1e-7 and abs(shore_ll[1] - inland_ll[1]) < 1e-7:
            raise SystemExit("the land events share one position; set land_path: [[lon, lat], [lon, lat]] "
                             "(shoreline, inland) in config/transects.yml")
    else:
        shore_ll, inland_ll = tuple(map(float, land_path[0])), tuple(map(float, land_path[1]))
    inland_end = extend_point(shore_ll, inland_ll, land_extension_m)
    points = [
        {"name": "inland end", "lon": inland_end[0], "lat": inland_end[1], "kind": "inland"},
        {"name": "shoreline", "lon": shore_ll[0], "lat": shore_ll[1], "kind": "shoreline"},
    ]
    waypoints = {w["after"]: w["points"] for w in tcfg.get("waypoints") or []}
    for lon, lat in waypoints.get("shoreline", []):
        points.append({"name": "waypoint", "lon": float(lon), "lat": float(lat), "kind": "waypoint"})
    for _, st in stations.iterrows():
        points.append({"name": st.alias, "lon": float(st.lon), "lat": float(st.lat), "kind": "station",
                       "site_id": st.site_id})
        for lon, lat in waypoints.get(st.site_id, []):
            points.append({"name": "waypoint", "lon": float(lon), "lat": float(lat), "kind": "waypoint"})
    if len(stations) == 0 and sea_extension_m > 0:
        sea_end = extend_point(inland_ll, shore_ll, sea_extension_m)
        points.append({"name": "offshore end", "lon": sea_end[0], "lat": sea_end[1], "kind": "sea"})
    return points


def load_path(params: dict, tcfg: dict) -> tuple[list[dict], pd.DataFrame, pd.DataFrame]:
    """(path points, land positions, stations) for a transect."""
    land = load_land_positions(params, tcfg["land_site"])
    stations = load_marine_stations(params, list(tcfg.get("marine_sites") or []))
    points = path_points(tcfg, land, stations, params["profile"]["land_extension_m"],
                         float(params["profile"].get("sea_extension_m", 0) or 0))
    return points, land, stations


def load_shoreline_water(params: dict, site_id: str) -> dict | None:
    """The land site's own shoreline water samples (marine export, Sample type 'shoreline'): count, depth range, position."""
    df = read_export(params, "marine_xlsx")
    rows = df[(df["Site ID"] == site_id)]
    if "Sample type" in rows.columns:
        rows = rows[rows["Sample type"].astype(str).str.strip().str.lower() == "shoreline"]
    if rows.empty:
        return None
    dmin = _numeric(rows, "Depth range min (m)").median()
    dmax = _numeric(rows, "Depth range max (m)").median()
    lat = _numeric(rows, "Event latitude start (DD)").fillna(_numeric(rows, "Site latitude (DD)")).median()
    lon = _numeric(rows, "Event longitude start (DD)").fillna(_numeric(rows, "Site longitude (DD)")).median()
    return {"n_samples": int(len(rows)), "depth_min_m": 0.0 if math.isnan(dmin) else float(dmin),
            "depth_max_m": 3.0 if math.isnan(dmax) else float(dmax), "lat": float(lat), "lon": float(lon)}
