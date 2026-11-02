"""Sites, land sample positions and marine stations from the TREC data exports.

The xlsx exports carry a title in row 1 and the header in row 2. The attribute
columns ("Specific attributes", "Event specific attributes", ...) hold JSON lists of
{"attribute", "value", "unit"} records.
"""
from __future__ import annotations

import json
import math

import pandas as pd
from pyproj import Geod

GEOD = Geod(ellps="WGS84")


def read_trec_xlsx(path) -> pd.DataFrame:
    df = pd.read_excel(path, header=1)
    df.columns = [str(c).strip() for c in df.columns]
    if "Site ID" in df.columns:
        df["Site ID"] = df["Site ID"].astype(str).str.strip()
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


def load_sites(params: dict) -> pd.DataFrame:
    """One row per TREC site, indexed by site ID: alias, locality, type, lat, lon, elevation, bathymetry, station name."""
    df = read_trec_xlsx(params["paths"]["sites_xlsx"])
    attrs = df["Specific attributes"].map(parse_attributes) if "Specific attributes" in df.columns else [{}] * len(df)
    attrs = pd.Series(list(attrs), index=df.index)
    sites = pd.DataFrame({
        "site_id": df["Site ID"],
        "alias": df.get("Site alias", pd.Series("", index=df.index)).astype(str).str.strip(),
        "locality": df.get("Locality", pd.Series("", index=df.index)),
        "site_type": df.get("Site type", pd.Series("", index=df.index)),
        "lat": _numeric(df, "Latitude (DD)"),
        "lon": _numeric(df, "Longitude (DD)"),
        "elevation_m": attrs.map(lambda a: to_float(a.get("elevation"))),
        "bathymetry_m": attrs.map(lambda a: to_float(a.get("bathymetry"))),
        "station_name": attrs.map(lambda a: a.get("station_name")),
    })
    return sites.set_index("site_id", drop=False)


def load_land_positions(params: dict, site_id: str) -> pd.DataFrame:
    """Distinct sampling positions of a land transect, sorted by distance from the shoreline.

    Columns: context, zone (first word of the context), shore_distance_m, lat, lon,
    n_samples, medium. Event coordinates are used when present, else the site's.
    """
    df = read_trec_xlsx(params["paths"]["samples_xlsx"])
    rows = df[df["Site ID"] == site_id].copy()
    if rows.empty:
        raise SystemExit(f"no soil/sediment samples for site {site_id} in {params['paths']['samples_xlsx'].name}")
    attrs = rows["Event specific attributes"].map(parse_attributes)
    rows["shore_distance_m"] = attrs.map(lambda a: to_float(a.get("shore_distance")))
    rows["context"] = rows["Local environmental context"].astype(str).str.strip()
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
    df = read_trec_xlsx(params["paths"]["marine_xlsx"])
    records = []
    for sid in site_ids:
        if sid not in sites.index:
            raise SystemExit(f"site {sid} is not in {params['paths']['sites_xlsx'].name}")
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
    return pd.DataFrame(records)


def extend_point(from_lonlat, to_lonlat, distance_m: float) -> tuple[float, float]:
    """The point `distance_m` beyond `to` on the geodesic from `from` through `to`."""
    az, _, _ = GEOD.inv(from_lonlat[0], from_lonlat[1], to_lonlat[0], to_lonlat[1])
    lon, lat, _ = GEOD.fwd(to_lonlat[0], to_lonlat[1], az, distance_m)
    return float(lon), float(lat)


def path_points(tcfg: dict, land: pd.DataFrame, stations: pd.DataFrame, land_extension_m: float) -> list[dict]:
    """Ordered path: inland end -> shoreline -> [waypoints] -> stations, each {name, lon, lat, kind}."""
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
    return points


def load_path(params: dict, tcfg: dict) -> tuple[list[dict], pd.DataFrame, pd.DataFrame]:
    """(path points, land positions, stations) for a transect."""
    land = load_land_positions(params, tcfg["land_site"])
    stations = load_marine_stations(params, list(tcfg.get("marine_sites") or []))
    points = path_points(tcfg, land, stations, params["profile"]["land_extension_m"])
    return points, land, stations
