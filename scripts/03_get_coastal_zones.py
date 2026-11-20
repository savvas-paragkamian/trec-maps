#!/usr/bin/env python
"""Save the Copernicus Coastal Zones 2018 polygons over a transect's grid (network step; the EEA's public map service)."""
import time

from _common import C, parser, setup
from trecmaps import fetch, sites, terrain


def main() -> None:
    ap = parser(__doc__)
    ap.add_argument("--dry-run", action="store_true", help="count the polygons without reading them")
    args = ap.parse_args()
    params, tcfg = setup(args)
    S = params["sources"]["coastal_zones"]
    log = lambda m: C.log("coastal zones", m)
    points, _, _ = sites.load_path(params, tcfg)
    _, lonlat = terrain.grid_box(params["grid"], points)
    dest = terrain.coastal_zones_path(params, tcfg["name"])
    if dest.exists():
        if terrain.polygons_cover(dest, lonlat):
            log(f"have {dest.name} and it covers the grid, skipped")
            return
        log(f"{dest.name} does not cover the current grid box; fetching it again")
    # the ids of the polygons that touch the box (an id list has no feature limit), then the polygons
    # `chunk` ids per request; the few huge ones (the sea) one by one, kept in _big/ for the next transect
    query = S["service"].rstrip("/") + "/query"
    box = {"geometry": ",".join(f"{v:.6f}" for v in lonlat), "geometryType": "esriGeometryEnvelope", "inSR": 4326,
           "spatialRel": "esriSpatialRelIntersects", "returnIdsOnly": "true", "f": "json"}
    ids = sorted(fetch.get_json(query, {**box, "where": "1=1"}).get("objectIds") or [])
    big = sorted(fetch.get_json(query, {**box, "where": f"AREA_HA > {float(S['big_ha']):g}"}).get("objectIds") or []) if ids else []
    log(f"{len(ids)} polygons touch the grid box" + (f", {len(big)} of them over {float(S['big_ha']):g} ha" if ids else
                                                     " (outside the layer's coverage)"))
    if args.dry_run:
        return
    geometry = {"outFields": S["field"], "outSR": 4326, "maxAllowableOffset": S["generalise_deg"],
                "geometryPrecision": 6, "f": "geojson"}
    small, chunk, features = [i for i in ids if i not in set(big)], int(S.get("chunk", 300)), []
    for k in range(0, len(small), chunk):
        features += fetch.get_json(query, {**geometry, "objectIds": ",".join(map(str, small[k:k + chunk]))})["features"]
        time.sleep(1)
    for i in big:
        cache = dest.parent / "_big" / f"cz_2018_{i}.geojson.gz"
        if not cache.exists():
            C.write_json(fetch.get_json(query, {**geometry, "objectIds": i}), cache)
            log(f"saved {cache.name} ({cache.stat().st_size / 2**20:.1f} MiB)")
            time.sleep(1)
        features += C.read_json(cache)["features"]
    if len(features) != len(ids):
        raise SystemExit(f"the service listed {len(ids)} polygons and returned {len(features)}")
    clipped = terrain.clip_polygons(features, lonlat, S["field"])
    C.write_json({"type": "FeatureCollection", "bbox": list(lonlat), "source": S["service"],
                  "fetched": time.strftime("%Y-%m-%d %H:%M"), "features": clipped}, dest)
    log(f"saved {dest.name} ({len(clipped)} polygons inside the box, {len({f['properties']['code'] for f in clipped})} classes, "
        f"{dest.stat().st_size / 2**10:.0f} KiB)")


if __name__ == "__main__":
    main()
