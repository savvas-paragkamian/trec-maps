# AGENTS.md — project memory for trec-maps

Read this first. README.md is the user guide; this file holds the rules, decisions and
gotchas, and a dated log of the work (requests quoted verbatim).

## Purpose

Python map scripts for TREC transects: land-sea profile maps from Copernicus GLO-30 (land)
and EMODnet Bathymetry (sea), as 2D profiles and forge3d 3D blocks. Pilot: Anavyssos (041-03)
with Tara stations TE_176 and TE_177. Owner: Savvas Paragkamian (EMBL CPP fellow, TREC soil
metagenomes); the sibling repository `trec-coastal-microbiome` (R, targets) holds the
analysis pipeline and the per-transect 2D maps.

## Hard rules

- **TREC data never go to git**: `data/*` is ignored except `data/README.md`; no coordinates or
  sample tables in `config/` (transects are defined by site IDs, read from the xlsx at run time).
  Publication of TREC data needs the consortium's consent.
- **Network only in `scripts/0*_get_*.py`** (through `trecmaps/fetch.py`): `--dry-run`, write
  `<file>.part` then rename, skip files that exist. `trecmaps/` never downloads anything.
- **Conda from free channels only**: `conda-forge`, `bioconda`, `nodefaults`, strict priority,
  as in `trec-coastal-microbiome`. The environment lives in conda's envs dir, never inside
  this Dropbox folder. forge3d comes from PyPI (pip section), 1.44.0 is the verified version.
- **No `Co-Authored-By: Claude` lines** in commits or PRs. The user commits and merges.
- Keep it small ("less is more"): six modules, five scripts, two YAML files. Nothing hard-coded
  in `trecmaps/` that belongs in `config/params.yml`.
- Do not name a data source that was not verified; the citations in `params.yml` are the ones to use.

## Layout

```
config/params.yml, transects.yml     every tunable; transects by site ID
trecmaps/config.py                   YAML loading, --set overrides, hex colours
trecmaps/sites.py                    xlsx -> sites, land positions (shore_distance), stations (cable_out), path points
trecmaps/fetch.py                    the download helper (HEAD size, .part, 404 = ocean tile)
trecmaps/terrain.py                  grid box, GLO-30 tile names, WCS request, reprojection, land/sea merge, build_grid
trecmaps/profile.py                  polyline, sampling, sample/station placement, checks, the 2D figure
trecmaps/render.py                   forge3d tiled render (from crete3d), albedo ramps, camera -> pixel projection
trecmaps/compose.py                  keyed background, vignette, Page (matplotlib in pixel space), furniture
scripts/00..30                       thin CLIs; scripts/_common.py adds the repo root to sys.path
data/raw, data/derived/<t>, output/<t>
```

## Decisions

- **Transect path** = inland end (inland-most sample + 50 m) -> shoreline sample -> stations in
  the order of `marine_sites`, with optional waypoints. For Anavyssos the TE_176 -> TE_177 leg
  would cross the HCMR cape and the cape west of it, hence two waypoints south of them; the
  path is a V (TE_176 lies east of the shoreline, TE_177 10 km west). `s` is the distance along
  that path, 0 at the shoreline, negative inland; land samples sit at `s = -shore_distance`.
- **Sampling depths** of the Tara stations come from `cable_out` in the event attributes; the
  sheet's `Depth range (m)` is "0-3" for nearly every marine sample and is not used. Water depth
  per station: the site's `bathymetry` attribute (sounder values as fallback).
- **Grid**: EPSG:2100 (the user's default for Greek maps; EPSG:3035 is what the R pipeline
  uses), 10 m cells, margins 2.5 km north / 1.5 km elsewhere around the path's box. Anavyssos:
  1337 x 638 cells, 0.85 M, far under forge3d's limit (about 7 M cells).
- **Land/sea merge**: land = GLO-30 > 0.5 m minus components < 8 cells; sea floor = EMODnet
  (cubic, land set to nodata first) capped at -0.5 m, ramped to the cap over 6 cells from the
  shoreline. The grid agrees with the recorded station depths within 1.5 m (TE_176 -9.2 vs
  -8.6, TE_177 -129.6 vs -131).
- **3D style**: one oblique block, land hypsometric, sea floor depth-tinted, no water surface;
  path and stations are drawn after the render by projecting world points with the camera
  (vector overlay in the PDF). Camera 170°/38°, fov 21°, sun 125°/35°, exaggeration 3.
- **2D figure**: two panels with a broken x axis at 1.5 km (metres, then kilometres), shared
  y, a zoom inset of the land transect in the right panel, zone labels only where they fit.

## Gotchas

- forge3d (1.44.0, Metal): `sun_azimuth_deg` is counter-clockwise (pass 360 - bearing); hard
  512 MiB per call; tile sizes must be frame / 2^k; `min_frames >= 2`; `earth_model="flat"`
  needs `refraction_model="none"`. The depth AOV is the **ray distance** (verified: median
  rendered / projected = 1.0003 with `render.world_to_pixel`), NaN on misses.
- The plate is centred on the origin, column -> +X, row -> +Z, height -> +Y, in metres
  (`spacing = cell`); `render.grid_to_world` relies on this.
- EMODnet `emodnet:mean` holds elevations: negative at sea, positive on land (up to 469 m in
  the Anavyssos box). `terrain.read_bathy` drops `>= 0` and flips the sign only when a grid has
  no negative value at all. The WCS answers EPSG:4326 only; request on the native 1/960° grid.
- The bathymetry extract is cut to the grid box, so a changed path (waypoints) can leave it too
  small: `01_get_bathymetry.py` re-fetches when the saved file does not cover the box and
  `build_grid` refuses to run with a short one.
- GLO-30 is a surface model (buildings, trees); `grid.land_smooth_sigma_cells` can soften it.
  The DEM at the 041-03 site is 1.9 m where the site export says 7 m (the site point is on the beach).
- xlsx exports: title in row 1, header in row 2 (`pd.read_excel(header=1)`); attribute columns
  are JSON lists of `{attribute, value, unit}`.
- Dropbox: `PYTHONDONTWRITEBYTECODE=1` in the Makefile, `*conflicted copy*` ignored, no venv here.
- Timings (M1 Max): preview 1600x900 x 6 frames about 1 s; final 3840x2160 x 64 frames about 12 s (plus a few seconds for the page).

## Verify

```bash
make check                                   # forge3d 1.44.0 gpu True
python scripts/10_build_grid.py --force      # prints z range, land fraction, station depths vs recorded
python scripts/20_profile_2d.py              # "marine legs stay at sea" and the station depth checks
python scripts/30_render_3d.py               # "median rendered/projected = 1.000x" and "path: 100% visible"
```

## Log

### 2026-11-02 — repository created, Anavyssos pilot

Request: "based on ~/projects/crete3d, trec-coastal-microbiome/, /Users/talos/Dropbox/10_projects/BIAS-II/spilaio-angiti/map-angitis.
I want this repo to be a map scripts repository with conda (miniforge4, only with channels like the
trec-coastal-microbiome) for 3d maps with forge3d. This will be a python based repo. I want profile
maps of athens along of DEM and then bathymetry data of EMODNET. Do a pilot first with one transect
for 3d. check the other repos and read their memory and agent and readme."

Choices confirmed with the user: pilot Anavyssos + TE_176 + TE_177; 3D as a terrain block with the
transect overlaid (not a ribbon); the environment created with the existing miniconda3 from the
Miniforge-style `environment.yml` (no Miniforge installed; the current release is still named
Miniforge3, there is no Miniforge4); Copernicus GLO-30 from AWS as the DEM.

Done: environment, config, the six modules and five scripts, downloads for Anavyssos (one GLO-30
tile, a 161 x 70 EMODnet extract), the merged grid, `profile.png/.pdf`, `block_preview` and
`block_final` with the overlay. Not done: the other four Athens transects (041-01, -02, -04, -05;
add entries to `transects.yml`), a ribbon/section 3D style, GitHub remote and first commit (the user's).
