# AGENTS.md — project memory for trec-maps

Read this first. README.md is the user guide; this file holds the rules, decisions and
gotchas, and a dated log of the work (requests quoted verbatim).

## Purpose

Python map scripts for TREC transects: land-sea profile maps from Copernicus GLO-30 (land)
and EMODnet Bathymetry (sea), as 2D profiles and forge3d 3D blocks, for all 114 LSI sites of
the data exports (115 transects with the Anavyssos zoom). Pilot: Anavyssos (041-03)
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
- Keep it small ("less is more"): eight modules, seven scripts, two YAML files. Nothing hard-coded
  in `trecmaps/` that belongs in `config/params.yml`.
- Do not name a data source that was not verified; the citations in `params.yml` are the ones to use.

## Layout

```
config/params.yml, transects.yml     every tunable; transects by site ID
trecmaps/config.py                   YAML loading, --set overrides, hex colours
trecmaps/sites.py                    exports (pickle-cached) -> sites, land positions (shore_distance), stations (cable_out), auto transects, path points
trecmaps/fetch.py                    the download helper (HEAD size, .part, 404 = ocean tile) and the COG window reader
trecmaps/schematic.py                the sampling-concept figure (symlog distance axis, two vertical scales, tubes, cones, boats)
trecmaps/terrain.py                  grid box (UTM auto, cell cap), tile names, DEM/WorldCover windows, WCS request, land/sea merge, build_grid
trecmaps/profile.py                  polyline, sea routing of marine legs, sampling, sample/station placement, checks, the 2D figure
trecmaps/render.py                   forge3d tiled render (from crete3d), albedo ramps, camera -> pixel projection
trecmaps/compose.py                  keyed background, vignette, Page (matplotlib in pixel space), furniture
scripts/00..30                       thin CLIs; scripts/_common.py adds the repo root to sys.path, applies transect overrides, lists transects
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
- **Grid**: first EPSG:2100 (the user's default for Greek maps; EPSG:3035 is what the R pipeline
  uses), 10 m cells, margins 2.5 km north / 1.5 km elsewhere around the path's box. Anavyssos:
  1337 x 638 cells, 0.85 M, far under forge3d's limit (about 7 M cells). Since the all-transects
  round the CRS is the path's UTM zone (`grid.crs: auto`, Athens = EPSG:32634/32635) and
  `grid.max_cells` (6 M) coarsens the cell for very wide boxes (none needed it so far).
- **Land/sea merge**: land = GLO-30 > 0.5 m minus components < 8 cells; sea floor = EMODnet
  (cubic, land set to nodata first) capped at -0.5 m, ramped to the cap over 6 cells from the
  shoreline. The grid agrees with the recorded station depths within 1.5 m (TE_176 -9.2 vs
  -8.6, TE_177 -129.6 vs -131).
- **3D style**: one block seen from the south at a low angle (camera azimuth 180°, elevation
  20°, fov 21°; the zoom 175°/22°), so it reads as a section through the coast rather than a
  map: the cut front wall shows the sea floor profile, the exaggerated depth is seen side-on.
  The first version at 170°/38° looked like a 2D map, which the user rejected. Sun 125°/35°,
  exaggeration 3 on land. Path and stations are drawn after the render by projecting world
  points with the camera (vector overlay in the PDF). The land started hypsometric and is
  WorldCover now; the sea started as a flat depth tint and is a post-processed water surface.
- **2D figure**: three panels on one distance line: the land transect in detail with its own
  vertical scale (the inset it replaced overlapped the sea floor), then the nearshore (metres)
  and the offshore (kilometres) sharing one elevation/depth axis; zone labels only where they fit.
- **3D block as a solid**: the grid is rimmed with one cell at a base depth (`render.block`), so
  the heightfield drops into cut walls; everything below sea level is tinted by its true depth in
  post-processing (`render.water`, from the depth AOV and the camera), which paints the water
  column on the walls as a cross-section; a sea-level line runs along the front wall and a
  vertical axis in metres stands outside the south-west corner. Depths are exaggerated
  `exaggeration * sea_exaggeration_factor` (3 x 2 = 6), land x3, both stated in the legend;
  without the factor the shelf break read as flat. The nearshore zoom uses a shorter water
  depth scale (8 m) because its bay is under 10 m deep.
- **Sea colours**: the floor albedo is a sand-to-dark-brown ramp (like sea-floor relief models), the
  water body a depth ramp from turquoise (0-5 m) to deep blue (160 m+, `render.water.color_stops`),
  and the surface veil a pale blue; the floor shows through the shallows and fades under the blue
  offshore (`max_alpha` 0.78 keeps some relief visible). Waves are long gentle swells (520-90 m,
  slopes 0.03-0.01) with a soft ripple field; the first, shorter and steeper set was called harsh.
- **Water surface** (`render.water_surface`, `render.water.surface`): per pixel, the ray's crossing
  of sea level inside the block is the surface point and the hit below it the floor; the floor is
  displaced a little by the wave slope (refraction), faded to the water colour by its true depth,
  then veiled by a semi-transparent surface (`base_alpha` plus Fresnel) whose brightness follows
  the wave tilt towards the sun, amplified (`wave_shading`, a mirror reflects a sky gradient, so
  a few degrees of tilt are a lot) with white added on crests. Waves: four directional swells
  (260-40 m, slopes 0.05-0.02, from about 200°) plus a smoothed random ripple field; the nearshore
  zoom scales the wavelengths by 0.3. Sun glints were tried and dropped: with the sun and the
  camera both near 35° elevation the mirror direction never meets the sun. Cut walls below sea
  level keep the depth tint only (the water column in cross-section).
- **Schematic** (`schematic.py`, after the TREC "sampling concept" figure in the design PDF): one
  section on a symmetric-log distance axis (linear within 20 m of the shoreline), so the
  20-255 m soil positions and the 1.2 / 12 km stations sit side by side; two stacked axes with
  their own vertical scales (elevation above, depth below); WorldCover classes as a ribbon on
  the land surface; soil tubes (green) and sediment tubes (orange) on the ground; water samples
  as cones at `cable_out` depth plus the land site's 0-3 m shoreline water samples; a boat per
  Tara station; the Inland / Coastal transition zone / Offshore band.
- **Land cover**: ESA WorldCover 2021 (the layer of the R repository), read as a window of the
  3-degree AWS COG, nearest-neighbour onto the grid, legend colours softened in HSV (0.72/0.96
  for 2D as in R, 0.55/0.92 under the path tracer), water class in the sea colour. The 3D block
  uses it by default (`render.albedo.land_mode: landcover`); `hypsometric` is the fallback.
- **The other Athens transects** (2026-11-02): `psatha` (041-01 + TE_172), `salamina` (041-02 +
  TE_173), `thoriko` (041-04 + TE_178, one waypoint south of the headland east of the marina),
  `schinias` (041-05, no station within 40 km: the path continues `profile.sea_extension_m` =
  2 km offshore along the transect bearing). The camera azimuth is `auto`: the seaward bearing of
  the land transect (inland -> shoreline), snapped to 90° so the front wall is a plate edge
  (Anavyssos, Thoriko, Schinias 180°; Salamina 0°; Psatha 270° by snapping, overridden to 295°
  with exaggeration 1.5 and a higher camera because 576 m peaks wall the bay). The sun follows
  the camera (`offset_from_camera_deg` -55). The camera distance also fits the plate's depth and
  the block's height (`fill_vertical`), which deep plates like Salamina's needed.
- **Side view** (`--view side`, `block_<quality>_side`): the plate-edge bearing that lays the path
  (first point -> last point) out along the viewer's right, chosen by `render.side_view_azimuth`:
  score = |left-to-right fraction| + 0.5 if it runs land -> sea + a penalty for a camera on the
  land side of the coast (its front wall would cut the land, not the sea floor). A plain "turn
  90° from the seaward bearing" sent the long axis of the east-west plates (Anavyssos, Thoriko)
  into the depth, a tall narrow strip. Result: Psatha 0°, Salamina 90°, Schinias 270° (turned,
  land left / sea right); Anavyssos and Thoriko 180°, their front view, which already runs
  TE_177 -> shore and shore -> TE_178 left to right, so `block_final_side` repeats it there.
  `render.side_overrides` and the beside-the-pin labels apply only when the view is turned.
  `camera.side_azimuth_deg` in a transect's `render` overrides the choice.
- **TE_178 depth**: the site export says 22 m (sounder 22-35 m) but EMODnet has 76 m at the site
  position; the labels show both (`profile.floor_label`) whenever they differ by more than 25%.
  `check_profile` ignores the first 50 m off the shoreline (beach cells in a 30 m DEM).
- **Nearshore zoom** (`anavyssos_nearshore`): 041-03 + TE_176 only, 5 m cells, smaller margins,
  camera 175°/22°, exaggeration 3, as per-transect overrides in `transects.yml`.
- **Sea floor where EMODnet is silent** (the shallow bay, about 12 000 cells in the zoom):
  a distance-weighted slope from the shoreline cap to the nearest known depth, then a 40 m
  Gaussian over the sea floor (`grid.bathy_smooth_sigma_m`). Before that the nearest-neighbour
  fill produced flat terraces that the 3D render showed as shelves. This floor is modelled.

- **All transects** (2026-11-02): the two exports the user named hold every LSI site (114, 20
  countries) and every Tara station (188); the site table is derived from their `Site ...`
  columns (no separate site export any more). `transects.yml` keeps the hand-made Athens
  entries and an `auto` rule: every other LSI site becomes `<alias lower-cased>`, with its
  *related* Tara stations (`Site related sites`) within 20 km, nearest first, at most two; if
  none is that close, the nearest station within 10 km (`max_unrelated_km`; 20 km pulled
  La Tremblade 1 to a station behind a peninsula); else no station and the 2 km offshore
  extension (13 sites: Baie de Veys 2, Texel 2, Sylt 1, Klaipeda 1-2, Tallinn 2, Southampton
  1-2, Ferrol 2, Porto 1, Cadiz 3, Malaga 1, Barcelona 3, La Tremblade 1). Distance beats the
  relation because the relation is wrong in places (Barcelona 1/2 and Banyuls 1/2 list each
  other's stations, 100 km away). `auto.overrides` tweaks a generated transect in place (wider
  margins for 12 whose sea route left the box); an explicit entry with the same `land_site`
  replaces it. `build_grid` rebuilds when the box changed and the get-scripts re-fetch a
  window that falls short of it. The grid CRS is the path's UTM zone (`grid.crs: auto`); the DEM is read as a
  window of the GLO-30 COG like WorldCover (full tiles would be gigabytes in Dropbox); the
  WorldCover classes other than water count as land (`grid.landcover_as_land`, for polders and
  flats under the 0.5 m threshold); samples without `shore_distance` (Galway 1 shoreline,
  Faro 1 inland) take the zone's nominal distance; the block's exaggerations shrink when it
  would be taller than 0.45 x the plate's shorter side (`render.max_block_height_frac`:
  Bergen's fjord at x6 was a tower); a station without a `station_name` is titled by its
  alias only (`profile.station_title`); the soil label moves beside its markers when a
  station pin is close (Bergen, Roscoff 2). `make all Q=<quality>` runs every transect in
  sequence (one GPU render at a time), `make maps T=<name>` one of them with both views.
- **Sea routing** (`profile.route_around_land`, `profile.route_around_land` in `profile.py`): the
  first batch had 55 transects whose straight marine legs crossed land (peninsulas, islands,
  flats). A leg with land cells on it (ignoring the first 50 m off the shoreline) is replaced by
  the shortest route through sea cells: Dijkstra (scipy csgraph) on the land mask reduced to
  50 m cells and dilated by 100 m, retried without the dilation when a channel is too narrow
  ("hugging the shore"), endpoints snapped to the nearest sea cell within 1.5 km, the route
  simplified to a few waypoints, given up when it is longer than 5 x the straight leg or
  leaves the grid box (then the straight leg stays and `check_profile` warns; the fix is a
  waypoint or wider margins in `auto.overrides`). A route that follows the grid edge is
  logged as such (Tallinn 1 got 5 km margins for it). A station-less transect whose 2 km
  offshore leg hits land turns its bearing by up to 45°, or, with land all round (an estuary,
  a ria: La Tremblade 1, Ferrol 1), ends on the nearest sea cell. `s` along the profile is
  then the sailed distance; the subtitle keeps the straight distance to the station. The
  routed path is `derived/<t>/path.json`, read by the schematic for the shoreline-water
  position. Outcome of the batch: 45 legs routed in the default boxes, 7 more with 5 km
  margins; still crossing land (no sea route in the box, or a channel the land masks close):
  Baie de Veys 1 (TE_007 behind the Carentan flats), Banyuls 2, Bilbao 2 (TE_086 in the
  Abra), Helsingor 1, Thrace 4 (a lagoon site), La Tremblade 1 (the Seudre estuary, under
  100 m wide). Their profiles show land on the marine leg; they need hand-set waypoints.

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
  are JSON lists of `{attribute, value, unit}`. The marine export (29 118 rows) takes 11 s to
  parse, so `sites.read_export` caches each export as a pickle in `data/derived/_cache`
  (invalidated by the xlsx mtime). `Site type` is `LSI` or `TARA`; the marine export's "site"
  name is misleading, it is a sample export too. Helsingor 2 has a terrestrial sample 3 km
  inland, so its land path is 3 km long.
- pandas refuses `rows.loc[mask, float_col] = <empty string series>`: fill a float column with
  `fillna(series.astype(float))`, never with a masked assignment of a mapped series.
- `stations.csv` stores the cast depths as `d1;d2`; a one-cast station comes back from pandas as a
  float, so `profile.casts_of` accepts numbers too (the first renders had no depth markers).
- The shoreline water samples (0-3 m) sit where the filled floor is only 0.5 m deep; the schematic
  draws that marker no deeper than half the local floor depth.
- Symlog axes: label positions near the panel edges are decided in display space after
  `fig.canvas.draw()`; zone-band labels sit at the visual midpoint, not the arithmetic one.
- Water tint: a first version scaled by the ray's path length through water was invisible at these
  exaggerations; depth under the hit point is what reads. `render.pixel_rays` is the inverse of
  `world_to_pixel`, so hit points come from `origin + dir * depth`.
- Wave shading: a Lambert term on 3° wave slopes changes brightness by under 1% (invisible), a
  gain of 8 was a blotchy mess with blown crests; 3.5 with `base_alpha` 0.22 reads as water. The
  ripple field is sampled with cubic splines on a grid of a quarter of its scale; linear sampling
  on a coarser grid showed as a dotted moiré. Check the sea on 1:1 crops of `block_final.png`,
  the 1600 px preview hides all of this.
- Symmetric padding keeps the real cells' world X/Z unchanged (plate centred on the origin), so the
  overlay projection needs no offset; only Y changes with the sea factor (`render.world_y`).
- Page furniture: credits and the symbol key sit top right (sky); the legend row at the bottom left
  must clear the front wall, hence `fill` 0.84 and `target_shift_north_frac` -0.10 for the main
  block (negative lifts the block in the frame). At a low camera the front wall is tall in the
  frame; check the axis labels at the left edge and the legend row after any camera change.
- `make_camera` target shifts are camera-relative ("north" = away from the camera, "east" = to its
  right); a sign slip here dropped the block onto the legend once. The vertical axis stands at the
  block corner nearest the camera, moved towards it and to the viewer's left (`front_left_corner`),
  so its labels stay on the background for diagonal views too; the sea-level line follows the edge
  nearest the camera (`front_edge`).
- Shallow sandy floors mixed with turquoise read grey-beige until the water tint was made stronger
  in the first metres (`depth_scale_m` 8, `max_alpha` 0.85).
- Dropbox: `PYTHONDONTWRITEBYTECODE=1` in the Makefile, `*conflicted copy*` ignored, no venv here.
- Timings (M1 Max): preview 1600x900 x 6 frames about 1 s; final 3840x2160 x 64 frames about 12 s (plus a few seconds for the page).

## Verify

```bash
make check                                   # forge3d 1.44.0 gpu True
make list                                    # 115 transects (6 explicit + 109 generated); 15 s the first time (xlsx -> pickle cache)
python scripts/10_build_grid.py --force      # prints z range, land fraction, station depths vs recorded
python scripts/20_profile_2d.py              # "leg ... routed through N waypoints", "marine legs stay at sea", station depth checks
python scripts/30_render_3d.py               # "median rendered/projected = 1.000x" and "path: 100% visible"
make all Q=preview                           # every transect, about 40 min; failures land in output/failed.txt
```

## Status and next steps (2026-11-09)

- Every transect has final figures (`make all Q=final`, 2026-11-02/03; 230 renders at 9 s). Code
  and config are committed ("make maps for all transects of trec"); the figures and data are not
  (ignored). The grids all match the current config (checked 2026-11-09).
- Ten transects still have marine legs over land and need hand-set waypoints in
  `auto.overrides` (see Sea routing): Baie de Veys 1, Banyuls 2, Helsingor 1, Rostock 1,
  Sylt 2, Roscoff 3, Thrace 4, Kristineberg 3, Bilbao 2, La Tremblade 1. After a waypoint:
  `make maps T=<name> Q=final` (the bathymetry and the windows re-fetch if the box grew).
- Only about a dozen of the 109 generated figures were reviewed by eye; a review pass would look
  at label collisions on short legs (Naples 2), the station labels over hills (Kotor 1), the
  Wadden Sea flats counted as land (Sylt 2, Roscoff 3), and whether the 2 km offshore extension
  is wanted for the 13 station-less transects.
- Candidates for deletion in `data/`: the two older exports (`TREC_sample_data_export.xlsx`,
  `TREC_sample_data_athens_marine_metagenome.xlsx`) and the four full GLO-30 tiles in
  `data/raw/dem` (57 MB); nothing reads them.
- The Mac moved from EET to CET between the batch and the check, so file times read one hour
  earlier than the batch logs say.
- Possible follow-ups the user has not asked for: a contact-sheet of all 115 blocks, an index
  page listing transects with their stations, EPSG:3035 for a pan-European overview map.

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

### 2026-11-02 — schematic figure, land cover, sample depths in 3D, nearshore zoom

Request: "now take a look at screenshot. Create a similar one with real data. you can crop/log
transform the x axis of the distance in order to move closer the LSI with the ocean samples. Also
add the land cover data as in trec-coastal repo. And the 3d also show the depth of the sea samples
all in one map. Recap: make another 2-3 figures based on instuctions"

The screenshot is the TREC "sampling concept" figure (Box 2 of the design document). Done:
`scripts/02_get_landcover.py` (WorldCover window), land cover in the grid, the profile tables and
the 3D albedo with a class legend; `scripts/21_profile_schematic.py` and `trecmaps/schematic.py`;
pins with sea level, sampling depth and sea floor per Tara station in the 3D block; the
`anavyssos_nearshore` transect with per-transect overrides; the nearshore sea-floor blend and
smoothing. Figures: `output/anavyssos/schematic`, `output/anavyssos/block_final` (land cover +
pins), `output/anavyssos_nearshore/block_final` and `schematic`. Final renders take about 12 s.

### 2026-11-02 — profile overlap, depth in the 3D block

Request: "fix the profile.png ovelapping maps. Again the 3d doesnt have a z axis for depth fo the sea
like the terrain with dem in the block_preview."

Done: the land-transect zoom became the first panel of `profile.png` (own vertical scale, no
overlap); the 3D block became a solid with cut walls, a depth-tinted water body (also on the walls),
a sea-level line, a vertical axis in metres and twice the vertical exaggeration below sea level;
credits and the symbol key moved to the top right. Both transects re-rendered (`block_final`).

### 2026-11-02 — water surface

Request: "the water is not visible. Can you add something like atransparency and wavy texture to
visualise depth beter?" then "is it ready?"

Done: `render.water_surface` (semi-transparent wavy surface at sea level with depth attenuation,
refraction and amplified wave shading; settings under `render.water.surface`, nearshore override
`wave_scale: 0.3`), the sea legend bar shows the veiled colour, both `block_final` recomposed from
the cached renders. Tuned on 1:1 crops: the first two settings were invisible and then far too harsh.

### 2026-11-02 — section view

Request: "the angle and azimuth is also an issue: make the block view closer to profile view because
now it is more like a 2d map."

Done: camera from the south at 20° (zoom 22°) instead of 38°, fill and target shift retuned so the
axis and the legend clear the taller front wall, the soil-transect label moved above the inland-most
position. Both `block_final` re-rendered.

### 2026-11-02 — sea colours

Request: "the waves are harch, make them a lot more smouth and add a slight blueish color to the sea.
make the seacolor from tirquaz to deep blue based on depth and the seaflor dark brown like the 3d
models of seafloor."

Done: brown sea-floor albedo, depth-ramped water colour (`color_stops`), bluish surface veil, long
smooth swells with a softer ripple field; both `block_final` re-rendered (the floor colour is baked
into the render).

### 2026-11-02 — the other Athens transects

Request: "move to the other Athens transects like anavyssos"

Done: `psatha`, `salamina`, `thoriko`, `schinias` in `transects.yml`; `sea_extension_m` for a
transect without a station; auto camera azimuth from the land transect with a tied sun and a
distance that fits the plate's depth; the axis and the sea-level line placed by azimuth; the
recorded-vs-EMODnet depth note in the labels; the Thoriko waypoint; stronger shallow-water tint.
Downloads: GLO-30 tiles N38E023, N37E024, N38E024; EMODnet and WorldCover windows per transect.
All six transects have profile, schematic and block figures.

### 2026-11-02 — side view

Request: "make also a version of each 3d map that has the transect from left to right in all transects."

Done: `--view side` in `30_render_3d.py` (Makefile `V=side`), `render.side_view_azimuth`,
`render.side_overrides`, label placement beside the pin for turned views; `block_final_side`
rendered for all six transects (Anavyssos and Thoriko: the front view, see Decisions).

### 2026-11-02 — all transects

Request: "based on TREC_sample_data_export-all-MetaBGT and
TREC_site_data_export_all-sites_marine_matagenome . do the maps for all transects. one by one to
limit gpu usage"

Done: the two exports replace the three old ones (`paths`), sites derived from them, the `auto`
transect rule (109 generated + 6 explicit), UTM grids, DEM windows, WorldCover in the land
mask, nominal shore distances, the block-height cap, station titles, the soil-label rule,
`make list | maps | all`. Checked by eye on Roscoff 2 (Atlantic, two stations, 90° camera) and
Bergen 1 (fjord, 250 m deep at 600 m). Batches (2026-11-02/03): `make all Q=preview` for every
transect (115, none failed; 55 had marine legs over land, which led to the sea routing, the
wider-margin overrides and the 10 km rule for unrelated stations), then `make all Q=final`:
230 final renders at 9 s each, every transect has profile, schematic, `block_final` and
`block_final_side` (PNG + PDF). Still crossing land (see Sea routing): Baie de Veys 1,
Banyuls 2, Helsingor 1, Rostock 1 (its offshore end), Sylt 2 and Roscoff 3 (tidal flats
counted as land on the last stretch), Thrace 4, Kristineberg 3, Bilbao 2, La Tremblade 1.
Not reviewed: most of the 109 generated figures were only spot-checked (Roscoff 1-2, Bergen 1,
Tallinn 1, Texel 2, Naples 2, Barcelona 2, Kotor 1, Aarhus 1, Sylt 2, Helsingor 2); label
collisions near short legs (Naples 2) are left for per-transect overrides.
