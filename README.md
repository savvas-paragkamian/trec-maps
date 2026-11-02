# trec-maps

Map scripts for the TREC expedition (EMBL, *Traversing European Coastlines*): land-sea
**profile maps** of the sampling transects, built from a DEM on land and EMODnet bathymetry
at sea, as 2D profiles (matplotlib) and 3D path-traced blocks (forge3d). Python only.

Every TREC land-sea transect is a map: the five Athens transects are configured by hand
(`config/transects.yml`), the other 109 LSI sites of the data exports become transects
automatically, paired with their Tara stations (`make list` prints all 115). The Athens entries:

| transect | land site | Tara station(s) |
|---|---|---|
| `anavyssos` (+ `anavyssos_nearshore`, a 5 m zoom) | 041-03, HCMR supersite | TE_176 (1.2 km, 9 m), TE_177 (10 km, 131 m) |
| `psatha` | 041-01, Gulf of Corinth | TE_172 (1.6 km, 77 m) |
| `salamina` | 041-02, Saronic Gulf | TE_173 (0.6 km, 18 m) |
| `thoriko` | 041-04, Thoriko marina | TE_178 (6 km; 22 m recorded, 76 m in EMODnet at that position) |
| `schinias` | 041-05, Marathon | none nearby; the profile continues 2 km offshore along the transect bearing |

| output | what |
|---|---|
| `output/<transect>/profile.png`, `.pdf` | elevation / depth along the transect path in three panels: the land transect in detail (own vertical scale), the nearshore in metres and the offshore in kilometres, with TREC zones, soil positions and the stations' sampling depths |
| `output/<transect>/schematic.png`, `.pdf` | the TREC sampling-concept figure with real data: symmetric-log distance axis around the shoreline, ESA WorldCover along the land surface, soil and sediment tubes, water samples at their sampling depths, the inland / coastal transition / offshore band |
| `output/<transect>/block_<quality>.png`, `.pdf` | path-traced solid block seen from the south at a low angle (a section through the coast rather than a map): land (WorldCover classes) and sea floor under a semi-transparent, wavy sea surface (the floor fades to the water colour with depth and shows through the waves), cut walls down to a base depth with the water column in cross-section, a sea-level line and a vertical axis in metres at the front corner, the path, the soil positions and a pin per Tara station (sea level, sampling depth, sea floor); vector overlay in the PDF. Depths are exaggerated twice as much as land (stated in the legend) |

## Setup

A conda environment from free community channels only (`conda-forge`, `bioconda`,
`nodefaults`), the convention of the TREC repositories, so it works with Miniforge or with a
Miniconda whose `.condarc` avoids the Anaconda `defaults` channel. The environment lives in
conda's envs directory, never inside this Dropbox folder.

```bash
CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yml   # or: make env
conda activate trec-maps
make check          # forge3d 1.44.0 gpu True
```

forge3d renders on the local GPU (Metal on Apple silicon); the 2D profile does not need it.

## Data

Put the two TREC data-portal exports in `data/` (they are not committed, see `data/README.md`):
`TREC_sample_data_export-all-MetaBGT.xlsx` (the soil and sediment samples of every LSI site)
and `TREC_site_data_export_all-sites_marine_matagenome.xlsx` (the seawater samples: Tara
stations and the LSI shoreline water). Sites, positions, shore distances, sampling depths and
the site-to-station relations are read from them (each export is cached as a pickle under
`data/derived/_cache`, the marine one takes 11 s to parse). Everything else is downloaded by
the three get-scripts:

- **Land:** Copernicus DEM GLO-30; only the window over the grid is read from the 1° AWS
  cloud-optimised tiles (`scripts/00_get_dem.py`).
- **Sea floor:** EMODnet Bathymetry DTM (2024 release, 1/16 arc minute, about 115 m) through
  the WCS service, requested on its native grid for the transect's box (`scripts/01_get_bathymetry.py`).
- **Land cover:** ESA WorldCover 2021 v200, 10 m; only the window over the grid is read from the
  3-degree AWS tiles through GDAL range requests (`scripts/02_get_landcover.py`).

## Run

```bash
make list                  # every transect: name, land site, stations, subtitle
make pilot                 # dem, bathy, landcover, grid, profile, schematic, render (preview) for T=anavyssos
make pilot T=roscoff_lsi_2
make maps T=bergen_lsi_1 Q=final   # everything for one transect, both 3D views
make all Q=preview         # every transect one after the other (one render at a time); failures in output/failed.txt
make final                 # 3840 x 2160 block, 64 frames
make final V=side          # the side view: the transect running left to right, land left where possible (block_final_side)
make render Q=draft        # 2560 x 1440
python scripts/30_render_3d.py --quality preview --tag _try --set camera.azimuth_deg=200 exaggeration=4
python scripts/30_render_3d.py --quality final --recompose    # reuse the render, redo the page
```

Steps, each a thin script over `trecmaps/`:

| step | script | writes |
|---|---|---|
| 0 | `00_get_dem.py` | `data/raw/dem/glo30_<transect>_<tile>.tif` (the only DEM network code; `--dry-run`, `.part` then rename, skips existing files) |
| 1 | `01_get_bathymetry.py` | `data/raw/bathy/emodnet_mean_<transect>.tif` + `.json` (re-fetched when the grid box grows) |
| 2 | `02_get_landcover.py` | `data/raw/landcover/worldcover_2021_<transect>_<tile>.tif` |
| 10 | `10_build_grid.py` | `data/derived/<transect>/grid.npz`, `z.tif`, `grid.json`: DEM, bathymetry and land cover on one metric grid, merged with a land mask |
| 20 | `20_profile_2d.py` | `path.json` (the path, marine legs routed around land), `profile.csv`, `land_samples.csv`, `stations.csv`, `output/<transect>/profile.png/.pdf` |
| 21 | `21_profile_schematic.py` | `output/<transect>/schematic.png/.pdf` (from the step-20 tables) |
| 30 | `30_render_3d.py` | `output/<transect>/block_<quality>[_side].png/.pdf` and a render cache for `--recompose`; `--view front` (from the sea) or `side` (transect left to right) |

## Configuration

- `config/params.yml`: paths, sources (with the WorldCover legend), grid (CRS, cell size,
  margins, land/sea merge), profile and schematic (step, zones, figures), render (camera, sun,
  exaggeration, quality presets, colours; `albedo.land_mode` is `landcover` or `hypsometric`).
- `config/transects.yml`: explicit entries by TREC site IDs only (`land_site`, `marine_sites`
  in path order, optional `waypoints` to keep a marine leg off the beach), plus optional
  `grid`, `profile`, `schematic` and `render` overrides merged over `params.yml`; and the
  `auto` rule that turns every other LSI site of the samples export into a transect named
  after its alias (`roscoff_lsi_2`): its related Tara stations within `max_station_km`
  (nearest first, at most `max_stations`), else the nearest station within
  `max_unrelated_km`, else none and the profile continues `profile.sea_extension_m`
  offshore. `auto.overrides` holds per-transect tweaks (wider `grid.margin_m` where the sea
  route needs room). Coordinates come from the xlsx exports at run time.

To tune a generated transect further (waypoints, camera), add it to `auto.overrides` or copy
it as an explicit entry with the same `land_site`. A marine leg that would cross land follows the
shortest sea route instead (`profile.route_around_land`: Dijkstra over the land mask on a 50 m
grid, 100 m off the shore, retried along the shore when a channel is narrower); the routed
path is saved as `data/derived/<transect>/path.json`. `20_profile_2d.py` warns when a leg still
crosses land, which happens when the route would leave the grid box (a station behind a
peninsula): add a waypoint or wider `grid.margin_m` in an explicit entry and rerun from
`01_get_bathymetry.py`. The 3D camera
sits over the sea on the transect's seaward bearing by default (`camera.azimuth_deg: auto`,
snapped to a plate edge), with the sun tied to it; a transect can override camera, exaggeration
and water settings, as `psatha` does for its mountain-walled bay.

## How the grid is built

The path (inland end, shoreline, waypoints, stations) is boxed with margins in the path's
UTM zone (`grid.crs: auto`; the Athens figures were first made in EPSG:2100) at 10 m cells,
coarsened when the box would exceed `grid.max_cells`. GLO-30 is resampled bilinearly, EMODnet
cubically after its land cells are set to nodata, WorldCover by nearest neighbour. Land is
DEM > 0.5 m (speckle removed) or any WorldCover class but water, which keeps polders and
flats below the threshold on land. Sea cells take the EMODnet floor where it exists; where the
115 m grid stops short of the coastline the floor slopes from −0.5 m at the shoreline down to
the nearest known depth (distance-weighted), and a 40 m Gaussian keeps the coarse source
cells from reading as terraces. That nearshore floor is modelled, not measured. GLO-30
heights are above EGM2008, EMODnet depths below LAT; the Saronic tidal range is about
0.1 m, so the offset is ignored. In the 3D block the exaggerations (x3 land, x6 sea) shrink
when the block would be taller than `render.max_block_height_frac` of the plate's shorter
side, as in the fjords (Bergen: x1.25 / x2.5).

## Sources and citations

- Copernicus DEM GLO-30, ESA (produced by Airbus/DLR), https://spacedata.copernicus.eu/collections/copernicus-digital-elevation-model
- EMODnet Bathymetry Consortium (2024). *EMODnet Digital Bathymetry (DTM 2024)*. https://doi.org/10.12770/cf51df64-56f9-4a99-b1aa-36b8d7b743a1 (CC BY)
- Zanaga, D. et al. (2022). *ESA WorldCover 10 m 2021 v200*. https://doi.org/10.5281/zenodo.7254221 (CC BY 4.0)
- Sample and site metadata: TREC data portal (EMBL). Publication of TREC data needs the consortium's consent.
- Rendering: forge3d (https://github.com/milos-agathon/forge3d).
