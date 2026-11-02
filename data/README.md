# data/

Nothing in this folder is committed except this file (TREC data are shared only with
the consortium's consent; the downloads are rebuilt by the scripts).

- `TREC_site_data_export_all-sites.xlsx`, `TREC_sample_data_export.xlsx`,
  `TREC_sample_data_athens_marine_metagenome.xlsx`: exports from the TREC data portal
  (sites, soil/sediment samples, Athens marine samples). Title in row 1, header in row 2.
- `raw/dem/`: Copernicus GLO-30 tiles, `scripts/00_get_dem.py`.
- `raw/bathy/`: EMODnet Bathymetry WCS extracts, `scripts/01_get_bathymetry.py`.
- `derived/<transect>/`: merged land-sea grid (`z.tif`, `grid.npz`, `grid.json`),
  profile tables (`profile.csv`, `land_samples.csv`, `stations.csv`) and render caches.
