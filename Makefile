# Run inside the conda environment:  conda activate trec-maps
T  ?= anavyssos
Q  ?= preview
V  ?= front
PY ?= python
export PYTHONDONTWRITEBYTECODE = 1   # no __pycache__ churn inside Dropbox

help:
	@echo "make env | check | list | dem | bathy | landcover | grid | profile | schematic | render | final | pilot | maps | all   (T=$(T) Q=$(Q) V=$(V): front or side)"
env:
	CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yml
env-update:
	CONDA_CHANNEL_PRIORITY=strict conda env update -n trec-maps -f environment.yml --prune
check:
	$(PY) -c "import forge3d as f3d, rasterio, shapely, pyproj; print('forge3d', f3d.__version__, 'gpu', f3d.has_gpu())"
dem:
	$(PY) scripts/00_get_dem.py --transect $(T)
bathy:
	$(PY) scripts/01_get_bathymetry.py --transect $(T)
landcover:
	$(PY) scripts/02_get_landcover.py --transect $(T)
grid:
	$(PY) scripts/10_build_grid.py --transect $(T)
profile:
	$(PY) scripts/20_profile_2d.py --transect $(T)
schematic:
	$(PY) scripts/21_profile_schematic.py --transect $(T)
render:
	$(PY) scripts/30_render_3d.py --transect $(T) --quality $(Q) --view $(V)
final:
	$(PY) scripts/30_render_3d.py --transect $(T) --quality final --view $(V)
pilot: dem bathy landcover grid profile schematic render
list:
	$(PY) scripts/_common.py
maps: dem bathy landcover grid profile schematic          # everything for one transect, both 3D views at Q
	$(PY) scripts/30_render_3d.py --transect $(T) --quality $(Q) --view front
	$(PY) scripts/30_render_3d.py --transect $(T) --quality $(Q) --view side
all:                                                       # every transect one after the other (one render at a time)
	@mkdir -p output; : > output/failed.txt
	@for t in $$($(PY) scripts/_common.py | cut -d' ' -f1); do echo "=== $$t"; \
	  $(MAKE) --no-print-directory maps T=$$t Q=$(Q) || echo "$$t" >> output/failed.txt; done
	@echo "failed: $$(cat output/failed.txt | tr '\n' ' ')"
.PHONY: help env env-update check dem bathy landcover grid profile schematic render final pilot list maps all
