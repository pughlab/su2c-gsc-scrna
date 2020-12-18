#!/bin/bash
set -e
# note that scvelo_conda environment in anaconda must be active
source activate scvelo_conda
jupyter nbconvert --execute scvelo_merge_loom.ipynb --to html --output-dir reports/scvelo --ExecutePreprocessor.timeout=-1
jupyter nbconvert --execute scvelo_subset_loom_GBM_GSCs.ipynb --to html --output-dir reports/scvelo --ExecutePreprocessor.timeout=-1
jupyter nbconvert --execute scvelo_subset_loom_GSCs.ipynb --to html --output-dir reports/scvelo --ExecutePreprocessor.timeout=-1
jupyter nbconvert --execute scvelo_diffusion_map_GBM_GSCs.ipynb --to html --output-dir reports/scvelo --ExecutePreprocessor.timeout=-1
jupyter nbconvert --execute scvelo_diffusion_map_GSCs.ipynb --to html --output-dir reports/scvelo --ExecutePreprocessor.timeout=-1
source deactivate
