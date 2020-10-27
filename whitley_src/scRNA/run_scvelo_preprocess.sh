#!/bin/bash
set -e
# note that scvelo_conda environment in anaconda must be active
source activate scvelo_conda
Rscript laura_prep_data_for_scvelo.R > laura_prep_data_for_scvelo_out.txt 2>laura_prep_data_for_scvelo_stderr.txt

source deactivate
