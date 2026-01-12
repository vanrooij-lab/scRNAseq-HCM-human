#!/bin/bash

# Set script dir
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/

# execute the scripts
Rscript --vanilla $script_dir/HCM_SCS_2021_08_SCENIC.R $commands
