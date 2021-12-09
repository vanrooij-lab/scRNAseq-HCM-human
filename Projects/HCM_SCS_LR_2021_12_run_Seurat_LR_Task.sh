#!/bin/bash

# Set script dir
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/

# execute the scripts
Rscript --vanilla $script_dir/HCM_SCS_LR_2021-12-07_HPC_LR_analysis.R $commands
