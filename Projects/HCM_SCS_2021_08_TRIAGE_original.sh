#!/bin/bash

# Script path
triage_path=/Users/m.wehrens/Documents/git_repos/TRIAGE/python

# Go to dir with data
cd /Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/TRIAGE/

# Run triage
# I modified the script, argument parsing had an issue
python ${triage_path}/disc.py -i exp_matrix_counts.tsv -o exp_matrix_triage_out.tsv -r ${triage_path}/human_rts.txt -f 1 -l 0 -p 0
  # Setting log-transform and pseudo+1 to false/0, since it was already performed in R by Seurat