#!/bin/bash

################################################################################

# Perform standard analysis on Rooij/TeichAllCells
processors=1
commands="run_LigRec_stAnalysis-CORES=${processors}"
mem=64G
run_days=1
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
script_name=HCM_SCS_LR_2021_12_run_Seurat_LR_Task.sh
#last_jobid=$(sbatch --dependency=afterany:${jobid_Teich_Cl} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=3-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=${run_days}-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})

# Now create overview heatmaps of receptor expression
processors=1
commands="LigRec_Expr-CORES=${processors}"
mem=64G
run_days=1
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
script_name=HCM_SCS_LR_2021_12_run_Seurat_LR_Task.sh
#last_jobid=$(sbatch --dependency=afterany:${jobid_Teich_Cl} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=3-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=${run_days}-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})


################################################################################