#!/bin/bash



# Run the full RaceID2 like analysis, without var. feature selection
commands=runf_all_RID2l
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
sbatch --output=slurm-${commands}-%x.%j.out --job-name=${commands} -c 1 --time=1-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh
  # --gres=tmpspace:500G
  # set c 8 when e.g. doing clustering
  
# Now the clustering of RID2l VAR
commands=runf_all_RID2l_VAR_cl
processors=20
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
sbatch --output=slurm-${commands}-%x.%j.out --job-name=${commands} -c ${processors} --time=1-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh

# split the pre-processed data into three sets
commands=split_datasets
processors=1
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
sbatch --output=slurm-${commands}-%x.%j.out --job-name=${commands} -c ${processors} --time=1-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh
