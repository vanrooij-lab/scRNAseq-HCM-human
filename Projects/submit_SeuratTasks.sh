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


# Run the default analysis for the pooled datasets
commands=runf_all_default
processors=1
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
sbatch --output=slurm-${commands}-%x.%j.out --job-name=${commands} -c ${processors} --time=1-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh
#
# Run the default analysis for the pooled datasets
# DEPENDENCY SET FOR PREVIOUS JOB (MANUALLY!)
commands=runf_all_default_cl
processors=20
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
sbatch --dependency=afterany:7707272 --output=slurm-${commands}-%x.%j.out --job-name=${commands} -c ${processors} --time=2-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh

#####

# little test (parsable makes output ID only)
last_jobid=$(sbatch --parsable test.sh)

################################################################################
# ANALYSIS OF TEICHMANN DATA SEPARATELY

# Note: Analysis "run_separate" not required!
# We go straight to ..

# Integration
commands="run_batch_corr_sep-TEICHMANNonly"
processors=1
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
mem=250G
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=3-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
# |
# V
# Analysis (plots + DE)
# Depend on above
commands="run_batch_corr_sep_nowplot_and_DE-TEICHMANNonly_Int1c-CORES=8"
processors=8
mem=150G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
last_jobid2=$(sbatch --dependency=afterany:${last_jobid} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)


# CONTINUE BELOW !!!!!

# STILL TO EXECUTE, RID2-like
# Now here, the analysis is required ..

commands="run_separate-TEICHMANN-SETTINGS_RID2l"
processors=1
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
# |
# V
commands="run_separate_nowclusterDE-TEICHMANNonly_RID2l-CORES=8"
processors=8
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
last_jobid2=$(sbatch --dependency=afterany:${last_jobid} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
#last_jobid2=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)

################################################################################
# Also do clustering DE for HU

commands="run_separate_nowclusterDE-HUonly_RID2l-CORES=8"
processors=8
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
mem=150G
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)


################################################################################
# Regulons for Hu and Teichmann
commands="run_regulon_step1-HUonly_RID2l-CORES=6"
processors=6
mem=128G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
script_name=run_SeuratRegulonTask.sh
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})

# Teichmann
commands="run_regulon_step1-TEICHMANNonly_RID2l-CORES=5"
processors=5
mem=255G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
script_name=run_SeuratRegulonTask.sh
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=3-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})

# Regulon analysis for Rooij (also did this locally, but convenient as a test)
commands="run_regulon_step1-ROOIJonly_RID2l-CORES=5"
processors=5
mem=64G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
script_name=run_SeuratRegulonTask.sh
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})




#####

# Testing new separation sign for commands (dash, -)
commands="fake1-fake2-fake3_hoi"
processors=1
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=0-00:05:00 --mem=1G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)







