#!/bin/bash



# Run the full RaceID2 like analysis, without var. feature selection
# This took to long, so wasn't performed
#commands=runf_all_RID2l
#script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
#sbatch --output=slurm-${commands}-%x.%j.out --job-name=${commands} -c 1 --time=1-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh
  # --gres=tmpspace:500G
  # set c 8 when e.g. doing clustering
  
################################################################################  

# This is not entirely up-to-date any more, the sections
# create_septal_all_dataset
# create_septal_all_btypSel_dataset  
# Are now important to first set up what we'll be splitting
  
  
################################################################################  
  
# Run the default analysis for the pooled datasets
commands=runf_all_default
processors=1
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name=${commands} -c ${processors} --time=1-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
# |
# V
# CLUSTER+DE DEFAULT
# Run the default analysis for the pooled datasets
# Then the clustering of the default analysis (if before not ran, use dependency)
processors=20
commands="runf_all_default_cl-cores=${processors}"
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
# With dependency:
sbatch --dependency=afterany:${last_jobid} --output=slurm-${commands}-%x.%j.out --job-name=${commands} -c ${processors} --time=2-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh
# Without dependency
# sbatch --output=slurm-${commands}-%x.%j.out --job-name=${commands} -c ${processors} --time=2-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh
  
################################################################################    

  
# Again analysis, but now for RID2l
commands=runf_all_RID2l_VAR  
processors=1
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name=${commands} -c ${processors} --time=1-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
# |
# V
# CLUSTER+DE RID2l VAR
# Now the clustering of RID2l VAR
processors=20
commands="runf_all_RID2l_VAR_cl-cores=${processors}"
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
# w/o dependency
# sbatch --output=slurm-${commands}-%x.%j.out --job-name=${commands} -c ${processors} --time=1-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh
# w dependency
sbatch --dependency=afterany:${last_jobid} --output=slurm-${commands}-%x.%j.out --job-name=${commands} -c ${processors} --time=1-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh


# SPLIT DATASETS
# split the pre-processed data into three sets
commands=split_datasets
processors=1
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
split_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name=${commands} -c ${processors} --time=1-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)


#####

# little test (parsable makes output ID only)
last_jobid=$(sbatch --parsable test.sh)






################################################################################
# Run all, but with filter for septal cells in Teichmann

# ALL, septal sel Teichmann
commands="run_separate-dataset=ALL.SP-settings=SETTINGS_RID2l"
processors=1
memory=128G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${memory} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
# |
# V
processors=20
commands="run_separate_nowclusterDE-dataset=ALL.SP_RID2l-cores=${processors}"
memory=128G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
jobid_Teich_Cl=$(sbatch --dependency=afterany:${last_jobid} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${memory} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
#last_jobid2=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)


################################################################################
# Now, execute the analyses of the separate datasets
# RACEID2 like settings

# TEICHMANN
commands="run_separate-dataset=TEICHMANNonly-settings=SETTINGS_RID2l"
processors=1
memory=128G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${memory} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
# |
# V
processors=20
commands="run_separate_nowclusterDE-dataset=TEICHMANNonly_RID2l-cores=${processors}"
memory=128G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
jobid_Teich_Cl=$(sbatch --dependency=afterany:${last_jobid} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${memory} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
#last_jobid2=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)

# TEICHMANN !SEPTAL-ONLY!
commands="run_separate-dataset=TEICHMANN.SP.only-settings=SETTINGS_RID2l"
processors=1
memory=128G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${memory} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
# |
# V
processors=20
commands="run_separate_nowclusterDE-dataset=TEICHMANN.SP.only_RID2l-cores=${processors}"
memory=128G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
jobid_Teich_Cl=$(sbatch --dependency=afterany:${last_jobid} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${memory} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
#last_jobid2=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)


# HU
commands="run_separate-dataset=HUonly-settings=SETTINGS_RID2l"
processors=1
memory=64G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${memory} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
# |
# V
processors=20
commands="run_separate_nowclusterDE-dataset=HUonly_RID2l-cores=${processors}"
memory=128G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
jobid_Hu_Cl=$(sbatch --dependency=afterany:${last_jobid} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${memory} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
#last_jobid2=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)

# ROOIJ
commands="run_separate-dataset=ROOIJonly-settings=SETTINGS_RID2l"
processors=1
memory=32G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${memory} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
# |
# V
processors=10
commands="run_separate_nowclusterDE-dataset=ROOIJonly_RID2l-cores=${processors}"
memory=32G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
jobid_Rooij_Cl=$(sbatch --dependency=afterany:${last_jobid} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${memory} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
#last_jobid2=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)

################################################################################
# For supp figures, and for SCENIC input, each patient is also run
# separately, this is done manually through an R file, 
# HCM_SCS_2021_08_Split_Seurat_per_patient.R

################################################################################
# ROOIJ analysis, default settings

# ROOIJ
commands="run_separate-dataset=ROOIJonly-settings=SETTINGS_default"
processors=1
memory=32G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${memory} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
# |
# V
processors=10
commands="run_separate_nowclusterDE-dataset=ROOIJonly_default-cores=${processors}"
memory=32G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
jobid_Rooij_Cl=$(sbatch --dependency=afterany:${last_jobid} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${memory} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
#last_jobid2=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=150G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)


################################################################################
# ANALYSIS OF TEICHMANN DATA, INTEGRATED [ F A I L E D ]

# Note: Analysis "run_separate" not required!
# We go straight to ..

# Integration (This fails because the data is so large.)
commands="run_batch_corr_sep-dataset=TEICHMANNonly"
processors=1
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
mem=250G
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=3-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
# |
# V
# Analysis (plots + DE)
# Depend on above
commands="run_batch_corr_sep_nowplot_and_DE-dataset=TEICHMANNonly_Int1c-cores=8"
processors=8
mem=150G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
last_jobid2=$(sbatch --dependency=afterany:${last_jobid} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)

################################################################################
# ANALYSIS OF ROOIJ DATA, INTEGRATED 
# (We didn't use this -- I'm highly suspicious of these integration methods
# as they also flawlessly integrate and mix different cell types, which obviously
# is undesired behavior.)

# Integration 
commands="run_batch_corr_sep-dataset=ROOIJonly"
processors=1
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
mem=10G
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=1-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)
# |
# V
# Analysis (plots + DE)
# Depend on above
commands="run_batch_corr_sep_nowplot_and_DE-dataset=ROOIJonly_Int1c-cores=8"
processors=8
mem=32G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
last_jobid2=$(sbatch --dependency=afterany:${last_jobid} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=1-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)

################################################################################
# This can be skipped, but used if you wanna re-run some plots

for DATASET in TEICHMANN HU ROOIJ
do
  processors=1
  commands="rerun_plotsOnly-seuratObject_name=${DATASET}only_RID2l-cores=${processors}"
  memory=32G
  script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
  sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=1-00:00:00 --mem=${memory} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh
done

processors=1
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
memory=64G
commands="rerun_plotsMerged-seuratObject_name=all_RID2l_VAR-cores=${processors}"
sbatch --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=1-00:00:00 --mem=${memory} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh
commands="rerun_plotsMerged-seuratObject_name=default-cores=${processors}"
sbatch --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=1-00:00:00 --mem=${memory} --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh

################################################################################
# Regulons analyses for all datasets
# Note: amount of processors is set to amount of patients

# Regulon analysis for Rooij 
processors=5
commands="run_regulon_step1-dataset=ROOIJonly_RID2l-CORES=${processors}"
mem=64G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
script_name=run_SeuratRegulonTask.sh
#last_jobid=$(sbatch --dependency=afterany:${jobid_Rooij_Cl} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})

# Regulons for HU
processors=7 # 4
commands="run_regulon_step1-dataset=HUonly_RID2l-CORES=${processors}"
mem=200G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
script_name=run_SeuratRegulonTask.sh
#last_jobid=$(sbatch --dependency=afterany:${jobid_Hu_Cl} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=2-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})

# (INCLUDING NON-SEPTAL, DON'T USE)
# Regulons for Teichmann 
processors=14
commands="run_regulon_step1-dataset=TEICHMANNonly_RID2l-CORES=${processors}"
mem=255G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
script_name=run_SeuratRegulonTask.sh
#last_jobid=$(sbatch --dependency=afterany:${jobid_Teich_Cl} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=3-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=3-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})

# Regulons for Teichmann !SEPTAL ONLY!
processors=14
commands="run_regulon_step1-dataset=TEICHMANN.SP.only_RID2l-CORES=${processors}"
mem=255G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
script_name=run_SeuratRegulonTask.sh
#last_jobid=$(sbatch --dependency=afterany:${jobid_Teich_Cl} --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=3-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=3-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})


################################################################################
# Correlation analysis

processors=1
commands="correlations_of_interest-CORES=${processors}"
mem=156G
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
script_name=run_SeuratCorrTask.sh
sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=3-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name}

################################################################################
# SCENIC analysis [ D O   N O T   R U N ]
# Don't run this; in the end I used the command line interface;
# see HCM_SCS_2021_08_SCENIC-py-ICL.sh


# R.P1
# R.P2 R.P3 R.P4 R.P5

for patient in R.P1 R.P2 R.P3 R.P4 R.P5
do
    
  processors=10
  commands="prep_corr_genie-patient=${patient}-cores=${processors}"
  mem=32G
  script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
  script_name=run_SeuratScenicTask.sh
  last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=3-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name})

  processors=20
  commands="run_SCENIC-patient=${patient}-cores=${processors}"
  mem=200G
  script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
  script_name=run_SeuratScenicTask.sh
  sbatch --parsable --dependency=afterany:${last_jobid} --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=3-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name}
  # sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=3-00:00:00 --mem=${mem} --export=ALL,commands="${commands}" ${script_dir}/${script_name}

done
################################################################################
# Additionally, for Homer and Lisa, see:

# HCM_SCS_2021_08_LISA.sh
# HCM_SCS_2021_08_HOMER.sh

################################################################################

#####

# Testing new separation sign for commands (dash, -)
commands="fake1-fake2-fake3_hoi"
processors=1
script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/
last_jobid=$(sbatch --parsable --output=slurm-${commands}-%x.%j.out --job-name="${commands}" -c ${processors} --time=0-00:05:00 --mem=1G --export=ALL,commands="${commands}" ${script_dir}/run_SeuratTask.sh)


################################################################################













