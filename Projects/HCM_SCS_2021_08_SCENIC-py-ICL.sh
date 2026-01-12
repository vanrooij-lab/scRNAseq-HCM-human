#!/bin/bash

# ==============================================================================
# Note that loom files are generated in the R files that serve as input
# for this file.

# See R file with pre-processing for comments with more info.

# ==============================================================================
# This code is not run
# 
# It's only here as a reminder how to submit the jobs.

if [[ 0 -eq 1 ]]; then
  
  # Execute this code manually on the cluster to submit the job:
  script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis
  processors=20
  mem=50G
  for patient in R.P1 R.P2 R.P3 R.P4 R.P5 
  # for patient in R.P1
  do
    sbatch --output=slurm-scenic-${patient}-%x.%j.out --job-name=SCN_${patient} -c ${processors} --time=1-00:00:00 --mem=${mem} --export=ALL,patient="${patient}" ${script_dir}/HCM_SCS_2021_08_SCENIC-py-ICL.sh
    # --parsable 
  done
  
  # Now also start the analysis for Hu and Teichmann data
  # First Hu data
  script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis
  processors=20
  mem=50G
  for patient in H.N1 H.N2 H.N3 H.N4 H.N5 H.N13 
  # H.N14 was also not present
  do
    sbatch --output=slurm-scenic-${patient}-%x.%j.out --job-name=SCN_${patient} -c ${processors} --time=1-00:00:00 --mem=${mem} --export=ALL,patient="${patient}" ${script_dir}/HCM_SCS_2021_08_SCENIC-py-ICL.sh
    # --parsable 
  done
  
  # Teichmann probably needs a bit more memory
  script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis
  processors=20
  mem=120G  
  for patient in T.H5 T.H6 T.H3 T.H2 T.H7 T.H4 T.D2 T.D3 T.D4 T.D5 T.D6 T.D7 T.D11
  # T.D1 not present in sep data set
  do
    sbatch --output=slurm-scenic-${patient}-%x.%j.out --job-name=SCN_${patient} -c ${processors} --time=1-00:00:00 --mem=${mem} --export=ALL,patient="${patient}" ${script_dir}/HCM_SCS_2021_08_SCENIC-py-ICL.sh
    # --parsable 
  done
  
fi

# ==============================================================================
# Note added later on how to install scenic

# conda create -y -n pyscenic python=3.10
# conda activate pyscenic
#
# pip install pyscenic

# See also: https://pyscenic.readthedocs.io/en/latest/installation.html

# ==============================================================================
# Actual script

# prepare script

# activate correct environment where SCENIC is installed
# Note, use "conda env list" to see all env
# First make conda accessible in this shell (https://github.com/conda/conda/issues/7980)
condapath=$(conda info --base)
source ${condapath}/etc/profile.d/conda.sh
# Activate
conda activate pyscenic

# paths
scenicdatadir=/hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3b/SCENIC/DATABASES/

datadir=/hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3b/
patient=${patient} # Give this via command line "R.P1" 
datasuffix="RID2l"

currentworkdir="${datadir}SCENIC/${patient}/"

# go to appropriate SCENIC/patient directory
mkdir -p ${currentworkdir}
cd  ${currentworkdir}

# ========= GRN (gene reaction network) ========================================
echo "GRN step started"

# pyscenic grn ${datadir}Ldata/${patient}${datasuffix}.loom ${loomdatadir}hs_hgnc_curated_tfs.txt -o adj.csv --num_workers 20
NUMBA_THREADING_LAYER='omp' pyscenic grn ${datadir}Ldata/${patient}${datasuffix}.loom ${scenicdatadir}hs_hgnc_curated_tfs.txt -o adj.csv --num_workers 20
  # Note: tbb warning rose from numba py lib, see also https://github.com/aertslab/pySCENIC/issues/182 for fix
  
# ========= CTX (includes motif stuff I guess) ==================================
echo "CTX step started"

# See also https://resources.aertslab.org/cistarget/ for a list of resources with similar databases
f_db_names="${scenicdatadir}hg19-tss-centered-10kb-7species.mc9nr.feather  ${scenicdatadir}hg19-500bp-upstream-7species.mc9nr.feather"
# motif downloaded via link @ https://rdrr.io/github/aertslab/RcisTarget/man/importAnnotations.html
f_motif_path=motifs-v9-nr.hgnc-m0.001-o0.0.tbl

NUMBA_THREADING_LAYER='omp' pyscenic ctx ${currentworkdir}adj.csv \
    ${f_db_names} \
    --annotations_fname ${scenicdatadir}${f_motif_path} \
    --expression_mtx_fname ${datadir}Ldata/${patient}${datasuffix}.loom \
    --output ${currentworkdir}reg.csv \
    --num_workers 20    
    # I don't think I want to mask dropouts
    # --mask_dropouts \
    
# ========= AUC, area under curve measurements for "expression" regulon ========
echo "AUC step started"

mkdir -p ${currentworkdir}output/    
NUMBA_THREADING_LAYER='omp' pyscenic aucell \
        ${datadir}Ldata/${patient}${datasuffix}.loom \
        ${currentworkdir}reg.csv \
        --output ${currentworkdir}output/SCENIC_out_${patient}${datasuffix}.loom \
        --num_workers 20    
    
echo "SCENIC bash done .."
    
    
    
    
    
    
