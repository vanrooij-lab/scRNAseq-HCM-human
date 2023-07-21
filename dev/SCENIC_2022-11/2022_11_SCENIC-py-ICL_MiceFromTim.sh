#!/bin/bash

# ==============================================================================
# Note that loom files are generated in the R files that serve as input
# for this file.

# See R file with pre-processing for comments with more info.

# ==============================================================================
# This code is not run
# 
# It's only here as a reminder how to submit the jobs.

# This isn't executed
if [[ 0 -eq 1 ]]; then
  
  # Execute this code manually on the cluster to submit the job:
  script_dir=/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis
  processors=20
  mem=50G # 50G

  loomfilename=dataTim_exprMat

  # sbatch this very file
  sbatch --output=slurm-scenic-${loomfilename}-%x.%j.out --job-name=SCN_TIM -c ${processors} --time=3-00:00:00 --mem=${mem} --export=ALL,loomfilename="${loomfilename}" ${script_dir}/2022_Tim/2022_11_SCENIC-py-ICL_MiceFromTim.sh
    
fi

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
scenicdatadir=/hpc/hub_oudenaarden/mwehrens/data/2022_Tim/DATABASES_mice/

datadir=/hpc/hub_oudenaarden/mwehrens/data/2022_Tim/

currentworkdir="${datadir}SCENIC/${loomfilename}/"

# go to appropriate SCENIC/patient directory
mkdir -p ${currentworkdir}
cd  ${currentworkdir}

# ========= GRN (gene reaction network) ========================================
echo "GRN step started"

tflist=allTFs_mm.txt

# pyscenic grn ${datadir}Ldata/${patient}${datasuffix}.loom ${loomdatadir}hs_hgnc_curated_tfs.txt -o adj.csv --num_workers 20
NUMBA_THREADING_LAYER='omp' pyscenic grn ${datadir}Ldata/${loomfilename}.loom ${scenicdatadir}${tflist} -o adj.csv --num_workers 20
  # Note: tbb warning rose from numba py lib, see also https://github.com/aertslab/pySCENIC/issues/182 for fix
  
# ========= CTX (includes motif stuff I guess) ==================================
echo "CTX step started"

# See also https://resources.aertslab.org/cistarget/ for a list of resources with similar databases
# Mouse resource downloaded from: https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/
f_db_names="${scenicdatadir}mm9-tss-centered-10kb-7species.mc9nr.feather ${scenicdatadir}mm9-500bp-upstream-7species.mc9nr.feather"
# motif downloaded via link @ https://rdrr.io/github/aertslab/RcisTarget/man/importAnnotations.html
f_motif_path=motifs-v9-nr.mgi-m0.001-o0.0.tbl

NUMBA_THREADING_LAYER='omp' pyscenic ctx ${currentworkdir}adj.csv \
    ${f_db_names} \
    --annotations_fname ${scenicdatadir}${f_motif_path} \
    --expression_mtx_fname ${datadir}Ldata/${loomfilename}.loom \
    --output ${currentworkdir}reg.csv \
    --num_workers 20    
    # I don't think I want to mask dropouts
    # --mask_dropouts \
    
# ========= AUC, area under curve measurements for "expression" regulon ========
echo "AUC step started"

mkdir -p ${currentworkdir}output/    
NUMBA_THREADING_LAYER='omp' pyscenic aucell \
        ${datadir}Ldata/${loomfilename}.loom \
        ${currentworkdir}reg.csv \
        --output ${currentworkdir}output/SCENIC_out_${loomfilename}.loom \
        --num_workers 20    
    
echo "SCENIC bash done .."
    
    
    
    
    
    
