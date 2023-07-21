#!/bin/bash

#####
# Misc notes

#so this the command that found on the github discussion in your comment: 

#Just fyi, I installed pyscenic just now using

#conda create -y -n pyscenic python=3.7
#conda activate pyscenic
#conda install -y numpy
#conda install -y -c anaconda cytoolz
#pip install pyscenic

#Now aucell ran after running pip install loompy==2.0.17."

### Set up some information

#directory to the scenic folder where the three databses files are and where each patient has a sub folder and 
#scenicdatadir="/home/e.mili_cbs-niob.local/Desktop/Single_nucleus_Chaffin_2020/run/new_harmony/CM/Final_ALL_patients_res0.23_algorithm3/scenic/"
scenicdatadir="/hpc/hub_oudenaarden/mwehrens/data/2023_SCENIC_Eloi/"

f_db_names="${scenicdatadir}hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather  ${scenicdatadir}hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather"

# This code is to be manually executed in the command line of the HPC
# It will just call this very script; the code below is then *not*
# executed because of the if statement below; instead, the applicable
# code below is executed to perform the SCENIC code for each respective patient.
if [[ 0 -eq 1 ]]; then
  
  # Execute this code manually on the cluster to submit the job:
  patients=("P1622" "P1422" "P1722" "P1462" "P1558" "P1540" "P1718" "P1602" "P1430" "P1630" "P1516" "P1515" "P1603" "P1447" "P1678" "P1547" "P1610" "P1510" "P1371" "P1479" "P1300" "P1617" "P1290" "P1600" "P1437" "P1685" "P1582" "P1702" "P1561" "P1358" "P1549" "P1735" "P1304" "P1707" "P1425" "P1508" "P1631" "P1539" "P1726" "P1504" "P1472" "P1606")
  script_dir=/hpc/hub_oudenaarden/mwehrens/data/2023_SCENIC_Eloi
  processors=20
  mem=15G # note: use the "seff <jobid>" command to check memory usage
  for i in {0..41}; do
    
    patient=${patients[$i]}
    echo "${patient}"

    echo "AUC step started"
    sbatch --output=slurm-scenic-${patient}-%x.%j.out --job-name=SCN_${patient} -c ${processors} --time=1-00:00:00 --mem=${mem} --export=ALL,patient="${patient}" ${script_dir}/last_step_scenicanalysis_bashscripte_v2-slurm.sh
    # --parsable 
  done
  
  # just a little line to execute after jobs are done to check all files are there
  for i in {0..41}; do
      patient=${patients[$i]}
      ls -lhat $patient
  done
  
  # collect output to one directory
  mkdir collected_output
  for i in {0..41}; do
      patient=${patients[$i]}
      mv ${scenicdatadir}/${patient}/SCENIC_out_${patient}.loom ${scenicdatadir}/collected_output/SCENIC_out_${patient}.loom      
  done
  
fi

# ==============================================================================
# prepare script

# activate correct environment where SCENIC is installed
# Note, use "conda env list" to see all env
# First make conda accessible in this shell (https://github.com/conda/conda/issues/7980)
condapath=$(conda info --base)
source ${condapath}/etc/profile.d/conda.sh
# Activate
conda activate pyscenic

# ==============================================================================
# Start of the actual script
mkdir -p ${currentworkdir}output/    

echo "AUC step started"   
NUMBA_THREADING_LAYER='omp' pyscenic aucell \
     ${scenicdatadir}${patient}/${patient}.loom \
     ${scenicdatadir}/${patient}/reg.csv \
     --output ${scenicdatadir}/${patient}/SCENIC_out_${patient}.loom \
     --num_workers 20    

echo "SCENIC bash done .."


