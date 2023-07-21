#so this the command that found on the github discussion in your comment: 

"Just fyi, I installed pyscenic just now using

conda create -y -n pyscenic python=3.7
conda activate pyscenic
conda install -y numpy
conda install -y -c anaconda cytoolz
pip install pyscenic

Now aucell ran after running pip install loompy==2.0.17."


### this is the start of the script
patients=("P1622" "P1422" "P1722" "P1462" "P1558" "P1540" "P1718" "P1602" "P1430" "P1630" "P1516" "P1515" "P1603" "P1447" "P1678" "P1547" "P1610" "P1510" "P1371" "P1479" "P1300" "P1617" "P1290" "P1600" "P1437" "P1685" "P1582" "P1702" "P1561" "P1358" "P1549" "P1735" "P1304" "P1707" "P1425" "P1508" "P1631" "P1539" "P1726" "P1504" "P1472" "P1606")

#directory to the scenic folder where the three databses files are and where each patient has a sub folder and 
#scenicdatadir="/home/e.mili_cbs-niob.local/Desktop/Single_nucleus_Chaffin_2020/run/new_harmony/CM/Final_ALL_patients_res0.23_algorithm3/scenic/"
scenicdatadir="/hpc/hub_oudenaarden/mwehrens/data/2023_SCENIC_Eloi/"

f_db_names="${scenicdatadir}hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather  ${scenicdatadir}hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather"

# Loop over each patient
for i in {1..42}; do
    patient_id=${patients[@]}

    echo "AUC step started"   
    NUMBA_THREADING_LAYER='omp' pyscenic aucell \
         ${scenicdatadir}${patient_id}/${patient_id}.loom \
         ${scenicdatadir}/${patient_id}reg.csv \
         --output ${scenicdatadir}/${patient_id}SCENIC_out_${patient_id}.loom \
         --num_workers 20    
    
    echo "SCENIC bash done .."

done

