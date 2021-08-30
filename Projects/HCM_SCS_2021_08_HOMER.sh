#!/bin/bash

# analysis of the data in the following dir:
dir="/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/Homer"

for i in {1..6}
do
   { findMotifs.pl $dir/ROOIJonly_RID2l_HOMER_core_regulons_s.R.$i.txt human $dir/output_s.R.$i -bg $dir/ROOIJonly_RID2l_HOMER_backgroundGenes.txt ;
    echo "Regulon $i done"} &
done

wait

echo "All Homer regulon analyses done .."

################################################################################

dir="/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/Homer"
cd $dir

# Now the clusters

echo "Start @ $(date)"

for i in {1..5}
do
   { findMotifs.pl $dir/ClusterHits_ROOIJonly_RID2l_clExtended_table_cl${i}.txt human $dir/output_cl.$i -bg $dir/ROOIJonly_RID2l_HOMER_backgroundGenes.txt ;
    echo "Regulon $i done"} &
done

wait

echo "All Homer cluster analyses done .."
echo "End @ $(date)"

################################################################################
