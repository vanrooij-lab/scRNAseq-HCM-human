#!/bin/bash

# Lisa

################################################################################
# To start these jobs

# To start these jobs
if [[ "" == "never" ]]; then
  cd /hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3/LISA
  for cl_idx in {1..5}; do
    Down=''
    sbatch --job-name=LISA_${cl_idx}${Down} -c 1 --time=1-00:00:00 --mem=32G --export=ALL,cl_idx="${cl_idx}" --output=slurm-%x.%j.out /hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/HCM_SCS_2021_08_LISA.sh
    Down='Down'
    sbatch --job-name=LISA_${cl_idx}${Down} -c 1 --time=1-00:00:00 --mem=32G --export=ALL,cl_idx="${cl_idx}",Down="${Down}" --output=slurm-%x.%j.out /hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/HCM_SCS_2021_08_LISA.sh
  done
fi


################################################################################

################################################################################
# Installation of lisa

if [[ "" == "install" ]]; then
  conda create -n LISA
  conda install -c liulab-dfci lisa2
  
  curl http://cistrome.org/~alynch/data/lisa_data/hg38_1000_2.0.h5 -o hg38_1000_2.0.h5
  lisa install hg38 oneshot ./hg38_1000_2.0.h5 --remove
  
fi

################################################################################
# Running for Clusters

cd /hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3/LISA

# activate correct environment
# Note, use "conda env list" to see all env
# First make conda accessible in this shell (https://github.com/conda/conda/issues/7980)
condapath=$(conda info --base)
source ${condapath}/etc/profile.d/conda.sh
# Activate
conda activate LISA

# multi doesn't appear to allow setting a background
# lisa multi hg38 clusters/ClusterHits_ROOIJonly_RID2l_clExtended_table_cl*.txt --rp_map enhanced_10K -o clusters/ --save_metadata

echo "Starting LISA .."
echo "Start @ $(date)"

lisa oneshot hg38 clusters/ClusterHits${Down}_ROOIJonly_RID2l_clExtended_table_cl${cl_idx}.txt \
    --rp_map enhanced_10K -o clusters/oneshot_${Down}cl.${cl_idx}_out --save_metadata --background_strategy provided \
    --background_list clusters/ROOIJonly_RID2l_HOMER_backgroundGenes.txt
    #--rp_map enhanced_10K -o clusters/strictbg_oneshot_${Down}cl.${cl_idx}_out --save_metadata --background_strategy provided \
    #--background_list clusters/ROOIJonly_RID2l_clExtended_background_table_symbol_strict_0.2.txt
    
echo "LISA cluster analysis done .."
echo "End @ $(date)"

echo "END"


# ClusterHits_ROOIJonly_RID2l_clExtended_table_cl1.txt	GATA4, KLF10, MEF2B
# ClusterHits_ROOIJonly_RID2l_clExtended_table_cl2.txt	JUN, GATA4, CEBPB, AR, MEF2B
# ClusterHits_ROOIJonly_RID2l_clExtended_table_cl3.txt	MYOD1, KLF10, GATA4, MED1, TEAD1
# ClusterHits_ROOIJonly_RID2l_clExtended_table_cl4.txt	TP53, CDK9, FOXA1, TAF1, SUPT5H, CEBPB, MED1
# ClusterHits_ROOIJonly_RID2l_clExtended_table_cl5.txt	GATA2, TAL1, JUN, RELA, AR, FOS

