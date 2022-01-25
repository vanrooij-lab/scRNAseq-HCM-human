
################################################################################

# Load new dataset
DATASET_NAME='ROOIJonly_RID2l'
current_analysis[['ROOIJonly_RID2l']] =
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))

# Load the Rooij dataset to check whether we still have consistency after bug fix
DATASET_NAME='ROOIJonly_RID2l'
base_dir_old = "/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3a/"
current_analysis[['ROOIJonly_RID2l_3a']] =
    LoadH5Seurat(file = paste0(base_dir_old,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))

# All gene names in the data are the same, except that the name is sometimes slightly modified in the hcgn part
r1=rownames(current_analysis[['ROOIJonly_RID2l_3a']]@assays$RNA@counts)
r2=rownames(current_analysis[['ROOIJonly_RID2l']]@assays$RNA@counts)
df_r1r2 = data.frame(r1, r2)
View(df_r1r2[df_r1r2$r1!=df_r1r2$r2,])

# Sanity check, checks out
count_comparison = 
    current_analysis[['ROOIJonly_RID2l_3a']]@assays$RNA@counts == current_analysis[['ROOIJonly_RID2l']]@assays$RNA@counts
all(count_comparison)


# CONCLUSION
# We have consistency, yay!

################################################################################

# Now @HPC check also the complete dataset

DATASET_NAME='ALL.SP'

base_dir_old = "/hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3a/"

# Load new set
current_analysis = list()
current_analysis[[DATASET_NAME]] =
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))
# Load old set
current_analysis[[paste0(DATASET_NAME,'_old')]] =
        LoadH5Seurat(file = paste0(base_dir_old,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))
    

# Check consistency between count tables
dim(current_analysis[['ALL.SP']]@assays$RNA@counts)
dim(current_analysis[['ALL.SP_old']]@assays$RNA@counts)

# So this is noticably different, as expected.



