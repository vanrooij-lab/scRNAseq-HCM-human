
# Let's check whether the matrices are indeed perfectly identical

DATASET_NAME='ROOIJonly.sp.bt_RID2l_clExtended'
if (!exists('current_analysis')) {current_analysis = list()}
current_analysis[[DATASET_NAME]] =
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))
    # current_analysis$ROOIJonly.sp.bt_RID2l_clExtended

DATASET_NAME_old='ROOIJonly_RID2l_clExtended'
if (!exists('current_analysis')) {current_analysis = list()}
current_analysis[[paste0(DATASET_NAME_old,'_3a')]] =
    LoadH5Seurat(file = paste0("/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3a/",'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME_old,'.h5seurat'))
    
# Length is of course the same
length(rownames(current_analysis$ROOIJonly_RID2l_3a))
length(rownames(current_analysis$ROOIJonly.sp.bt_RID2l_clExtended))

# Names are not entirely though
all(rownames(current_analysis$ROOIJonly_RID2l_3a) == rownames(current_analysis$ROOIJonly.sp.bt_RID2l_clExtended))
n_old = rownames(current_analysis$ROOIJonly_RID2l_3a)
n_new = rownames(current_analysis$ROOIJonly.sp.bt_RID2l_clExtended)
n_new[n_old!=n_new]

# What about the matrices?
dim(current_analysis$ROOIJonly_RID2l_3a@assays$RNA@counts) == dim(current_analysis$ROOIJonly.sp.bt_RID2l_clExtended@assays$RNA@counts)
all(current_analysis$ROOIJonly_RID2l_3a@assays$RNA@counts == current_analysis$ROOIJonly.sp.bt_RID2l_clExtended@assays$RNA@counts)
littletest=current_analysis$ROOIJonly_RID2l_3a@assays$RNA@counts == current_analysis$ROOIJonly.sp.bt_RID2l_clExtended@assays$RNA@counts
    # dim(littletest)
all(littletest) # TRUE 
rm('littletest')
    # Matrices are indeed the same

# So if we want to use 3a analysis output as input for something else, we better check
# whether the gene names are indeed identical
#
# For Lisa/Homer input of cluster enrichments, this is indeed the case
for (cl_idx in 1:6) {

    tcla=read.table(paste0('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3a/GeneLists/ClusterHits_ROOIJonly_RID2l_clExtended_table_cl',cl_idx,'.txt'))
    tcla
    tclb=read.table(paste0('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3b/GeneLists/ClusterHits_ROOIJonly.sp.bt_RID2l_clExtended_table_cl',cl_idx,'.txt'))
    tclb
    
    print(all(tc1a==tc1b))

}

# Note that TRIAGE anyways works on proteins, so also can use old analysis output
# ---> Re-ran this script anyways
#
# For custom regulons/modules, check whether gene lists are the same
#
# Turns out that NIBAN1 is now known as FAM129A (seems OK)
#
for (r_idx in 1:5) {

    # r_idx=4
    
    trega=read.table(paste0('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3a/Homer/REGULONS/ROOIJonly_RID2l_HOMER_core_regulons_s.R.',r_idx,'.txt'))
    trega_v=trega[,1]
    tregb=read.table(paste0('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3b/Homer/REGULONS/ROOIJonly.sp.bt_RID2l_HOMER_core_regulons_s.R.',r_idx,'.txt'))
    tregb_v=tregb[,1]
    
    print(paste0(r_idx,': ',all(trega==tregb)))

    # trega_v[!(trega_v %in% tregb_v)]
    # tregb_v[!(tregb_v %in% trega_v)]
    
}

# Note that backgrounds do change! --> So re-run is required
#
bg_3a = read.table('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3a/Homer/REGULONS/ROOIJonly_RID2l_HOMER_backgroundGenes_0.05.txt')
bg_3b = read.table('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3b/Homer/REGULONS/ROOIJonly.sp.bt_RID2l_HOMER_backgroundGenes_0.05.txt')
bg_3a_v = bg_3a[,1]
bg_3b_v = bg_3b[,1]
bg_3a_v[!(bg_3a_v %in% bg_3b_v)]
bg_3b_v[!(bg_3b_v %in% bg_3a_v)]

bg_changes_df = data.frame(bg_3a_v, bg_3b_v)
View(bg_changes_df[bg_changes_df$bg_3a_v!=bg_changes_df$bg_3b_v,])
dim(bg_changes_df[bg_changes_df$bg_3a_v!=bg_changes_df$bg_3b_v,])
dim(bg_changes_df)
145/6504



