

# base_dir = /Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis

# Load Seurat datasets of interest
for (seuratObject_name in c('ROOIJonly_RID2l','ROOIJonly_default_Int1c')) {
    current_analysis[[seuratObject_name]] =
            LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',seuratObject_name,'.h5seurat'))
}
    
# Loading grouped SCS
load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/Previous_analysis_for_reference/export_groupedSCS_patientAllMod_only.Rdata')
# Load required RaceID2 class
source(paste0('/Users/m.wehrens/Documents/git_repos/SCS_Joep/',"/Functions/RaceID2_StemID_class.R"))

# get cell names joep
cell_names_Joep = colnames(hcm_scs@fdata)
cluster_assignments_Joep = hcm_scs@cluster$kpart

# General well plate layout
wellplate=data.frame(colnr = rep(1:24, times=16), 
                     rowsym=rep(c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'), each=24),
                     wellcount=1:384)
# now convert
cell_names_Joep_inSeuratstyle =
    sapply(cell_names_Joep, function(x) {
        s=strsplit(x,split = '_')[[1]]; 
        paste0(switch(s[1],JE1='JE01',JE2='JE02',JE3='JE03',JE4='JE04'),'_',
              sprintf("X%03d", wellplate[wellplate$colnr==s[3]&wellplate$rowsym==s[2],]$wellcount))
        })

names(cluster_assignments_Joep) = cell_names_Joep_inSeuratstyle

# Now add 
current_analysis$ROOIJonly_default_Int1c[['RaceID2_original_clusters']] = 
    as.factor(cluster_assignments_Joep[colnames(current_analysis$ROOIJonly_default_Int1c)])
current_analysis$ROOIJonly_default_Int1c[['RaceID2_included']] = 
    colnames(current_analysis$ROOIJonly_default_Int1c) %in% cell_names_Joep_inSeuratstyle

# Make a little plot
p=DimPlot(current_analysis$ROOIJonly_default_Int1c, pt.size = 3, group.by = 'RaceID2_original_clusters', cols= col_vector_60)#, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/ROOIJonly_default_Int1c_custom_comparisonpreviousClustering_all.png'), plot = p, height=5, width=5)
p
p=DimPlot(subset(current_analysis$ROOIJonly_default_Int1c, RaceID2_included==T), pt.size = 3, group.by = 'RaceID2_original_clusters', cols= col_vector_60)#, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/ROOIJonly_default_Int1c_custom_comparisonpreviousClustering.png'), plot = p, height=5, width=5)
p

# Now do same for RID2l analysis
current_analysis$ROOIJonly_RID2l[['RaceID2_original_clusters']] = 
    as.factor(cluster_assignments_Joep[colnames(current_analysis$ROOIJonly_RID2l)])
current_analysis$ROOIJonly_RID2l[['RaceID2_included']] = 
    colnames(current_analysis$ROOIJonly_RID2l) %in% cell_names_Joep_inSeuratstyle

p=DimPlot(current_analysis$ROOIJonly_RID2l, pt.size = 2, group.by = 'RaceID2_original_clusters', cols= col_vector_60)#, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/ROOIJonly_RID2l_custom_comparisonpreviousClustering_all.png'), plot = p, height=5, width=5)
p
p=DimPlot(subset(current_analysis$ROOIJonly_RID2l, RaceID2_included==T), pt.size = 2, group.by = 'RaceID2_original_clusters', cols= col_vector_60)#, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/ROOIJonly_RID2l_custom_comparisonpreviousClustering.png'), plot = p, height=5, width=5)
p

# ggsave(filename = paste0(base_dir,'Rplots/',XXXX,'.png'), plot = p_cm, height=5, width=5)
        
    


