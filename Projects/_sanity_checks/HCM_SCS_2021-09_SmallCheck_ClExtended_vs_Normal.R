

# load the files

if (!exists('current_analysis')) {current_analysis = list()}

for (DATASET_NAME in c('ROOIJonly_RID2l_clExtended','ROOIJonly_default')) {
    # DATASET_NAME = 'ROOIJonly_RID2l_clExtended'
    current_analysis[[DATASET_NAME]] =
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))
}

p1=DimPlot(current_analysis$ROOIJonly_default, group.by='annotation_patient_fct')+ggtitle('Default Seurat analysis')+theme_void()
p2=DimPlot(current_analysis$ROOIJonly_RID2l_clExtended, group.by='annotation_patient_fct')+ggtitle('Adjusted Seurat')+theme_void()
p1+p2
