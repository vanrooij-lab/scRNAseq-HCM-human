
# New batch correction, C1b, now including all genes, but correcting only for the patients

################################################################################

# Load analysis object (again)
current_analysis=list()
DEFAULT_SET="ROOIJonly_default"
current_analysis[[DEFAULT_SET]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DEFAULT_SET,'.h5seurat'))
RIDL_SET = 'ROOIJonly_RID2l'
current_analysis[[RIDL_SET]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',RIDL_SET,'.h5seurat'))

################################################################################

currentSeuratObject_list3 <- SplitObject(current_analysis[[DEFAULT_SET]], split.by = "annotation_patient_str")
currentSeuratObject_list3 <- lapply(X = currentSeuratObject_list3, FUN = SCTransform); beepr::beep()

all_genenames = Reduce(intersect, lapply(currentSeuratObject_list3, rownames))

#currentFeatures3 <- SelectIntegrationFeatures(object.list = currentSeuratObject_list3, nfeatures = 3000)
currentFeatures3 <- all_genenames

currentSeuratObject_list3 <- PrepSCTIntegration(object.list = currentSeuratObject_list3, anchor.features = currentFeatures3); beepr::beep()
currentAnchors3 <- FindIntegrationAnchors(object.list = currentSeuratObject_list3, 
    normalization.method = 'SCT', anchor.features = currentFeatures3, dims=1:30); beepr::beep() 
    # dims=1:10, k.anchor = ..

currentSeuratObject_recombined3 <- IntegrateData(anchorset = currentAnchors3, normalization.method = 'SCT'); beepr::beep()
#currentSeuratObject_recombined3 <- IntegrateData(anchorset = currentAnchors3, normalization.method = 'SCT', k.weight=10); beepr::beep()
 
# SaveH5Seurat(object = currentSeuratObject_recombined3, overwrite = T,
#         filename = paste0(base_dir,'Rdata/',DEFAULT_SET,'_Integrated1c.h5seurat'))
SaveH5Seurat(object = currentSeuratObject_recombined3,
            filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',paste0(DEFAULT_SET, '_Int1c'),'.h5seurat'))
# currentSeuratObject_recombined3 = LoadH5Seurat(file = paste0(base_dir,'Rdata/',DEFAULT_SET,'_Integrated1c.h5seurat'))
    

# Run analysis
currentSeuratObject_recombined3 = mySeuratAnalysis_verybasic_part2only(mySeuratObject = currentSeuratObject_recombined3)

DimPlot(currentSeuratObject_recombined3, reduction = "umap", group.by = "ident", label = F) 
DimPlot(currentSeuratObject_recombined3, reduction = "umap", group.by = "ident", label = T) 

mySeuratCommonPlots(mySeuratObject = currentSeuratObject_recombined3, run_name = 'ROOIJonly_default_batch.C1c')

FeaturePlot(currentSeuratObject_recombined3, features = c('MYH7'), cols = rainbow_colors)
FeaturePlot(currentSeuratObject_recombined3, features = c('TTN'), cols = rainbow_colors)

# Now show the RIDL clustering projected on the current umap
currentSeuratObject_recombined3[['ridl_clusters']] = Idents(current_analysis[[RIDL_SET]])
DimPlot(currentSeuratObject_recombined3, reduction = "umap", group.by = "ridl_clusters", label = T) 

################################################################################
# Now also run cluster DE

MYMCCORES=1

DE_clusters=list()
DE_clusters[['ROOIJonly_default_batch.C1c']] = diff_express_clusters(currentSeuratObject_recombined3, mc.cores = MYMCCORES)
diff_express_clusters_save_results(all_markers = DE_clusters[['ROOIJonly_default_batch.C1c']], 
    run_name = 'ROOIJonly_default_batch.C1c', base_dir = base_dir, topX = 30)
beepr::beep()





