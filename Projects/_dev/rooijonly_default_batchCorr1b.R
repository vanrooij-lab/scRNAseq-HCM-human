


# Load analysis object (again)
current_analysis=list()
DEFAULT_SET="ROOIJonly_default"
current_analysis[[DEFAULT_SET]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DEFAULT_SET,'.h5seurat'))
RIDL_SET = 'ROOIJonly_RID2l'
current_analysis[[RIDL_SET]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',RIDL_SET,'.h5seurat'))

# New batch correction, C1b, now including all genes

currentSeuratObject_list2 <- SplitObject(current_analysis[[DEFAULT_SET]], split.by = "annotation_sample_str")
currentSeuratObject_list2 <- lapply(X = currentSeuratObject_list2, FUN = SCTransform); beepr::beep()

all_genenames = Reduce(intersect, lapply(currentSeuratObject_list2, rownames))

#currentFeatures2 <- SelectIntegrationFeatures(object.list = currentSeuratObject_list2, nfeatures = 3000)
currentFeatures2 <- all_genenames

currentSeuratObject_list2 <- PrepSCTIntegration(object.list = currentSeuratObject_list2, anchor.features = currentFeatures2); beepr::beep()
currentAnchors2 <- FindIntegrationAnchors(object.list = currentSeuratObject_list2, 
    normalization.method = 'SCT', anchor.features = currentFeatures2, dims=1:30); beepr::beep() 
    # dims=1:10, k.anchor = ..
currentSeuratObject_recombined2 <- IntegrateData(anchorset = currentAnchors2, normalization.method = 'SCT', k.weight=10); beepr::beep()
    #currentSeuratObject_recombined2 <- IntegrateData(anchorset = currentAnchors2, normalization.method = 'SCT'); beepr::beep()

# Run analysis
currentSeuratObject_recombined2 = mySeuratAnalysis_verybasic_part2only(mySeuratObject = currentSeuratObject_recombined2)

shorthand_someplots(currentSeuratObject_recombined2)
DimPlot(currentSeuratObject_recombined2, reduction = "umap", group.by = "ident", label = F) 
DimPlot(currentSeuratObject_recombined2, reduction = "umap", group.by = "ident", label = T) 

mySeuratCommonPlots(mySeuratObject = currentSeuratObject_recombined2, run_name = 'ROOIJonly_default_batch.C1b')
FeaturePlot(currentSeuratObject_recombined2, features = c('MYH7'), cols = rainbow_colors)
FeaturePlot(currentSeuratObject_recombined2, features = c('TTN'), cols = rainbow_colors)







