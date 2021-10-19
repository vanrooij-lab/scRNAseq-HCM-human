
# Apply sevaral batch corrections to separated datasets to see what happens

################################################################################

# Load analysis object (again)
current_analysis=list()
DEFAULT_SET="ROOIJonly_default"
current_analysis[[DEFAULT_SET]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DEFAULT_SET,'.h5seurat'))
RIDL_SET = 'ROOIJonly_RID2l'
current_analysis[[RIDL_SET]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',RIDL_SET,'.h5seurat'))

################################################################################
# Run combat

library(sva)

#perform the batch correction
groups = current_analysis[[DEFAULT_SET]]$annotation_patient_str # sapply(as.character(conditions), switch, "UHR" = 1, "HBR" = 2, USE.NAMES = F)
batches = current_analysis[[DEFAULT_SET]]$annotation_sample_str # sapply(as.character(library_methods), switch, "Ribo" = 1, "Poly" = 2, USE.NAMES = F)
corrected_data = ComBat_seq(counts = as.matrix(current_analysis[[DEFAULT_SET]]@assays$RNA@counts), batch = batches, group = groups)
# confounding error

################################################################################
# Run "default" batch correction

# # split the dataset into a list of two seurat objects
currentSeuratObject.b1.list = SplitObject(current_analysis[[DEFAULT_SET]], split.by = "annotation_sample_str")

# normalize and identify variable features for each dataset independently
currentSeuratObject.b1.list <- lapply(X = currentSeuratObject.b1.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
b1.list.features <- SelectIntegrationFeatures(object.list = currentSeuratObject.b1.list)

# Now integrate 
b1.object.anchors <- FindIntegrationAnchors(object.list = currentSeuratObject.b1.list, anchor.features = b1.list.features); beepr::beep()

# this command creates an 'integrated' data assay
currentSeuratObject.b1.combined <- IntegrateData(anchorset = b1.object.anchors, k.weight = 20)

# Run analysis
currentSeuratObject.b1.combined = mySeuratAnalysis_verybasic_part2only(mySeuratObject = currentSeuratObject.b1.combined)

################################################################################
# Run Seurat batch correction, more conservative version

# Perform batch-corretion

# In case there's too few cells in a sample
# WRL_object_sel <- subset(current_analysis[[DEFAULT_SET]], subset = annotation_batch_n>50)

currentSeuratObject_list <- SplitObject(current_analysis[[DEFAULT_SET]], split.by = "annotation_sample_str")
currentSeuratObject_list <- lapply(X = currentSeuratObject_list, FUN = SCTransform); beepr::beep()
currentFeatures <- SelectIntegrationFeatures(object.list = currentSeuratObject_list, nfeatures = 3000)
currentSeuratObject_list <- PrepSCTIntegration(object.list = currentSeuratObject_list, anchor.features = currentFeatures); beepr::beep()
currentAnchors <- FindIntegrationAnchors(object.list = currentSeuratObject_list, 
    normalization.method = 'SCT', anchor.features = currentFeatures, dims=1:30); beepr::beep() 
    # dims=1:10, k.anchor = ..
currentSeuratObject_recombined <- IntegrateData(anchorset = currentAnchors, normalization.method = 'SCT', k.weight=10); beepr::beep()
    #currentSeuratObject_recombined <- IntegrateData(anchorset = currentAnchors, normalization.method = 'SCT'); beepr::beep()
# Run analysis
currentSeuratObject_recombined = mySeuratAnalysis_verybasic_part2only(mySeuratObject = currentSeuratObject_recombined)

shorthand_someplots(currentSeuratObject_recombined)
DimPlot(currentSeuratObject_recombined, reduction = "umap", group.by = "ident", label = F) 
DimPlot(currentSeuratObject_recombined, reduction = "umap", group.by = "ident", label = T) 

mySeuratCommonPlots(mySeuratObject = currentSeuratObject_recombined, run_name = 'ROOIJonly_default_batch.C1')
FeaturePlot(currentSeuratObject_recombined, features = c('MYH7'), cols = rainbow_colors)
FeaturePlot(currentSeuratObject_recombined, features = c('TTN'), cols = rainbow_colors)

rownames_integrated = rownames(currentSeuratObject_recombined@assays$integrated@scale.data)
rownames_integrated[grepl('TTN',rownames_integrated)]
rownames_integrated[grepl('MALAT1',rownames_integrated)]
rownames_integrated[grepl('NPPA',rownames_integrated)]
rownames_integrated[grepl('MYH7',rownames_integrated)]
rownames_integrated[grepl('MYH6',rownames_integrated)]

dim(currentSeuratObject_recombined@assays$RNA@counts)
dim(currentSeuratObject_recombined@assays$RNA@data)
rownames(currentSeuratObject_recombined@assays$RNA@data)[grepl('MYH7',rownames(currentSeuratObject_recombined@assays$RNA@data))]
rownames(currentSeuratObject_recombined@assays$RNA@data)[grepl('TTN',rownames(currentSeuratObject_recombined@assays$RNA@data))]

currentSeuratObject_recombined@assays$RNA@data['MYH7',]


# Also perform cluster enrichment analysis
DE_clusters=list()
DE_clusters[['ROOIJonly_default_batch.C1']] = diff_express_clusters(currentSeuratObject_recombined, mc.cores = 1)
diff_express_clusters_save_results(all_markers = DE_clusters[['ROOIJonly_default_batch.C1']], 
    run_name = 'ROOIJonly_default_batch.C1', base_dir = base_dir, topX = 30)

################################################################################
# Some custom plots on those analyses

shorthand_someplots = function(seuratobject) {
    p=DimPlot(seuratobject, reduction = "umap", group.by = "annotation_sample_str")
    print(p)
    p=DimPlot(seuratobject, reduction = "umap", group.by = "annotation_patient_str")
    print(p)
    
    # ident from default analysis
    p=DimPlot(seuratobject, reduction = "umap", group.by = "ident") 
    print(p)
    # now show ident from RID2l analysis
    #current_analysis[[DEFAULT_SET]]$Idents_RIDL = Idents(current_analysis[[RIDL_SET]])
    seuratobject$Idents_RIDL = Idents(current_analysis[[RIDL_SET]])
    p=DimPlot(seuratobject, reduction = "umap", group.by = "Idents_RIDL") 
    print(p)
    
    # Some features of interest
    p=FeaturePlot(seuratobject, features = c("TTN","NPPA","NPPB",'KCNQ1OT1'), cols = rainbow_colors)
    print(p)
}

shorthand_someplots(currentSeuratObject_recombined)
shorthand_someplots(currentSeuratObject.b1.combined.test)

# Some extra plots for 2nd correction
DimPlot(currentSeuratObject_recombined, reduction = "umap", group.by = "ident", label = T) 

# Some extra for 1st correction
FeaturePlot(currentSeuratObject.b1.combined.test, features = c("MALAT1"), cols = rainbow_colors)





