
# This part of the script follows the tutorial
# at https://satijalab.org/seurat/articles/integration_introduction.html

################################################
# My data is in:
# Teichmann_Sampled



library(Seurat)
library(patchwork)

# MW
library(RColorBrewer)

################################################

# Let's first perform a standard Seurat analysis on the data as is
# Run the standard workflow for visualization and clustering
# See:
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# Before starting, let's deduce some annotation
annotation_source=unlist(lapply(  strsplit(colnames(Teichmann_Sampled), split='_'), 
    function(splitted_string) {splitted_string[2]}))
annotation_sample=unlist(lapply(  strsplit(colnames(Teichmann_Sampled), split='-'), 
    function(splitted_string) {splitted_string[3]}))

# Initialize object
Teichmann_object = CreateSeuratObject(Teichmann_Sampled, project = 'Teichmann_s')

# Add annotation
Teichmann_object[["source_tissue"]]     <- as.factor(annotation_source)
Teichmann_object[["annotation_sample"]] <- as.factor(annotation_sample)

VlnPlot(Teichmann_object, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


# Remove some cells
Teichmann_object <- subset(Teichmann_object, subset = nFeature_RNA > 200 & nCount_RNA > 3000)
# Normalize data
Teichmann_object <- NormalizeData(Teichmann_object, normalization.method = "LogNormalize", scale.factor = 10000) # values are defaults
# Find variable features
Teichmann_object <- FindVariableFeatures(Teichmann_object, selection.method = "vst", nfeatures = 2000)

# ### some plots var features
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Teichmann_object), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Teichmann_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

all.genes <- rownames(Teichmann_object)
Teichmann_object <- ScaleData(Teichmann_object, features = all.genes)

Teichmann_object <- RunPCA(Teichmann_object, npcs = 30, verbose = FALSE)
Teichmann_object <- RunUMAP(Teichmann_object, reduction = "pca", dims = 1:30)
Teichmann_object <- FindNeighbors(Teichmann_object, reduction = "pca", dims = 1:30)
Teichmann_object <- FindClusters(Teichmann_object, resolution = 0.5)


# Analysis done, show some plots

print(Teichmann_object[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Teichmann_object, dims = 1:5, reduction = "pca", ncol = 5)

DimPlot(Teichmann_object, reduction = "pca")

JackStrawPlot(Teichmann_object, dims = 1:5) # doesn't work for some reason?
ElbowPlot(Teichmann_object)

DimPlot(Teichmann_object, reduction = "umap")

FeaturePlot(Teichmann_object, features = c("MALAT1"))
FeaturePlot(Teichmann_object, features = c("TTN"))

FeaturePlot(Teichmann_object, features = c('MSRB2', 'MRPS21', 'DNAJA4', 'COX6B1', 'MTRNR2L1', 'TRDN'))

DimPlot(Teichmann_object, group.by = 'source_tissue')

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector_60 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
DimPlot(Teichmann_object, group.by = 'annotation_sample', cols = col_vector_60,
            label = T, repel = T, label.size = 3)+theme(legend.position = 'none')



################################################



# For now, we don't have different samples, so we can just run it on the one thing
Teichmann_object <- NormalizeData(Teichmann_object)
Teichmann_object <- FindVariableFeatures(Teichmann_object, selection.method = "vst", nfeatures = 2000)
Teichmann_anchors = FindIntegrationAnchors(Teichman_object)






