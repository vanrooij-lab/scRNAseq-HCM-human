# Enter commands in R (or R studio, if installed)
#install.packages('BiocManager')
#BiocManager::install('multtest')
#install.packages('Seurat')
library(Seurat)
library(multtest)
library(dplyr)
library(cowplot)
library(ggplot2)
library(patchwork)

# devtools::install_github(repo = "satijalab/seurat", ref = "develop")

options(future.globals.maxSize = 4000 * 1024^2)


# Load environment post integration 2-4-2020
load("/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/Environment_Final_script_all_cells_Seurat_analysis_Wang_et_al_post_integration_2-4-2020.RData")

# Load environment post LV CM selection 2-4-2020
load("/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/Environment_Final_script_all_cells_Seurat_analysis_Wang_et_al_post_LVCM_selection_2-4-2020.RData")

# Load environment post LV CM DEG 17-4-2020
load("/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/Environment_Final_script_all_cells_Seurat_analysis_Wang_et_al_post_LVCM_DEG_17-4-2020.RData")

# Load data

metadata_wang_N <- 
  read.csv('/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/Data_Wang_et_al/editMW_GSE109816_normal_heart_cell_info.txt', sep = '\t')
count_table_wang_N <- 
  read.csv('/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/Data_Wang_et_al/GSE109816_normal_heart_umi_matrix.csv')

metadata_wang_HF = 
  read.csv('/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/Data_Wang_et_al/editMW_GSE121893_human_heart_sc_info.txt', sep = '\t')
count_table_wang_HF = 
  read.csv('/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/Data_Wang_et_al/human_heart_sc_umi.csv')

# make genes rownames
rownames(count_table_wang_N) = count_table_wang_N[,1]
count_table_wang_N = count_table_wang_N[,-1]

rownames(count_table_wang_HF) = count_table_wang_HF[,1]
count_table_wang_HF = count_table_wang_HF[,-1]

# To construct a reference, we will identify ‘anchors’ between the individual datasets.
# First we merge the two different raw count tables of Wang et al into a Seurat Object by first creating 
# two different objects (from each count table) and then merging these. 
# Second, we split the combined object into a list, with each dataset as an element.

Human_heart_N <- CreateSeuratObject(counts = count_table_wang_N)
Human_heart_HF <- CreateSeuratObject(counts = count_table_wang_HF)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Human_heart_N[["percent.mt"]] <- PercentageFeatureSet(Human_heart_N, pattern = "^MT-")
Human_heart_HF[["percent.mt"]] <- PercentageFeatureSet(Human_heart_HF, pattern = "^MT-")

any(grepl("^MT-",Human_heart_N@assays[["RNA"]]@counts@Dimnames[[1]])) 
any(grepl("^MT-",Human_heart_HF@assays[["RNA"]]@counts@Dimnames[[1]])) 
# The above argument indicates no Mt reads in HF dataset but the N data set contains the reads still.
# Check if this data is already filtered

min(Human_heart_N@meta.data[["nCount_RNA"]]) # 1
min(Human_heart_N@meta.data[["nFeature_RNA"]]) # 1
min(Human_heart_HF@meta.data[["nCount_RNA"]]) # 916
min(Human_heart_HF@meta.data[["nFeature_RNA"]]) # 472

median(Human_heart_N@meta.data[["nCount_RNA"]]) # 
median(Human_heart_N@meta.data[["nFeature_RNA"]]) #  
median(Human_heart_HF@meta.data[["nCount_RNA"]]) # 
median(Human_heart_HF@meta.data[["nFeature_RNA"]]) #

VlnPlot(Human_heart_N, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Human_heart_HF, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# This indicates that one of the datasets (HF) is already filtered, but the N data set is not filtered yet

# Calculate the maximal amount of UMI's accepted as described in Wang et al. (UMIs limited to 2 s.d. 
# from the mean of log10)
meanlog10_umi <- mean(log10(Human_heart_N@meta.data$nCount_RNA))
sdlog10_umi <- sd(log10(Human_heart_N@meta.data$nCount_RNA))
meanlog10_umi + (2*sdlog10_umi) # 6.004741
max_umi <- 10^6.004741 # 1010976

# Filter Human_heart_N with at least 500 genes per cell (from Wang et al.) 
Human_heart_N <- CreateSeuratObject(counts = count_table_wang_N, min.features = 500)
Human_heart_N[["percent.mt"]] <- PercentageFeatureSet(Human_heart_N, pattern = "^MT-")

VlnPlot(Human_heart_N, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(Human_heart_N, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Human_heart_N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Wang et al. also apply a filter step whereby cells containing <50% unique reads are removed. 
metadata_wang_HF$unique_aligned <- metadata_wang_HF$Ambiguity/metadata_wang_HF$Aligned
metadata_wang_N$unique_aligned <- metadata_wang_N$Ambiguity/metadata_wang_N$Aligned # largest is 0.5, inidicates this is already performed

plot1 <- FeatureScatter(Human_heart_HF, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Human_heart_HF, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Although there does not seem to be a correlation between the amount of total reads and % mitochondrial reads,
# filter out cells of >72% mito reads because this was also done by Wang et al in the HF data set. 
Human_heart_N <- subset(Human_heart_N, percent.mt < 72 & nCount_RNA < max_umi)

# Remove MT genes of N data set
count_table_wang_N_post_filter <- as.matrix(GetAssayData(Human_heart_N, slot = "counts")) # pull count table from object
count_table_wang_N_post_filter_MT <- count_table_wang_N_post_filter[grepl("^MT-", rownames(count_table_wang_N_post_filter)),] # Check if we can grepl the MT
str(count_table_wang_N_post_filter_MT)
count_table_wang_N_post_filter_no_MT <- count_table_wang_N_post_filter[!grepl("^MT-", rownames(count_table_wang_N_post_filter)),] # Remove MT
str(count_table_wang_N_post_filter_no_MT)

Human_heart_N <- CreateSeuratObject(counts = count_table_wang_N_post_filter_no_MT)

min(Human_heart_N@meta.data[["nCount_RNA"]]) # 764
min(Human_heart_N@meta.data[["nFeature_RNA"]]) # 478

# To later split based on individual data we first add the individual data as metadata
# Since the N dataset had to be filterd, the metadata is not of same size (rows and collumns) anymore.
# therefore this has first to be fixed by selecting only the remaining cells that are in the filtered set. 

metadata_wang_N_filtered <- NULL

for (i in 1:length(colnames(count_table_wang_N_post_filter_no_MT))) {
  x <- metadata_wang_N[grepl(paste0('^',colnames(count_table_wang_N_post_filter_no_MT)[i],'$'),metadata_wang_N$ID),] 
  metadata_wang_N_filtered <- rbind(metadata_wang_N_filtered,x)
}

# Add cell ID to metadata
CellsMeta = Human_heart_N@meta.data
CellsMeta["ID"] <- metadata_wang_N_filtered$ID
head(CellsMeta)
Human_heart_N <- AddMetaData(Human_heart_N, CellsMeta$ID, col.name = "ID")
head(Human_heart_N@meta.data)

CellsMeta = Human_heart_HF@meta.data
CellsMeta["ID"] <- metadata_wang_HF$ID
head(CellsMeta)
Human_heart_HF <- AddMetaData(Human_heart_HF, CellsMeta$ID, col.name = "ID")
head(Human_heart_HF@meta.data)

# Remove cells from CM digest that have <10.000 distinct UMI's (filter from Wang et al.), again only necessary in the N dataset.

CellsMeta = Human_heart_N@meta.data
CellsMeta["Distinct.UMIs"] <- metadata_wang_N_filtered$Distinct.UMIs
head(CellsMeta)
Human_heart_N <- AddMetaData(Human_heart_N, CellsMeta$Distinct.UMIs, col.name = "Distinct.UMIs")
head(Human_heart_N@meta.data)

less_than_10000_UMIs <- metadata_wang_N_filtered[metadata_wang_N_filtered$Distinct.UMIs<10000,]
less_than_10000_UMIs_CM <- less_than_10000_UMIs[grepl("_CM", less_than_10000_UMIs$Type),]

Human_heart_N <- subset(Human_heart_N, cells = c(as.character(less_than_10000_UMIs_CM$ID)), invert = T)
metadata_wang_N_filtered <- metadata_wang_N_filtered[(grepl("_CM", metadata_wang_N_filtered$Type) & metadata_wang_N_filtered$Distinct.UMIs>10000) | 
                                                        grepl("_NCM", metadata_wang_N_filtered$Type) ,]

# To later split based on individual data we first add the individual data as metadata

CellsMeta = Human_heart_N@meta.data
CellsMeta["indvidual"] <- metadata_wang_N_filtered$Individual
head(CellsMeta)
Human_heart_N <- AddMetaData(Human_heart_N, CellsMeta$indvidual, col.name = "individual")
head(Human_heart_N@meta.data)

CellsMeta = Human_heart_HF@meta.data
CellsMeta["indvidual"] <- metadata_wang_HF$Individual
head(CellsMeta)
Human_heart_HF <- AddMetaData(Human_heart_HF, CellsMeta$indvidual, col.name = "individual")
head(Human_heart_HF@meta.data)

# Also addd cell metadata condition

CellsMeta = Human_heart_N@meta.data
CellsMeta["condition"] <- NA
CellsMeta$condition[grepl("HF", metadata_wang_N_filtered$Type )] <- "HF"
CellsMeta$condition[grepl("N_", metadata_wang_N_filtered$Type )] <- "N"
Human_heart_N <- AddMetaData(Human_heart_N, CellsMeta$condition, col.name = "condition")
head(Human_heart_N@meta.data)

CellsMeta = Human_heart_HF@meta.data
CellsMeta["condition"] <- NA
CellsMeta$condition[grepl("HF", metadata_wang_HF$Type )] <- "HF"
CellsMeta$condition[grepl("N_", metadata_wang_HF$Type )] <- "N"
Human_heart_HF <- AddMetaData(Human_heart_HF, CellsMeta$condition, col.name = "condition")
head(Human_heart_HF@meta.data)

# Also addd cell metadata Ventricle/Atria

CellsMeta = Human_heart_N@meta.data
CellsMeta["location"] <- NA
CellsMeta$condition[grepl("LV", metadata_wang_N_filtered$Type )] <- "LV"
CellsMeta$condition[grepl("LA", metadata_wang_N_filtered$Type )] <- "LA"
Human_heart_N <- AddMetaData(Human_heart_N, CellsMeta$condition, col.name = "location")
head(Human_heart_N@meta.data)

CellsMeta = Human_heart_HF@meta.data
CellsMeta["location"] <- NA
CellsMeta$condition[grepl("LV", metadata_wang_HF$Type )] <- "LV"
CellsMeta$condition[grepl("LA", metadata_wang_HF$Type )] <- "LA"
Human_heart_HF <- AddMetaData(Human_heart_HF, CellsMeta$condition, col.name = "location")
head(Human_heart_HF@meta.data)

# Merge all cells
Human_heart_combined <- merge(Human_heart_N, y = Human_heart_HF, add.cell.ids = c("Healthy set", "HF set"))

# Check how it looks like
Human_heart_N
Human_heart_HF
Human_heart_combined

# Split data
Human_heart_combined.list <- SplitObject(Human_heart_combined, split.by = "individual")

# Use the SCTransform function to normalize each dataset
for (i in 1:length(Human_heart_combined.list)) {
  Human_heart_combined.list[[i]] <- SCTransform(Human_heart_combined.list[[i]], verbose = FALSE)
}

# Next, select features for downstream integration, and run PrepSCTIntegration, 
# which ensures that all necessary Pearson residuals have been calculated.

Human_heart_combined.features <- SelectIntegrationFeatures(object.list = Human_heart_combined.list, nfeatures = 3000)
Human_heart_combined.list <- PrepSCTIntegration(object.list = Human_heart_combined.list, anchor.features = Human_heart_combined.features, 
                                    verbose = FALSE)

# Next, identify anchors and integrate the datasets. Commands are identical to the standard workflow, 
# but make sure to set normalization.method = 'SCT':
  
Human_heart_combined.anchors <- FindIntegrationAnchors(object.list = Human_heart_combined.list, normalization.method = "SCT", 
                                             anchor.features = Human_heart_combined.features, verbose = FALSE)

Human_heart_combined.integrated <- IntegrateData(anchorset = Human_heart_combined.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

# Now proceed with downstream analysis (i.e. visualization, clustering) on the integrated dataset. 
# Commands are identical to the standard workflow, but do not run the ScaleData function after integration. 
# You can see that after integration, cells group by their biological cell type (which has been pre-annotated), 
# instead of by their underlying technology.

Human_heart_combined.integrated <- RunPCA(Human_heart_combined.integrated, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(Human_heart_combined.integrated[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Human_heart_combined.integrated, dims = 1:2, reduction = "pca")
DimPlot(Human_heart_combined.integrated, reduction = "pca")
DimHeatmap(Human_heart_combined.integrated, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Human_heart_combined.integrated, dims = 1:15, cells = 500, balanced = TRUE)

# Elbow plot indicates to include at least the first 6 PC's for further analysis
ElbowPlot(Human_heart_combined.integrated) 

# JackStrawPlot shows that the first 7 PC's are significant
Human_heart_combined.integrated <- JackStraw(Human_heart_combined.integrated, num.replicate = 100)
Human_heart_combined.integrated <- ScoreJackStraw(Human_heart_combined.integrated, dims = 1:20)
JackStrawPlot(Human_heart_combined.integrated, dims = 1:20) 

# Dimension for UMAP and clustering will be set at 7
Human_heart_combined.integrated <- RunUMAP(Human_heart_combined.integrated, dims = 1:7, min.dist = 0.2)

# Cluster cells

# FindClusters performs graph-based clustering on the neighbor graph that is constructed with the 
# FindNeighbors function call. This neighbor graph is constructed using PCA space when you specifiy reduction = "pca". 
# You shouldn't add reduction = "pca" to FindClusters.

Human_heart_combined.integrated <- FindNeighbors(Human_heart_combined.integrated, dims = 1:7)
Human_heart_combined.integrated <- FindClusters(Human_heart_combined.integrated, resolution = 0.6)

DimPlot(object = Human_heart_combined.integrated, reduction = "umap", label = TRUE, repel = TRUE) 
DimPlot(object = Human_heart_combined.integrated, reduction = "umap", group.by = "condition")
DimPlot(object = Human_heart_combined.integrated, reduction = "umap", group.by = "individual")
DimPlot(object = Human_heart_combined.integrated, reduction = "umap", group.by = "location")

outputDir <- "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/Images all cells"

ggsave(paste0(outputDir,'/All cells clusters heatmap','.tiff'),
       units='mm',width=200,height=200,dpi=300)

# Look at cluster IDs of the first 5 cells
head(Idents(Human_heart_combined.integrated), 5)

# Switch to normalized data for DEG
DefaultAssay(object = Human_heart_combined.integrated) <- "RNA"
Human_heart_combined.integrated <- NormalizeData(Human_heart_combined.integrated)

# find markers for every cluster compared to all remaining cells, report only the positive ones
Human_heart_combined.integrated.markers <- FindAllMarkers(Human_heart_combined.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Human_heart_combined.integrated.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top10 <- Human_heart_combined.integrated.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DefaultAssay(object = Human_heart_combined.integrated) <- "RNA"
all.genes <- rownames(Human_heart_combined.integrated)
Human_heart_combined.integrated_rescaled <- ScaleData(Human_heart_combined.integrated, features = all.genes)
DoHeatmap(Human_heart_combined.integrated_rescaled, features = top10$gene, size = 3) +
  theme(text = element_text(size = 4))

# Plot UMAP for figures
DimPlot(Human_heart_combined.integrated, reduction = "umap"  , label = F, label.size = 6) + NoAxes() 
FeaturePlot(Human_heart_combined.integrated, "SORBS2",reduction = "umap"  , label = F, label.size = 6) + NoAxes() 

markergenes <- c("RYR2","ATP2A2","TNNT2","TNNI3")
markergenes <- c("PECAM1","FLT1","DCN","COL1A2")
markergenes <- c("RGS5","PDGFRB","CALD1","CNN1")
markergenes <- c("MS4A6A","TRBC2","CCL4","CD163")

p <- FeaturePlot(Human_heart_combined.integrated, markergenes, combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + NoAxes()
}

cowplot::plot_grid(plotlist = p)

# Saving images
outputDir <- "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/Images all cells"

ggsave(paste0(outputDir,'/All cells Sorbs2','.tiff'),
       units='mm',width=250,height=150,dpi=300)

# CM = cluster 0, 1, 2, 3, 4, 13.
FeaturePlot(Human_heart_combined.integrated, features = c("NPPB"))
FeaturePlot(Human_heart_combined.integrated, features = c("NPPA"))
FeaturePlot(Human_heart_combined.integrated, features = c("MYH6"))
FeaturePlot(Human_heart_combined.integrated, features = c("MYH7"))
FeaturePlot(Human_heart_combined.integrated, features = c("TNNT2"))      
FeaturePlot(Human_heart_combined.integrated, features = c("ATP2A2"))
FeaturePlot(Human_heart_combined.integrated, features = c("ACTC1"))
FeaturePlot(Human_heart_combined.integrated, features = c("RYR2"))
FeaturePlot(Human_heart_combined.integrated, features = c("TNNI3"))

# cluster 5, 6 = periostin - EC
FeaturePlot(Human_heart_combined.integrated, features = c("PECAM1"))
FeaturePlot(Human_heart_combined.integrated, features = c("FLT1"))
FeaturePlot(Human_heart_combined.integrated, features = c("PLVAP"))

# cluster 7 = periostin/plvap + EC
FeaturePlot(Human_heart_combined.integrated, features = c("POSTN"))
FeaturePlot(Human_heart_combined.integrated, features = c("PECAM1"))

FeaturePlot(Human_heart_combined.integrated, features = c("CCDC80"))

# FB = cl 8, 10
FeaturePlot(Human_heart_combined.integrated, features = c("PDGFRA")) 
FeaturePlot(Human_heart_combined.integrated, features = c("COL1A2")) 
FeaturePlot(Human_heart_combined.integrated, features = c("COL6A3")) 
FeaturePlot(Human_heart_combined.integrated, features = c("DCN"))

# cluster 9 pericytes/stromal cells
FeaturePlot(Human_heart_combined.integrated, features = c("RGS5")) # pericyte marker
FeaturePlot(Human_heart_combined.integrated, features = c("SPARC")) # stromal marker
FeaturePlot(Human_heart_combined.integrated, features = c("PDGFRB")) # pericyte marker

# cluster 11 SMC
FeaturePlot(Human_heart_combined.integrated, features = c("ACTA2")) # smc marker
FeaturePlot(Human_heart_combined.integrated, features = c("CALD1")) # smc marker
FeaturePlot(Human_heart_combined.integrated, features = c("CNN1")) # smc marker

# Monocytes/lymphocytes = cluster 12
FeaturePlot(Human_heart_combined.integrated, features = c("MRC1")) # monocyte
FeaturePlot(Human_heart_combined.integrated, features = c("CCL3")) # monocyte
FeaturePlot(Human_heart_combined.integrated, features = c("CCL4")) # monocyte
FeaturePlot(Human_heart_combined.integrated, features = c("TRBC2")) # T cell receptor

FeaturePlot(Human_heart_combined.integrated, features = c("MS4A6A"))

# EPC = part of cluster 8
FeaturePlot(Human_heart_combined.integrated, features = c("ANXA8")) 
FeaturePlot(Human_heart_combined.integrated, features = c("UPK1B")) 
FeaturePlot(Human_heart_combined.integrated, features = c("MSLN")) 

FeaturePlot(Human_heart_combined.integrated, features = c("LOXL2")) 

head(Human_heart_combined.integrated.markers[Human_heart_combined.integrated.markers$cluster == 12,])
Human_heart_combined.integrated.markers[Human_heart_combined.integrated.markers$cluster == 8,]

# assigning cell types based on marker gene expression

new.cluster.ids <- c("CM", "CM", "CM", "CM","CM", "EC", "EC", "EC", "FB",
                      "Pericytes","FB", "SMC", "Monocytes/Lymphocytes", "CM")
names(new.cluster.ids) <- levels(Human_heart_combined.integrated)
Human_heart_combined.integrated <- RenameIdents(Human_heart_combined.integrated, new.cluster.ids)
DimPlot(Human_heart_combined.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

DefaultAssay(object = Human_heart_combined.integrated) <- "RNA"
Human_heart_combined.integrated <- NormalizeData(Human_heart_combined.integrated)
CM_vs_all <- FindMarkers(Human_heart_combined.integrated, ident.1 = "CM", logfc.threshold = 0, only.pos = T, min.pct = 0)

CM_vs_all$avg_FC <- exp(CM_vs_all$avg_logFC)
CM_vs_all$avg_log2FC <- log2(CM_vs_all$avg_FC)
write.table(CM_vs_all, "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/DEG/CM_vs_all_fc0_pct0.txt", sep = "\t")


######### Get the DEG between CM vs non-CM in normal controls

Human_heart_combined.integrated <- subset(Human_heart_combined.integrated, subset = condition == "N") 
DefaultAssay(object = Human_heart_combined.integrated) <- "RNA"
Human_heart_combined.integrated <- NormalizeData(Human_heart_combined.integrated)
CM_vs_all <- FindMarkers(Human_heart_combined.integrated, ident.1 = "CM", logfc.threshold = 0, only.pos = T, min.pct = 0)

CM_vs_all$avg_FC <- exp(CM_vs_all$avg_logFC)
CM_vs_all$avg_log2FC <- log2(CM_vs_all$avg_FC)
write.table(CM_vs_all, "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/DEG/N_CM_vs_all_fc0_pct0.txt", sep = "\t")

# select all LV CM's

# Return A Subset Of The Seurat Object
# Creates a Seurat object containing only a subset of the cells in the original object. 
# Takes either a list of cells to use as a subset, or a parameter (for example, a gene), to subset on.

Human_heart_combined.integrated_CM <- subset(Human_heart_combined.integrated, idents = "CM") # left are 6783 cells

# Now select the LV
Human_heart_combined.integrated_CM <- subset(Human_heart_combined.integrated_CM, subset = location == "LV") # left are 3253 cells 

# If you want to work on a subset (e.g. dHF and control), select the dilated HF and control subject
# Idents(Human_heart_combined.integrated_CM) <- "individual"
# Human_heart_combined.integrated_CM_d <- subset(Human_heart_combined.integrated_CM, idents = c("D1","D2","D4","D5","N1","N2","N3","N4","N5","N13","N14"))

################################################################################################################
# If you want to calculate correlations of genes use the code below, otherwise skip the lines below
################################################################################################################
# If testing correlated genes in only HF or N failure set, or NPPA positive CM:
Human_heart_combined.integrated_CM <- subset(Human_heart_combined.integrated_CM, subset = condition == "HF") 
Human_heart_combined.integrated_CM <- subset(Human_heart_combined.integrated_CM, subset = NPPA > 1) 

DefaultAssay(object = Human_heart_combined.integrated_CM) <- "RNA"
Human_heart_combined.integrated_CM <- NormalizeData(Human_heart_combined.integrated_CM)
library("Hmisc")

matrix<- GetAssayData(object = Human_heart_combined.integrated_CM, slot = "data")
matrix_mod <- as.matrix(matrix)
gene <- as.numeric(matrix_mod["NPPA",])
correlations <- apply(matrix_mod,1,function(x){cor.test(gene,x)}) # the function cor only returns correlation coefficient and no p-value

correlation_data_frame <- data.frame("Gene" = character(), "Pearsons R" = numeric(), "p-value" = numeric())
correlation_data_frame <- rbind(correlation_data_frame, list(NA,NA,NA))   

for (i in 1:length(correlations)){
  correlation_data_frame <- rbind(correlation_data_frame, c(correlation_data_frame$Gene <- names(correlations[i]),
                                                            correlation_data_frame$Pearsons.R <- as.numeric(correlations[[i]]$estimate),
                                                            correlation_data_frame$p.value <- as.numeric(correlations[[i]]$p.value)))
}
  
colnames(correlation_data_frame) <- c("Gene","Pearsons R","p-value")
correlation_data_frame <- correlation_data_frame[-1,]
rownames(correlation_data_frame) <- correlation_data_frame[,1]

correlation_data_frame$`Pearsons R` <- as.numeric(correlation_data_frame$`Pearsons R`)
correlation_data_frame$`p-value` <- as.numeric(correlation_data_frame$`p-value`)

write.table(correlation_data_frame, "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/Correlated genes/NPPB_LVCM.txt", sep = "\t")

BiocManager::install("EnhancedVolcano")
library(ggrepel)
library(EnhancedVolcano)
library(magrittr)

EnhancedVolcano(correlation_data_frame,
                lab = rownames(correlation_data_frame),
                x = 'Pearsons R',
                y = 'p-value',
                subtitle = "",
                xlim = c(-0.7,0.7),
                title = 'SORBS2 correlated genes',
                xlab = "Pearsons R",
                gridlines.minor = F,
                FCcutoff = 0.1,
                pointSize = 1,
                pCutoff = 0.05,
                labSize = 4,
                col = c("black","black","black","firebrick"),
                selectLab = c("NPPB","SORBS2"),
                drawConnectors = T,
                #widthConnectors = 1.0,
                #colConnectors = 'black',
                boxedLabels = T,
                legendPosition = "none",
                legendLabels = c('', "",
                                 "", "")
)

outputDir <- "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/Correlated genes"

ggsave(paste0(outputDir,'/SORBS2 correlated genes_R0.1','.tiff'),
       units='mm',width=200,height=150,dpi=600)

###############################################################################################################################
################################################################################################################################
################################################################################################################################

# If you want to compare cells based on gene expression (i.e. NPPA negative vs positive cells) use code below, otherwise skip these lines

# If testing NPPA positive vs negative in only HF or N failure set:
Human_heart_combined.integrated_CM <- subset(Human_heart_combined.integrated_CM, subset = condition == "N") # left are 3253 cells 

DefaultAssay(Human_heart_combined.integrated_CM) <- "RNA"
Human_heart_combined.integrated_CM <- NormalizeData(Human_heart_combined.integrated_CM)

FeaturePlot(Human_heart_combined.integrated_CM, features = "NPPA")

NPPA_pos <- subset(Human_heart_combined.integrated_CM, NPPA> 2) # left are 578 cells (if working with HF dataset)
NPPA_neg <- subset(Human_heart_combined.integrated_CM, NPPA< 2) # left are 2569 cells (if working with whole dataset)
Human_heart_combined.integrated_CM <- SetIdent(Human_heart_combined.integrated_CM, cells = Cells(NPPA_pos), value = "NPPA pos")
Human_heart_combined.integrated_CM <- SetIdent(Human_heart_combined.integrated_CM, cells = Cells(NPPA_neg), value = "NPPA neg")

DefaultAssay(object = Human_heart_combined.integrated_CM) <- "RNA"
Human_heart_combined.integrated_CM <- NormalizeData(Human_heart_combined.integrated_CM)

NPPA_DEG <- FindMarkers(Human_heart_combined.integrated_CM, ident.1 = "NPPA pos", ident.2 = "NPPA neg", min.pct = 0, logfc.threshold = 0.1)
NPPA_DEG$avg_FC <- exp(NPPA_DEG$avg_logFC)
NPPA_DEG$avg_log2FC <- log2(NPPA_DEG$avg_FC)

write.table(NPPA_DEG, "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/DEG/Analysis from all cells followed by subsetting LV/NPPA_pos_vs_neg_both_HF+N_LVCM.txt", sep = "\t")
write.table(NPPA_DEG, "/Users/l.timmer/Desktop/listRayan_NPPA_DEG_LVCM.txt", sep = "\t")
Human_heart_combined.integrated_CM <- NPPA_pos
# Select RNA and perform either SCTransform or again FindVariableFeatures and Scale data. 
# If too much influence of indivduals, we could continue with the integrated data set (set "integrated" in stead of "RNA")

# Notes form Seurat team:
# Would not recommend repeating the integration on a subset of the cells, using the integrated assay computed on the full dataset should be sufficient for subclustering.

####################################################################################################
# Conclusion integrated data: No bias due to individuals, no seperation due to condition,         ##
# but gradients due to condition which are confirmed by stress marker genes as Nppa, Nppb, Tmsb4x ##
####################################################################################################

# if working on a specific subset; e.g. dHF + N CM set: 
# DefaultAssay(object = Human_heart_combined.integrated_CM_d) <- "integrated"
# Human_heart_combined.integrated_CM <- Human_heart_combined.integrated_CM_d

# if working on all LV CM:
DefaultAssay(object = Human_heart_combined.integrated_CM) <- "integrated"

# note that if you wish to perform additional rounds of clustering after subsetting we recommend
# re-running FindVariableFeatures() and ScaleData()

Human_heart_combined.integrated_CM <- FindVariableFeatures(Human_heart_combined.integrated_CM, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10<- head(VariableFeatures(Human_heart_combined.integrated_CM), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Human_heart_combined.integrated_CM)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + NoLegend()
plot2

# Scale data 
# Scaling data is used to give lowly and highly expressed genes equal weight such that 
# these can be used for dimensionality reduction. Scaled data should therefore never be used for DEG
all.genes <- rownames(Human_heart_combined.integrated_CM)
Human_heart_combined.integrated_CM <- ScaleData(Human_heart_combined.integrated_CM, features = all.genes)

Human_heart_combined.integrated_CM <- RunPCA(Human_heart_combined.integrated_CM, features = VariableFeatures(object = Human_heart_combined.integrated_CM))

print(Human_heart_combined.integrated_CM[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Human_heart_combined.integrated_CM, dims = 1:2, reduction = "pca")
DimHeatmap(Human_heart_combined.integrated_CM, dims = 1:10, cells = 500, balanced = TRUE)
ElbowPlot(Human_heart_combined.integrated_CM)

Human_heart_combined.integrated_CM <- JackStraw(Human_heart_combined.integrated_CM, num.replicate = 100)
Human_heart_combined.integrated_CM <- ScoreJackStraw(Human_heart_combined.integrated_CM, dims = 1:20)
JackStrawPlot(Human_heart_combined.integrated_CM, dims = 1:20)

# based on the elbowplot we set dimensions at 4
Human_heart_combined.integrated_CM <- RunUMAP(Human_heart_combined.integrated_CM, dims = 1:4)

# Cluster cells
Human_heart_combined.integrated_CM <- FindNeighbors(Human_heart_combined.integrated_CM, dims = 1:4)
Human_heart_combined.integrated_CM <- FindClusters(Human_heart_combined.integrated_CM, resolution = 0.3)

DimPlot(Human_heart_combined.integrated_CM, reduction = "umap", label = TRUE)
DimPlot(Human_heart_combined.integrated_CM, reduction = "umap", group.by = "condition")
DimPlot(Human_heart_combined.integrated_CM, reduction = "umap", group.by = "individual")
DimPlot(Human_heart_combined.integrated_CM, reduction = "umap", group.by = "location")

# Switch to normalized data for DEG
DefaultAssay(object = Human_heart_combined.integrated_CM) <- "RNA"
Human_heart_combined.integrated_CM <- NormalizeData(Human_heart_combined.integrated_CM)

FeaturePlot(Human_heart_combined.integrated_CM, features = "TTN")
FeaturePlot(Human_heart_combined.integrated_CM, features = "MYH6")
FeaturePlot(Human_heart_combined.integrated_CM, features = "MYL4")
FeaturePlot(Human_heart_combined.integrated_CM, features = "TNNT2")
FeaturePlot(Human_heart_combined.integrated_CM, features = "IFITM3")
FeaturePlot(Human_heart_combined.integrated_CM, features = "PECAM1")
FeaturePlot(Human_heart_combined.integrated_CM, features = "NPPA")
FeaturePlot(Human_heart_combined.integrated_CM, features = "TMSB4X")
FeaturePlot(Human_heart_combined.integrated_CM, features = "HSPA1B")

# Can we see a gradient of the genes below? 
FeaturePlot(Human_heart_combined.integrated_CM, features = "ACTA1") # YES
FeaturePlot(Human_heart_combined.integrated_CM, features = "SORBS2") # NO
FeaturePlot(Human_heart_combined.integrated_CM, features = "TNFRSF12A") # YES
FeaturePlot(Human_heart_combined.integrated_CM, features = "MYH6") # YES
FeaturePlot(Human_heart_combined.integrated_CM, features = "MYL4") # YES
FeaturePlot(Human_heart_combined.integrated_CM, features = "DKK3") # YES
FeaturePlot(Human_heart_combined.integrated_CM, features = "SP1") # TOO FEW CELLS
FeaturePlot(Human_heart_combined.integrated_CM, features = "LRP2") # TOO FEW CELLS
FeaturePlot(Human_heart_combined.integrated_CM, features = "LRP5") # TOO FEW CELLS
FeaturePlot(Human_heart_combined.integrated_CM, features = "AP5Z1") # TOO FEW CELLS
FeaturePlot(Human_heart_combined.integrated_CM, features = "CSTF3") # TOO FEW CELLS
FeaturePlot(Human_heart_combined.integrated_CM, features = "SPP1") # YES
FeaturePlot(Human_heart_combined.integrated_CM, features = "VAV2") # TOO FEW CELLS
FeaturePlot(Human_heart_combined.integrated_CM, features = "TMEM101") # NO
FeaturePlot(Human_heart_combined.integrated_CM, features = "ANKRD1") # NO
FeaturePlot(Human_heart_combined.integrated_CM, features = "XIRP2") # REVERSED
FeaturePlot(Human_heart_combined.integrated_CM, features = "TTN") #  REVERSED
FeaturePlot(Human_heart_combined.integrated_CM, features = "TNNT2") # NO
FeaturePlot(Human_heart_combined.integrated_CM, features = "CTGF") # YES (BUT VERY FEW CELLS)
FeaturePlot(Human_heart_combined.integrated_CM, features = "LGALS3") # NO
FeaturePlot(Human_heart_combined.integrated_CM, features = "CLU") # NO
FeaturePlot(Human_heart_combined.integrated_CM, features = "BGN") # YES (BUT VERY FEW CELLS)
FeaturePlot(Human_heart_combined.integrated_CM, features = "GPX3") # NO
FeaturePlot(Human_heart_combined.integrated_CM, features = "CSRP3") # NO
FeaturePlot(Human_heart_combined.integrated_CM, features = "SPRR1A") # TOO FEW CELLS
FeaturePlot(Human_heart_combined.integrated_CM, features = "MMP9") # TOO FEW CELLS
FeaturePlot(Human_heart_combined.integrated_CM, features = "SYNPO2L") # NO
FeaturePlot(Human_heart_combined.integrated_CM, features = "PRG4") # TOO FEW CELLS
FeaturePlot(Human_heart_combined.integrated_CM, features = "SERPINA3") # TOO FEW CELLS
FeaturePlot(Human_heart_combined.integrated_CM, features = "ANKRD23") # TOO FEW CELLS
FeaturePlot(Human_heart_combined.integrated_CM, features = "COL8A1") # TOO FEW CELLS
FeaturePlot(Human_heart_combined.integrated_CM, features = "LOX") # TOO FEW CELLS
FeaturePlot(Human_heart_combined.integrated_CM, features = "COMP") # TOO FEW CELLS
FeaturePlot(Human_heart_combined.integrated_CM, features = "MB") # NO
FeaturePlot(Human_heart_combined.integrated_CM, features = "LTBP2") # TOO FEW CELLS
FeaturePlot(Human_heart_combined.integrated_CM, features = "POSTN") # NO

Human_heart_combined.integrated_CM.markers <- FindAllMarkers(Human_heart_combined.integrated_CM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Human_heart_combined.integrated_CM.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top10 <- Human_heart_combined.integrated_CM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Human_heart_combined.integrated_CM, features = top10$gene) + NoLegend()

head(Human_heart_combined.integrated_CM.markers[Human_heart_combined.integrated_CM.markers$cluster == 4,])
Human_heart_combined.integrated_CM.markers[Human_heart_combined.integrated_CM.markers$cluster ==4,]

# Make a heatmap plot for the N and HF cardiomyocytes
Idents(object = Human_heart_combined.integrated_CM) <- "condition"
Human_heart_combined.integrated_CM.markers <- FindAllMarkers(Human_heart_combined.integrated_CM, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

all.genes <- rownames(Human_heart_combined.integrated_CM)
Human_heart_combined.integrated_CM <- ScaleData(Human_heart_combined.integrated_CM, features = all.genes)

top10 <- Human_heart_combined.integrated_CM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
Human_heart_combined.integrated_CM.markers$avg_FC <- exp(Human_heart_combined.integrated_CM.markers$avg_logFC)
top10 <- Human_heart_combined.integrated_CM.markers[Human_heart_combined.integrated_CM.markers$p_val_adj<0.05 & Human_heart_combined.integrated_CM.markers$avg_FC>1.5,]

DoHeatmap(Human_heart_combined.integrated_CM, features = top10$gene, size = 3) +
  theme(text = element_text(size = 4))

# Saving images
outputDir <- "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/Images LV CM/LV CM integrated assay"

ggsave(paste0(outputDir,'/Violin LVCM SORBS2','.tiff'),
       units='mm',width=150,height=120,dpi=300)

# Make a violin plot for GOI between HF and N.
VlnPlot(Human_heart_combined.integrated_CM, features = "SORBS2")

cl0_vs_all <- Human_heart_combined.integrated_CM.markers[Human_heart_combined.integrated_CM.markers$cluster == 0,]
cl1_vs_all <- Human_heart_combined.integrated_CM.markers[Human_heart_combined.integrated_CM.markers$cluster == 1,]
cl2_vs_all <- Human_heart_combined.integrated_CM.markers[Human_heart_combined.integrated_CM.markers$cluster == 2,]
cl3_vs_all <- Human_heart_combined.integrated_CM.markers[Human_heart_combined.integrated_CM.markers$cluster == 3,]
cl4_vs_all <- Human_heart_combined.integrated_CM.markers[Human_heart_combined.integrated_CM.markers$cluster == 4,]
cl5_vs_all <- Human_heart_combined.integrated_CM.markers[Human_heart_combined.integrated_CM.markers$cluster == 5,]

# Add FC and Log2FC
cl0_vs_all$avg_FC <- exp(cl0_vs_all$avg_logFC)
cl0_vs_all$avg_log2FC <- log2(cl0_vs_all$avg_FC)

cl1_vs_all$avg_FC <- exp(cl1_vs_all$avg_logFC)
cl1_vs_all$avg_log2FC <- log2(cl1_vs_all$avg_FC)

cl2_vs_all$avg_FC <- exp(cl2_vs_all$avg_logFC)
cl2_vs_all$avg_log2FC <- log2(cl2_vs_all$avg_FC)

cl3_vs_all$avg_FC <- exp(cl3_vs_all$avg_logFC)
cl3_vs_all$avg_log2FC <- log2(cl3_vs_all$avg_FC)

cl4_vs_all$avg_FC <- exp(cl4_vs_all$avg_logFC)
cl4_vs_all$avg_log2FC <- log2(cl4_vs_all$avg_FC)

cl5_vs_all$avg_FC <- exp(cl5_vs_all$avg_logFC)
cl5_vs_all$avg_log2FC <- log2(cl5_vs_all$avg_FC)

# Get files of DEG per LV CM cluster
write.table(cl0_vs_all, "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/DEG/Analysis from all cells followed by subsetting LV/Enriched_genes_per_CM_cluster/cl0_vs_all.txt", sep = "\t")
write.table(cl1_vs_all, "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/DEG/Analysis from all cells followed by subsetting LV/Enriched_genes_per_CM_cluster/cl1_vs_all.txt", sep = "\t")
write.table(cl2_vs_all, "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/DEG/Analysis from all cells followed by subsetting LV/Enriched_genes_per_CM_cluster/cl2_vs_all.txt", sep = "\t")
write.table(cl3_vs_all, "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/DEG/Analysis from all cells followed by subsetting LV/Enriched_genes_per_CM_cluster/cl3_vs_all.txt", sep = "\t")
write.table(cl4_vs_all, "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/DEG/Analysis from all cells followed by subsetting LV/Enriched_genes_per_CM_cluster/cl4_vs_all.txt", sep = "\t")
write.table(cl5_vs_all, "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/DEG/Analysis from all cells followed by subsetting LV/Enriched_genes_per_CM_cluster/cl5_vs_all.txt", sep = "\t")

# Plot UMAP for figures
DimPlot(Human_heart_combined.integrated_CM, reduction = "umap", group.by = "individual" ,label = F, label.size = 6) + NoAxes() 
FeaturePlot(Human_heart_combined.integrated_CM, reduction = "umap", "SORBS2", label = F, label.size = 6) + NoAxes() 

markergenes <- c("TNNI3","RYR2", "PECAM1","FLT1","PDGFRA","COL1A2","ACTA2")
markergenes <- c("NPPA", "NPPB", "ACTA1","MYH6","TNFRSF12A","SORBS2","MYL4","DKK3")
p <- FeaturePlot(Human_heart_combined.integrated_CM, markergenes, combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]]  + NoAxes()
}

cowplot::plot_grid(plotlist = p)

# Saving images
outputDir <- "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/Images LV CM/LV CM integrated assay"

ggsave(paste0(outputDir,'/LV CM SORBS2','.tiff'),
       units='mm',width=250,height=150,dpi=600)

# DEG between conditions
DefaultAssay(object = Human_heart_combined.integrated_CM) <- "RNA"
Human_heart_combined.integrated_CM <- NormalizeData(Human_heart_combined.integrated_CM)
Idents(object = Human_heart_combined.integrated_CM) <- Human_heart_combined.integrated_CM@meta.data$individual
cHFvsN <- FindMarkers(Human_heart_combined.integrated_CM, ident.1 = c("C1","C2"), ident.2 = c("N1","N2","N3","N4","N5","N13","N14"), logfc.threshold = 0)
dHFvsN <- FindMarkers(Human_heart_combined.integrated_CM, ident.1 = c("D1","D2","D4","D5"), ident.2 = c("N1","N2","N3","N4","N5","N13","N14"), logfc.threshold = 0)

Idents(object = Human_heart_combined.integrated_CM) <- Human_heart_combined.integrated_CM@meta.data$condition
HFvsN <- FindMarkers(Human_heart_combined.integrated_CM, ident.1 = "HF", ident.2 = "N", logfc.threshold = 0)

HFvsN$avg_FC <- exp(HFvsN$avg_logFC)
HFvsN$avg_log2FC <- log2(HFvsN$avg_FC)

dHFvsN$avg_FC <- exp(dHFvsN$avg_logFC)
dHFvsN$avg_log2FC <- log2(dHFvsN$avg_FC)

cHFvsN$avg_FC <- exp(cHFvsN$avg_logFC)
cHFvsN$avg_log2FC <- log2(cHFvsN$avg_FC)

write.table(HFvsN, "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/DEG/Analysis from all cells followed by subsetting LV/HFvsN.txt", sep = "\t")
write.table(cHFvsN, "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/DEG/Analysis from all cells followed by subsetting LV/cHFvsN.txt", sep = "\t")
write.table(dHFvsN, "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/DEG/Analysis from all cells followed by subsetting LV/dHFvsN.txt", sep = "\t")

# Plot UMAP for figures
DimPlot(Human_heart_combined.integrated_CM, reduction = "pca", group.by = "individual"  ,label = F,repel = T, label.size = 6) + NoAxes() 

FeaturePlot(Human_heart_combined.integrated_CM, "SORBS2",reduction = "umap", label = F,repel = T, label.size = 6) + NoAxes() 

markergenes <- c("NPPA","NPPB", "MALAT1","TMSB4X","IFI27","RGS5")
p <- FeaturePlot(Human_heart_combined.integrated_CM, markergenes, combine = FALSE, reduction = "pca")

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

cowplot::plot_grid(plotlist = p)

# Saving images
outputDir <- "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/Images LV CM/"

ggsave(paste0(outputDir,'vulcano LV CM cHFvsN','.tiff'),
       units='mm',width=150,height=120,dpi=600)

# Plot Vulcanopolot

library(ggrepel)
library(EnhancedVolcano)
library(magrittr)

EnhancedVolcano(HFvsN,
                lab = rownames(HFvsN),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlim = c(-3,5),
                title = 'Heart Failure vs Normal',
                subtitle = "",
                FCcutoff = log2(1.5),
                pCutoff = 0.05,
                pointSize = 2,
                drawConnectors = T,
                gridlines.minor = F,
                boxedLabels = T,
                labhjust = -0.5,
                labvjust = 1,
                selectLab = c("NPPA","ACTA1", "MYL2", "MALAT1","SORBS2","NPPB",
                              "LOXL1", "HSPA1B")
)

EnhancedVolcano(dHFvsN,
                lab = rownames(dHFvsN),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlim = c(-3,5),
                title = 'dilated Heart Failure vs Normal',
                subtitle = "",
                FCcutoff = log2(1.5),
                pCutoff = 0.05,
                pointSize = 2,
                drawConnectors = T,
                gridlines.minor = F,
                boxedLabels = T,
                labhjust = -0.5,
                labvjust = 1,
                selectLab = c("NPPA","ACTA1", "MYL2", "MALAT1","SORBS2","NPPB",
                              "LOXL1", "HSPA1B")
)


EnhancedVolcano(cHFvsN,
                lab = rownames(cHFvsN),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlim = c(-3,5),
                title = 'coronary Heart Failure vs Normal',
                subtitle = "",
                FCcutoff = log2(1.5),
                pCutoff = 0.05,
                pointSize = 2,
                drawConnectors = T,
                gridlines.minor = F,
                boxedLabels = T,
                labhjust = -0.5,
                labvjust = 1,
                selectLab = c("NPPA","ACTA1", "MYL2", "MALAT1","SORBS2","NPPB",
                              "LOXL1", "HSPA1B")
)

outputDir <- "/Users/l.timmer/Documents/Seq analysis/Sequence_Analysis_LT/Heart failure project/Study_Wang_et_al/Images LV CM/"

ggsave(paste0(outputDir,'vulcano LV CM cHFvsN_beautified','.jpeg'),
       units='mm',width=150,height=120,dpi=300)

# Calculate cellular distribution per cluster of LV CM

# Total, hence this is the expected distribution if all cells would be equally divided
ncol(Human_heart_combined.integrated_CM) # 3253
ncol(subset(Human_heart_combined.integrated_CM, subset = condition == "HF")) # 827 (25%)
ncol(subset(Human_heart_combined.integrated_CM, subset = condition == "N")) # 2426 (75%)

# Cluster 0: 
ncol(subset(Human_heart_combined.integrated_CM, subset = seurat_clusters == 0)) # 1031
ncol(subset(Human_heart_combined.integrated_CM, subset = seurat_clusters == 0 & condition == "HF")) # 379 (37%%)
ncol(subset(Human_heart_combined.integrated_CM, subset = seurat_clusters == 0 & condition == "N")) # 652 (63%%)

# Cluster 1: 
ncol(subset(Human_heart_combined.integrated_CM, subset = seurat_clusters == 1)) # 836
ncol(subset(Human_heart_combined.integrated_CM, subset = seurat_clusters == 1 & condition == "HF")) # 142 (17%)
ncol(subset(Human_heart_combined.integrated_CM, subset = seurat_clusters == 1 & condition == "N")) # 694 (83%)

# Cluster 2: 
ncol(subset(Human_heart_combined.integrated_CM, subset = seurat_clusters == 2)) # 623
ncol(subset(Human_heart_combined.integrated_CM, subset = seurat_clusters == 2 & condition == "HF")) # 67 (11%)
ncol(subset(Human_heart_combined.integrated_CM, subset = seurat_clusters == 2 & condition == "N")) # 556 (89%)

# Cluster 3: 
ncol(subset(Human_heart_combined.integrated_CM, subset = seurat_clusters == 3)) # 500
ncol(subset(Human_heart_combined.integrated_CM, subset = seurat_clusters == 3 & condition == "HF")) # 168 (34%%)
ncol(subset(Human_heart_combined.integrated_CM, subset = seurat_clusters == 3 & condition == "N")) # 332 (66%)

# Cluster 4: 
ncol(subset(Human_heart_combined.integrated_CM, subset = seurat_clusters == 4)) # 263
ncol(subset(Human_heart_combined.integrated_CM, subset = seurat_clusters == 4 & condition == "HF")) # 71 (27%%)
ncol(subset(Human_heart_combined.integrated_CM, subset = seurat_clusters == 4 & condition == "N")) # 192 (73%)

