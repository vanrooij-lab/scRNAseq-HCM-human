
# OK, so Seurat didn't give the best results for this dataset in terms of correction, so 
# let's try the "Combat" package.

# See: 
# Zhang, Y., Parmigiani, G., & Johnson, W. E. (2020). ComBat-seq: batch effect adjustment for RNA-seq count data. NAR Genomics and Bioinformatics, 2(3), 1–10. https://doi.org/10.1093/nargab/lqaa078
# See also e.g.:
# Tutorials:
# https://rnabio.org/module-03-expression/0003/05/01/Batch-Correction/
# https://towardsdatascience.com/how-to-batch-correct-single-cell-7bad210c7ae1
# https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html
# Andrews, T. S., Kiselev, V. Y., McCarthy, D., & Hemberg, M. (2021). Tutorial: guidelines for the computational analysis of single-cell RNA sequencing data. Nature Protocols, 16(1). https://doi.org/10.1038/s41596-020-00409-w

################################################

# Previous analysis:
# load("/Volumes/workdrive_m.wehrens_hubrecht/temporary_files/WRL_data_seurat-2021-06.RData")

# custom libs
source('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Functions_2021/Load-Pool-Scale_Simple_MW.R')

# some libs
library(stringr)
library(Seurat)
library(ggplot2)

# Combat-seq related
library(BiocManager)
library("sva") # BiocManager::install('sva')

# Some more libs
library("gridExtra")
library("edgeR")
library("UpSetR")

library(reshape)

# checking script performance/progress
library(pryr)
library(beepr)

# Some colors
library(RColorBrewer)
rainbow_colors = rev( brewer.pal(n = 11, name = "Spectral") )
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector_60 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Again some markers
markers=list()
markers$CM=c('TTN','MYH7','TNNT2','MYH6')
markers$EC=c('VWF','IFI27')
markers$FB=c('DCN','COL1A2')
markers$SMC=c('ACTA2','CALD1','MYH11')

################################################

# Let's first load the data here again, now also including the selection of Teichmann cells
data_dir1='/Volumes/fastq_m.wehrens/Mapping/HCM_SCS/mapping.93.may25/counttables/'
data_dir2='/Volumes/fastq_m.wehrens/Mapping/WANG2/counttables/'
dataset_list_paths=list('AL1'=paste0(data_dir1, 'HUB-AL-s001_HG25TBGXF_S5_cat_pT_uniaggGenes_spliced.UFICounts.tsv'),
                        'JE5'=paste0(data_dir1, 'JE5_AHFL77BGX5_S6_cat_pT_uniaggGenes_spliced.UFICounts.tsv'),
                        'WANG13'=paste0(data_dir2, 'GSM3449619_N13_cat_nc_uniaggGenes_spliced.UFICounts.tsv'),
                        'WANG14'=paste0(data_dir2, 'GSM3449620_N14_cat_nc_uniaggGenes_spliced.UFICounts.tsv'))
# Load data
SCS_df_list_data = loadData_MW(dataset_list_paths, toPool = NULL)
# Teichmann is for now in saved workspace (small subselection)
load('/Volumes/workdrive_m.wehrens_hubrecht/data/2020_12_Litvinukova2020/HPC_sync/Teichmann_subset.Rdata')

# To look at some old names:
# test =loadData_MW(dataset_list_paths[1], toPool = NULL)
# rm('test')

# Note however, that gene names are now not compatible..
rownames(Teichmann_Sampled)
rownames(SCS_df_list_data$AL1)

SCS_df_list_data = lapply(SCS_df_list_data, preprocess_convertAAnames_toSymbol)

# Now extend the list with Teichmann data
SCS_df_list_data$LITVI_samp = as.data.frame(Teichmann_Sampled) # this requires Seurat lib to be loaded
colnames(SCS_df_list_data$LITVI_samp) = paste0('LITVI.',colnames(SCS_df_list_data$LITVI_samp))
rownames(SCS_df_list_data$LITVI_samp) = gsub(pattern = '-',replacement = '\\.',x = rownames(SCS_df_list_data$LITVI_samp))

# And now pool those (will take a while)
# This might be something that needs to be done @ the HPC, 
# Now takes ±3 mins and object size is ±2GB
start_time <- Sys.time()
SCS_df_pooled_data = pool_df_mw(SCS_df_list_data)
end_time <- Sys.time()
print(end_time - start_time)
beepr::beep()
object_size(SCS_df_pooled_data)
object_size(SCS_df_list_data) 

# End of custom code
########################################################################

########################################################################
# Start of Seurat / Combat analysis

# Now convert this object to Seurat 

# Before starting, let's deduce some annotation
WRL_annotation_samples=unlist(lapply(  strsplit(colnames(SCS_df_pooled_data), split='\\.'), 
    function(splitted_string) {splitted_string[1]}))

# Initialize object
WRL_object = CreateSeuratObject(SCS_df_pooled_data, project = 'WRL')
object_size(WRL_object) 
    # clearly better mem management
    # perhaps eventually I should use Seurat merging to merge stuff?

# Add annotation
WRL_object[["annotation_sample_fct"]]     <- as.factor(WRL_annotation_samples)
WRL_object[["annotation_sample_str"]]     <- WRL_annotation_samples

# Count mitochondrial reads
rownames(WRL_object)[grepl("^MT\\.",rownames(WRL_object))]
WRL_object[["percent.mt"]] <- PercentageFeatureSet(WRL_object, pattern = "^MT\\.")

VlnPlot(WRL_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'annotation_sample')

# Remove some cells
# Note that percent.mt < 80 is not very stringent
WRL_object <- subset(WRL_object, subset = nFeature_RNA > 200 & nCount_RNA > 3000 & percent.mt < 80)
# Normalize data
WRL_object <- NormalizeData(WRL_object, normalization.method = "LogNormalize", scale.factor = 10000) # values are defaults
# Find variable features
WRL_object <- FindVariableFeatures(WRL_object, selection.method = "vst", nfeatures = 2000)

# ### some plots var features
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(WRL_object), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(WRL_object)
print(plot1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

all.genes <- rownames(WRL_object)
WRL_object <- ScaleData(WRL_object, features = all.genes)

WRL_object <- RunPCA(WRL_object, npcs = 30, verbose = FALSE)
WRL_object <- RunUMAP(WRL_object, reduction = "pca", dims = 1:30)
WRL_object <- FindNeighbors(WRL_object, reduction = "pca", dims = 1:30)
WRL_object <- FindClusters(WRL_object, resolution = 0.5)


# Analysis done, show some plots

print(WRL_object[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(WRL_object, dims = 1:5, reduction = "pca", ncol = 5)

DimPlot(WRL_object, reduction = "pca")
DimPlot(WRL_object, reduction = "pca", group.by = 'annotation_sample')

JackStrawPlot(WRL_object, dims = 1:5) # doesn't work for some reason?
ElbowPlot(WRL_object)



# Some plots relevant to our scientific questions

# First sample origin
DimPlot(WRL_object, group.by = 'annotation_sample', #cols = col_vector_60,
            label = T, repel = T, label.size = 5)+theme(legend.position = 'none')
DimPlot(WRL_object, reduction = "umap")
#DimPlot(WRL_object, group.by = 'annotation_sample')

FeaturePlot(WRL_object, features = c("MALAT1"), cols = rainbow_colors)
FeaturePlot(WRL_object, features = c("TTN"), cols = rainbow_colors)
FeaturePlot(WRL_object, features = c("MYH7"), cols = rainbow_colors)
FeaturePlot(WRL_object, features = c("MYH6"), cols = rainbow_colors)
FeaturePlot(WRL_object, features = c("NPPA"), cols = rainbow_colors)
FeaturePlot(WRL_object, features = c("NPPB"), cols = rainbow_colors)

p=FeaturePlot(WRL_object, features = markers$CM, cols = rainbow_colors)
FeaturePlot(WRL_object, features = markers$FB, cols = rainbow_colors)
FeaturePlot(WRL_object, features = markers$SMC, cols = rainbow_colors)
FeaturePlot(WRL_object, features = markers$EC, cols = rainbow_colors)

########################################################################
########################################################################
# OK, so this is more or less as expected, the cells are not really overlapping ..
# Let's try to use COMBAT now ..

# quickly determine batches present in the litvi data ..
litvi_batch_names = colnames(WRL_object@assays$RNA@scale.data)[
                        WRL_object[["annotation_sample_str"]][[1]]=='LITVI'
                    ]   
litvi_batch_names = sapply(str_split(litvi_batch_names,'-'), function(x) {x[3]})
unique(litvi_batch_names)

#perform the batch correction
#groups = sapply(as.character(conditions), switch, "UHR" = 1, "HBR" = 2, USE.NAMES = F)
#batches = sapply(as.character(library_methods), switch, "Ribo" = 1, "Poly" = 2, USE.NAMES = F)
batches_simple = (WRL_object[["annotation_sample_str"]][[1]])
group_map = c(AL1=1, JE5=1, WANG13=2, WANG14=2, LITVI=3)
groups = group_map[batches_simple]
# However, Litvi data also consists of batches, which
# were already pooled before ..
batches = batches_simple
batches[batches_simple=='LITVI']=litvi_batch_names

df_annotation = data.frame(groups, batches, n=1:length(batches))
df_annotation_sorted=df_annotation[order(df_annotation$groups, df_annotation$batches),]
df_annotation_sorted$ns=1:length(batches)
df_annotation_melted=melt(df_annotation_sorted, id.vars = c('n','ns'), measure.vars = c('groups','batches'))

ggplot(df_annotation_melted)+
    geom_tile(aes(x=ns,y=variable,fill=value))+theme_minimal()+theme(legend.position='none')+
     scale_fill_manual(values = col_vector_60)+give_better_textsize_plot(15)+
    xlab('Sample')+ylab(element_blank())

# Combat-seq can't handle samples with only 1 entry;
# won't matter later, for now make selection n>1
table_batch_names=table(batches)
sel_cells_combat_n2 = table_batch_names[batches]>1

corrected_data = ComBat_seq(counts = as.matrix(WRL_object@assays$RNA@data)[,sel_cells_combat_n2], batch = batches[sel_cells_combat_n2], group = groups[sel_cells_combat_n2])
    # this gives a confounding error (doh)

WRL_object[['annotation_batch']] = batches
WRL_object[['annotation_batch_n']] = as.double(table_batch_names[batches])

# corrected_data = cbind(uncorrected_data[,c("Gene","Chr")], corrected_data)
# 
# #compare dimensions of corrected and uncorrected data sets
# dim(uncorrected_data)
# dim(corrected_data)
# 
# #visually compare values of corrected and uncorrected data sets
# head(uncorrected_data)
# head(corrected_data)

# End of Combat part
########################################################################


########################################################################
# Perhaps the conservative Seurat stuff has something to offer after all?

WRL_object_sel <- subset(WRL_object, subset = annotation_batch_n>50)

WRL_object.list3 <- SplitObject(WRL_object_sel, split.by = "annotation_batch")
WRL_object.list3 <- lapply(X = WRL_object.list3, FUN = SCTransform); beepr::beep()
WRL.features3 <- SelectIntegrationFeatures(object.list = WRL_object.list3, nfeatures = 3000)
WRL_object.list3 <- PrepSCTIntegration(object.list = WRL_object.list3, anchor.features = WRL.features3); beepr::beep()
WRL.anchors3 <- FindIntegrationAnchors(object.list = WRL_object.list3, normalization.method = 'SCT', anchor.features = WRL.features3, dims=1:10, k.anchor = 5); beepr::beep()
WRL.combined.sct <- IntegrateData(anchorset = WRL.anchors3, normalization.method = 'SCT', k.weight = 10); beepr::beep()
WRL.combined.sct <- RunPCA(WRL.combined.sct, verbose = FALSE); beepr::beep()
WRL.combined.sct <- RunUMAP(WRL.combined.sct, reduction = "pca", dims = 1:30); beepr::beep()
WRL.combined.sct <- RunTSNE(WRL.combined.sct, reduction = "pca", dims = 1:30); beepr::beep()

# Plots
DimPlot(WRL.combined.sct, reduction = "umap", group.by = "annotation_sample")+
    ggtitle('UMAP of batch-corrected sample data (Seurat)')
DimPlot(WRL.combined.sct, reduction = "umap", group.by = "annotation_batch")
#p2 <- DimPlot(WRL.combined.sct, reduction = "umap", group.by = 'seurat_annotations',label = TRUE, repel = TRUE)
#p1 + p2

# tsne
DimPlot(WRL.combined.sct, reduction = "tsne", group.by = "annotation_sample")
DimPlot(WRL.combined.sct, reduction = "tsne", group.by = "annotation_batch")

# multi-plot
DimPlot(WRL.combined.sct, reduction = "umap", split.by = "annotation_sample")
DimPlot(WRL.combined.sct, reduction = "tsne", split.by = "annotation_sample")

# Non-corrected plot for reference
DimPlot(WRL_object, reduction = "umap", split.by = "annotation_sample")
DimPlot(WRL_object, reduction = "umap")

p1+FeaturePlot(WRL.combined.sct, features = c("ENSG00000251562-MALAT1-lincRNA"), cols = rainbow_colors)


DimPlot(WRL_object, reduction = "umap")

FeaturePlot(WRL_object, features = c("NPPA"), cols = rainbow_colors)




# Let's create some additional annotation to recognize cells
barcodes = sapply(str_split(pattern = '\\.',string = colnames(WRL_object)), function(x) {x[4]})
conversion_table = c(WANG13='N13', WANG14='N14')
individuals=conversion_table[WRL_object$annotation_sample_str]
WRL_object[['cell_name_mw']]=paste0(individuals,'-',barcodes)
WRL_object[['vCM']] = WRL_object$cell_name_mw %in% desired_cells_mwName

WRL_object_sel_vCM = subset(WRL_object, subset = vCM == T)

p1=DimPlot(WRL_object, reduction = "umap", group.by = 'annotation_sample') # full
print(p1)
p2=DimPlot(WRL_object_sel_vCM, reduction = "umap", group.by = 'annotation_sample') # selection
print(p2)
(p1+xlim(c(-25,15))+ylim(c(-10,15)))+(p2+xlim(c(-25,15))+ylim(c(-10,15)))


# For reference, individual annotation again
DimPlot(WRL_object, reduction = "umap", group.by = 'annotation_sample')

# show whether in selection or not, plus some markers
DimPlot(WRL_object, reduction = "umap", group.by = 'vCM')+
FeaturePlot(WRL_object, features = 'TNNT2', cols = rainbow_colors)+
FeaturePlot(WRL_object, features = 'MYH6', cols = rainbow_colors)+
FeaturePlot(WRL_object, features = 'MYH7', cols = rainbow_colors)

# highlight vCM markers suggested by Anne
markers$vCM = c('RYR2', 'TTN', 'ATP2A', 'ACTC1') # ryr2, ttn, atp2a, actc1
FeaturePlot(WRL_object, features = markers$vCM, cols = rainbow_colors)
# slightly different styled plot, same markers
DimPlot(WRL_object, reduction = "umap", group.by = 'vCM')+
FeaturePlot(WRL_object, features = 'RYR2', cols = rainbow_colors)+
FeaturePlot(WRL_object, features = 'TTN', cols = rainbow_colors)+
FeaturePlot(WRL_object, features = 'ACTC1', cols = rainbow_colors)
# 
(DimPlot(WRL_object, reduction = "umap", group.by = 'annotation_sample')+ggtitle('Source individual'))+
(DimPlot(WRL_object, reduction = "umap", group.by = 'vCM')+ggtitle('vCM selection'))+
(FeaturePlot(WRL_object, features = 'MYL2', cols = rainbow_colors)+ggtitle('MYL2, ventricular'))+
(FeaturePlot(WRL_object, features = 'MYL7', cols = rainbow_colors)+ggtitle('MYL7, atrial'))


# yet another plot
DimPlot(WRL_object_sel_vCM, reduction = "umap", group.by = 'annotation_sample')
FeaturePlot(WRL_object, features = markers$CM, cols = rainbow_colors)
FeaturePlot(WRL_object_sel_vCM, features = markers$CM, cols = rainbow_colors)

# NPPA
(FeaturePlot(WRL_object, features = 'NPPB', cols = rainbow_colors)+xlim(c(-25,15))+ylim(c(-10,15)))+
(FeaturePlot(WRL_object_sel_vCM, features = 'NPPB', cols = rainbow_colors)+xlim(c(-25,15))+ylim(c(-10,15)))



