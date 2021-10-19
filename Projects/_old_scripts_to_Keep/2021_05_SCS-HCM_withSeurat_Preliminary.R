
# This part of the script follows the tutorials at 
# - https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# - https://satijalab.org/seurat/articles/integration_introduction.html

################################################
# My data is in:

library(Seurat)
library(patchwork)

# MW
library(RColorBrewer)
library(stringr)

rainbow_colors = rev( brewer.pal(n = 11, name = "Spectral") )
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector_60 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

################################################
# Some marker lists

markers=list()
markers$CM=c('TTN','MYH7','TNNT2','MYH6')
markers$EC=c('VWF','IFI27')
markers$FB=c('DCN','COL1A2')
markers$SMC=c('ACTA2','CALD1','MYH11')
    
markers_ext=list()
markers_ext$CM=c('ENSG00000155657-TTN-ProteinCoding','ENSG00000092054-MYH7-ProteinCoding','ENSG00000118194-TNNT2-ProteinCoding','ENSG00000197616-MYH6-ProteinCoding')
markers_ext$EC=c('ENSG00000110799-VWF-ProteinCoding','ENSG00000165949-IFI27-ProteinCoding')
markers_ext$FB=c('ENSG00000011465-DCN-ProteinCoding','ENSG00000164692-COL1A2-ProteinCoding')
markers_ext$SMC=c('ENSG00000107796-ACTA2-ProteinCoding','ENSG00000122786-CALD1-ProteinCoding','ENSG00000133392-MYH11-ProteinCoding')


################################################

# Let's load the data here 
data_dir1='/Volumes/fastq_m.wehrens/Mapping/HCM_SCS/mapping.93.may25/counttables/'
data_dir2='/Volumes/fastq_m.wehrens/Mapping/WANG2/counttables/'
dataset_list_paths=list('AL1'=paste0(data_dir1, 'HUB-AL-s001_HG25TBGXF_S5_cat_pT_uniaggGenes_spliced.TranscriptCounts.tsv'),
                        'JE5'=paste0(data_dir1, 'JE5_AHFL77BGX5_S6_cat_pT_uniaggGenes_spliced.TranscriptCounts.tsv'),
                        'WANG13'=paste0(data_dir2, 'GSM3449619_N13_cat_nc_uniaggGenes_spliced.TranscriptCounts.tsv'),
                        'WANG14'=paste0(data_dir2, 'GSM3449620_N14_cat_nc_uniaggGenes_spliced.TranscriptCounts.tsv'))
SCS_df_list_data = loadData_MW(dataset_list_paths, toPool = NULL)
SCS_df_pooled_data = pool_df_mw(multiple_dfs)
# We'll do scaling with Seurat ..
# manual_scale_table(countTable) 

################################################

# Let's first put these sample sets through the standard Seurat
# pipeline
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# Before starting, let's deduce some annotation
# ("WR" = Wang & v. Rooij)
WR_annotation_samples=unlist(lapply(  strsplit(colnames(SCS_df_pooled_data), split='\\.'), 
    function(splitted_string) {splitted_string[1]}))

# Initialize object
WR_object = CreateSeuratObject(SCS_df_pooled_data, project = 'Wang_and_vRooij')

# Add annotation
WR_object[["annotation_sample"]]     <- as.factor(WR_annotation_samples)

VlnPlot(WR_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'annotation_sample')

# Count mito reads
# TODO: check whether this is workign
rownames(WR_object)[grepl("-MT\\.",rownames(WR_object))]
WR_object[["percent.mt"]] <- PercentageFeatureSet(WR_object, pattern = "-MT\\.")

# Remove some cells
# Note that percent.mt < 80 is not very stringent
WR_object <- subset(WR_object, subset = nFeature_RNA > 200 & nCount_RNA > 3000 & percent.mt < 80)
# Normalize data
WR_object <- NormalizeData(WR_object, normalization.method = "LogNormalize", scale.factor = 10000) # values are defaults
# Find variable features
WR_object <- FindVariableFeatures(WR_object, selection.method = "vst", nfeatures = 2000)

# ### some plots var features
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(WR_object), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(WR_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

all.genes <- rownames(WR_object)
WR_object <- ScaleData(WR_object, features = all.genes)

WR_object <- RunPCA(WR_object, npcs = 30, verbose = FALSE)
WR_object <- RunUMAP(WR_object, reduction = "pca", dims = 1:30)
WR_object <- FindNeighbors(WR_object, reduction = "pca", dims = 1:30)
WR_object <- FindClusters(WR_object, resolution = 0.5)


# Analysis done, show some plots

print(WR_object[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(WR_object, dims = 1:5, reduction = "pca", ncol = 5)

DimPlot(WR_object, reduction = "pca")

JackStrawPlot(WR_object, dims = 1:5) # doesn't work for some reason?
ElbowPlot(WR_object)

DimPlot(WR_object, reduction = "umap")

# Some plots relevant to our scientific questions

# First sample origin
DimPlot(WR_object, group.by = 'annotation_sample', cols = col_vector_60,
            label = T, repel = T, label.size = 5)+theme(legend.position = 'none')
#DimPlot(WR_object, group.by = 'annotation_sample')

# short_gene_name = 'MALAT1'
shorthand_expand_gene_name = function(the_scs_object, short_gene_name, showcounts=NULL) {
    gene_names = 
        rownames(the_scs_object)[grepl(paste0('-',short_gene_name,'-'), rownames(the_scs_object))]
    if(length(gene_names)>1) {
        print('Multiple genes found witht this identifier, unique chosen from:')
        print(paste0(gene_names,';\n'))
        the_duplicity=(sapply(gene_names, str_count, pattern='-')+1)/3
        if (!is.null(showcounts)) {
            p=ggplot(data.frame(expr=apply(the_scs_object@assays$RNA[gene_names,],1,sum),
                              gene=gene_names,
                              duplicity=the_duplicity,
                              duplicity_f=as.factor((sapply(gene_names, str_count, pattern='-')+1)/3),
                              duplicity_xlab=paste0(
                                  sprintf("%03d",       (sapply(gene_names, str_count, pattern='-')+1)/3     ),
                                  '_(',     sprintf("%03d",     1:length(gene_names)     )    ,')'
                                  )
                                ))+
                #geom_bar(aes(x=duplicity_xlab, y=expr), stat='identity')+
                geom_boxplot(aes(x=duplicity_f, y=expr))+
                geom_jitter(aes(x=duplicity_f, y=expr))+theme_bw()+
                    theme(axis.text.x = element_text(angle = 90))
            #+coord_flip()
            print(paste0('Unique one: ',gene_names[the_duplicity==1]))
            print(p)
        } else {
            return(gene_names[the_duplicity==1])    
        }
        return(NULL)
    } else {
        return(gene_names)   
    }
}


shorthand_expand_gene_name(WR_object,'NPPB')
shorthand_expand_gene_name(WR_object,'NPPA', showcounts=T)
shorthand_expand_gene_name(WR_object,'TTN', showcounts=T)

FeaturePlot(WR_object, features = c("ENSG00000251562-MALAT1-lincRNA"), cols = rainbow_colors)
FeaturePlot(WR_object, features = c("ENSG00000155657-TTN-ProteinCoding"), cols = rainbow_colors)
FeaturePlot(WR_object, features = c("ENSG00000092054-MYH7-ProteinCoding"), cols = rainbow_colors)
FeaturePlot(WR_object, features = c("ENSG00000175206-NPPA-ProteinCoding"), cols = rainbow_colors)
FeaturePlot(WR_object, features = c("ENSG00000120937-NPPB-ProteinCoding"), cols = rainbow_colors)

FeaturePlot(WR_object, features = markers_ext$CM, cols = rainbow_colors)
FeaturePlot(WR_object, features = markers_ext$FB, cols = rainbow_colors)
FeaturePlot(WR_object, features = markers_ext$SMC, cols = rainbow_colors)
FeaturePlot(WR_object, features = markers_ext$EC, cols = rainbow_colors)

#####
# Just checking how many double-features there are
gene_names_var_all = WR_object@assays$RNA@var.features
doublefeature_df = data.frame(gene=gene_names_var_all, duplicity=(sapply(gene_names_var_all, str_count, pattern='-')+1)/3)
View(doublefeature_df)
ggplot(doublefeature_df)+
    geom_histogram(aes(x=duplicity))+theme_bw()
# counts
genes_single=doublefeature_df[doublefeature_df$duplicity==1,]$gene
genes_multi=doublefeature_df[doublefeature_df$duplicity>1,]$gene
# percentage reads from multi-named genes is only 1.5%
sum(WR_object@assays$RNA[genes_multi,])/sum(WR_object@assays$RNA[genes_single,])
#####

shorthand_expand_gene_name(WR_object,'TRDN', showcounts = T)
FeaturePlot(WR_object, features = c('ENSG00000148450-MSRB2-ProteinCoding', 
                                    'ENSG00000266472-MRPS21-ProteinCoding', 
                                    'ENSG00000140403-DNAJA4-ProteinCoding', 
                                    'ENSG00000126267-COX6B1-ProteinCoding', 
                                    'ENSG00000256618-MTRNR2L1-ProteinCoding', 
                                    'ENSG00000186439-TRDN-ProteinCoding'))





################################################
# So let's now try to integrate these datasets
# Let's follow https://satijalab.org/seurat/articles/integration_introduction.html first
# (Note: according to tutorial, this might be prone to over-correction -- 
# let's see what comes out.)

# # split the dataset into a list of two seurat objects
WR_object.list = SplitObject(WR_object, split.by = "annotation_sample")

# normalize and identify variable features for each dataset independently
WR_object.list <- lapply(X = WR_object.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
WR.list.features <- SelectIntegrationFeatures(object.list = WR_object.list)

# Now integrate 
WR_object.anchors <- FindIntegrationAnchors(object.list = WR_object.list, anchor.features = WR.list.features)

# this command creates an 'integrated' data assay
WR_object.combined <- IntegrateData(anchorset = WR_object.anchors, k.weight = 99)
# https://github.com/satijalab/seurat/issues/3784
# see my own reaction there also ..

DefaultAssay(WR_object.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
WR_object.combined <- ScaleData(WR_object.combined, verbose = FALSE)
WR_object.combined <- RunPCA(WR_object.combined, npcs = 30, verbose = FALSE)
WR_object.combined <- RunUMAP(WR_object.combined, reduction = "pca", dims = 1:30)
WR_object.combined <- FindNeighbors(WR_object.combined, reduction = "pca", dims = 1:30)
WR_object.combined <- FindClusters(WR_object.combined, resolution = 0.5)

DimPlot(WR_object.combined, reduction = "umap", group.by = "annotation_sample")+ggtitle('Samples')
DimPlot(WR_object.combined, reduction = "umap", split.by = "annotation_sample")
DimPlot(WR_object.combined, reduction = "umap", label = TRUE, repel = TRUE)+ggtitle('Clusters')

# I need to read the paper, but it seems highly unlikely that 
# all the different cell types in the Wang paper mix so nicely
# with the cells from our datasets ..
################################################################################

################################################################################
# Let's try their other method ..
# https://satijalab.org/seurat/articles/integration_rpca.html

# split the dataset into a list of two seurat objects (stim and CTRL)
WR.list2 <- SplitObject(WR_object, split.by = "annotation_sample")

# normalize and identify variable features for each dataset independently
WR.list2 <- lapply(X = WR.list2, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
WR.features2 <- SelectIntegrationFeatures(object.list = WR.list2)
WR.list2 <- lapply(X = WR.list2, FUN = function(x) {
    x <- ScaleData(x, features = WR.features2, verbose = FALSE)
    x <- RunPCA(x, features = WR.features2, verbose = FALSE)
})

WR.anchors2 <- FindIntegrationAnchors(object.list = WR.list2, anchor.features = WR.features2, reduction = "rpca")

WR.combined2 <- IntegrateData(anchorset = WR.anchors2, k.weight = 10)
WR.combined2 <- IntegrateData(anchorset = WR.anchors2, k.weight = 99)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(WR.combined2) <- "integrated"

# Run the standard workflow for visualization and clustering
WR.combined2 <- ScaleData(WR.combined2, verbose = FALSE)
WR.combined2 <- RunPCA(WR.combined2, npcs = 30, verbose = FALSE)
WR.combined2 <- RunUMAP(WR.combined2, reduction = "pca", dims = 1:30)
WR.combined2 <- FindNeighbors(WR.combined2, reduction = "pca", dims = 1:30)
WR.combined2 <- FindClusters(WR.combined2, resolution = 0.5)

# Visualization
p1 <- DimPlot(WR.combined2, reduction = "umap", group.by = "annotation_sample")+ggtitle('Samples')
print(p1)
p2 <- DimPlot(WR.combined2, reduction = "umap", label = TRUE, repel = TRUE)+ggtitle('Clusters')
print(p2)

p1 + p2

p1+FeaturePlot(WR.combined2, features = c("ENSG00000251562-MALAT1-lincRNA"), cols = rainbow_colors)
p1+FeaturePlot(WR.combined2, features = c("ENSG00000155657-TTN-ProteinCoding"), cols = rainbow_colors)
p1+FeaturePlot(WR.combined2, features = c("ENSG00000092054-MYH7-ProteinCoding"), cols = rainbow_colors)
p1+FeaturePlot(WR.combined2, features = c("ENSG00000175206-NPPA-ProteinCoding"), cols = rainbow_colors)
p1+FeaturePlot(WR.combined2, features = c("ENSG00000120937-NPPB-ProteinCoding"), cols = rainbow_colors)

################################################################################

WR_object.list3 <- SplitObject(WR_object, split.by = "annotation_sample")
WR_object.list3 <- lapply(X = WR_object.list3, FUN = SCTransform) 
WR.features3 <- SelectIntegrationFeatures(object.list = WR_object.list3, nfeatures = 3000)
WR_object.list3 <- PrepSCTIntegration(object.list = WR_object.list3, anchor.features = WR.features3)

WR.anchors3 <- FindIntegrationAnchors(object.list = WR_object.list3, normalization.method = 'SCT', anchor.features = WR.features3)
WR.combined.sct <- IntegrateData(anchorset = WR.anchors3, normalization.method = 'SCT', k.weight = 30)
WR.combined.sct <- IntegrateData(anchorset = WR.anchors3, normalization.method = 'SCT', k.weight = 60)
WR.combined.sct <- RunPCA(WR.combined.sct, verbose = FALSE)
WR.combined.sct <- RunUMAP(WR.combined.sct, reduction = "pca", dims = 1:30)
WR.combined.sct <- RunTSNE(WR.combined.sct, reduction = "pca", dims = 1:30)

# Plots
p1 <- DimPlot(WR.combined.sct, reduction = "umap", group.by = "annotation_sample")
print(p1)
#p2 <- DimPlot(WR.combined.sct, reduction = "umap", group.by = 'seurat_annotations',label = TRUE, repel = TRUE)
#p1 + p2

DimPlot(WR.combined.sct, reduction = "umap", split.by = "annotation_sample")

p1+FeaturePlot(WR.combined.sct, features = c("ENSG00000251562-MALAT1-lincRNA"), cols = rainbow_colors)

# CM markers
p1+FeaturePlot(WR.combined.sct, features = c("ENSG00000155657-TTN-ProteinCoding"), cols = rainbow_colors)
p1+FeaturePlot(WR.combined.sct, features = c("ENSG00000092054-MYH7-ProteinCoding"), cols = rainbow_colors)
p1+FeaturePlot(WR.combined.sct, features = c("ENSG00000118194-TNNT2-ProteinCoding"), cols = rainbow_colors)
p1+FeaturePlot(WR.combined.sct, features = c("ENSG00000197616-MYH6-ProteinCoding"), cols = rainbow_colors)

# EC markers
p1+FeaturePlot(WR.combined.sct, features = c("ENSG00000110799-VWF-ProteinCoding"), cols = rainbow_colors)
p1+FeaturePlot(WR.combined.sct, features = c("ENSG00000165949-IFI27-ProteinCoding"), cols = rainbow_colors)
#p1+FeaturePlot(WR.combined.sct, features = c("ENSG00000261371-PECAM1-ProteinCoding"), cols = rainbow_colors)
    # PECAM1 not found

shortshort_gname3 = function(gnames) {
    sapply(gnames, function(x) {shorthand_expand_gene_name(WR.combined.sct,x)}) }

shortshort_gname3(c('DCN','COL1A2'))

# Fibro markers
shortshort_gname3(c('DCN','COL1A2'))
FeaturePlot(WR.combined.sct, features = shortshort_gname3(c('DCN','COL1A2')), cols = rainbow_colors)

# MP markers
# (none found)
#shortshort_gname3(c('CD163','CCL4','CXCL8'))
#p1+FeaturePlot(WR.combined.sct, features = shortshort_gname3(c('CD163','CCL4','CXCL8')), cols = rainbow_colors)

# SMC
shortshort_gname3(c('ACTA2','CALD1','MYH11'))
FeaturePlot(WR.combined.sct, features = shortshort_gname3(c('ACTA2','CALD1','MYH11')), cols = rainbow_colors)

DimPlot(WR.combined.sct, reduction = "umap")
DimPlot(WR.combined.sct, reduction = "tsne")

# Using tSNE representation
p1=DimPlot(WR.combined.sct, reduction = "tsne", group.by = "annotation_sample")
print(p1)
FeaturePlot(WR.combined.sct, features = markers_ext$CM, cols = rainbow_colors, reduction = 'tsne')
FeaturePlot(WR.combined.sct, features = markers_ext$FB, cols = rainbow_colors, reduction = 'tsne')
FeaturePlot(WR.combined.sct, features = markers_ext$SMC, cols = rainbow_colors, reduction = 'tsne')
FeaturePlot(WR.combined.sct, features = markers_ext$EC, cols = rainbow_colors, reduction = 'tsne')

# Let's combine all these corrections
VlnPlot(object = WR_object, features = 'ENSG00000092054-MYH7-ProteinCoding', group.by = 'annotation_sample')+theme(legend.position='none')
VlnPlot(object = WR_object.combined, features = 'ENSG00000092054-MYH7-ProteinCoding', group.by = 'annotation_sample')+theme(legend.position='none')
VlnPlot(object = WR.combined2, features = 'ENSG00000092054-MYH7-ProteinCoding', group.by = 'annotation_sample')+theme(legend.position='none')
VlnPlot(object = WR.combined.sct, features = 'ENSG00000092054-MYH7-ProteinCoding', group.by = 'annotation_sample')+theme(legend.position='none')


# Now let's see "batch correction" for some other important markers
# Fibroblast
VlnPlot(object = WR_object, features = markers_ext$FB, group.by = 'annotation_sample', ncol=1)+
    theme(legend.position='none')+plot_annotation(title = 'No Correction')
VlnPlot(object = WR_object.combined, features = markers_ext$FB, group.by = 'annotation_sample', ncol=1)+
    theme(legend.position='none')+plot_annotation(title = 'Correction 1')
VlnPlot(object = WR.combined2, features = markers_ext$FB, group.by = 'annotation_sample', ncol=1)+
    theme(legend.position='none')+plot_annotation(title = 'Correction 2')
VlnPlot(object = WR.combined.sct, features = markers_ext$FB, group.by = 'annotation_sample', ncol=1)+
    theme(legend.position='none')+plot_annotation(title = 'Correction 3')

# Smooth muscle cells
VlnPlot(object = WR_object, features = markers_ext$SMC, group.by = 'annotation_sample', ncol=1)+
    theme(legend.position='none')+plot_annotation(title = 'No Correction')
VlnPlot(object = WR_object.combined, features = markers_ext$SMC, group.by = 'annotation_sample', ncol=1)+
    theme(legend.position='none')+plot_annotation(title = 'Correction 1')
VlnPlot(object = WR.combined2, features = markers_ext$SMC, group.by = 'annotation_sample', ncol=1)+
    theme(legend.position='none')+plot_annotation(title = 'Correction 2')
VlnPlot(object = WR.combined.sct, features = markers_ext$SMC, group.by = 'annotation_sample', ncol=1)+
    theme(legend.position='none')+plot_annotation(title = 'Correction 3')







