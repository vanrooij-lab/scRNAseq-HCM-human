
# Basically, we'll try to follow 
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
#
# We'll be loading all the data sets separately,
# and then wel'll merge the samples using 
# seurat
# https://satijalab.org/seurat/archive/v3.1/merge_vignette.html

# To determine which cells are vCM cells, use
# Wang_2021_05_looking_processed_data.R
# (See /Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Projects/)
# base_dir = '/Users/m.wehrens/Data/_2020_03_Wang/'
# base_dir = '/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/'
# load(paste0(base_dir,'Rdata/desired_cells_mwName.Rdata'))

# For the HPC, load the following file:
# base_dir = '/hpc/hub_oudenaarden/mwehrens/data/Wang/'
# load(file = paste0(base_dir,'Rdata/WRL_whole_analysis.Rdata'))

# Some notes about the Seurat object (provided by Andrew Butler, author*)
# Yes it expected that both the counts and data slot contain the raw counts immediately after converting based on the commands you ran. This way of doing things is fine.
# Yes, after normalizing in Seurat, the data slot should contain the normalized data (and the counts slot still contains the raw data).
# Yes, ScaleData works off of the normalized data (data slot). The source code for ScaleData is here
# NormalizeData always works off of the counts slots and will overwrite the data slot.
# *) source: https://github.com/satijalab/seurat/issues/2362

########################################################################
# Load some libraries, custom scripts, set some directories

library(Seurat)
library(pryr)

library(ggplot2)
library(stringr)

source('/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/Load-Pool-Scale_Simple_MW.R')

base_dir = '/hpc/hub_oudenaarden/mwehrens/data/Wang/'

data_dir1='/hpc/hub_oudenaarden/mwehrens/fastq/HCM_SCS/mapping.93.may25/counttables/'
data_dir2='/hpc/hub_oudenaarden/mwehrens/fastq/WANG2/mapping_may25/counttables/'

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

########################################################################
# Define where data is

HCM_SCS_ids = 
    c(  JE10='HUB-JE-010_HGVN3BGX9_S1_cat', 
        AL01='HUB-AL-s001_HG25TBGXF_S5_cat', 
        AL02='HUB-AL-s002_HG25TBGXF_S6_cat', 
        JE11='HUB-JE-011_HGVN3BGX9_S2_cat', 
        JE01='rJE1_AHFL7NBGX5_S3_cat', 
        JE02='JE2_AHY3WGBGX3_S1_cat', 
        JE03='JE4_AHFL7NBGX5_S4_cat', 
        JE04='JE3_AHY3WGBGX3_S2_cat', 
        JE05='JE5_AHFL77BGX5_S6_cat', 
        JE06='JE6_AHFL77BGX5_S7_cat', 
        JE07='JE7_AHFL7NBGX5_S16_cat', 
        JE08='JE8_AHFL7NBGX5_S17_cat', 
        MW05='HUB-MW-005_AH32W2BGX9_S5_cat', 
        MW06='HUB-MW-006_AH32W2BGX9_S6_cat', 
        MW07='HUB-MW-007_HC3GFBGX9_S6_cat', 
        MW08='HUB-MW-008_HC3GFBGX9_S7_cat')
dataset_list_paths_HCM = paste0(data_dir1, HCM_SCS_ids, '_pT_uniaggGenes_spliced.UFICounts.tsv')
names(dataset_list_paths_HCM) = names(HCM_SCS_ids)


WANG_ids = c(N1='GSM2970361_N1_LV_cat', 
             N3='GSM2970362_N3_LV_cat', 
             N4='GSM2970366_N4_LV_cat', 
             N5='GSM2970360_N5_LV_cat',
             N13='GSM3449619_N13_cat',
             N14='GSM3449620_N14_cat')
dataset_list_paths_WANG = paste0(data_dir2, WANG_ids, '_nc_uniaggGenes_spliced.UFICounts.tsv')
names(dataset_list_paths_WANG) = names(WANG_ids)

dataset_list_paths=c(dataset_list_paths_HCM,dataset_list_paths_WANG)

# For annotation purposes
samples_patient1 = c('JE01', 'JE02', 'JE03', 'JE04')
samples_patient2 = c('JE05', 'JE06', 'JE07', 'JE08')
samples_patient3 = c('MW05', 'MW06', 'MW07', 'MW08')
samples_patient4 = c('JE10', 'JE11')
samples_patient5 = c('AL01', 'AL02')
samples_Rooij = c(samples_patient1,samples_patient2,samples_patient3,samples_patient4,samples_patient5)

########################################################################                        
# Now actually load the data
# This takes rather long, perhaps would be nice to convert
# the loadData_MW function such that it uses lapply,
# or even a multi-core implementation of lapply

SCS_df_list_data_raw = loadData_MW(dataset_list_paths, toPool = NULL)
pryr::object_size(SCS_df_list_data_raw)
save(list = c('SCS_df_list_data_raw'),file = paste0(base_dir,'Rdata/SCS_df_list_data_raw.Rdata'))
# load(paste0(base_dir,'Rdata/SCS_df_list_data_raw.Rdata'))

SCS_df_list_data = lapply(SCS_df_list_data_raw, preprocess_convertAAnames_toSymbol)
save(list = c('SCS_df_list_data'),file = paste0(base_dir,'Rdata/SCS_df_list_data.Rdata'))

########################################################################                        
# Start the Seurat analysis

WRL_object_list = lapply(1:length(SCS_df_list_data), 
        function(idx) {
            object=CreateSeuratObject(counts = SCS_df_list_data[[idx]], project = names(SCS_df_list_data)[idx])
            print(paste0(names(SCS_df_list_data)[idx],' done .'))
            return(object)
            })
names(WRL_object_list) = names(SCS_df_list_data)
object_size(WRL_object_list) 
save(list = c('WRL_object_list'),file = paste0(base_dir,'Rdata/WRL_object_list.Rdata'))
# load(paste0(base_dir,'Rdata/WRL_object_list.Rdata'))

##########
# Now, we also load the Teichmann data
Teichmann <- LoadH5Seurat('/hpc/hub_oudenaarden/mwehrens/data/Teichmann/hca_heart_ventricular_CM_raw.h5seurat')

# However, now we also need create consistent gene names
# Seurat doesn't support renaming, so we create a new object (https://github.com/satijalab/seurat/issues/2617)
Teichmann_rawcounts = Teichmann@assays$RNA@counts
Teichmann_new_rownames = gsub(pattern = '-',replacement = '\\.',x = rownames(Teichmann_rawcounts))
rownames(Teichmann_rawcounts) = Teichmann_new_rownames
Teichmann_new = CreateSeuratObject(counts = Teichmann_rawcounts, project = 'TEICH')

##########

# Now, we can combine the data
WRL_object_big <- merge(WRL_object_list[[1]], y = unlist(WRL_object_list)[2:length(WRL_object_list)], add.cell.ids = names(WRL_object_list), project = "WRL")
WRL_object_big
object_size(WRL_object_big) 
save(list = c('WRL_object_big'),file = paste0(base_dir,'Rdata/WRL_object_big.Rdata'))
dim(WRL_object_big@assays$RNA@counts)

# Also merge Teichmann in
WRL_object_list$TEICH = Teichmann_new
WRL_object_bigbig <- merge(WRL_object_list[[1]], y = unlist(WRL_object_list)[2:length(WRL_object_list)], add.cell.ids = names(WRL_object_list), project = "WRL")
object_size(WRL_object_big) 
save(list = c('WRL_object_bigbig'),file = paste0(base_dir,'Rdata/WRL_object_bigbig.Rdata'))
# load(paste0(base_dir,'Rdata/WRL_object_bigbig.Rdata'))
dim(WRL_object_bigbig@assays$RNA@counts)
    # note that it might have been more convenient to save it as a H5Seurat
    # object also ..
    # but maybe can do that later when analysis is done ..

# Alternative way
# WRL_object_bigbig = merge(WRL_object_big, y = Teichmann_new, add.cell.ids = c("", "TEICH"), project = "WRL")

head(colnames(WRL_object_bigbig))
head(rownames(WRL_object_bigbig))

# add some annotation

WRL_annotation_samples=unlist(lapply(  strsplit(colnames(WRL_object_bigbig), split='_'), 
    function(splitted_string) {splitted_string[1]}))
WRL_object_bigbig[["annotation_sample_fct"]]     <- as.factor(WRL_annotation_samples)
WRL_object_bigbig[["annotation_sample_str"]]     <- WRL_annotation_samples

# Count mitochondrial reads
rownames(WRL_object_bigbig)[grepl("^MT\\.",rownames(WRL_object_bigbig))]
WRL_object_bigbig[["percent.mt"]] <- PercentageFeatureSet(WRL_object_bigbig, pattern = "^MT\\.")

################################################################################
# Now add vCM selection
# This is only aimed towards filtering the proper Wang cells
# (So could have been done earlier also, I guess)

# base_dir_local = '/Volumes/workdrive_m.wehrens_hubrecht/R-sessions/HCM_SCS_combining-at-HPC/'
# load(paste0(base_dir_local,'Rdata/desired_cells_mwName.Rdata'))
load(paste0(base_dir,'Rdata/desired_cells_mwName.Rdata'))

barcodes_forWang = sapply(str_split(pattern = '\\.',string = colnames(WRL_object_bigbig)), function(x) {x[4]})

individuals=WRL_object_bigbig$annotation_sample_str
# now generate cell names
cell_name_mw=paste0(individuals,'-',barcodes_forWang) # Start with Wang
cell_name_mw[individuals %in% names(HCM_SCS_ids)] = 
     colnames(WRL_object_bigbig)[individuals %in% names(HCM_SCS_ids)] # Then also add v. Rooij
cell_name_mw[individuals == 'TEICH'] = 
     colnames(WRL_object_bigbig)[individuals == 'TEICH'] # Then also add Teichmann

WRL_object_bigbig[['cell_name_mw']] = cell_name_mw 
# Now create vCM selection
vCM_sel = rep(F, dim(WRL_object_bigbig)[2])
vCM_sel[WRL_object_bigbig$cell_name_mw %in% desired_cells_mwName] = T
vCM_sel[individuals %in% c(names(HCM_SCS_ids),'TEICH')] = T
WRL_object_bigbig[['vCM']] = vCM_sel

################################################################################

# Select cells (TO DO: determine proper cutoff to also keep enough of our data)
WRL_object_bigbig_sel <- subset(WRL_object_bigbig, subset = nFeature_RNA > 50 & nCount_RNA > 1000 & percent.mt < 80 & vCM == T)
save(list = c('WRL_object_bigbig_sel'),file = paste0(base_dir,'Rdata/WRL_object_bigbig_sel.Rdata'))
# load(paste0(base_dir,'Rdata/WRL_object_bigbig_sel.Rdata'))
object_size(WRL_object_bigbig_sel)

# Some minor stats
table(WRL_object_bigbig_sel$orig.ident)
    # With current criteria, we get
    # AL01   AL02   JE01   JE02   JE03   JE04   JE05   JE06   JE07   JE08   JE10 
    #    265    285    169    145    179    140    213    257    190    322     97 
    #   JE11   MW05   MW06   MW07   MW08     N1    N13    N14     N3     N4     N5 
    #    116    254    246    221    165    259     73     28    360     68    337 
    #  TEICH 
    # 121831 

################################################################################
# First run analysis without vCM selection

# Collect some gene info
all.genes <- rownames(WRL_object_bigbig_sel)
all_mito_features = all.genes[grepl(pattern = '^MT\\.',x=all.genes)]
# determine non-mito gene set
all_genes_except_mito = all.genes[!(all.genes %in% all_mito_features)]

# Now remove them
dim(WRL_object_bigbig_sel)
WRL_object_bigbig_sel = subset(WRL_object_bigbig_sel, features = all_genes_except_mito)
dim(WRL_object_bigbig_sel)

# Normalize data
WRL_object_bigbig_sel <- NormalizeData(WRL_object_bigbig_sel, normalization.method = "LogNormalize", scale.factor = 10000) # values are defaults
# Find variable features
WRL_object_bigbig_sel <- FindVariableFeatures(WRL_object_bigbig_sel, selection.method = "vst", nfeatures = 2000)

# ### some plots var features
all_variable_features = VariableFeatures(WRL_object_bigbig_sel)

# Identify the 10 most highly variable genes
top30 <- head(VariableFeatures(WRL_object_bigbig_sel), 30)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(WRL_object_bigbig_sel)
ggsave(filename = paste0(base_dir,'Rplots/1_VariableFeatures.png'), plot = plot1)
#print(plot1)
plot2 <- LabelPoints(plot = plot1, points = top30, repel = TRUE, xnudge=0, ynudge=0)
ggsave(filename = paste0(base_dir,'Rplots/1_VariableFeatures_labeled.png'), plot = plot2, width=15, height=15)
#print(plot2)
#plot3 <- LabelPoints(plot = plot1, points = all_mito_features, repel = TRUE)
#ggsave(filename = paste0(base_dir,'Rplots/1_VariableFeatures_mitochondrial.png'), plot = plot3)

# not necessary any more
#all_variable_features_except_mito = all_variable_features[!(all_variable_features %in% all_mito_features)]
WRL_object_bigbig_sel <- ScaleData(WRL_object_bigbig_sel, features = all_variable_features) # all_variable_features is default

save(list = c('WRL_object_bigbig_sel'), file = paste0(base_dir,'Rdata/WRL_object_bigbig_sel_scaled.Rdata'))

WRL_object_bigbig_sel <- RunPCA(WRL_object_bigbig_sel, npcs = 30, verbose = FALSE)
WRL_object_bigbig_sel <- RunUMAP(WRL_object_bigbig_sel, reduction = "pca", dims = 1:30)
WRL_object_bigbig_sel <- FindNeighbors(WRL_object_bigbig_sel, reduction = "pca", dims = 1:30)
WRL_object_bigbig_sel <- FindClusters(WRL_object_bigbig_sel, resolution = 0.5)

save(list = c('WRL_object_bigbig_sel'), file = paste0(base_dir,'Rdata/WRL_object_bigbig_sel_AnaDone.Rdata'))

# Add annotation for source, only indicating Rooij, Wang or Litvi
individuals2=WRL_object_bigbig_sel$annotation_sample_str
cell_paper_source = rep(NA, length(individuals2))
cell_paper_source[individuals2 %in% names(HCM_SCS_ids)] = 'vRooij'
cell_paper_source[individuals2 %in% names(WANG_ids)] = 'Hu'
cell_paper_source[individuals2 %in% c('TEICH')] = 'Teichmann'
WRL_object_bigbig_sel[["from_paper"]] = cell_paper_source

################################################################################
# Now we can look at some plots 
# load('/Volumes/workdrive_m.wehrens_hubrecht/R-sessions/HCM_SCS_combining-at-HPC/WRL_object_big_sel_AnaDone.Rdata')

# Analysis done, show some plots

print(WRL_object_bigbig_sel[["pca"]], dims = 1:5, nfeatures = 5)

plot_pca = VizDimLoadings(WRL_object_bigbig_sel, dims = 1:5, reduction = "pca", ncol = 5)
ggsave(filename = paste0(base_dir,'Rplots/1_dimred_pca.png'), plot = plot_pca, height=10, width=15)

plot_pca2=DimPlot(WRL_object_bigbig_sel, reduction = "pca")
ggsave(filename = paste0(base_dir,'Rplots/1_dimred_pca_scatter.png'), plot = plot_pca2, height=10, width=15)

plot_pca3=DimPlot(WRL_object_bigbig_sel, reduction = "pca", group.by = 'annotation_sample_str')
ggsave(filename = paste0(base_dir,'Rplots/1_dimred_pca_scatter_annotation.png'), plot = plot_pca3, height=10, width=15)


JackStrawPlot(WRL_object_bigbig_sel, dims = 1:5) # doesn't work for some reason?
p_elbow = ElbowPlot(WRL_object_bigbig_sel)
ggsave(filename = paste0(base_dir,'Rplots/1_dimred_elbow_pca.png'), plot = p_elbow, height=7, width=15)

# Some plots relevant to our scientific questions

# First sample origin
plot_umap = DimPlot(WRL_object_bigbig_sel, group.by = 'annotation_sample_fct', #cols = col_vector_60,
                        label = T, repel = T, label.size = 5)+theme(legend.position = 'none')
ggsave(filename = paste0(base_dir,'Rplots/1_dimred_umap_annotated.png'), plot = plot_umap, height=10, width=15)
# Then also Seurat clustering
plot_umap = DimPlot(WRL_object_bigbig_sel, #cols = col_vector_60,
                        label = T, repel = T, label.size = 5)+theme(legend.position = 'none')
ggsave(filename = paste0(base_dir,'Rplots/1_dimred_umap_SeuratClustering.png'), plot = plot_umap, height=10, width=15)

# Now also get some information about the cluster of interest
#cl9_cells=colnames(WRL_object_bigbig_sel)[WRL_object_bigbig_sel@meta.data$seurat_clusters==9]
#other_cells=colnames(WRL_object_bigbig_sel)[WRL_object_bigbig_sel@meta.data$seurat_clusters!=9]
#WRL_cl9_markers <- FindMarkers(WRL_object_bigbig_sel, cells.1 = cl9_cells, cells.2 = other_cells, min.pct = 0.25)
WRL_cl9_markers <- FindMarkers(WRL_object_bigbig_sel, ident.1 = 9, min.pct = 0.25)
head(cluster2.markers, n = 25)

###
DimPlot(WRL_object_bigbig_sel, reduction = "umap")
#DimPlot(WRL_object_bigbig_sel, group.by = 'annotation_sample')

markers$vCM = c('RYR2', 'TTN', 'ATP2A', 'ACTC1') # ryr2, ttn, atp2a, actc1
FeaturePlot(WRL_object_bigbig_sel, features = markers$vCM, cols = rainbow_colors)

for (marker in c('MALAT1', 'TTN', 'MYH7', 'MYH6', 'NPPA', 'NPPB', 'ACTA1','MYL2','SORBS2','CSRP3','NDUFA4','CRYAB','HSPB1')) {
    p_cm = FeaturePlot(WRL_object_bigbig_sel, features = marker, cols = rainbow_colors)
    ggsave(filename = paste0(base_dir,'Rplots/2_umap_markers_',marker,'.png'), plot = p_cm, height=10, width=10)
}
# (Some stress and other markers)    
for (marker in c('MYH7', 'MYH6', 'NPPA', 'NPPB','ACTA1','MYL2','SORBS2','CSRP3','NDUFA4','CRYAB','HSPB1','TTN','MALAT1')) {    
    # Showing same on Violin
    pViol_m = VlnPlot(object = WRL_object_bigbig_sel, features = marker)
    ggsave(filename = paste0(base_dir,'Rplots/4_Violin_markers_perCluster_',marker,'.png'), plot = pViol_m, height=7.5, width=7.5)
}


p=FeaturePlot(WRL_object_bigbig_sel, features = markers$CM, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/2_umap_markers_CM-overview.png'), plot = p, height=10, width=10)
p=FeaturePlot(WRL_object_bigbig_sel, features = markers$FB, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/2_umap_markers_FB-overview.png'), plot = p, height=10, width=10)
p=FeaturePlot(WRL_object_bigbig_sel, features = markers$SMC, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/2_umap_markers_SMC-overview.png'), plot = p, height=10, width=10)
p=FeaturePlot(WRL_object_bigbig_sel, features = markers$EC, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/2_umap_markers_EC-overview.png'), plot = p, height=10, width=10)


p=DimPlot(WRL_object_bigbig_sel, reduction = "umap", group.by = 'vCM')+
FeaturePlot(WRL_object_bigbig_sel, features = 'RYR2', cols = rainbow_colors)+
FeaturePlot(WRL_object_bigbig_sel, features = 'TTN', cols = rainbow_colors)+
FeaturePlot(WRL_object_bigbig_sel, features = 'ACTC1', cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/2_umap_markers_vCM-custom-overview.png'), plot = p, height=15, width=15)

# Now some violin plots showing same marker stuff
pViol_m = VlnPlot(object = WRL_object_bigbig_sel, features = c('MYH6','MYH7'), group.by = 'from_paper')
ggsave(filename = paste0(base_dir,'Rplots/4_Violin_markers_MYH.png'), plot = pViol_m, height=7.5, width=7.5)

pViol_m = VlnPlot(object = WRL_object_bigbig_sel, features = c('NPPA','NPPB'), group.by = 'from_paper')
ggsave(filename = paste0(base_dir,'Rplots/4_Violin_markers_NPP.png'), plot = pViol_m, height=7.5, width=7.5)


# Show umap with annotation source by paper
p_source=DimPlot(WRL_object_bigbig_sel, group.by = 'from_paper', cols = col_vector_60,
    label = T, repel = T, label.size = 5)
ggsave(filename = paste0(base_dir,'Rplots/2_umap_bypaper.png'), plot = p_source, height=10, width=10)

################################################################################
################################################################################
# Wang vs. Rooij only
################################################################################
################################################################################

# Now compare only Wang and us, becuase I've seen some previous plots that looked much better
# Just subselect from the large thing
# 
# Note that now mitochondrial gene filtering and cell selection are already performed

# First create the object
sel_WRonly = WRL_object_bigbig_sel$annotation_sample_str %in% c(names(HCM_SCS_ids), names(WANG_ids))
WRL_object_bigbig_sel[['sel_WRonly']] = sel_WRonly

WRL_object_bigbig_WRonly = subset(WRL_object_bigbig_sel, sel_WRonly == T)

# Then perform the pipeline again
# ===

# Normalize data
WRL_object_bigbig_WRonly <- NormalizeData(WRL_object_bigbig_WRonly, normalization.method = "LogNormalize", scale.factor = 10000) # values are defaults
# Find variable features
WRL_object_bigbig_WRonly <- FindVariableFeatures(WRL_object_bigbig_WRonly, selection.method = "vst", nfeatures = 2000)

# ### some plots var features
all_variable_features = VariableFeatures(WRL_object_bigbig_WRonly)

# Identify the 10 most highly variable genes
top30 <- head(VariableFeatures(WRL_object_bigbig_WRonly), 30)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(WRL_object_bigbig_WRonly)
ggsave(filename = paste0(base_dir,'Rplots/WRonly_1_VariableFeatures.png'), plot = plot1)
#print(plot1)
plot2 <- LabelPoints(plot = plot1, points = top30, repel = TRUE, xnudge=0, ynudge=0)
ggsave(filename = paste0(base_dir,'Rplots/WRonly_1_VariableFeatures_labeled.png'), plot = plot2, width=15, height=15)
#print(plot2)
#plot3 <- LabelPoints(plot = plot1, points = all_mito_features, repel = TRUE)
#ggsave(filename = paste0(base_dir,'Rplots/WRonly_1_VariableFeatures_mitochondrial.png'), plot = plot3)

# not necessary any more
#all_variable_features_except_mito = all_variable_features[!(all_variable_features %in% all_mito_features)]
WRL_object_bigbig_WRonly <- ScaleData(WRL_object_bigbig_WRonly, features = all_variable_features) # all_variable_features is default

save(list = c('WRL_object_bigbig_WRonly'), file = paste0(base_dir,'Rdata/WRL_object_bigbig_WRonly_scaled.Rdata'))

WRL_object_bigbig_WRonly <- RunPCA(WRL_object_bigbig_WRonly, npcs = 30, verbose = FALSE)
WRL_object_bigbig_WRonly <- RunUMAP(WRL_object_bigbig_WRonly, reduction = "pca", dims = 1:30)
WRL_object_bigbig_WRonly <- FindNeighbors(WRL_object_bigbig_WRonly, reduction = "pca", dims = 1:30)
WRL_object_bigbig_WRonly <- FindClusters(WRL_object_bigbig_WRonly, resolution = 0.5)

save(list = c('WRL_object_bigbig_WRonly'), file = paste0(base_dir,'Rdata/WRL_object_bigbig_WRonly_AnaDone.Rdata'))

# Add annotation for source, only indicating Rooij, Wang or Litvi
# (This is redundant, for historic reasons)
individuals2=WRL_object_bigbig_WRonly$annotation_sample_str
cell_paper_source = rep(NA, length(individuals2))
cell_paper_source[individuals2 %in% names(HCM_SCS_ids)] = 'vRooij'
cell_paper_source[individuals2 %in% names(WANG_ids)] = 'Hu'
cell_paper_source[individuals2 %in% c('TEICH')] = 'Teichmann'
WRL_object_bigbig_WRonly[["from_paper2"]] = cell_paper_source

# Show umap with annotation source by paper
p_source=DimPlot(WRL_object_bigbig_WRonly, group.by = 'from_paper2', cols = col_vector_60,
    label = T, repel = T, label.size = 5)
ggsave(filename = paste0(base_dir,'Rplots/WRonly_2_umap_bypaper.png'), plot = p_source, height=10, width=10)


###################################
# Now again make create some plots

print(WRL_object_bigbig_WRonly[["pca"]], dims = 1:5, nfeatures = 5)

plot_pca = VizDimLoadings(WRL_object_bigbig_WRonly, dims = 1:5, reduction = "pca", ncol = 5)
ggsave(filename = paste0(base_dir,'Rplots/WRonly_1_dimred_pca.png'), plot = plot_pca, height=10, width=15)

plot_pca2=DimPlot(WRL_object_bigbig_WRonly, reduction = "pca")
ggsave(filename = paste0(base_dir,'Rplots/WRonly_1_dimred_pca_scatter.png'), plot = plot_pca2, height=10, width=15)

plot_pca3=DimPlot(WRL_object_bigbig_WRonly, reduction = "pca", group.by = 'annotation_sample_str')
ggsave(filename = paste0(base_dir,'Rplots/WRonly_1_dimred_pca_scatter_annotation.png'), plot = plot_pca3, height=10, width=15)


#JackStrawPlot(WRL_object_bigbig_WRonly, dims = 1:5) # doesn't work for some reason?
p_elbow = ElbowPlot(WRL_object_bigbig_WRonly)
ggsave(filename = paste0(base_dir,'Rplots/WRonly_1_dimred_elbow_pca.png'), plot = p_elbow, height=7, width=15)

# Some plots relevant to our scientific questions

# First sample origin
plot_umap = DimPlot(WRL_object_bigbig_WRonly, group.by = 'annotation_sample_fct', #cols = col_vector_60,
                        label = T, repel = T, label.size = 3)+theme(legend.position = 'none')
ggsave(filename = paste0(base_dir,'Rplots/WRonly_1_dimred_umap_annotated.png'), plot = plot_umap, height=20, width=15)
plot_umap2 = DimPlot(WRL_object_bigbig_WRonly, group.by = 'annotation_sample_fct', cols = col_vector_60,
    label = T, repel = T, label.size = 3)
ggsave(filename = paste0(base_dir,'Rplots/WRonly_1_dimred_umap_srclegend.png'), plot = plot_umap2, height=10, width=10)

for (marker in c('MALAT1', 'TTN', 'MYH7', 'MYH6', 'NPPA', 'NPPB')) {
    p_cm = FeaturePlot(WRL_object_bigbig_WRonly, features = marker, cols = rainbow_colors)
    ggsave(filename = paste0(base_dir,'Rplots/WRonly_2_umap_markers_',marker,'.png'), plot = p_cm, height=10, width=10)
}

p=FeaturePlot(WRL_object_bigbig_WRonly, features = markers$CM, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/WRonly_2_umap_markers_CM-overview.png'), plot = p, height=10, width=10)
p=FeaturePlot(WRL_object_bigbig_WRonly, features = markers$FB, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/WRonly_2_umap_markers_FB-overview.png'), plot = p, height=10, width=10)
p=FeaturePlot(WRL_object_bigbig_WRonly, features = markers$SMC, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/WRonly_2_umap_markers_SMC-overview.png'), plot = p, height=10, width=10)
p=FeaturePlot(WRL_object_bigbig_WRonly, features = markers$EC, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/WRonly_2_umap_markers_EC-overview.png'), plot = p, height=10, width=10)


p=DimPlot(WRL_object_bigbig_WRonly, reduction = "umap", group.by = 'vCM')+
FeaturePlot(WRL_object_bigbig_WRonly, features = 'RYR2', cols = rainbow_colors)+
FeaturePlot(WRL_object_bigbig_WRonly, features = 'TTN', cols = rainbow_colors)+
FeaturePlot(WRL_object_bigbig_WRonly, features = 'ACTC1', cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/WRonly_2_umap_markers_vCM-custom-overview.png'), plot = p, height=15, width=15)



################################################################################
################################################################################
# Rooij only
################################################################################
################################################################################

# Now compare only Wang and us, becuase I've seen some previous plots that looked much better
# Just subselect from the large thing
# 
# Note that now mitochondrial gene filtering and cell selection are already performed

# First create the object
sel_Ronly = WRL_object_bigbig_sel$annotation_sample_str %in% c(names(HCM_SCS_ids))
WRL_object_bigbig_sel[['sel_Ronly']] = sel_Ronly

WRL_object_bigbig_Ronly = subset(WRL_object_bigbig_sel, sel_Ronly == T)

# Then perform the pipeline again
# ===

# Normalize data
WRL_object_bigbig_Ronly <- NormalizeData(WRL_object_bigbig_Ronly, normalization.method = "LogNormalize", scale.factor = 10000) # values are defaults
# Find variable features
WRL_object_bigbig_Ronly <- FindVariableFeatures(WRL_object_bigbig_Ronly, selection.method = "vst", nfeatures = 2000)

    # got some singularity warnings

# ### some plots var features
all_variable_features = VariableFeatures(WRL_object_bigbig_Ronly)

# Identify the 10 most highly variable genes
top30 <- head(VariableFeatures(WRL_object_bigbig_Ronly), 30)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(WRL_object_bigbig_Ronly)
ggsave(filename = paste0(base_dir,'Rplots/Ronly_1_VariableFeatures.png'), plot = plot1)
#print(plot1)
plot2 <- LabelPoints(plot = plot1, points = top30, repel = TRUE, xnudge=0, ynudge=0)
ggsave(filename = paste0(base_dir,'Rplots/Ronly_1_VariableFeatures_labeled.png'), plot = plot2, width=15, height=15)
#print(plot2)
#plot3 <- LabelPoints(plot = plot1, points = all_mito_features, repel = TRUE)
#ggsave(filename = paste0(base_dir,'Rplots/Ronly_1_VariableFeatures_mitochondrial.png'), plot = plot3)

# not necessary any more
#all_variable_features_except_mito = all_variable_features[!(all_variable_features %in% all_mito_features)]
WRL_object_bigbig_Ronly <- ScaleData(WRL_object_bigbig_Ronly, features = all_variable_features) # all_variable_features is default

save(list = c('WRL_object_bigbig_Ronly'), file = paste0(base_dir,'Rdata/WRL_object_bigbig_Ronly_scaled.Rdata'))

WRL_object_bigbig_Ronly <- RunPCA(WRL_object_bigbig_Ronly, npcs = 30, verbose = FALSE)
WRL_object_bigbig_Ronly <- RunUMAP(WRL_object_bigbig_Ronly, reduction = "pca", dims = 1:30)
WRL_object_bigbig_Ronly <- FindNeighbors(WRL_object_bigbig_Ronly, reduction = "pca", dims = 1:30)
WRL_object_bigbig_Ronly <- FindClusters(WRL_object_bigbig_Ronly, resolution = 0.5)

save(list = c('WRL_object_bigbig_Ronly'), file = paste0(base_dir,'Rdata/WRL_object_bigbig_Ronly_AnaDone.Rdata'))

# Add annotation for source, only indicating Rooij, Wang or Litvi
# (This is redundant, for historic reasons)
individuals2=WRL_object_bigbig_Ronly$annotation_sample_str
cell_paper_source = rep(NA, length(individuals2))
cell_paper_source[individuals2 %in% names(HCM_SCS_ids)] = 'vRooij'
cell_paper_source[individuals2 %in% names(WANG_ids)] = 'Hu'
cell_paper_source[individuals2 %in% c('TEICH')] = 'Teichmann'
WRL_object_bigbig_Ronly[["from_paper2"]] = cell_paper_source

# Now also annotate our patients
individuals2=WRL_object_bigbig_Ronly$annotation_sample_str
patient_annotation = rep(NA, length(individuals2))
patient_annotation[individuals2 %in% samples_patient1] = 'P.1'
patient_annotation[individuals2 %in% samples_patient2] = 'P.2'
patient_annotation[individuals2 %in% samples_patient3] = 'P.3'
patient_annotation[individuals2 %in% samples_patient4] = 'P.4'
patient_annotation[individuals2 %in% samples_patient5] = 'P.5'
WRL_object_bigbig_Ronly[["patient_annotation"]] = patient_annotation

# Show umap with annotation source by paper
p_source=DimPlot(WRL_object_bigbig_Ronly, group.by = 'from_paper2', cols = col_vector_60,
    label = T, repel = T, label.size = 5)
ggsave(filename = paste0(base_dir,'Rplots/Ronly_2_umap_bypaper.png'), plot = p_source, height=10, width=10)

# Show umap with annotation patients
p_source=DimPlot(WRL_object_bigbig_Ronly, group.by = 'patient_annotation', cols = col_vector_60,
    label = T, repel = T, label.size = 7)
ggsave(filename = paste0(base_dir,'Rplots/Ronly_2_umap_bypatient.png'), plot = p_source, height=7, width=7)


###################################
# Now again make create some plots

print(WRL_object_bigbig_Ronly[["pca"]], dims = 1:5, nfeatures = 5)

plot_pca = VizDimLoadings(WRL_object_bigbig_Ronly, dims = 1:5, reduction = "pca", ncol = 5)
ggsave(filename = paste0(base_dir,'Rplots/Ronly_1_dimred_pca.png'), plot = plot_pca, height=10, width=15)

plot_pca2=DimPlot(WRL_object_bigbig_Ronly, reduction = "pca")
ggsave(filename = paste0(base_dir,'Rplots/Ronly_1_dimred_pca_scatter.png'), plot = plot_pca2, height=10, width=15)

plot_pca3=DimPlot(WRL_object_bigbig_Ronly, reduction = "pca", group.by = 'annotation_sample_str')
ggsave(filename = paste0(base_dir,'Rplots/Ronly_1_dimred_pca_scatter_annotation.png'), plot = plot_pca3, height=10, width=15)


#JackStrawPlot(WRL_object_bigbig_Ronly, dims = 1:5) # doesn't work for some reason?
p_elbow = ElbowPlot(WRL_object_bigbig_Ronly)
ggsave(filename = paste0(base_dir,'Rplots/Ronly_1_dimred_elbow_pca.png'), plot = p_elbow, height=7, width=15)

# Some plots relevant to our scientific questions

# First sample origin
plot_umap = DimPlot(WRL_object_bigbig_Ronly, group.by = 'annotation_sample_fct', #cols = col_vector_60,
                        label = T, repel = T, label.size = 3)+theme(legend.position = 'none')
ggsave(filename = paste0(base_dir,'Rplots/Ronly_1_dimred_umap_annotated.png'), plot = plot_umap, height=20, width=15)
plot_umap2 = DimPlot(WRL_object_bigbig_Ronly, group.by = 'annotation_sample_fct', cols = col_vector_60,
    label = T, repel = T, label.size = 3)
ggsave(filename = paste0(base_dir,'Rplots/Ronly_1_dimred_umap_srclegend.png'), plot = plot_umap2, height=10, width=10)

for (marker in c('MALAT1', 'TTN', 'MYH7', 'MYH6', 'NPPA', 'NPPB', 'ACTC1', 'RYR2', 'CRYAB')) {
    p_cm = FeaturePlot(WRL_object_bigbig_Ronly, features = marker, cols = rainbow_colors)
    ggsave(filename = paste0(base_dir,'Rplots/Ronly_2_umap_markers_',marker,'.png'), plot = p_cm, height=5, width=5)
}

p=FeaturePlot(WRL_object_bigbig_Ronly, features = markers$CM, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/Ronly_2_umap_markers_CM-overview.png'), plot = p, height=10, width=10)
p=FeaturePlot(WRL_object_bigbig_Ronly, features = markers$FB, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/Ronly_2_umap_markers_FB-overview.png'), plot = p, height=10, width=10)
p=FeaturePlot(WRL_object_bigbig_Ronly, features = markers$SMC, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/Ronly_2_umap_markers_SMC-overview.png'), plot = p, height=10, width=10)
p=FeaturePlot(WRL_object_bigbig_Ronly, features = markers$EC, cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/Ronly_2_umap_markers_EC-overview.png'), plot = p, height=10, width=10)


p=DimPlot(WRL_object_bigbig_Ronly, reduction = "umap", group.by = 'vCM')+
FeaturePlot(WRL_object_bigbig_Ronly, features = 'RYR2', cols = rainbow_colors)+
FeaturePlot(WRL_object_bigbig_Ronly, features = 'TTN', cols = rainbow_colors)+
FeaturePlot(WRL_object_bigbig_Ronly, features = 'ACTC1', cols = rainbow_colors)
ggsave(filename = paste0(base_dir,'Rplots/Ronly_2_umap_markers_vCM-custom-overview.png'), plot = p, height=15, width=15)

################################################################################

gene_in_cell_count_Ronly = 
    apply(WRL_object_bigbig_Ronly@assays$RNA@counts>0, 1, sum)

p_genes = ggplot(data.frame(nr_cells_detected_in=gene_in_cell_count_Ronly))+
    geom_freqpoly(aes(x=nr_cells_detected_in),bins=50)
ggsave(filename = paste0(base_dir,'Rplots/Ronly_0_histrogram_QC_geneInCells.png'), plot = p_genes, height=7, width=10)
    
table(gene_in_cell_count_Ronly)

sel_genes_Ronly = names(gene_in_cell_count_Ronly)[gene_in_cell_count_Ronly>4]

WRL_object_bigbig_Ronly_selG = subset(WRL_object_bigbig_Ronly, features= sel_genes_Ronly)

WRL_object_bigbig_Ronly_selG <- NormalizeData(WRL_object_bigbig_Ronly_selG, normalization.method = "LogNormalize", scale.factor = 10000) # values are defaults
WRL_object_bigbig_Ronly_selG <- FindVariableFeatures(WRL_object_bigbig_Ronly_selG, selection.method = "vst", nfeatures = 2000)
WRL_object_bigbig_Ronly_selG <- ScaleData(WRL_object_bigbig_Ronly_selG, features = sel_genes_Ronly) # this is now basically all that were selected
WRL_object_bigbig_Ronly_selG <- RunPCA(WRL_object_bigbig_Ronly_selG, npcs = 30, verbose = FALSE)
WRL_object_bigbig_Ronly_selG <- RunUMAP(WRL_object_bigbig_Ronly_selG, reduction = "pca", dims = 1:30)
WRL_object_bigbig_Ronly_selG <- FindNeighbors(WRL_object_bigbig_Ronly_selG, reduction = "pca", dims = 1:30)
WRL_object_bigbig_Ronly_selG <- FindClusters(WRL_object_bigbig_Ronly_selG, resolution = 0.5)

# Show umap with annotation patients
p_source=DimPlot(WRL_object_bigbig_Ronly_selG, group.by = 'patient_annotation', cols = col_vector_60,
    label = T, repel = T, label.size = 7)
ggsave(filename = paste0(base_dir,'Rplots/RonlySG_2_umap_bypatient.png'), plot = p_source, height=7, width=7)

for (marker in c('MALAT1', 'TTN', 'MYH7', 'MYH6', 'NPPA', 'NPPB', 'ACTC1', 'CRYAB', 'RYR2')) {
    p_cm = FeaturePlot(WRL_object_bigbig_Ronly_selG, features = marker, cols = rainbow_colors)
    ggsave(filename = paste0(base_dir,'Rplots/RonlySG_2_umap_markers_',marker,'.png'), plot = p_cm, height=5, width=5)
}

######################################################################
# repeat, but don't scale the data

# WRL_object_bigbig_Ronly_selG_NS <- FindVariableFeatures(WRL_object_bigbig_Ronly_selG_NS, selection.method = "vst", nfeatures = 2000)
WRL_object_bigbig_Ronly_selG_NS <- ScaleData(WRL_object_bigbig_Ronly_selG, features = sel_genes_Ronly, do.scale = F, do.center = F) # this is now basically all that were selected
WRL_object_bigbig_Ronly_selG_NS <- RunPCA(WRL_object_bigbig_Ronly_selG_NS, npcs = 30, verbose = FALSE)
WRL_object_bigbig_Ronly_selG_NS <- RunUMAP(WRL_object_bigbig_Ronly_selG_NS, reduction = "pca", dims = 1:30)
WRL_object_bigbig_Ronly_selG_NS <- FindNeighbors(WRL_object_bigbig_Ronly_selG_NS, reduction = "pca", dims = 1:30)
WRL_object_bigbig_Ronly_selG_NS <- FindClusters(WRL_object_bigbig_Ronly_selG_NS, resolution = 0.5)

# Show umap with annotation patients
p_source=DimPlot(WRL_object_bigbig_Ronly_selG_NS, group.by = 'patient_annotation', cols = col_vector_60,
    label = T, repel = T, label.size = 7)
ggsave(filename = paste0(base_dir,'Rplots/RonlySGNS_2_umap_bypatient.png'), plot = p_source, height=7, width=7)

# rownames_WRL_sG_NS = rownames(WRL_object_bigbig_Ronly_selG_NS@assays$RNA@scale.data)
# rownames_WRL_sG_NS[grepl('KCNQ',rownames_WRL_sG_NS)]
for (marker in c('MALAT1', 'TTN', 'MYH7', 'MYH6', 'NPPA', 'NPPB', 'ACTC1', 'CRYAB', 'RYR2')) {
    # marker = KCNQ1OT1
    p_cm = FeaturePlot(WRL_object_bigbig_Ronly_selG_NS, features = marker, cols = rainbow_colors)
    ggsave(filename = paste0(base_dir,'Rplots/RonlySGNS_2_umap_markers_',marker,'.png'), plot = p_cm, height=5, width=5)
}

######################################################################
# Now repeat this at single patient level

WRL_object_bigbig_Ronly_selG.split = SplitObject(WRL_object_bigbig_Ronly_selG, split.by = "patient_annotation")

WRL_object_bigbig_Ronly_selG.split = lapply(WRL_object_bigbig_Ronly_selG.split, 
    function(mySeuratObject) {
        mySeuratObject <- NormalizeData(mySeuratObject, normalization.method = "LogNormalize", scale.factor = 10000) # values are defaults
        mySeuratObject <- ScaleData(mySeuratObject, features = sel_genes_Ronly) # this is now basically all that were selected
        mySeuratObject <- RunPCA(mySeuratObject, npcs = 30, verbose = FALSE)
        mySeuratObject <- RunUMAP(mySeuratObject, reduction = "pca", dims = 1:30)
        mySeuratObject <- FindNeighbors(mySeuratObject, reduction = "pca", dims = 1:30)
        mySeuratObject <- FindClusters(mySeuratObject, resolution = 0.5)
    })

lapply(names(WRL_object_bigbig_Ronly_selG.split), 
    function(current_obj_name) {
        
        mySeuratObject=WRL_object_bigbig_Ronly_selG.split[[current_obj_name]]
        
        # Show umap with annotation patients
        p_source=DimPlot(mySeuratObject, cols = col_vector_60, label = T, repel = T, label.size = 7)
        ggsave(filename = paste0(base_dir,'Rplots/RonlySG_PATIENT_',current_obj_name,'_2_umap.png'), plot = p_source, height=7, width=7)
        
        for (marker in c('MALAT1', 'TTN', 'MYH7', 'MYH6', 'NPPA', 'NPPB')) {
            p_cm = FeaturePlot(mySeuratObject, features = marker, cols = rainbow_colors)+ggtitle(paste0(current_obj_name,': ', marker))
            ggsave(filename = paste0(base_dir,'Rplots/RonlySG_PATIENT_',current_obj_name,'_2_umap_markers_',marker,'.png'), plot = p_cm, height=5, width=5)
        }

})

################################################################################
# Repeat above, but now also don't use log-normalization of data
#
# Note: sel_genes_Ronly should be set properly (done above)

# To use RC scaling, we might need to re-create the Seurat object

# Now, we can combine the data
to_combine = which(names(WRL_object_list) %in% samples_Rooij)
WRL_object_Rooij <- merge(WRL_object_list[[to_combine[1]]], y = unlist(WRL_object_list[to_combine[2:length(to_combine)]]), add.cell.ids = names(WRL_object_list)[to_combine], project = "WRL")
object_size(WRL_object_Rooij) 
save(list = c('WRL_object_Rooij'),file = paste0(base_dir,'Rdata/WRL_object_Rooij.Rdata'))

# add some annotation
Rooij_annotation_samples=unlist(lapply(  strsplit(colnames(WRL_object_Rooij), split='_'), 
    function(splitted_string) {splitted_string[1]}))
WRL_object_Rooij[["annotation_sample_fct"]]     <- as.factor(Rooij_annotation_samples)
WRL_object_Rooij[["annotation_sample_str"]]     <- Rooij_annotation_samples

# Count mitochondrial reads
rownames(WRL_object_Rooij)[grepl("^MT\\.",rownames(WRL_object_Rooij))]
WRL_object_Rooij[["percent.mt"]] <- PercentageFeatureSet(WRL_object_Rooij, pattern = "^MT\\.")
mito_genes_rooij = rownames(WRL_object_Rooij)[grepl("^MT\\.",rownames(WRL_object_Rooij))]

# Remove mitochondrial genes
all_genes_except_mito_Rooij = rownames(WRL_object_Rooij)
all_genes_except_mito_Rooij = all_genes_except_mito_Rooij[!(all_genes_except_mito_Rooij %in% mito_genes_rooij)]
WRL_object_Rooij = subset(WRL_object_Rooij, features= all_genes_except_mito_Rooij)

# filter cells
WRL_object_Rooij_sel <- subset(WRL_object_Rooij, subset = nFeature_RNA > 50 & nCount_RNA > 1000 & percent.mt < 80)

# normalize
WRL_object_Rooij_sel_NSRC <- NormalizeData(WRL_object_Rooij_sel, normalization.method = "RC", scale.factor = 1e6) # RC: relative counts, 1e6 -> reads per million
# cellTotals_Rooij = apply(WRL_object_Rooij_sel@assays$RNA@counts, 2, sum); print(cellTotals_Rooij)
# #cellTotals_Rooij = apply(WRL_object_Rooij@assays$RNA@counts, 2, sum); print(cellTotals_Rooij)

# Remove few-found genes now also
gene_in_cell_count_Rooij = 
    apply(WRL_object_Rooij_sel_NSRC@assays$RNA@counts>0, 1, sum)
sel_genes_Ronly = names(gene_in_cell_count_Rooij)[gene_in_cell_count_Rooij>4]
# table(gene_in_cell_count_Rooij) # just to show 
WRL_object_Rooij_sel_NSRC = subset(WRL_object_Rooij_sel_NSRC, features= sel_genes_Ronly)

#WRL_object_Rooij_sel_NSRC <- FindVariableFeatures(WRL_object_Rooij_sel_NSRC, selection.method = "vst", nfeatures = 2000)
WRL_object_Rooij_sel_NSRC <- ScaleData(WRL_object_Rooij_sel_NSRC, features = sel_genes_Ronly, do.scale = F, do.center = F, scale.max=Inf) # this is now basically all that were selected
WRL_object_Rooij_sel_NSRC <- RunPCA(WRL_object_Rooij_sel_NSRC, npcs = 30, verbose = FALSE, features = sel_genes_Ronly)
WRL_object_Rooij_sel_NSRC <- RunUMAP(WRL_object_Rooij_sel_NSRC, reduction = "pca", dims = 1:30)
WRL_object_Rooij_sel_NSRC <- FindNeighbors(WRL_object_Rooij_sel_NSRC, reduction = "pca", dims = 1:30)
WRL_object_Rooij_sel_NSRC <- FindClusters(WRL_object_Rooij_sel_NSRC, resolution = 0.5)

individuals2=WRL_object_Rooij_sel_NSRC$annotation_sample_str
patient_annotation = rep(NA, length(individuals2))
patient_annotation[individuals2 %in% samples_patient1] = 'P.1'
patient_annotation[individuals2 %in% samples_patient2] = 'P.2'
patient_annotation[individuals2 %in% samples_patient3] = 'P.3'
patient_annotation[individuals2 %in% samples_patient4] = 'P.4'
patient_annotation[individuals2 %in% samples_patient5] = 'P.5'
WRL_object_Rooij_sel_NSRC[["patient_annotation"]] = patient_annotation

# Show umap with annotation patients
p_source=DimPlot(WRL_object_Rooij_sel_NSRC, group.by = 'patient_annotation', cols = col_vector_60,
    label = T, repel = T, label.size = 7)
ggsave(filename = paste0(base_dir,'Rplots/RonlySGNSRC_2_umap_bypatient.png'), plot = p_source, height=7, width=7)

# rownames_WRL_sG_NS = rownames(WRL_object_bigbig_Ronly_selG_NS@assays$RNA@scale.data)
# rownames_WRL_sG_NS[grepl('KCNQ',rownames_WRL_sG_NS)]
for (marker in c('MALAT1', 'TTN', 'MYH7', 'MYH6', 'NPPA', 'NPPB', 'ACTC1', 'CRYAB', 'RYR2')) {
    # marker = KCNQ1OT1
    p_cm = FeaturePlot(WRL_object_Rooij_sel_NSRC, features = marker, cols = rainbow_colors)
    ggsave(filename = paste0(base_dir,'Rplots/RonlySGNSRC_2_umap_markers_',marker,'.png'), plot = p_cm, height=5, width=5)
}

mean_gene_expression_Rooij = apply(WRL_object_Rooij_sel_NSRC@assays$RNA@scale.data,1,mean)
total_cell_expression_Rooij_sc = apply(WRL_object_Rooij_sel@assays$RNA@counts,2,sum) # should be 1e6
total_cell_expression_Rooij_sc = apply(WRL_object_Rooij_sel_NSRC@assays$RNA@data,2,sum); total_cell_expression_Rooij_sc[1:10]  # should be 1e6
total_cell_expression_Rooij_sc = apply(WRL_object_Rooij_sel_NSRC@assays$RNA@scale.data,2,sum); total_cell_expression_Rooij_sc[1:10] # should be 1e6

total_cell_expression_Rooij_sc[1:10] 

pViol_m = VlnPlot(object = WRL_object_Rooij_sel_NSRC, features = c('TTN')) #, group.by = 'from_paper')
ggsave(filename = paste0(base_dir,'Rplots/RonlySGNSRC_4_Violin_markers_TTN.png'), plot = pViol_m, height=7.5, width=7.5)


#####

# Save the whole analysis
save.image(file = paste0(base_dir,'Rdata/WRL_whole_analysis_v2.Rdata'))










