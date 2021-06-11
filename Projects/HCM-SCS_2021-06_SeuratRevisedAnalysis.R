
########################################################################

# Revised version of earlier script to analyze Van Rooij single cell data of
# HCM tissue samples from myectomies.
#
# Previously, we used RaceID2 to analyze data (see other Van Rooij repository), 
# but since we now also incorporate published data of healthy hearts, 
# this scripts repeats earlier analyses using the Seurat package.
#
# Original RaceID2 analysis by Joep Eding
# Seurat analysis including healhty hearst by m.wehrens@hubrecht.eu
# 2021-06
# 
# Outline of script:
# First I load the raw data, and create Seurat objects for each sample,
# then I define the general analysis that is performed as a function;
# this since I want to perform several versions of the analysis. 
# (The Seurat default versions emphasize different features than
# the RaceID2 analysis due to normalization choices.)
#
# I will then
# 1. Create UMAP to compare the different data sets
# (..)

########################################################################

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
# load(file = paste0(base_dir,'Rdata/RHL_whole_analysis.Rdata'))

# Some notes about the Seurat object (provided by Andrew Butler, author*)
# Yes it expected that both the counts and data slot contain the raw counts immediately after converting based on the commands you ran. This way of doing things is fine.
# Yes, after normalizing in Seurat, the data slot should contain the normalized data (and the counts slot still contains the raw data).
# Yes, ScaleData works off of the normalized data (data slot). The source code for ScaleData is here
# NormalizeData always works off of the counts slots and will overwrite the data slot.
# *) source: https://github.com/satijalab/seurat/issues/2362

########################################################################
# Load some libraries, custom scripts, set some directories

MYMCCORES = 4
script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'
base_dir = '/hpc/hub_oudenaarden/mwehrens/data/Wang/'
# script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/'

data_dir1='/hpc/hub_oudenaarden/mwehrens/fastq/HCM_SCS/mapping.93.may25/counttables/'
data_dir2='/hpc/hub_oudenaarden/mwehrens/fastq/WANG2/mapping_may25/counttables/'

library(Seurat)
library(SeuratDisk)
library(pryr)

library(ggplot2)
library(stringr)

source(paste0(script_dir,'Load-Pool-Scale_Simple_MW.R'))

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

# Now actually load the data we've mapped ourselves
# ===
# We'll be loading all data here, an we'll make selections later.
# (Note: Wang et al. was mapped again after downloading fastq files,
# according to same pipelin as our data.)

SCS_df_list_data_raw = loadData_MW_parallel(dataset_list_paths, mc.cores=MYMCCORES, prefix=F)
pryr::object_size(SCS_df_list_data_raw)
#save(list = c('SCS_df_list_data_raw'),file = paste0(base_dir,'Rdata/SCS_df_list_data_raw.Rdata'))
# load(paste0(base_dir,'Rdata/SCS_df_list_data_raw.Rdata'))

SCS_df_list_data = lapply(SCS_df_list_data_raw, preprocess_convertAAnames_toSymbol, revert_to_hgnc=T, script_dir=script_dir)
#save(list = c('SCS_df_list_data'),file = paste0(base_dir,'Rdata/SCS_df_list_data.Rdata'))

########################################################################
# Start working with Seurat objects

RHL_SeuratObject_list = mclapply(1:length(SCS_df_list_data), 
        function(idx) {
            object=CreateSeuratObject(counts = SCS_df_list_data[[idx]], project = names(SCS_df_list_data)[idx])
            print(paste0(names(SCS_df_list_data)[idx],' done .'))
            return(object)
            })
names(RHL_SeuratObject_list) = names(SCS_df_list_data)
object_size(RHL_SeuratObject_list) 
#save(list = c('RHL_SeuratObject_list'),file = paste0(base_dir,'Rdata/RHL_SeuratObject_list.Rdata'))
# load(paste0(base_dir,'Rdata/RHL_object_list.Rdata'))

########################################################################
# Now, we also load the Teichmann data (already a Seurat object)

# Note that I have previously converted the Teichmann data format to 
# the Seurat HD5 format; see script 
# /SCS_More_analyses/Projects/HCM_SCS_2021_05_HU_litvinukova_HPC.R

RHL_SeuratObject_list$TEICH <- LoadH5Seurat('/hpc/hub_oudenaarden/mwehrens/data/Teichmann/hca_heart_ventricular_CM_raw.h5seurat')
save(list = c('RHL_SeuratObject_list'),file = paste0(base_dir,'Rdata/RHL_SeuratObject_list.Rdata'))
    # Alt. way
    # Teichmann <- LoadH5Seurat('/hpc/hub_oudenaarden/mwehrens/data/Teichmann/hca_heart_ventricular_CM_raw.h5seurat')
    # RHL_SeuratObject_list$TEICH = Teichmann
    # rm('Teichmann')

# Previously, I converted Teichmann names to "Anna names" (dots instead of dashes), 
# but I now have converted all names to hgnc symbols.
#
    # # However, now we also need create consistent gene names with the other datasets
    # # Seurat doesn't support renaming, so we create a new object (https://github.com/satijalab/seurat/issues/2617)
    # Teichmann_rawcounts = Teichmann@assays$RNA@counts
    # Teichmann_new_rownames = gsub(pattern = '-',replacement = '\\.',x = rownames(Teichmann_rawcounts))
    # rownames(Teichmann_rawcounts) = Teichmann_new_rownames
    # Teichmann_new = CreateSeuratObject(counts = Teichmann_rawcounts, project = 'TEICH')
    # rm(c('Teichmann','Teichmann_rawcounts'))

########################################################################
# Merge the Seurat datasets

# Previously, I made a separate object for the Hu/Rooij but this seems unnecessary
    # # Now, we can combine the data
    # WRL_object_big <- merge(RHL_SeuratObject_list[[1]], y = unlist(RHL_SeuratObject_list)[2:length(RHL_SeuratObject_list)], add.cell.ids = names(RHL_SeuratObject_list), project = "WRL")
    # WRL_object_big
    # object_size(WRL_object_big) 
    # save(list = c('WRL_object_big'),file = paste0(base_dir,'Rdata/WRL_object_big.Rdata'))
    # dim(WRL_object_big@assays$RNA@counts)

# Merge Seurat list into one object
RHL_SeuratObject_merged <- merge(RHL_SeuratObject_list[[1]], y = unlist(RHL_SeuratObject_list)[2:length(RHL_SeuratObject_list)], add.cell.ids = names(RHL_SeuratObject_list), project = "HLR")
object_size(RHL_SeuratObject_merged)
save(list = c('RHL_SeuratObject_merged'),file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged.Rdata'))
# load(paste0(base_dir,'Rdata/RHL_SeuratObject_merged.Rdata'))
dim(RHL_SeuratObject_merged@assays$RNA@counts)
    # note that it might have been more convenient to save it as a H5Seurat
    # object also ..
    # but maybe can do that later when analysis is done ..

########################################################################
# add some annotation
#
# Note: extra metadata information for ..
# .. Wang, Hu et al can be found in metadata files, 
#   which are processed in files 
#   /SCS_More_analyses/Projects/HCM-SCS_2021_05_WANG_looking_processed_data.R
# .. Litvinukova, Teichmann et al can be found in: 
    #RHL_SeuratObject_list$TEICH@meta.data$donor
    #RHL_SeuratObject_list$TEICH@meta.data$cell_source
    #RHL_SeuratObject_list$TEICH@meta.data$sample

# Add sample information (patients in the case of Hu)
# This just extracts the prefix added by the merge function 
RHL_annotation_samples=unlist(lapply(  strsplit(colnames(RHL_SeuratObject_merged), split='_'), 
    function(splitted_string) {splitted_string[1]}))
RHL_SeuratObject_merged[["annotation_sample_fct"]]     <- as.factor(RHL_annotation_samples)
RHL_SeuratObject_merged[["annotation_sample_str"]]     <- RHL_annotation_samples

# Annotation for Teichmann?
TEICH_colnames = colnames(RHL_SeuratObject_merged)[RHL_annotation_samples=='TEICH']
    # Convenient meta-data information for Teichmann
    #RHL_SeuratObject_list$TEICH@meta.data$donor
    #RHL_SeuratObject_list$TEICH@meta.data$cell_source
    #RHL_SeuratObject_list$TEICH@meta.data$sample

# Count mitochondrial reads
rownames(WRL_object_bigbig)[grepl("^MT\\.",rownames(WRL_object_bigbig))]
WRL_object_bigbig[["percent.mt"]] <- PercentageFeatureSet(WRL_object_bigbig, pattern = "^MT\\.")

################################################################################

# To do:
# Select the vCMs


################################################################################

save.image(file = paste0(base_dir,'Rdata/RHL_whole_analysis.Rdata'))
# load(file = paste0(base_dir,'Rdata/RHL_whole_analysis.Rdata'))
