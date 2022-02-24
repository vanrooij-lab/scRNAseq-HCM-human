

# V2b! (newest version)




########################################################################

# Revised version of earlier script to analyze Van Rooij single cell data of
# HCM tissue samples from myectomies.
#
# Previously, we used RaceID2 to analyze data (see other Van Rooij repository), 
# but since we now also incorporate published data of healthy hearts, 
# this scripts repeats earlier analyses using the Seurat package.
#
# Original RaceID2 analysis by Joep Eding
# Seurat analysis including healthy hearts by m.wehrens@hubrecht.eu
# 2021-06
# 
# Outline of script:
# First I load the raw data, and create Seurat objects for each sample,
# then I define the general analysis that is performed as a function;
# this since I want to perform several versions of the analysis. 


########################################################################

# For convenience, on HPC, this script can be sourced by:
# script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')
#
# Local:
# LOCAL=1; script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Projects/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

########################################################################

# Basically, we'll try to follow 
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
#
# We'll be loading all the data sets separately,
# and then wel'll merge the samples using 
# seurat
# https://satijalab.org/seurat/archive/v3.1/merge_vignette.html

# Datasets used here:
# 
# - Our own dataset
# - Wang et al. dataset
#       Wang, L., Yu, P., Zhou, B., Song, J., Li, Z., Zhang, M., … Hu, S. (2020). Single-cell reconstruction of the adult human heart during heart failure and recovery reveals the cellular landscape underlying cardiac function. Nature Cell Biology, 22(1), 108–119. https://doi.org/10.1038/s41556-019-0446-7
#       Metadata was obtained via 
#       https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121893 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109816 
#       Raw data was downloaded using fasterq-dump and mapped with custom pipeline to make it more comparable to our own data
# - Litviňuková et al. dataset
#       Litviňuková, M., Talavera-López, C., Maatz, H., Reichart, D., Worth, C. L., Lindberg, E. L., … Teichmann, S. A. (2020). Cells of the adult human heart. Nature, (September). https://doi.org/10.1038/s41586-020-2797-4
#       Data was downloaded from https://www.heartcellatlas.org/#DataSources, where "Heart Ventricular Cardiomyocytes" where downloaded, 
#       and this file was converted to an H5Seurat file.

########################################################################

# Some small personal notes
# ==
# About the Seurat object (provided by Andrew Butler, author*):
# Yes it expected that both the counts and data slot contain the raw counts immediately after converting based on the commands you ran. This way of doing things is fine.
# Yes, after normalizing in Seurat, the data slot should contain the normalized data (and the counts slot still contains the raw data).
# Yes, ScaleData works off of the normalized data (data slot). The source code for ScaleData is here
# NormalizeData always works off of the counts slots and will overwrite the data slot.
# *) source: https://github.com/satijalab/seurat/issues/2362

########################################################################
# To be able to run this script from the command line on the cluster,
# let's make sure we can pass arguments to it
#
# Note that some sections are executed together.

print('Firing main module.')

if (!exists('myargs')) {
    myargs = commandArgs(trailingOnly = T)
    print(paste0('(main) myargs=',myargs))
} else { print(paste0('(main) myargs=',myargs)) }

library(stringr)
command_interpreter_mw = function(myargs) {
    
    if (length(myargs)==1) {
        
        print('Interpreting command line input')
        
        # split input argument in vector
        desired_command = unlist(str_split(string = myargs, pattern = '-'))
        
        # quickly parse parameters (there might be a lib for this, but ok)
        config=list()
        for (command in desired_command[grepl(pattern = '=',x = desired_command)]) {
            print(command)
            spl=str_split(string = command,pattern = '=')[[1]]
            config[[spl[1]]]=spl[2]
            print(paste0('Setting config$',spl[1],' to ' ,spl[2]))
        }
        # convert numerical strings to numeric
        config[grepl(config, pattern = '^[0-9]')]=as.numeric(config[grepl(config, pattern = '^[0-9]')])
           
        print(paste0('Sections to execute: ',paste0(desired_command,collapse = ', ')))
    } else {stop('Please pass 1 string, in the form \'arg1-arg2-arg3=x\' to give what section(s) to execute.')}

    return(list(desired_command=desired_command, config=config))
}

if (!exists('desired_command')) {
    cmd=command_interpreter_mw(myargs)
    desired_command=cmd$desired_command
    config=cmd$config
}

########################################################################
# Load some libraries, custom scripts, set some directories

# Note: always execute this library loading

if (exists('config')) {
    if ('cores' %in% names(config)) {
        MYMCCORES=config$cores
        print(paste0('Cores set to: ',MYMCCORES))
    } else { 
        MYMCCORES = 8 # not always used properly, currently
    }
} else {if (exists('LOCAL')) {MYMCCORES=1} else {MYMCCORES = 8}}


# Local  use
if (exists('LOCAL')) {
    script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/'
    base_dir = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3b/'
    
    data_dir1 = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_UmiTools_CountTables_Wang_HCMSCS/ROOIJ/' # NULL
    data_dir2 = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_UmiTools_CountTables_Wang_HCMSCS/HU/'
        # '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_Wang_Counttables/'
} else {
    script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'
    base_dir = '/hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3b/'
    
    data_dir1 = '/hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3b/counttables_tsv_exclmulti_filtered/'
        # '/hpc/hub_oudenaarden/mwehrens/fastq/HCM_SCS/mapping.93.may25/'
    data_dir2 = '/hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3b/counttables_tsv_exclmulti_filtered/'
        # '/hpc/hub_oudenaarden/mwehrens/fastq/WANG4/mapping/counttables/'
}

library(Seurat) 
    # conda install -c bioconda r-seurat
library(SeuratDisk) 
    # conda install hdf5 
        # Now make sure that Rstudio uses the python version
        # linked to this installation
    # install.packages('hdf5r')
    # install.packages("remotes")
    # remotes::install_github("mojaveazure/seurat-disk")
library(pryr)

library(ggplot2)
library(stringr)

library(patchwork)
library(scales)

library(ggrepel)

library(reshape2) 

library(pheatmap)

# also required:
# package manager: mutoss, openxlsx
# biocmanager: biomart, multtest, limma, GOstats, org.Mm.eg.db, org.Hs.eg.db

#### CONTINUE LOADING STUFF HERE 
# (install biomart next)

source(paste0(script_dir,'Functions/Load-Pool-Scale_Simple_MW.R'))
source(paste0(script_dir,'Functions/HCM-SCS_2021-06_SeuratAnalysisPipeline.R'))
source(paste0(script_dir,'Functions/MW_general_functions.R'))
source(paste0(script_dir,'Functions/MW_customized_plot_functions.R'))

# Some colors
library(RColorBrewer)
rainbow_colors = rev( brewer.pal(n = 11, name = "Spectral") )
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector_60 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors_distinguishable = c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff', '#000000')
colors_distinguishable_subset = c('#ffe119', '#4363d8', '#f58231', '#dcbeff', '#800000', '#000075', '#a9a9a9', '#ffffff', '#000000')

# Again some markers
markers=list()
markers$CM=c('TTN','MYH7','TNNT2','MYH6')
markers$EC=c('VWF','IFI27')
markers$FB=c('DCN','COL1A2')
markers$SMC=c('ACTA2','CALD1','MYH11')

PANEL_HEIGHT=47.597
PANEL_WIDTH=172/3

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
        MW05='HUB-MW-005_AH32W2BGX9_S5_cat', # MW 5-8 are all patient 3
        MW06='HUB-MW-006_AH32W2BGX9_S6_cat', 
        MW07='HUB-MW-007_HC3GFBGX9_S6_cat', 
        MW08='HUB-MW-008_HC3GFBGX9_S7_cat')
#dataset_list_paths_HCM = paste0(data_dir1, HCM_SCS_ids, '_pT_uniaggGenes_spliced.UFICounts.tsv')
dataset_list_paths_HCM = paste0(data_dir1, HCM_SCS_ids, '_pT.nonRibo_E99_Aligned.out.counts.table.tsv')
names(dataset_list_paths_HCM) = names(HCM_SCS_ids)

WANG_ids = c(   N1.97474.a = 'p.N1.plate.97474.part.1_cat',
                N1.97474.b = 'p.N1.plate.97474.part.2_cat',
    
                N2.97452.a = 'p.N2.plate.97452.part.1_cat',
                N2.97452.b = 'p.N2.plate.97452.part.2_cat',
                N2.97493.a = 'p.N2.plate.97493.part.1_cat',
                N2.97493.b = 'p.N2.plate.97493.part.2_cat',
    
                N3.97438.a = 'p.N3.plate.97438.part.1_cat',
                N3.97438.b = 'p.N3.plate.97438.part.2_cat',
    
                N4.97461.a = 'p.N4.plate.97461.part.1_cat',
                N4.97461.b = 'p.N4.plate.97461.part.2_cat',
    
                N5.97458.a = 'p.N5.plate.97458.part.1_cat',
                N5.97458.b = 'p.N5.plate.97458.part.2_cat',
    
                N13.100355.a = 'p.N13.plate.100355.part.1_cat',
                N13.100355.b = 'p.N13.plate.100355.part.2_cat',
    
                N14.104720.a = 'p.N14.plate.104720.part.1_cat',
                N14.104720.b = 'p.N14.plate.104720.part.2_cat')
# Annotation to patients
WANG_conversion_patients = sapply(str_split(WANG_ids, pattern = '\\.'), function(x) {x[[2]]})
names(WANG_conversion_patients) = names(WANG_ids)
# Annotation to samples
WANG_conversion_samples = sapply(str_split(WANG_ids, pattern = '\\.'), function(x) {x[[4]]})
names(WANG_conversion_samples) = names(WANG_ids)
# Generate full data paths
#dataset_list_paths_WANG = paste0(data_dir2, WANG_ids, '_nc_uniaggGenes_spliced.UFICounts.tsv')
dataset_list_paths_WANG = paste0(data_dir2, WANG_ids, '_nc.nonRibo_E99_Aligned.out.counts.table.tsv')
names(dataset_list_paths_WANG) = names(WANG_ids)

dataset_list_paths=c(dataset_list_paths_HCM,dataset_list_paths_WANG)

# For annotation purposes
samples_patient1 = c('JE01', 'JE02', 'JE03', 'JE04')
samples_patient2 = c('JE05', 'JE06', 'JE07', 'JE08')
samples_patient3 = c('MW05', 'MW06', 'MW07', 'MW08')
samples_patient4 = c('JE10', 'JE11')
samples_patient5 = c('AL01', 'AL02')
samples_Rooij = c(samples_patient1,samples_patient2,samples_patient3,samples_patient4,samples_patient5)
Rooij_conversion_pts = c('JE01'='R.P1', 'JE02'='R.P1', 'JE03'='R.P1', 'JE04'='R.P1',
                        'JE05'='R.P2', 'JE06'='R.P2', 'JE07'='R.P2', 'JE08'='R.P2',
                        'MW05'='R.P3', 'MW06'='R.P3', 'MW07'='R.P3', 'MW08'='R.P3',
                        'JE10'='R.P4', 'JE11'='R.P4',
                        'AL01'='R.P5', 'AL02'='R.P5')

########################################################################

if ('loadraw' %in% desired_command) {
    
    # Now actually load the data we've mapped ourselves
    # ===
    # We'll be loading all data here, an we'll make selections later.
    # (Note: Wang et al. was mapped again after downloading fastq files,
    # according to same pipelin as our data.)
    
    SCS_df_list_data_raw = loadData_MW_parallel(dataset_list_paths, mc.cores=MYMCCORES, prefix=F)
    pryr::object_size(SCS_df_list_data_raw)
    
    # Seurat can't handle underscores
    SCS_df_list_data = lapply(SCS_df_list_data_raw, 
        function(x) {rownames(x) = gsub(x = rownames(x), pattern = '_', replacement = ':'); return(x)})
    
    # Old code
    #save(list = c('SCS_df_list_data_raw'),file = paste0(base_dir,'Rdata/SCS_df_list_data_raw.Rdata'))
    # load(paste0(base_dir,'Rdata/SCS_df_list_data_raw.Rdata'))
    # This is necessary when using Anna's count tables
    # SCS_df_list_data = lapply(SCS_df_list_data_raw, preprocess_convertAAnames_toSymbol, revert_to_hgnc=T, script_dir=script_dir)
    
    if (!dir.exists(paste0(base_dir,'Rdata/SCS_df_list_data.Rdata'))) {dir.create(paste0(base_dir,'Rdata/'))} 
    save(list = c('SCS_df_list_data'),file = paste0(base_dir,'Rdata/SCS_df_list_data.Rdata'))
        # load(file = paste0(base_dir,'Rdata/SCS_df_list_data.Rdata'))

}

########################################################################
# Start working with Seurat objects

if ('generate_seurat' %in% desired_command) {
    
    load(file = paste0(base_dir,'Rdata/SCS_df_list_data.Rdata'))
    
    RHL_SeuratObject_list = mclapply(1:length(SCS_df_list_data), 
            function(idx) {
                object=CreateSeuratObject(counts = SCS_df_list_data[[idx]], project = names(SCS_df_list_data)[idx])
                print(paste0(names(SCS_df_list_data)[idx],' done .'))
                return(object)
                })
    names(RHL_SeuratObject_list) = names(SCS_df_list_data)
    object_size(RHL_SeuratObject_list) 
    
    save(list = c('RHL_SeuratObject_list'),file = paste0(base_dir,'Rdata/RHL_SeuratObject_list.Rdata'))
    # load(paste0(base_dir,'Rdata/RHL_SeuratObject_list.Rdata'))
    
}

########################################################################
# Now, we also load the Teichmann data (already a Seurat object)

# Note that I have previously converted the Teichmann data format to 
# the Seurat HD5 format; see script 
# /SCS_More_analyses/Projects/HCM_SCS_2021_05_HU_litvinukova_HPC.R

if ('add_teichmann' %in% desired_command) {
    
    load(paste0(base_dir,'Rdata/RHL_object_list.Rdata'))
    
    RHL_SeuratObject_list$TEICH <- LoadH5Seurat('/hpc/hub_oudenaarden/mwehrens/data/Teichmann/hca_heart_ventricular_CM_raw.h5seurat')
        # Alt. way
        # Teichmann <- LoadH5Seurat('/hpc/hub_oudenaarden/mwehrens/data/Teichmann/hca_heart_ventricular_CM_raw.h5seurat')
        # RHL_SeuratObject_list$TEICH = Teichmann
        # rm('Teichmann')
    
    # For some reason, the Teichmann data doesn't have the standard nCount and nFeature fields, so I'll just
    # recalculate those:
    CalcN_out = Seurat:::CalcN(RHL_SeuratObject_list$TEICH)
    RHL_SeuratObject_list$TEICH[['nCount_RNA']]   = CalcN_out$nCount
    RHL_SeuratObject_list$TEICH[['nFeature_RNA']] = CalcN_out$nFeature
    
    # Retreive ensembl names from the object
    ensembl_IDs = RHL_SeuratObject_list$TEICH@assays$RNA@meta.features$'gene_ids-Harvard-Nuclei'
    names(ensembl_IDs) = rownames(RHL_SeuratObject_list$TEICH@assays$RNA@meta.features)
        # Note that all(RHL_SeuratObject_list$TEICH@assays$RNA@meta.features$'gene_ids-Harvard-Nuclei' == RHL_SeuratObject_list$TEICH@assays$RNA@meta.features$'gene_ids-Sanger-Nuclei')
    
    ensembl_IDs_check = ensembl_IDs[rownames(RHL_SeuratObject_list$TEICH)]
        # this is redundant, since naming is consistent, i.e.
        # all(ensembl_IDs_check == ensembl_IDs)
        # still nice to do sanity check
    
    # Load our conversion table
    load(file = paste0(base_dir,'Rdata/ens_to_sym_conv_table__','93','_v2.Rdata'))  # ens_to_sym_conv_table__XX_v2
        # load(file = paste0("/hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3b/",'Rdata/ens_to_sym_conv_table__','93','_v2.Rdata'))  # ens_to_sym_conv_table__XX_v2
    
    # Additional checks & notes
    # hcgn_names  = rownames(RHL_SeuratObject_list$TEICH)
    # hcgn_names2 = names(ensembl_IDs)
    
    # hcgn_check2 = ens_to_sym_conv_table__XX_v2[ensembl_IDs]
    # all(hcgn_names[hcgn_names!=''&hcgn_check2!='']==hcgn_check2[hcgn_names!=''&hcgn_check2!=''])
    # check2_result_=hcgn_check2[hcgn_names!=hcgn_check2]
    # check2_result_[check2_result_!='']
    # 
    # Naming is not entirely consistent
    # Even though 23353 genes are correctly converted to hcgn symbol
    # sum(hcgn_names==hcgn_check2)
    # Issue is that 10166 genes do not get any annotation through the ens_to_sym_conv_table__XX_v2 table.
    # sum(hcgn_check2=="") # 10166
    # But on the other hand, some names are different in the Teichmann data vs. my conversion table, 
    # but this is just a very tiny list (TBCE, LINC01238, PRSS50, ATXN7, TXNRD3NB, CCDC39, MATR3, POLR2J3, ABCF2, PINX1, LINC01505, IGF2, HSPA14, LINC00856, EMG1, DIABLO, LINC02203, CCL3L3, H2BFS)
    # check2_result_[check2_result_!='']
    # 
    # Since we're using ens_to_sym_conv_table__XX_v2 to convert our own names, consistency can only be achieved
    # by using our own conversion table. This means that we view ENS IDs as primary in determining the genes.
    # We will store the teichmann hcgn annotation in a parameter
    # BUT we can store also the teichmann conversion to add some missing hcgn names later to the total dataset
    
    # Save the Teichmann annotation separately
    # Inversion of ensembl_IDs table
    hcgn_conversion_Teichmann = names(ensembl_IDs)
    names(hcgn_conversion_Teichmann) = ensembl_IDs
    save(list = c('hcgn_conversion_Teichmann'),file = paste0(base_dir,'Rdata/hcgn_conversion_Teichmann.Rdata'))
    
    # actually update names    
    # hgnc_names_Teich=rownames(RHL_SeuratObject_list$TEICH@assays$RNA@data)
    rownames(RHL_SeuratObject_list$TEICH@assays$RNA@data)   = paste0(ensembl_IDs_check,':',ens_to_sym_conv_table__XX_v2[ensembl_IDs_check])
    rownames(RHL_SeuratObject_list$TEICH@assays$RNA@counts) = paste0(ensembl_IDs_check,':',ens_to_sym_conv_table__XX_v2[ensembl_IDs_check])
        # can now also be queried by
        # rownames(RHL_SeuratObject_list$TEICH)
        # 
        # Potential other data containers for which to rename gene names
        # see also https://github.com/satijalab/seurat/issues/1049
        
    # Code/comments from older version script
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
    
    save(list = c('RHL_SeuratObject_list'),file = paste0(base_dir,'Rdata/RHL_SeuratObject_list.Rdata'))
    # load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_list.Rdata'))
    
}

########################################################################
# Merge the Seurat datasets

if ('create_merged_seurat' %in% desired_command) {
    
    load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_list.Rdata'))
    
    # Previously, I made a separate object for the Hu/Rooij but this seems unnecessary
        # # Now, we can combine the data
        # WRL_object_big <- merge(RHL_SeuratObject_list[[1]], y = unlist(RHL_SeuratObject_list)[2:length(RHL_SeuratObject_list)], add.cell.ids = names(RHL_SeuratObject_list), project = "WRL")
        # WRL_object_big
        # object_size(WRL_object_big) 
        # save(list = c('WRL_object_big'),file = paste0(base_dir,'Rdata/WRL_object_big.Rdata'))
        # dim(WRL_object_big@assays$RNA@counts)
    
    # Sanity check re. gene names
    # all(rownames(RHL_SeuratObject_list[[2]]) %in% rownames(RHL_SeuratObject_list$TEICH)) # TRUE
    
    # Merge Seurat list into one object
    RHL_SeuratObject_merged <- merge(RHL_SeuratObject_list[[1]], y = unlist(RHL_SeuratObject_list)[2:length(RHL_SeuratObject_list)], 
                                        add.cell.ids = names(RHL_SeuratObject_list), project = "RHL")
    object_size(RHL_SeuratObject_merged)
    
    dim(RHL_SeuratObject_merged@assays$RNA@counts)
        # note that it might have been more convenient to save it as a H5Seurat
        # object also ..
        # but maybe can do that later when analysis is done ..
    
}

########################################################################
# add some annotation
#
# Note: extra metadata information for ..
# 
# .. Wang, Hu et al can be found in metadata files, 
#   which are processed in files 
#   /SCS_More_analyses/Projects/HCM-SCS_2021_05_WANG_looking_processed_data.R
#   as far as I can see there isn't multiple samples per patient,
#   but since we're not using that information, i haven't looked into it more thoroughly.
#
# .. Litvinukova, Teichmann et al can be found in: 
    #RHL_SeuratObject_list$TEICH@meta.data$donor
    #RHL_SeuratObject_list$TEICH@meta.data$cell_source
    #RHL_SeuratObject_list$TEICH@meta.data$sample

if ('create_merged_seurat' %in% desired_command) {
    
    # This is required for e.g. LA/LV annotation in Wang sample
    # Metadata table for Wang..Hu paper
    load(file = paste0(base_dir,'Rdata/metadata_Wang_full_table_selection.Rdata'))
    
    # Add sample information (patients in the case of Hu)
    # This just extracts the prefix added by the merge function 
    RHL_annotation_prelim=unlist(lapply(  strsplit(colnames(RHL_SeuratObject_merged), split='_'), 
        function(splitted_string) {splitted_string[1]}))
    
    # Sanity check for Teichmann
    if (any(paste0('TEICH_',colnames(RHL_SeuratObject_list$TEICH))!=colnames(RHL_SeuratObject_merged)[RHL_annotation_prelim=='TEICH'])) {
        print('Annotation inconsistency.'); stop('Annotation inconsistency.')
        } else {print('Annotation Teichmann consistent.')}
    
    # Define sample annotation
    RHL_annotation_samples=RHL_annotation_prelim # note that our data is named after samples, as they map to patients
    RHL_annotation_samples[RHL_annotation_samples %in% names(WANG_ids)] = paste0('H.',WANG_conversion_samples[    RHL_annotation_samples[RHL_annotation_samples %in% names(WANG_ids)]    ])
    RHL_annotation_samples[RHL_annotation_samples=='TEICH']=paste0('T.',RHL_SeuratObject_list$TEICH@meta.data$sample)
    # unique(RHL_annotation_samples)
    # Determine annotation for from which paper a cell comes
    RHL_annotation_paper = rep(NA, length(RHL_annotation_prelim))
    RHL_annotation_paper[RHL_annotation_prelim %in% samples_Rooij]='vRooij'
    RHL_annotation_paper[RHL_annotation_prelim == 'TEICH']='Teichmann'
    RHL_annotation_paper[RHL_annotation_prelim %in% names(WANG_ids)]='Hu' 
    # Get patient information for our data and Teichmann's
    RHL_annotation_patients=RHL_annotation_prelim
    RHL_annotation_patients[RHL_annotation_patients=='TEICH']=paste0('T.',RHL_SeuratObject_list$TEICH@meta.data$donor)
    RHL_annotation_patients[RHL_annotation_patients %in% names(Rooij_conversion_pts)] = Rooij_conversion_pts[RHL_annotation_patients[RHL_annotation_patients %in% names(Rooij_conversion_pts)]]
    RHL_annotation_patients[RHL_annotation_patients %in% names(WANG_ids)] = paste0('H.',WANG_conversion_patients[    RHL_annotation_patients[RHL_annotation_patients %in% names(WANG_ids)]    ])
    # Now also retrieve the "region" data for the Teichmann data, which correspond for the vCMs to AX, SP, LV, RV (apex, septum, left ventricle, right v)
    # For our data, we only have septal cells
    RHL_annotation_region=RHL_annotation_prelim
    RHL_annotation_region[RHL_annotation_region %in% samples_Rooij] = 'SP'
    # For the Wang data, there's both LV and LA cells, but we selected for cells categorized as LV
    RHL_annotation_region[RHL_annotation_region %in% names(WANG_ids)] = 
        metadata_Wang_full_table_selection[gsub(pattern = '\\.a|\\.b',replacement='',x=colnames(RHL_SeuratObject_merged)[RHL_annotation_region %in% names(WANG_ids)]),]$condition
    RHL_annotation_region=gsub(RHL_annotation_region, pattern = 'N_',replacement = '')
    # For Teichmann, we need the region data
    RHL_annotation_region[RHL_annotation_region %in% 'TEICH'] = as.character(RHL_SeuratObject_list$TEICH@meta.data$region)
    
    # Add information to Seurat object
    RHL_SeuratObject_merged[["annotation_sample_fct"]]     <- as.factor(RHL_annotation_samples)
    RHL_SeuratObject_merged[["annotation_sample_str"]]     <- RHL_annotation_samples
    RHL_SeuratObject_merged[["annotation_patient_fct"]]     <- as.factor(RHL_annotation_patients)
    RHL_SeuratObject_merged[["annotation_patient_str"]]     <- RHL_annotation_patients
    RHL_SeuratObject_merged[["annotation_paper_fct"]]     <- as.factor(RHL_annotation_paper)
    RHL_SeuratObject_merged[["annotation_paper_str"]]     <- RHL_annotation_paper
    RHL_SeuratObject_merged[["annotation_region_fct"]]     <- as.factor(RHL_annotation_region)
    RHL_SeuratObject_merged[["annotation_region_str"]]     <- RHL_annotation_region
    
    # Some notes on Teichmann annotation
    #TEICH_colnames = colnames(RHL_SeuratObject_merged)[RHL_annotation_samples=='TEICH']
    # Convenient meta-data information for Teichmann
    #RHL_SeuratObject_list$TEICH@meta.data$donor
    #RHL_SeuratObject_list$TEICH@meta.data$cell_source
    #RHL_SeuratObject_list$TEICH@meta.data$sample
    
    # Count mitochondrial reads
    # checking whether we take the right genes:
    # rownames(RHL_SeuratObject_merged)[grepl("^MT-",rownames(RHL_SeuratObject_merged))]
    RHL_SeuratObject_merged[["percent.mt"]] <- PercentageFeatureSet(RHL_SeuratObject_merged, pattern = ":MT-")
    
    SaveH5Seurat(RHL_SeuratObject_merged,filename = paste0(base_dir,'Rdata/RHL_SeuratObject_merged.h5seurat'))
    # save(list = c('RHL_SeuratObject_merged'),file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged.Rdata'))
    # load(paste0(base_dir,'Rdata/RHL_SeuratObject_merged.Rdata'))
    
}

################################################################################

if ('create_separate_raw' %in% desired_command) {

    # load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged.Rdata'))
    RHL_SeuratObject_merged = LoadH5Seurat(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged.h5seurat'))
    
    # Create raw object for statistics
    RHL_SeuratObject_Rooij_raw = subset(RHL_SeuratObject_merged, annotation_paper_str=='vRooij')
    SaveH5Seurat(object = RHL_SeuratObject_Rooij_raw, overwrite = T,
            filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_Rooij_raw.h5seurat'))
    # Create raw object for statistics
    RHL_SeuratObject_Hu_raw = subset(RHL_SeuratObject_merged, annotation_paper_str=='Hu')
    SaveH5Seurat(object = RHL_SeuratObject_Hu_raw, overwrite = T,
            filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_Hu_raw.h5seurat'))

}

################################################################################

# Previously, I mapped all data from Hu, and it was necessary to select for 
# vCMs in that data.
# However, I now performed this selection before mapping, using their metadata, 
# and only mapped vCM cells. So this section has become obsolete.
# (Keeping here for later in commented form.)

# # Van Rooij data:   myectomy tissue, size-sorted to obtain CMs only;
# #                   we'll need some additional selection in downstream analysis
# #                   to filter out non-CMs.
# # Teichmann data:   according to their annotation, downloaded data only contains
# #                   vCMs
# # Hu data:          we'll need to take cells that finally obtained the annotation
# #                   of vCM, which are cells in the clusters LV1-LV7. So
# #                   my strategy here is to use the metadata, and select for cells
# #                   using '^LV[0-9]*$' grepl search.
# #
# # Hu data is already processed in HCM-SCS_2021_05_WANG_looking_processed_data.R,
# # where I identify cellnames that should be vCMs.
# 
# if ('select_vCM_genes_cells' %in% desired_command) {
#     
#     # So now we just need to select based on cellnames (Hu) / samples (Rooij/Teichmann)
#     # For Hu, see desired_cells_mwName determined in aforementioned script 
#     load(file = paste0(base_dir,'Rdata/desired_cells_mwName.Rdata'))
#     
#     # First need to simplify names of Hu
#     cellNames_Hu = colnames(RHL_SeuratObject_merged)[RHL_SeuratObject_merged$annotation_paper_str=='Hu']
#     Hu_CellBarcodes = sapply(str_split(pattern = '\\.',string = cellNames_Hu), function(x) {x[3]})
#     Hu_CellPatients =  sapply(str_split(pattern = '\\_',string = cellNames_Hu), function(x) {x[1]})
#     simplified_names = paste0(Hu_CellPatients,'-',Hu_CellBarcodes)
#     vCM_sel_Hu = simplified_names %in% desired_cells_mwName 
#     
#     # Now we can annotate which cells are vCM
#     vCM_sel = rep(F, dim(RHL_SeuratObject_merged)[2])
#     vCM_sel[RHL_SeuratObject_merged$annotation_paper_str == 'Hu'] = vCM_sel_Hu
#     vCM_sel[RHL_SeuratObject_merged$annotation_paper_str == 'vRooij'] = T
#     vCM_sel[RHL_SeuratObject_merged$annotation_paper_str == 'Teichmann'] = T
#     
#     RHL_SeuratObject_merged[['vCM']] = vCM_sel
#     
#     # Now save this object, as it's the basis from which we'll be working
#     save(list = c('RHL_SeuratObject_merged'),file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged.Rdata'))
#     # load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged.Rdata'))
# 
# }

################################################################################
# Now we have a data frame.
# 
# Note that for the Teichmann cells, we only obtained data from cells that were
# kept in the final analysis.
# For the Hu data, we also made a selection using their clustering, which implicitly
# also removed cells that were discarded in Hu pre-processing.
#
# For our data, no quality control is yet performed. So we'll first perform the 
# most basic level of quality control (removing cells with few reads etc).
# 
# Then, we'll run:
# > Analysis on all cells (Seurat default params)
# > Analysis on all cells (RaceID2 similar params)
# > We'll use those to investigate whether our cells are indeed disease-specific.
#   If so, we're in business, and can proceed with further analysis of 
# > our own cells separately
# > additionally, we can run separate analyses on the other datasets
#   separately to see whether our findings are specific to the diseased heart.

# First, we'll remove mitochondrial reads.

# Previously, I removed pseudo genes during post-processing.
# But now, I have filtered the GTF file beforehand (following the Cell Ranger procedure)
# So this is not necessary any more
# Keeping this as option since might be convenient later
DISCARD_GENE_OPTION=F
ENSEMBL_VERSION_BIOTYPES=93

if ('select_genes_cells' %in% desired_command) {
    
    # load(paste0(base_dir,'Rdata/RHL_SeuratObject_merged.Rdata'))
    RHL_SeuratObject_merged = LoadH5Seurat(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged.h5seurat'))
    
    # Let's remove mitochondrial and other undesired gene expression information
    
    # First, let's annotate gene biotypes
    # NOTE: strictly this code is not necessary any more for other stuff to function,
    # but it might still be nice as a sanity check to see that biotypes check out
    load(file = paste0(base_dir,'Rdata/gene_biotypes_',ENSEMBL_VERSION_BIOTYPES,'.Rdata')) 
    current_genes_ens = sapply(str_split(string = rownames(RHL_SeuratObject_merged), pattern = ':'), function(x) {x[[1]]})
    RHL_SeuratObject_merged@misc[[paste0('biotypes',ENSEMBL_VERSION_BIOTYPES)]] = gene_biotypes_XX[current_genes_ens]
        # Sanity checks
        # unique(RHL_SeuratObject_merged@misc[[paste0('biotypes',ENSEMBL_VERSION_BIOTYPES)]])
        # --> checks out
        # Also, in Teichmann, biotypes check out
        # genes_TEICH_ens = sapply(str_split(string = rownames(RHL_SeuratObject_list$TEICH), pattern = ':'), function(x) {x[[1]]})
        # gene_types_TEICH_ens = gene_biotypes_XX[genes_TEICH_ens]
        # unique(gene_types_TEICH_ens); table(gene_types_TEICH_ens)
        # For our data (testing one set), also checks out
        # genes_JE10_ens = sapply(str_split(string = rownames(RHL_SeuratObject_list$JE10), pattern = ':'), function(x) {x[[1]]})
        # gene_types_JE10_ens = gene_biotypes_XX[genes_JE10_ens]
        # unique(gene_types_JE10_ens); table(gene_types_JE10_ens)
    
        # XXXXX
    
        # Note that Rooij/Hu data only has protein_coding, lincRNA and antisense, while Teichmann has slightly more (due cellranger)
            #     gene_types_TEICH_ens
            #       antisense       IG_C_gene IG_C_pseudogene       IG_D_gene       IG_J_gene 
            #            5497              14               9              37              18 
            # IG_J_pseudogene       IG_V_gene IG_V_pseudogene         lincRNA  protein_coding 
            #               3             144             188            7484           19912 
            #       TR_C_gene       TR_D_gene       TR_J_gene TR_J_pseudogene       TR_V_gene 
            #               6               4              79               4             106 
            # TR_V_pseudogene 
            #              33 
    
    if (DISCARD_GENE_OPTION) {
        # Cell Ranger biotypes
        # Note that the biomart annotation is slightly different than what's listed on Cell Ranger website, below
        # is consistent with biomart; difference: should be lncRNA, not lincRNA
        desired_biotypes_cellranger = c('protein_coding', 'lncRNA', 'IG_LV_gene', 'IG_V_gene', 'IG_V_pseudogene', 'IG_D_gene', 'IG_J_gene', 'IG_J_pseudogene', 'IG_C_gene', 'IG_C_pseudogene', 'TR_V_gene', 'TR_V_pseudogene', 'TR_D_gene', 'TR_J_gene', 'TR_J_pseudogene', 'TR_C_gene')
        # desired_biotypes_cellranger = c('protein_coding', 'lincRNA', 'antisense', 'IG_LV_gene', 'IG_V_gene', 'IG_V_pseudogene', 'IG_D_gene', 'IG_J_gene', 'IG_J_pseudogene', 'IG_C_gene', 'IG_C_pseudogene', 'TR_V_gene', 'TR_V_pseudogene', 'TR_D_gene', 'TR_J_gene', 'TR_J_pseudogene', 'TR_C_gene')
        
        # Determine subsets of desired vs. undesired gene features
        RHL_SeuratObject_merged@misc$desired_genes = 
            rownames(RHL_SeuratObject_merged)[RHL_SeuratObject_merged@misc$biotypes %in% desired_biotypes_cellranger]
        RHL_SeuratObject_merged@misc$discarded_genes = 
            rownames(RHL_SeuratObject_merged)[!(RHL_SeuratObject_merged@misc$biotypes %in% desired_biotypes_cellranger)]
    
        # Then, let's (1) store non-desired data into different assay and (2) filter out non-desired data
        RHL_SeuratObject_merged[['rnaDiscardedCounts']] <- CreateAssayObject(counts = RHL_SeuratObject_merged@assays$RNA@counts[RHL_SeuratObject_merged@misc$discarded_genes,])
    }
    
    # Now let's also turn to mitochondrial reads
    mito_genes = rownames(RHL_SeuratObject_merged)[grepl(':MT-',rownames(RHL_SeuratObject_merged))]
    # Again save the non-desired data
    RHL_SeuratObject_merged[['rnaMitoCounts']] <- CreateAssayObject(counts = RHL_SeuratObject_merged@assays$RNA@counts[mito_genes,])
    # Create updated selection
    if (DISCARD_GENE_OPTION) {
        RHL_SeuratObject_merged@misc$desired_genes_excl_mito = 
            RHL_SeuratObject_merged@misc$desired_genes[!(RHL_SeuratObject_merged@misc$desired_genes %in% mito_genes)]
            # RHL_SeuratObject_merged@misc$desired_genes[(RHL_SeuratObject_merged@misc$desired_genes %in% mito_genes)]
    } else {
        RHL_SeuratObject_merged@misc$desired_genes_excl_mito = 
            rownames(RHL_SeuratObject_merged)[!(rownames(RHL_SeuratObject_merged) %in% mito_genes)]
    }
    
    # Then, subset the RNA assay 
    RHL_SeuratObject_merged_noMito =
        subset(RHL_SeuratObject_merged, features= RHL_SeuratObject_merged@misc$desired_genes_excl_mito)
    # And add the desired stored information again
    RHL_SeuratObject_merged_noMito[['rnaMitoCounts']] = RHL_SeuratObject_merged[['rnaMitoCounts']]
        # note: to look up; don't quite understand why subsetting can't be applied to all assays at once in Seurat
    if (DISCARD_GENE_OPTION) {
        RHL_SeuratObject_merged_noMito[['rnaDiscardedCounts']] = RHL_SeuratObject_merged[['rnaDiscardedCounts']]
    }
        # dim(RHL_SeuratObject_merged_noMito@assays$rnaMitoCounts@counts)
        # rownames(RHL_SeuratObject_merged_noMito@assays$rnaMitoCounts@counts)
        # colnames(RHL_SeuratObject_merged_noMito@assays$rnaMitoCounts@counts)[1:10]
        # dim(RHL_SeuratObject_merged_noMito@assays$rnaDiscardedCounts@counts)
        # rownames(RHL_SeuratObject_merged_noMito@assays$rnaDiscardedCounts@counts)[1:30]
        # colnames(RHL_SeuratObject_merged_noMito@assays$rnaDiscardedCounts@counts)[1:10]
    # Sanity check
    # DefaultAssay(object = RHL_SeuratObject_merged_noMito)
    #Key(object = RHL_SeuratObject_merged_noMito[["RNA"]])
    #Key(object = RHL_SeuratObject_merged_noMito[["rnaMitoCounts"]])
    #Key(object = RHL_SeuratObject_merged_noMito[["rnaDiscardedCounts"]])
    # all keys are same names but lowercase followed by underscore
    
    # Earlier code; I only removed selection of mitochondrial pseudogenes -- 
    # turns out Cell Ranger removes ALL pseudogenes, so it's more convenient that we also do that.
    #
    # # Now do the same for mitochondrial pseudogenes
    # mito_genes_sym = sapply(str_split(string = mito_genes, pattern = ':'), function(x) {x[[2]]})
    # pseudo_patterns = paste0( ':', gsub(x=mito_genes_sym, pattern = '-', replacement = '') , 'P')
    # mito_pseudogenes = rownames(RHL_SeuratObject_merged)[grepl(paste0(pseudo_patterns, collapse = '|'),rownames(RHL_SeuratObject_merged))]
    # # 
    # RHL_SeuratObject_merged@misc$mito.pseudocounts = RHL_SeuratObject_merged@assays$RNA@counts[mito_pseudogenes,]
    # rownames(RHL_SeuratObject_merged@misc$mito.pseudocounts) = mito_pseudogenes
    # colnames(RHL_SeuratObject_merged@misc$mito.pseudocounts) = colnames(RHL_SeuratObject_merged@assays$RNA@counts)
    # RHL_SeuratObject_merged@misc$mito.pseudocounts_colnames= colnames(RHL_SeuratObject_merged@assays$RNA@counts)
    # RHL_SeuratObject_merged@misc$mito.pseudocounts_rownames= mito_pseudogenes
    
    # Then determine the total well and detected feature counts again
    CalcN_out_merged = Seurat:::CalcN(RHL_SeuratObject_merged_noMito)
    RHL_SeuratObject_merged_noMito[['nCount.nMT']]   = CalcN_out_merged$nCount
    RHL_SeuratObject_merged_noMito[['nFeature.nMT']] = CalcN_out_merged$nFeature
    
    # Then perform cell selection
    # RHL_SeuratObject_merged_noMito_sel <- subset(RHL_SeuratObject_merged_noMito, subset = nFeature.nMT > 50 & nFeature.nMT > 1000 & percent.mt < 80)
    #
    # BUGFIX!! --> previously, I accidentically selected for nFeature.nMT > 1000
    RHL_SeuratObject_merged_noMito_sel <- subset(RHL_SeuratObject_merged_noMito, subset = nCount.nMT > 1000)
    
    # Little boilerplate, additional selection of Hu cells based on region
    # unique(RHL_SeuratObject_merged_noMito_sel[["annotation_region_str"]][RHL_SeuratObject_merged_noMito_sel[["annotation_paper_str"]]=='Hu'])
    RHL_SeuratObject_merged_noMito_sel = RHL_SeuratObject_merged_noMito_sel[,!(RHL_SeuratObject_merged_noMito_sel[["annotation_region_str"]]=='LA'&
               RHL_SeuratObject_merged_noMito_sel[["annotation_paper_str"]]=='Hu')]
    
    # Now also apply this to the stored "discarded" data
    RHL_SeuratObject_merged_noMito_sel[['rnaMitoCounts']] <- subset(RHL_SeuratObject_merged_noMito[['rnaMitoCounts']], cells = colnames(RHL_SeuratObject_merged_noMito_sel))
    if (DISCARD_GENE_OPTION) {
        RHL_SeuratObject_merged_noMito_sel[['rnaDiscardedCounts']] <- subset(RHL_SeuratObject_merged_noMito[['rnaDiscardedCounts']], cells = colnames(RHL_SeuratObject_merged_noMito_sel))
    }
    
    # table(RHL_SeuratObject_merged_sel$orig.ident) XXX
    cell_count_overview_tabl_paper    = as.data.frame(table(RHL_SeuratObject_merged_noMito_sel$annotation_paper_str))
    colnames(cell_count_overview_tabl_paper) = c('Source','Cell_count')
    cell_count_overview_table_patient = as.data.frame(table(RHL_SeuratObject_merged_noMito_sel$annotation_patient_str))
    colnames(cell_count_overview_table_patient) = c('Donor','Cell_count')
    openxlsx::write.xlsx(x = list(cellCounts_source=cell_count_overview_tabl_paper, 
                                  cellCounts_patients=cell_count_overview_table_patient) ,
                         file = paste0(base_dir,'Rplots/QC_general_statistics_cellcounts_alldata.xlsx'), overwrite = T)
    
    # Now also update names that have lacking information
    load(file = paste0(base_dir,'Rdata/hcgn_conversion_Teichmann.Rdata')) # hcgn_conversion_Teichmann
    # determine names lacking info
    full_names1 = rownames(RHL_SeuratObject_merged_noMito_sel@assays$RNA@counts)
    ens_names = shorthand_cutname( rownames(RHL_SeuratObject_merged_noMito_sel@assays$RNA@counts) , PART1OR2 = 1)
    hcgn_names = sapply(str_split(full_names1, pattern = ':'), function(Y) {Y[[2]]})
    
        # ens_names_data = shorthand_cutname( rownames(RHL_SeuratObject_merged_noMito_sel@assays$RNA@data) , PART1OR2 = 1)
        # all(ens_names==ens_names_data) # TRUE
        
    # Now create updated names
    newnames = hcgn_conversion_Teichmann[ens_names[hcgn_names=='']]
    hcgn_names_updated = hcgn_names
    hcgn_names_updated[hcgn_names==''] = newnames
    full_newnames = paste0(ens_names, ':', hcgn_names_updated)
        # check_df = data.frame(hcgn_names, hcgn_names_updated, full_newnames)
        # check_df[1:100,]
        # check_df[check_df$hcgn_names=='',][1:100,]
        
    # Actually update names
    rownames(RHL_SeuratObject_merged_noMito_sel@assays$RNA@counts) = full_newnames
    rownames(RHL_SeuratObject_merged_noMito_sel@assays$RNA@data) = full_newnames
    # rownames(RHL_SeuratObject_merged_noMito_sel@assays$RNA@scale.data) = full_newnames
        
    # Save it
    SaveH5Seurat(RHL_SeuratObject_merged_noMito_sel, filename = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_noMito_sel.Rdata'))
    # save(list = c('RHL_SeuratObject_merged_noMito_sel'), file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_noMito_sel.Rdata'))
    # load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_noMito_sel.Rdata'))
    
    # Remove old stuff
    rm('RHL_SeuratObject_merged')
    rm('RHL_SeuratObject_merged_noMito')
    
}
    
################################################################################
# Now we can start using the Seurat analysis pipeline as defined in
# HCM-SCS_2021-06_SeuratAnalysisPipeline.R

# Objtect that can be used as input for above function (see also below):

# load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_noMito_sel.Rdata'))

################################################################################
# Run the analysis using default Seurat settings

if ('runf_all_default' %in% desired_command) {
    CURRENT_RUNNAME='merged-s-def'
    
    load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_noMito_sel.Rdata'))
    
    # Run analysis
    if (!('plotsonly' %in% desired_command)) {
        RHL_SeuratObject_merged_noMito_sel.default = mySeuratAnalysis(RHL_SeuratObject_merged_noMito_sel, run_name = CURRENT_RUNNAME)
    }
    # Create plots
    mySeuratCommonPlots(RHL_SeuratObject_merged_noMito_sel.default, run_name = CURRENT_RUNNAME)
    # save(list = c('RHL_SeuratObject_merged_noMito_sel.default'), file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_noMito_sel.default.Rdata'))
    # load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_noMito_sel.default.Rdata'))
    
    SaveH5Seurat(object = RHL_SeuratObject_merged_noMito_sel.default, overwrite = T,
        filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_merged_noMito_sel.default.h5seurat'))
    #RHL_SeuratObject_merged_noMito_sel.default = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_merged_noMito_sel.default.h5seurat'))

}    

# desired_command='runf_all_default_cl'
if ('runf_all_default_cl' %in% desired_command) {
    CURRENT_RUNNAME='merged-s-def'
    
    H5_RHL_SeuratObject_merged_noMito_sel.default = 
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_merged_noMito_sel.default.h5seurat'))
    
    DE_cluster=list()
    DE_cluster[[CURRENT_RUNNAME]] = diff_express_clusters(H5_RHL_SeuratObject_merged_noMito_sel.default, mc.cores = MYMCCORES)
    save(list = c('DE_cluster'), file = paste0(base_dir,'Rdata/DE_cluster_',CURRENT_RUNNAME,'.Rdata'))
    
    diff_express_clusters_save_results(all_markers = DE_cluster[[CURRENT_RUNNAME]], run_name = CURRENT_RUNNAME,
                base_dir=base_dir, topX = 30)

    # all_default_clusters = diff_express_clusters(H5_RHL_SeuratObject_merged_noMito_sel.default, mc.cores = 8)
    # save(list = c('all_default_clusters'), file = paste0(base_dir,'Rdata/DE_cluster_',CURRENT_RUNNAME,'.Rdata'))
    
    
}
    
################################################################################
# Run the analysis for ALL sets, with RACEID2 like params

# NOTE THAT THIS WAS NEVER EXECUTED, BECAUSE FEATURES ALL IS VERY INTENSIVE AND NOT NECESSARY
# if ('runf_all_RID2l' %in% desired_command) {
#     CURRENT_RUNNAME='all_RID2l'
#     
#     load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_noMito_sel.Rdata'))
#     
#     RHL_SeuratObject_merged_noMito_sel.rid2like = mySeuratAnalysis(RHL_SeuratObject_merged_noMito_sel,
#         run_name=CURRENT_RUNNAME,
#         normalization.method='RC', scale.factor=median(RHL_SeuratObject_merged_noMito_sel$nCount.nMT),
#         do.scale=F,do.center=F,scale.max=Inf, features_to_use_choice = 'all')
#     mySeuratCommonPlots(RHL_SeuratObject_merged_noMito_sel.rid2like, run_name = CURRENT_RUNNAME)
#     # save(list = c('RHL_SeuratObject_merged_noMito_sel.rid2like'), file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_noMito_sel.rid2like.Rdata'))
#     # load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_noMito_sel.rid2like.Rdata'))
#     
#     SaveH5Seurat(object = RHL_SeuratObject_merged_noMito_sel.rid2like, overwrite = T,
#         filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_merged_noMito_sel.rid2like.h5seurat'))
#     #RHL_SeuratObject_merged_noMito_sel.rid2like = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_merged_noMito_sel.rid2like.h5seurat'))
#     
# }

################################################################################
# RID2-like but with var. features activated

if ('runf_all_RID2l_VAR' %in% desired_command) {
    CURRENT_RUNNAME='all_RID2l_VAR'
    
    load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_noMito_sel.Rdata'))
    
    # Run and make some standard plots
    RHL_SeuratObject_merged_noMito_sel.rid2like_VAR = mySeuratAnalysis(RHL_SeuratObject_merged_noMito_sel,
        run_name=CURRENT_RUNNAME,
        normalization.method='RC', scale.factor=median(RHL_SeuratObject_merged_noMito_sel$nCount.nMT),
        do.scale=F,do.center=F,scale.max=Inf, features_to_use_choice = 'variable')
    mySeuratCommonPlots(RHL_SeuratObject_merged_noMito_sel.rid2like_VAR, run_name = CURRENT_RUNNAME)
    
    SaveH5Seurat(object = RHL_SeuratObject_merged_noMito_sel.rid2like_VAR, overwrite = T,
        filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_merged_noMito_sel.',CURRENT_RUNNAME,'.h5seurat'))
    #RHL_SeuratObject_merged_noMito_sel.rid2like_VAR = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_merged_noMito_sel.',CURRENT_RUNNAME,'.h5seurat'))
    
}

# desired_command='runf_all_RID2l_VAR_cl'
if ('runf_all_RID2l_VAR_cl' %in% desired_command) {
    CURRENT_RUNNAME='all_RID2l_VAR'
    
    RHL_SeuratObject_merged_noMito_sel.rid2like_VAR = 
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_merged_noMito_sel.',CURRENT_RUNNAME,'.h5seurat'))
    
    # perhaps change var names?
    DE_cluster=list()
    DE_cluster[[CURRENT_RUNNAME]] = diff_express_clusters(RHL_SeuratObject_merged_noMito_sel.rid2like_VAR, mc.cores = MYMCCORES)
    save(list = c('DE_cluster'), file = paste0(base_dir,'Rdata/DE_cluster_',CURRENT_RUNNAME,'.Rdata'))
    #load(file = paste0(base_dir,'Rdata/DE_cluster_',CURRENT_RUNNAME,'.Rdata'))
    
    diff_express_clusters_save_results(all_markers = DE_cluster[[CURRENT_RUNNAME]], run_name = CURRENT_RUNNAME,
            base_dir=base_dir, topX = 30)

    print('Cluster analysis done')
    
}

################################################################################
# To facilitate a full comparison between Rooij and Teichmann, let's also 
# filter the pooled data sets for septal cells.
#
# Naming and procesisng after this selection filter will be consistent with
# the separate datasets I created.

if ('create_septal_all_dataset' %in% desired_command) {
    
    # Load all the data
    print('Loading')
    load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_noMito_sel.Rdata'))
    
    # Create "all" pooled set filtered for septal
    print('Filtering Teichmann for septal cells ..')
    # Note: RHL_SeuratObject_merged_noMito_sel_selSP would also have been a logical name, but 
    # I'll save it as H5Seurat below for further processing, so will save it with name
    # that's consistent with separate datasets. 
    RHL_SeuratObject_nM_sel_ALL.SP =
        subset(RHL_SeuratObject_merged_noMito_sel, annotation_region_str == 'SP' | annotation_paper_str == 'Hu')
    
    RHL_SeuratObject_nM_sel_ALL.SP[['rnaMitoCounts']] <- subset(RHL_SeuratObject_merged_noMito_sel[['rnaMitoCounts']], cells = colnames(RHL_SeuratObject_nM_sel_ALL.SP))
    if (DISCARD_GENE_OPTION) { RHL_SeuratObject_nM_sel_ALL.SP[['rnaDiscardedCounts']] <- subset(RHL_SeuratObject_merged_noMito_sel[['rnaDiscardedCounts']], cells = colnames(RHL_SeuratObject_nM_sel_ALL.SP)) }
    
    print('Saving RHL_SeuratObject_nM_sel_ALL.SP ..')
    
    # Now let's save this one as H5Seurat
    SaveH5Seurat(object = RHL_SeuratObject_nM_sel_ALL.SP, overwrite = T,
        filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_ALL.SP.h5seurat'))
    
    # NOTE: FIX: COMPROMISED H5_RHL_SeuratObject_nM_sel_HUonly.h5seurat
    
}

# Just for reference, what are the biotypes here?
if ('create_septal_all_btypSel_dataset' %in% desired_command) {
    
    seuratObject_name = paste0('RHL_SeuratObject_nM_sel_','ALL.SP')
    
    # load the file
    current_analysis = list()
    current_analysis[[seuratObject_name]] =
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_',seuratObject_name,'.h5seurat'))
    
    ENSEMBL_VERSION_BIOTYPES=93
    
    load(file = paste0(base_dir,'Rdata/gene_biotypes_',ENSEMBL_VERSION_BIOTYPES,'.Rdata'))
    
    gnames_spall = rownames(current_analysis[[seuratObject_name]])
    genes_ALL.SP_ens = sapply(str_split(string = gnames_spall, pattern = ':'), function(x) {x[[1]]})
    genes_types_ALL.SP_ens = gene_biotypes_XX[genes_ALL.SP_ens]
    IC_TR_X_genes_idx = !(genes_types_ALL.SP_ens %in% c('antisense','protein_coding','lincRNA'))
    gnames_spall[IC_TR_X_genes_idx]
    expressed_genes = rownames(current_analysis[[seuratObject_name]]@assays$RNA@counts)[rowSums(current_analysis[[seuratObject_name]]@assays$RNA@counts)>0]
    expressed_gene_biotypes = gene_biotypes_XX[sapply(str_split(string = expressed_genes, pattern = ':'), function(x) {x[[1]]})]
    expresesd_IC_TR_X_genes_idx = !(expressed_gene_biotypes %in% c('antisense','protein_coding','lincRNA'))
    expressed_genes[expresesd_IC_TR_X_genes_idx]
    
    # Pay attention to these genes, as they are now still in our analysis, and
    # might turn up as DE, but that 'd be artificial
    # expressed_genes[expresesd_IC_TR_X_genes_idx]
    #  [1] "ENSG00000233999:IGKV3OR2-268" "ENSG00000211592:IGKC"        
    #  [3] "ENSG00000239819:IGKV1D-8"     "ENSG00000228325:IGKV3D-7"    
    #  [5] "ENSG00000231292:IGKV1OR2-108" "ENSG00000227191:TRGC2"       
    #  [7] "ENSG00000211689:TRGC1"        "ENSG00000211695:TRGV9"       
    #  [9] "ENSG00000275743:TRBV14"       "ENSG00000211751:TRBC1"       
    # [11] "ENSG00000211772:TRBC2"        "ENSG00000255569:TRAV1-1"     
    # [13] "ENSG00000211799:TRAV19"       "ENSG00000211829:TRDC"        
    # [15] "ENSG00000277734:TRAC"         "ENSG00000211893:IGHG2"       
    # [17] "ENSG00000253755:IGHGP"        "ENSG00000211895:IGHA1"       
    # [19] "ENSG00000211896:IGHG1"        "ENSG00000211897:IGHG3"       
    # [21] "ENSG00000211898:IGHD"         "ENSG00000211899:IGHM"        
    # [23] "ENSG00000254174:IGHV1-12"     "ENSG00000211950:IGHV1-24"    
    # [25] "ENSG00000224650:IGHV3-74"     "ENSG00000211642:IGLV10-54"   
    # [27] "ENSG00000211643:IGLV5-52"     "ENSG00000211644:IGLV1-51"    
    # [29] "ENSG00000211670:IGLV3-9"      "ENSG00000211677:IGLC2"       
    # [31] "ENSG00000211679:IGLC3"       
    # Let's see how much counts they actually have
    counts_ICTRg = rowSums(current_analysis[[seuratObject_name]]@assays$RNA@counts[expresesd_IC_TR_X_genes_idx,])
    cellFraction_ICTRg = rowMeans(current_analysis[[seuratObject_name]]@assays$RNA@counts[expresesd_IC_TR_X_genes_idx,])
    cellFraction_ICTRg*100
    cellFraction_ICTRg[cellFraction_ICTRg>.001]
    # There are a few with relatively high expression
    # cellFraction_ICTRg[cellFraction_ICTRg>.001]
    # ENSG00000224363:FP700111.1  ENSG00000241767:LINC01324 
    #               0.061157975                0.002364519 
    # ENSG00000231519:AC007285.2     ENSG00000260097:SPDYE6 
    #                0.005591769                0.011087679 
    # ENSG00000272949:AC093668.2     ENSG00000205238:SPDYE2 
    #                0.030291411                0.380623722 
    # ENSG00000205236:AC105052.1 ENSG00000267645:AC105052.3 
    #                0.031793200                0.298983896 
    #    ENSG00000173678:SPDYE2B ENSG00000272896:AL355987.3 
    #                0.019746933                0.002779908 
    # ENSG00000279073:AL807752.6 ENSG00000229257:AL807752.2 
    #                0.002076943                0.002779908
    #
    # But when I remove them and perform this analysis, nothing seems to really change ..
    
    # Create a subset without this biotype
    current_analysis[[paste0(seuratObject_name,'_btypSel')]] = subset(current_analysis[[seuratObject_name]], features=gnames_spall[!IC_TR_X_genes_idx])
    current_analysis[[paste0(seuratObject_name,'_btypSel')]][['rnaMitoCounts']] <- subset(current_analysis[[seuratObject_name]][['rnaMitoCounts']], cells = colnames(current_analysis[[paste0(seuratObject_name,'_btypSel')]]))

    # Save it
    SaveH5Seurat(current_analysis[[paste0(seuratObject_name,'_btypSel')]], filename = paste0(base_dir,'Rdata/H5_',paste0(seuratObject_name,'_btypSel'),'.h5seurat'), overwrite = T)
        # SaveH5Seurat(current_analysis[[paste0(seuratObject_name,'_btypSel')]], filename = paste0(base_dir,'Rdata/H5_',paste0(seuratObject_name,'_btypSel'),'.h5seurat'))
        # current_analysis[[seuratObject_name]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_',seuratObject_name,'.h5seurat'))
    
    #
    #     table(expressed_gene_biotypes)
    # expressed_gene_biotypes
    #       antisense       IG_C_gene IG_C_pseudogene       IG_V_gene IG_V_pseudogene
    #            5153               9               1              10               1
    #         lincRNA  protein_coding       TR_C_gene       TR_V_gene
    #            6603           19262               6               4
    # 
    # For future analysis, might be nice to check consistency, as rooij and hu don't contain the IG_XX and TR_XX biotypes.
    
}
    
    
################################################################################
# Now perform above analysis on the separate datasets, 
# once with RaceID2-like settings, and once with default settings

if ('split_datasets' %in% desired_command) {
    
    # Load all the data
    print('Loading')
    
    seuratObject_name = paste0('RHL_SeuratObject_nM_sel_','ALL.SP','_btypSel')
    current_analysis[[seuratObject_name]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_',seuratObject_name,'.h5seurat'))
    
    # load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_noMito_sel.Rdata'))
    
    # Create separate subsets
    print('Splitting')
    RHL_SeuratObject_nM_sel_HUonly.sp.bt = subset(current_analysis[[seuratObject_name]], subset = annotation_paper_str == 'Hu')
    RHL_SeuratObject_nM_sel_HUonly.sp.bt[['rnaMitoCounts']] <- subset(RHL_SeuratObject_nM_sel_HUonly.sp.bt[['rnaMitoCounts']], cells = colnames(RHL_SeuratObject_nM_sel_HUonly.sp.bt))
    if (DISCARD_GENE_OPTION) { RHL_SeuratObject_nM_sel_HUonly.sp.bt[['rnaDiscardedCounts']] <- subset(RHL_SeuratObject_nM_sel_HUonly.sp.bt[['rnaDiscardedCounts']], cells = colnames(RHL_SeuratObject_nM_sel_HUonly.sp.bt)) }
    
    RHL_SeuratObject_nM_sel_TEICHMANNonly.sp.bt = subset(current_analysis[[seuratObject_name]], subset = annotation_paper_str == 'Teichmann')
    RHL_SeuratObject_nM_sel_TEICHMANNonly.sp.bt[['rnaMitoCounts']] <- subset(RHL_SeuratObject_nM_sel_TEICHMANNonly.sp.bt[['rnaMitoCounts']], cells = colnames(RHL_SeuratObject_nM_sel_TEICHMANNonly.sp.bt))
    if (DISCARD_GENE_OPTION) { RHL_SeuratObject_nM_sel_TEICHMANNonly.sp.bt[['rnaDiscardedCounts']] <- subset(RHL_SeuratObject_nM_sel_TEICHMANNonly.sp.bt[['rnaDiscardedCounts']], cells = colnames(RHL_SeuratObject_nM_sel_TEICHMANNonly.sp.bt)) }

    RHL_SeuratObject_nM_sel_ROOIJonly.sp.bt = subset(current_analysis[[seuratObject_name]], subset = annotation_paper_str == 'vRooij')
    RHL_SeuratObject_nM_sel_ROOIJonly.sp.bt[['rnaMitoCounts']] <- subset(RHL_SeuratObject_nM_sel_ROOIJonly.sp.bt[['rnaMitoCounts']], cells = colnames(RHL_SeuratObject_nM_sel_ROOIJonly.sp.bt))
    if (DISCARD_GENE_OPTION) { RHL_SeuratObject_nM_sel_ROOIJonly.sp.bt[['rnaDiscardedCounts']] <- subset(RHL_SeuratObject_nM_sel_ROOIJonly.sp.bt[['rnaDiscardedCounts']], cells = colnames(RHL_SeuratObject_nM_sel_ROOIJonly.sp.bt)) }

    # Save those
    print('Saving')
    SaveH5Seurat(object = RHL_SeuratObject_nM_sel_HUonly.sp.bt, overwrite = T,
        filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_HUonly.sp.bt.h5seurat'))
    SaveH5Seurat(object = RHL_SeuratObject_nM_sel_TEICHMANNonly.sp.bt, overwrite = T,
        filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_TEICHMANNonly.sp.bt.h5seurat'))
        # RHL_SeuratObject_nM_sel_TEICHMANNonly.sp.bt = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_TEICHMANNonly.sp.bt.h5seurat'))
    SaveH5Seurat(object = RHL_SeuratObject_nM_sel_ROOIJonly.sp.bt, overwrite = T,
        filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_ROOIJonly.sp.bt.h5seurat'))
    
    # I changed the order of things, so now this needn't be done separtely..
    # OLD CODE BELOW
    #
    # # Perhaps a more honest comparison between Rooij and Teichmann is to 
    # # select for septal cells.
    # # So let's also create a Teichmann version of the data which only contains its
    # # septal cells
    # # 
    # # Subset
    # RHL_SeuratObject_nM_sel_TEICHMANN.SP.only =
    #     subset(RHL_SeuratObject_nM_sel_TEICHMANNonly, annotation_region_str == 'SP')
    # # Save
    # SaveH5Seurat(object = RHL_SeuratObject_nM_sel_TEICHMANN.SP.only, overwrite = T,
    #     filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_TEICHMANN.SP.only.h5seurat'))
    
}


################################################################################



################################################################################

# Note: this can do one of six analysis
# required arg, dataset: 'ROOIJ','TEICHMANN','HU'
# required arg, settings: 'SETTINGS_RID2l' or 'SETTINGS_DEFAULT'

# Note to self: for local testing, e.g.:
# base_dir='/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/'
# desired_command='run_separate'
# config=list(); config$dataset='ROOIJ'; config$settings='SETTINGS_RID2l'

# Note to self, to run on the HPC, e.g. use:
# commands="run_separate-dataset=TEICHMANN-settings=SETTINGS_RID2l"

# current_dataset='ALL.SP_btypSel'
# current_dataset='ALL.SP'
# desired_settings='SETTINGS_RID2l'

if ('run_separate' %in% desired_command) {
    
    current_dataset=config$dataset
    desired_settings=config$settings
    
    seuratObject_name = paste0('RHL_SeuratObject_nM_sel_',current_dataset)
    
    # load the file
    current_analysis = list()
    current_analysis[[seuratObject_name]] =
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_',seuratObject_name,'.h5seurat'))
    
    # run the seurat analysis
    if (desired_settings == 'SETTINGS_RID2l') {
        CURRENT_RUNNAME = paste0(current_dataset,'_RID2l')
        # run RaceID2 like settings
        current_analysis[[CURRENT_RUNNAME]] = mySeuratAnalysis(current_analysis[[seuratObject_name]],
            run_name=CURRENT_RUNNAME,
            normalization.method='RC', scale.factor=median(current_analysis[[seuratObject_name]]$nCount.nMT),
            do.scale=F,do.center=F,scale.max=Inf, features_to_use_choice = 'variable') # variable because otherwise too large calculation ..
    } else {
        CURRENT_RUNNAME = paste0(current_dataset,'_default')
        # run default
        current_analysis[[CURRENT_RUNNAME]] = mySeuratAnalysis(current_analysis[[seuratObject_name]],
            run_name=CURRENT_RUNNAME)
    } 
    
    # create respective plots 
    mySeuratCommonPlots(current_analysis[[CURRENT_RUNNAME]], run_name = CURRENT_RUNNAME)
    
    # save the file
    SaveH5Seurat(object = current_analysis[[CURRENT_RUNNAME]], overwrite = T,
        filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',CURRENT_RUNNAME,'.h5seurat'))
    
    print('Plotting done')
    
}
    
# clustering    
# commands="run_separate_nowclusterDE-dataset=TEICHMANNonly_RID2l"
# note that settings name string now is suffix to dataset name
if ('run_separate_nowclusterDE' %in% desired_command) {
    
        CURRENT_RUNNAME=config$dataset
    
        print('Loading data')
    
        current_analysis=list()
        current_analysis[[CURRENT_RUNNAME]]=
            LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',CURRENT_RUNNAME,'.h5seurat'))
    
        print('Starting DE analysis of clusters')
        
        DE_cluster=list()
        DE_cluster[[CURRENT_RUNNAME]] = diff_express_clusters(current_analysis[[CURRENT_RUNNAME]], mc.cores = MYMCCORES)
        save(list = c('DE_cluster'), file = paste0(base_dir,'Rdata/DE_cluster_',CURRENT_RUNNAME,'.Rdata'))
        #load(file = paste0(base_dir,'Rdata/DE_cluster_',CURRENT_RUNNAME,'.Rdata'))
    
        diff_express_clusters_save_results(all_markers = DE_cluster[[CURRENT_RUNNAME]], run_name = CURRENT_RUNNAME,
                base_dir=base_dir, topX = 30)
        
        print('DE analysis of clusters done')
    
}


# commands="rerun_plotsOnly-seuratObject_name=ROOIJonly_RID2l"
# commands="rerun_plotsOnly-seuratObject_name=Huonly_RID2l"
# commands="rerun_plotsOnly-seuratObject_name=TEICHMANNonly_RID2l"
# Running manually: 
# config=list(); config$seuratObject_name = 'ROOIJonly_RID2l'
if ('rerun_plotsOnly' %in% desired_command) {
    
    seuratObject_name=config$seuratObject_name
    
    # load the file
    current_analysis = list()
    current_analysis[[seuratObject_name]] =
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',seuratObject_name,'.h5seurat'))
    
    # create respective plots 
    mySeuratCommonPlots(current_analysis[[seuratObject_name]], run_name = seuratObject_name)
    
    print('Plotting done')
    
}



# Already ran analysis
# CURRENT_RUNNAME = HUonly_RID2l, ROOIJonly_RID2l, ROOIJonly_default
# current_analysis=list();current_analysis[[CURRENT_RUNNAME]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',CURRENT_RUNNAME,'.h5seurat'))

################################################################################

# Run correlation analysis for all sets

# commands="correlations_of_interest"
if ('correlations_of_interest' %in% desired_command) {
 source(paste0(script_dir, 'HCM-SCS_2021-06_SeuratRevisedAnalysis_GeneCorrelationsOfInterest.R'))
}




################################################################################
# This code was not useful; not sure where I got the "weird" values, 
# but both methods produce ±similar values (except log2fold change
# is slightly off). So this code is not necessary.

# # Now let's do the enrichment analysis, but use log10 transformed values etc.
# # Note: Seurat enrichment analysis produced weird values when I used the 
# # RaceID2-similar settings; so to do proper DE analysis, I need to either
# # use another tool, or perform the transformations. Will do latter here.
# 
# # LOG10 clustering    
# # commands="run_separate_nowclusterDE-dataset=TEICHMANNonly_RID2l"
# # Manual:
# # config=list(); config$dataset = 'ROOIJonly_RID2l'
# # 'HUonly_RID2l'
# if ('run_separate_nowclusterDE' %in% desired_command) {
#     
#         CURRENT_RUNNAME=config$dataset
#     
#         print('Loading data')
#     
#         current_analysis=list()
#         current_analysis[[CURRENT_RUNNAME]]=
#             LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',CURRENT_RUNNAME,'.h5seurat'))
#     
#         
#         # Now first add data to this object that is log10 
#         current_analysis[[CURRENT_RUNNAME]][['RNAlog10']] = CreateAssayObject(counts = current_analysis[[CURRENT_RUNNAME]]@assays$RNA@counts)
#         # Now normalize & log scale according to standard seurat practices
#         current_analysis[[CURRENT_RUNNAME]] <- NormalizeData(current_analysis[[CURRENT_RUNNAME]], assay = 'RNAlog10')
#         current_analysis[[CURRENT_RUNNAME]] <- ScaleData(current_analysis[[CURRENT_RUNNAME]], assay = 'RNAlog10')
#             # The clustering should now still simply be present
#             # Idents(current_analysis[[CURRENT_RUNNAME]])
#         
#         
#         # Now perform the enrichment analysis
#         
#         THIS CODE IS NOT FINISHED YET ..
#
#         print('Starting DE analysis of clusters')
#         
#         DE_cluster=list()
#         DE_cluster[[CURRENT_RUNNAME]] = diff_express_clusters(current_analysis[[CURRENT_RUNNAME]], mc.cores = MYMCCORES)
#         save(list = c('DE_cluster'), file = paste0(base_dir,'Rdata/DE_cluster_',CURRENT_RUNNAME,'.Rdata'))
#         #load(file = paste0(base_dir,'Rdata/DE_cluster_',CURRENT_RUNNAME,'.Rdata'))
#     
#         diff_express_clusters_save_results(all_markers = DE_cluster[[CURRENT_RUNNAME]], run_name = CURRENT_RUNNAME,
#                 base_dir=base_dir, topX = 30)
#         
#         print('DE analysis of clusters done')
#     
# }


################################################################################

# commands="rerun_plotsMerged-seuratObject_name=all_RID2l_VAR"
# commands="rerun_plotsMerged-seuratObject_name=default"
if ('rerun_plotsMerged' %in% desired_command) {
    
    seuratObject_name=config$seuratObject_name
    
    # load the file
    current_analysis = list()
    current_analysis[[seuratObject_name]] =
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_merged_noMito_sel.',seuratObject_name,'.h5seurat'))
    
    # create respective plots 
    mySeuratCommonPlots(current_analysis[[seuratObject_name]], run_name = seuratObject_name)
    
    print('Plotting done')
    
}

################################################################################
# Apply batch corrections to generated datasets 
# (the ones with default settings)

# E.g. commands="run_batch_corr_sep-dataset=TEICHMANNonly"

if ('run_batch_corr_sep' %in% desired_command) {

    DATASET_TO_INTEGRATE=config$dataset
    
    Seurat_batch_correction_1c(set_name = DATASET_TO_INTEGRATE, base_dir = base_dir)
    
    corrected_set_name = paste0(set_name, '_Int1c')
    
}

# E.g. commands="run_batch_corr_sep_nowplot_and_DE-dataset=TEICHMANNonly_Int1c"

if ('run_batch_corr_sep_nowplot_and_DE' %in% desired_command) {
     
    corrected_set_name=config$dataset
    
    # Load corrected pre-processed analysis object 
    current_analysis=list()
    current_analysis[[corrected_set_name]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',corrected_set_name,'.h5seurat'))   

    # Plots
    mySeuratCommonPlots(current_analysis[[corrected_set_name]], run_name = corrected_set_name)
    
    # Differential cluster analysis
    DE_cluster=list()
    DE_cluster[[corrected_set_name]] = diff_express_clusters(current_analysis[[corrected_set_name]], mc.cores = MYMCCORES)
    # Save & export results 
    save(list = c('DE_cluster'), file = paste0(base_dir,'Rdata/DE_cluster_',corrected_set_name,'.Rdata'))
    diff_express_clusters_save_results(all_markers = DE_cluster[[corrected_set_name]], 
        run_name = corrected_set_name, base_dir = base_dir, topX = 30)
        
    
}

################################################################################
# Refine the cluster analysis on the pooled datasets to get a resolution that 
# has a desired level of coarse graining

if ('ALL.SP_btypSel_RID2l_ext_cluster_analysis' %in% desired_command) {
    
    # CURRENT_RUNNAME='ALL.SP_RID2l'
    CURRENT_RUNNAME='ALL.SP_btypSel_RID2l'
    CURRENT_RUNNAME_extended = paste0(CURRENT_RUNNAME, '_clExtended')
    
    current_analysis = list()
    current_analysis[[CURRENT_RUNNAME]] = 
        # LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_merged_noMito_sel.',CURRENT_RUNNAME,'.h5seurat'))
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',CURRENT_RUNNAME,'.h5seurat'))
    
    # Create the dataset
    current_analysis[[CURRENT_RUNNAME_extended]]=FindClusters(current_analysis[[CURRENT_RUNNAME]], resolution = seq(.1,1,.1))

    
    # Export plots
    for (resolution in seq(.1,1,.1)) {
        p=DimPlot(current_analysis[[CURRENT_RUNNAME_extended]], group.by = paste0('RNA_snn_res.',resolution), label = T, repel = T, label.size = 7/.pt, pt.size=.5)+
            theme_void()+ggtitle(element_blank())+theme(legend.position = 'none')
        ggsave(plot=p, filename = paste0(base_dir,'Rplots/',CURRENT_RUNNAME_extended,'_7_Clustering_multiLevels_','RNA_snn_res.',resolution,'.pdf'), 
                height=172/3-4, width=172/3-4, units='mm')
    }
      
    # Continue with resolution .1, which suits our purposes
    # Now run diff. expr.
    
    # First set identities to the custom ones
    # Note that this executes +1 action
    new_cls=as.numeric(current_analysis[[CURRENT_RUNNAME_extended]]$RNA_snn_res.0.1)
    current_analysis[[CURRENT_RUNNAME_extended]][['clusters_custom']]=factor(new_cls, levels=min(new_cls):max(new_cls))

    Idents(current_analysis[[CURRENT_RUNNAME_extended]]) <- current_analysis[[CURRENT_RUNNAME_extended]]$clusters_custom
    
    p=DimPlot(current_analysis[[CURRENT_RUNNAME_extended]], label = T, repel = T, label.size = 7/.pt, pt.size=.5)+
            theme_void()+ggtitle(element_blank())+theme(legend.position = 'none')
    ggsave(plot=p, filename = paste0(base_dir,'Rplots/',CURRENT_RUNNAME_extended,'_7_Clustering_multiLevels_','RNA_snn_res.',resolution,'_CHOSEN.pdf'), 
                height=172/3-4, width=172/3-4, units='mm')
    
    # Now perform clustering
    DE_cluster=list()
    DE_cluster[[CURRENT_RUNNAME_extended]] =
        diff_express_clusters(mySeuratObject = current_analysis[[CURRENT_RUNNAME_extended]], mc.cores = MYMCCORES)
    
    # Export results
    DE_out_ALL.SP_RID2l_clExtended = diff_express_clusters_save_results(
      all_markers = DE_cluster[[CURRENT_RUNNAME_extended]], run_name = CURRENT_RUNNAME_extended, base_dir = base_dir, topX = 30, extendedOutput = T, FC_cutoff = 1.1, pval_cutoff = .05)
    
    # Create convenient parameters for later use/export below
    table_topDE = DE_out_ALL.SP_RID2l_clExtended$topHitsPerCluster
    enriched_genes_lists_clusters_ALL.SP = DE_out_ALL.SP_RID2l_clExtended$enriched_genes_lists
    
    
    
    # Now save important data for later use
    
    # Save it
    SaveH5Seurat(object = current_analysis[[CURRENT_RUNNAME_extended]], overwrite = T, 
                 file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',CURRENT_RUNNAME_extended,'.h5seurat'))
        # current_analysis[[CURRENT_RUNNAME_extended]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',CURRENT_RUNNAME_extended,'.h5seurat'))

    
    # Save differential expression results
    save(list='DE_cluster', file = paste0(base_dir,'Rdata/DE_cluster__',CURRENT_RUNNAME_extended,'.Rdata'))
        # load(file = paste0(base_dir,'Rdata/DE_cluster__',CURRENT_RUNNAME_extended,'.Rdata')) # DE_cluster
    save(list='enriched_genes_lists_clusters_ALL.SP', file = paste0(base_dir,'Rdata/enriched_genes_lists_clusters_ALL.SP__',CURRENT_RUNNAME_extended,'.Rdata'))
        
    
}    
    
################################################################################
# More plots
# (Many plots used in paper are generated here)

# Code in:
# source(paste0(script_dir,'Functions/MW_customized_plot_functions.R'))

if ('more_custom_plots' %in% desired_command) {
    
    # --- load / set up data
    
    # CURRENT_RUNNAME='all_RID2l_VAR'
    # CURRENT_RUNNAME='ALL.SP_RID2l_clExtended'
    CURRENT_RUNNAME='ALL.SP_btypSel_RID2l_clExtended'

    current_analysis = list()
    current_analysis[[CURRENT_RUNNAME]] = 
        # LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_merged_noMito_sel.',CURRENT_RUNNAME,'.h5seurat'))
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',CURRENT_RUNNAME,'.h5seurat'))
    
    # Add one-letter shorthand for dataset
    name_change = c("vRooij"="R", "Hu"="H", "Teichmann"="T")
    current_analysis[[CURRENT_RUNNAME]]$annotation_paper_oneletter_fct = factor(
        name_change[current_analysis[[CURRENT_RUNNAME]]$annotation_paper_str], levels=c('R','H','T'))
    # Custom labeling
    name_change = c("vRooij"="HCM", "Hu"="Ctrl1", "Teichmann"="Ctrl2")
    current_analysis[[CURRENT_RUNNAME]]$annotation_paper_beatified = factor(
        name_change[current_analysis[[CURRENT_RUNNAME]]$annotation_paper_str], levels=c('HCM','Ctrl1','Ctrl2'))
    # Custom labeling also for patient
    current_analysis[[CURRENT_RUNNAME]]$annotation_patient_str_beauty = gsub("^R\\.", 'HCM\\.',current_analysis[[CURRENT_RUNNAME]]$annotation_patient_str)
    current_analysis[[CURRENT_RUNNAME]]$annotation_patient_str_beauty = gsub("^H\\.", 'Ctrl1\\.',current_analysis[[CURRENT_RUNNAME]]$annotation_patient_str_beauty)
    current_analysis[[CURRENT_RUNNAME]]$annotation_patient_str_beauty = gsub("^T\\.", 'Ctrl2\\.',current_analysis[[CURRENT_RUNNAME]]$annotation_patient_str_beauty)
    # Create custom order for annotation
    current_analysis[[CURRENT_RUNNAME]]$annotation_paper_fct = 
        factor(current_analysis[[CURRENT_RUNNAME]]$annotation_paper_str, levels=c("vRooij", "Hu", "Teichmann"))
    
    # --- end of load / set up data
    
    # Give overview of number of cells (2nd time I do this -- bit redundant -- though here .SP. dataset properly processed)
    table_per_paper = data.frame(table(current_analysis[[CURRENT_RUNNAME]]$annotation_paper_beatified))
    table_per_patient = data.frame(table(current_analysis[[CURRENT_RUNNAME]]$annotation_patient_str_beauty))
    colnames(table_per_paper) = c('Source', 'Cell_count')
    colnames(table_per_patient) = c('Donor', 'Cell_count')
    openxlsx::write.xlsx(x = list(cellCounts_source=table_per_paper, 
                                  cellCounts_patients=table_per_patient) ,
                         file = paste0(base_dir,'Rplots/QC_general_statistics_cellcounts_alldata_v2.xlsx'), overwrite = T)
    
    
    # Original plots, with a few extra tweaks
    mySeuratCommonPlots(current_analysis[[CURRENT_RUNNAME]], CURRENT_RUNNAME, add_custom_fields = 'annotation_paper_beatified')
    
    
    # Beautified UMAP with custom color
    # library(scales); library(wesanderson); show_col(wesanderson::wes_palettes$BottleRocket2)
    #
    custom_colors = c('#bd0020','#9d9d9c','#575756')
    p=DimPlot(current_analysis[[CURRENT_RUNNAME]], group.by = 'annotation_paper_beatified', cols = custom_colors, label = T, repel = T, label.size = 7/.pt, pt.size = NULL)+
        theme_void()+ggtitle(element_blank())+theme(legend.position = 'none')
    # test
    ggsave(plot = p, filename = paste0(base_dir,'Rplots/',CURRENT_RUNNAME,'_2_umapLabeled_by_','annotation_paper_beatified','_style4-customColors.pdf'), height=172/3-4, width=172/3-4, units='mm')
    
    
    # Stress gene plots
    LIST_NAME='stress'; CAT = 'annotation_paper_beatified'; CATNAME = 'Dataset'; PLOTTYPE='boxplot'
    genes_of_interest = c('NPPA', 'NPPB', 'ACTA1', 'MYH7')#, 'MYH6')
    plot_list = lapply(genes_of_interest, 
        function(x) {shorthand_plotViolinBox_custom(current_analysis,analysis_name=CURRENT_RUNNAME,cat_by = CAT, 
                                                    gene_of_interest=x,base_dir=base_dir, cat_name = CATNAME,
                                                    type = 'violin')})
    # Small 2x2 version
    p=wrap_plots(plot_list, nrow = 2)& theme(plot.margin = margin(t = .5, r = .5, b = .5, l = 2, unit = "mm"))
    ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',CURRENT_RUNNAME,'_6_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'.pdf'), width = 172/3-4, height= 172/3-4, units='mm', device=cairo_pdf)
    #p
    p=wrap_plots(plot_list, nrow = 1)& theme(plot.margin = margin(t = .5, r = .5, b = .5, l = 2, unit = "mm"))
    # ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',CURRENT_RUNNAME,'_6_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'-style2.pdf'), width = 172-4, height= (172-4)/4, units='mm', device=cairo_pdf)
    #ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',CURRENT_RUNNAME,'_6_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'-style2.pdf'), width = (172/3*2)-4, height= ((172/3*2)-4)/4*1.5, units='mm', device=cairo_pdf)
    ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',CURRENT_RUNNAME,'_6_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'-style2.pdf'), width = (172/3*2)-4, height= PANEL_HEIGHT, units='mm', device=cairo_pdf)
    
    
    
    # Distribution of source paper over clusters (bar plot)
    # library(scales); show_col(col_vector_60)
    # library(wesanderson); show_col(wesanderson::wes_palettes$BottleRocket2)
    # TO DO: Maybe put wes anderson here? Because color conflict?
    custom_colors = c('#bd0020','#9d9d9c','#575756')
    p=ggplot(data.frame( cluster = Idents(current_analysis[[CURRENT_RUNNAME]]),
                         Donor   = current_analysis[[CURRENT_RUNNAME]]$annotation_paper_beatified))+
        geom_bar(aes(x=cluster, fill=Donor))+theme_bw()+
        xlab('Cluster')+ylab('Number of cells')+
        give_better_textsize_plot(8)+
        theme(legend.position = 'right', legend.key.size = unit(3, "mm"))+
        scale_fill_manual(values = custom_colors)
    ggsave(filename = paste0(base_dir,'Rplots/',CURRENT_RUNNAME,'_5_Barplot_PatientCluster_distr_Datasets.pdf'), 
        plot = p, height=172/3-4, width=172/3-4, units='mm', device = cairo_pdf)
    
    # Now create MYH6 separately, also use to add legend
    TYPE='violin'
    p_=shorthand_plotViolinBox_custom(current_analysis,analysis_name=CURRENT_RUNNAME,cat_by = 'annotation_paper_fct', 
                                                    gene_of_interest='MYH6',base_dir=base_dir, cat_name = CATNAME,
                                                    type = TYPE)+theme(legend.position='right')
    ggsave(plot = p_,filename = paste0(base_dir, 'Rplots/',CURRENT_RUNNAME,'_6_',TYPE,'_FORLEGEND.pdf'), width = 7, height= 6, units='cm', device=cairo_pdf)

    
    # NPPA custom plot
    GENE='NPPA'
    p=shorthand_seurat_custom_expr(seuratObject = current_analysis[[CURRENT_RUNNAME]], 
                                         gene_of_interest = GENE,
                                         textsize = 8, pointsize = .1, custom_title = GENE, mymargin = .1) 
    # p
    ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_',GENE,'_small.pdf'), plot = p,
           height=(172/3-4), width=(172/3-4), units='mm', device = cairo_pdf)
    
    # Then look at how genes that came up during analysis are expressed
    # Load regulon file
    REGULON_DATASET = 'ROOIJonly.sp.bt_RID2l'
    if (!(  file.exists(paste0(base_dir,'Rdata/',REGULON_DATASET,'_core_regulons_sorted.Rdata'))  )) {
        print('Regulon file doesnt exist yet..') 
    } else { 
        load(file = paste0(base_dir,'Rdata/',REGULON_DATASET,'_core_regulons_sorted.Rdata')) # core_regulons_sorted 
    
        # Box plots
        shorthand_custom_boxplot(seuratObject_list=current_analysis, 
                                 gene_lists=core_regulons_sorted, 
                                 seuratObjectNameToTake=CURRENT_RUNNAME, 
                                 group.by='annotation_paper_beatified', 
                                 topX=10, mylimits=.01) 
        
        # Summary composite plots
        # seuratObject_list=current_analysis; gene_lists=core_regulons_sorted; seuratObjectNameToTake=CURRENT_RUNNAME; group.by='annotation_paper_beatified'; group.by2='annotation_patient_fct'; zscore=T
        p_lists=shorthand_custom_compositeplot(seuratObject_list=current_analysis, 
                                 gene_lists=core_regulons_sorted, 
                                 seuratObjectNameToTake=CURRENT_RUNNAME, 
                                 group.by='annotation_paper_beatified', 
                                 group.by2='annotation_patient_fct',
                                 zscore=T) 
            # To load data at later point:
            if (F) {
                REGULON_DATASET = 'ROOIJonly.sp.bt_RID2l'
                CURRENT_RUNNAME='ALL.SP_btypSel_RID2l_clExtended'
                group.by='annotation_paper_beatified'; group.by2='annotation_patient_fct'
                last_list_name = names(core_regulons_sorted)[length(core_regulons_sorted)]
                df_agr_combined = openxlsx::read.xlsx( paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customSummaryMean_ALL-', last_list_name,'-etc_for_',group.by,'_splt2',group.by2,'-data.xlsx'))
            }
        
        p=wrap_plots(p_lists$p_violin_list, nrow=1)
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customCOMPOSITE_REG.pdf'), plot = p,
               height=(PANEL_WIDTH*3-4)/length(core_regulons_sorted), width=(PANEL_WIDTH*3-4), units='mm', device=cairo_pdf)
        
        p=wrap_plots(p_lists$p_bar_list_g2, nrow=1)
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customCOMPOSITE_REG_g2.pdf'), plot = p,
               height=(PANEL_WIDTH*3-4)/length(core_regulons_sorted), width=(PANEL_WIDTH*3-4), units='mm', device=cairo_pdf)
        
        # core regulon expr on UMAPs
        p_list_modules = lapply(names(core_regulons_sorted), function(reg_name) {
                    shorthand_seurat_custom_expr(seuratObject = current_analysis[[CURRENT_RUNNAME]], 
                                     gene_of_interest = core_regulons_sorted[[reg_name]],
                                     textsize = 6, pointsize = .5, custom_title = reg_name, mymargin = .1, zscore = T) 
            })
        p_modules=wrap_plots(p_list_modules, nrow=1)
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_REG.pdf'), plot = p_modules,
               height=(PANEL_WIDTH*3-4)/5, width=(PANEL_WIDTH*3-4), units='mm', device=cairo_pdf)    
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_REG.png'), plot = p_modules,
               height=(PANEL_WIDTH*3-4)/5, width=(PANEL_WIDTH*3-4), units='mm', dpi=1200)      
        
        
        # Also do SCENIC regulons
        # === 
        ROOIJ_DATASET = 'ROOIJonly.sp.bt_RID2l'
        load(file=paste0(base_dir,'Rdata/',ROOIJ_DATASET,'__SCENIC_reg_top_genes_sorted_full.Rdata')) # SCENIC_reg_top_genes_sorted_full
        
        # Box plots
        shorthand_custom_boxplot(seuratObject_list=current_analysis, 
                                 gene_lists=SCENIC_reg_top_genes_sorted_full, 
                                 seuratObjectNameToTake=CURRENT_RUNNAME, 
                                 group.by='annotation_paper_oneletter_fct', 
                                 topX=10, mylimits=.01) 
        
        # Summary plots
        p_lists=shorthand_custom_compositeplot(seuratObject_list=current_analysis, 
                                 gene_lists=SCENIC_reg_top_genes_sorted_full, 
                                 seuratObjectNameToTake=CURRENT_RUNNAME, 
                                 group.by='annotation_paper_oneletter_fct', 
                                 group.by2='annotation_patient_fct',
                                 zscore=T) 
        
            # Comparison of which SCENIC regulons show most enrichement 
            # in the HCM data compared to our data.
            # 
            # Posisble to load data at later point
            if (F) {
                REGULON_DATASET = 'ROOIJonly.sp.bt_RID2l'
                CURRENT_RUNNAME='ALL.SP_btypSel_RID2l_clExtended'
                group.by='annotation_paper_oneletter_fct'; group.by2='annotation_patient_fct'
                last_list_name = names(SCENIC_reg_top_genes_sorted_full)[length(SCENIC_reg_top_genes_sorted_full)]
                df_agr_combined_SCENIC = openxlsx::read.xlsx( paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customSummaryMean_ALL-', last_list_name,'-etc_for_',group.by,'_splt2',group.by2,'-data.xlsx'))
                # Stats # label:scenic_HCM_enrich_stats
                agg2_scenic_out = aggregate(x = df_agr_combined_SCENIC$expression, by=list(df_agr_combined_SCENIC$annotation_paper_oneletter_fct, df_agr_combined_SCENIC$set_name), mean)
                R_values=agg2_scenic_out[agg2_scenic_out$Group.1=='R',]
                T_values=agg2_scenic_out[agg2_scenic_out$Group.1=='T',]
                H_values=agg2_scenic_out[agg2_scenic_out$Group.1=='H',]
                R_values$x_relative = R_values$x-((H_values$x+T_values$x)/2)
                R_values$x_compmax = R_values$x - apply(rbind(H_values$x,T_values$x), 2, max)
                
                toString(R_values[order(R_values$x_relative, decreasing = T),]$Group.2)
                toString(R_values[order(R_values$x_compmax, decreasing = T),]$Group.2)
                sum(R_values$x_relative>0)
                sum(R_values$x>0)
                ggplot(R_values, aes(x=Group.2, y=x_relative))+
                    geom_bar(stat='identity')+
                    coord_flip()
                ggplot(R_values, aes(x=Group.2, y=x_compmax))+
                    geom_bar(stat='identity')+
                    coord_flip()
            }
        
        
        p=wrap_plots(p_lists$p_violin_list, nrow=5)
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customCOMPOSITE_REG_SCENIC.pdf'), plot = p,
               height=(PANEL_WIDTH*3-4), width=(PANEL_WIDTH*3-4), units='mm', device=cairo_pdf)
        
        p=wrap_plots(p_lists$p_bar_list_g2, nrow=5)
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customCOMPOSITE_REG_SCENIC_g2.pdf'), plot = p,
               height=(PANEL_WIDTH*3-4), width=(PANEL_WIDTH*3-4),  units='mm', device=cairo_pdf)
        
        # UMAP
        p_list = lapply(names(SCENIC_reg_top_genes_sorted_full), function(reg_name) {
                    shorthand_seurat_custom_expr(seuratObject = current_analysis[[CURRENT_RUNNAME]], 
                                     gene_of_interest = SCENIC_reg_top_genes_sorted_full[[reg_name]],
                                     textsize = 6, pointsize = .5, custom_title = reg_name, mymargin = .1, zscore = T) 
            })
        p=wrap_plots(p_list, nrow=5)
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_REG_SCENIC.pdf'), plot = p,
               height=(PANEL_WIDTH*3-4), width=(PANEL_WIDTH*3-4), units='mm', device=cairo_pdf)      
        # Also in png to avoid difficult-to-load files
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_REG_SCENIC.png'), plot = p,
               height=(PANEL_WIDTH*3-4), width=(PANEL_WIDTH*3-4), units='mm', dpi=1200)
        # UMAP, smaller size
        p_list = lapply(names(SCENIC_reg_top_genes_sorted_full), function(reg_name) {
                    shorthand_seurat_custom_expr(seuratObject = current_analysis[[CURRENT_RUNNAME]], 
                                     gene_of_interest = SCENIC_reg_top_genes_sorted_full[[reg_name]],
                                     textsize = 6, pointsize = .1, custom_title = reg_name, mymargin = .1, zscore = T) 
            })
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_REG_SCENIC-smaller.pdf'), plot = p,
               height=(PANEL_WIDTH*3-4), width=(PANEL_WIDTH*3-4), units='mm', device=cairo_pdf)      
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_REG_SCENIC-smaller.png'), plot = p,
               height=(PANEL_WIDTH-4), width=(PANEL_WIDTH-4), units='mm', dpi=1200)    


        # Now do the SCENIC TFs
        shorthand_custom_boxplot_perpatient(seuratObject_list=current_analysis, 
                         gene_lists=list(SCENIC_TFs=gsub('\\(\\+\\)','',names(SCENIC_reg_top_genes_sorted_full))), 
                         seuratObjectNameToTake=CURRENT_RUNNAME, 
                         group.by='annotation_paper_beatified', 
                         aggr.by='annotation_patient_fct',
                         topX=10, mylimits=.01) 
        
    }
    
    if (!(  file.exists(paste0(base_dir,'Rplots/gene_lists_customcorrelated.Rdata'))  )) {
        print('Gene lists custom correlations file doesnt exist yet..') 
    } else { 
        load(paste0(base_dir,'Rdata/gene_lists_customcorrelated__Rooijbased.Rdata'))
        
        # Re-organize lists
        # Now added also "original" gene of interest itself to beginning of lists
        gene_lists_customcorrelated_reorganized=list()
        for (gene in names(gene_lists_customcorrelated)) {
            if (!is.na(gene_lists_customcorrelated[[gene]]$pos)) {
                gene_lists_customcorrelated_reorganized[[paste0('posCorrWith_',gene)]] = c(gene, gene_lists_customcorrelated[[gene]]$pos)
            }
            if (!is.na(gene_lists_customcorrelated[[gene]]$neg)) {
                gene_lists_customcorrelated_reorganized[[paste0('negCorrWith_',gene)]] = c(gene, gene_lists_customcorrelated[[gene]]$neg)
            }
        }
    
        # Box plots
        shorthand_custom_boxplot_perpatient(
                                 seuratObject_list=current_analysis, 
                                 gene_lists=gene_lists_customcorrelated_reorganized, 
                                 seuratObjectNameToTake=CURRENT_RUNNAME, 
                                 group.by='annotation_paper_beatified', 
                                 aggr.by='annotation_patient_fct',
                                 topX=10, mylimits=.01, mySize=c(PANEL_WIDTH*2-4,PANEL_HEIGHT), myFontSize=8) 
        
    }
    
    
    # ROOIJ CLUSTERS
    # Now also check whether the genes in the clusters are enriched in direct comparison of three samples
    # 
    ANALYSIS_NAME_clExtended = 'ROOIJonly.sp.bt_RID2l_clExtended'
    if (!(  file.exists(paste0(base_dir,'Rdata/enriched_genes_lists_clusters_ROOIJ__',ANALYSIS_NAME_clExtended,'.Rdata'))  )) {
        print('Gene lists custom correlations file doesnt exist yet..') 
    } else { 
        load(paste0(base_dir,'Rdata/enriched_genes_lists_clusters_ROOIJ__',ANALYSIS_NAME_clExtended,'.Rdata')) # enriched_genes_lists_clusters_ROOIJ
    
        enriched_genes_lists_clusters_ROOIJ_=enriched_genes_lists_clusters_ROOIJ
        names(enriched_genes_lists_clusters_ROOIJ_) = paste0('Cl.',names(enriched_genes_lists_clusters_ROOIJ))
        
        # per-patient plots
        # Summary plots
        # Note: these are not so useful, since the enrichment is really per cluster !
        p_lists=shorthand_custom_compositeplot(seuratObject_list=current_analysis, 
                                 gene_lists=enriched_genes_lists_clusters_ROOIJ_, 
                                 seuratObjectNameToTake=CURRENT_RUNNAME, 
                                 group.by='annotation_paper_oneletter_fct', 
                                 group.by2='annotation_patient_fct',
                                 zscore=T) 
        p=wrap_plots(p_lists$p_violin_list, nrow=1)
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customCOMPOSITE_CL.pdf'), plot = p,
               height=(PANEL_WIDTH*3-4)/length(enriched_genes_lists_clusters_ROOIJ_), width=(PANEL_WIDTH*3-4), units='mm', device=cairo_pdf)
        
        p=wrap_plots(p_lists$p_bar_list_g2, nrow=1)
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customCOMPOSITE_CL_g2.pdf'), plot = p,
               height=(PANEL_WIDTH*3-4)/length(enriched_genes_lists_clusters_ROOIJ_), width=(PANEL_WIDTH*3-4), units='mm', device=cairo_pdf)
        
        # enriched genes expr on UMAPs
        ZOOM_FACTOR=4
        p_list_clusters = lapply(names(enriched_genes_lists_clusters_ROOIJ_), function(cl_name) {
                    shorthand_seurat_custom_expr(seuratObject = current_analysis[[CURRENT_RUNNAME]], 
                                     gene_of_interest = enriched_genes_lists_clusters_ROOIJ_[[cl_name]],
                                     textsize = 6*ZOOM_FACTOR, pointsize = .5, custom_title = cl_name, mymargin = .5*ZOOM_FACTOR, zscore = T) 
                                        # note: text size twice as large, because i save at zoom 200%, as trick to reduce point size
            })
        # p_clusters=wrap_plots(p_list_clusters, nrow=1)
        # ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_CL.pdf'), plot = p_clusters,
        #        height=(PANEL_WIDTH*3-4)/5, width=(PANEL_WIDTH*3-4), units='mm', device=cairo_pdf)    
        # ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_CL.png'), plot = p_clusters,
        #        height=(PANEL_WIDTH*3-4)/5, width=(PANEL_WIDTH*3-4), units='mm', dpi=1200)    
        p_clusters=wrap_plots(p_list_clusters, nrow=2)
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_CL_v2.png'), plot = p_clusters,
               height=(PANEL_WIDTH-4)*2/3*ZOOM_FACTOR, width=(PANEL_WIDTH-4)*ZOOM_FACTOR, units='mm', dpi=1200)    
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_CL_v2.pdf'), plot = p_clusters,
               height=(PANEL_WIDTH-4)*2/3*ZOOM_FACTOR, width=(PANEL_WIDTH-4)*ZOOM_FACTOR, units='mm', device=cairo_pdf)
        # long aspect ratio
        p_clusters=wrap_plots(p_list_clusters, nrow=1)
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_CL_v3l.png'), plot = p_clusters,
               height=(3*PANEL_WIDTH-4)/6*ZOOM_FACTOR, width=(3*PANEL_WIDTH-4)*ZOOM_FACTOR, units='mm', dpi=1200)    
        ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_CL_v3l.pdf'), plot = p_clusters,
               height=(3*PANEL_WIDTH-4)/6*ZOOM_FACTOR, width=(3*PANEL_WIDTH-4)*ZOOM_FACTOR, units='mm', device=cairo_pdf)    
        
        # Show plot per cluster
        # current_analysis[[CURRENT_RUNNAME]]$
        # (..) unfinished
        
        # NOMURA OVERLAP LISTS
        if (file.exists(paste0(base_dir,'Rdata/zcustom__genelist_module2_nomura_overlap.Rdata'))&file.exists(paste0(base_dir,'Rdata/zcustom__genelist_SCENIC.CEBPZ_nomura_overlap.Rdata'))) {
            
            load(paste0(base_dir,'Rdata/zcustom__genelist_module2_nomura_overlap.Rdata')) # genelist_module2_nomura_overlap
            load(paste0(base_dir,'Rdata/zcustom__genelist_SCENIC.CEBPZ_nomura_overlap.Rdata')) # genelist_SCENIC.CEBPZ_nomura_overlap

            # module 2 expr on UMAPs
            ZOOM_FACTOR=4
            NR_GENES=10
            ACTUAL_NR_GENES = min(length(genelist_module2_nomura_overlap), NR_GENES)
            p_list_nomura1 = lapply(genelist_module2_nomura_overlap[1:ACTUAL_NR_GENES], function(current_gene) {
                        shorthand_seurat_custom_expr(seuratObject = current_analysis[[CURRENT_RUNNAME]], 
                                         gene_of_interest = current_gene,
                                         textsize = 6*ZOOM_FACTOR, pointsize = .5, custom_title = current_gene, mymargin = .5*ZOOM_FACTOR, zscore = T) 
                                            # note: text size twice as large, because i save at zoom 200%, as trick to reduce point size
                })
            p_nomura1=wrap_plots(p_list_nomura1, nrow=1)
            ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_Nomura-module2_v3l.png'), plot = p_nomura1,
               height=(3*PANEL_WIDTH-4)/NR_GENES*ZOOM_FACTOR, width=(3*PANEL_WIDTH-4)*ZOOM_FACTOR*ACTUAL_NR_GENES/10, units='mm', dpi=1200)    
        
            # SCENIC CEBPZ expr on UMAPs
            ZOOM_FACTOR=4
            NR_GENES=10
            ACTUAL_NR_GENES = min(length(genelist_SCENIC.CEBPZ_nomura_overlap), NR_GENES)
            p_list_nomura2 = lapply(genelist_SCENIC.CEBPZ_nomura_overlap[1:ACTUAL_NR_GENES], function(current_gene) {
                        shorthand_seurat_custom_expr(seuratObject = current_analysis[[CURRENT_RUNNAME]], 
                                         gene_of_interest = current_gene,
                                         textsize = 6*ZOOM_FACTOR, pointsize = .5, custom_title = current_gene, mymargin = .5*ZOOM_FACTOR, zscore = T) 
                                            # note: text size twice as large, because i save at zoom 200%, as trick to reduce point size
                })
            p_nomura2=wrap_plots(p_list_nomura2, nrow=1)
            ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_Nomura-scenicCEBPZ_v3l.png'), plot = p_nomura2,
               height=(3*PANEL_WIDTH-4)/NR_GENES*ZOOM_FACTOR, width=(3*PANEL_WIDTH-4)*ZOOM_FACTOR*ACTUAL_NR_GENES/10, units='mm', dpi=1200)    
        
                
        }
        

        
    }    
    
}

################################################################################
# Custom plots for our dataset 
# (Many plots used in paper are generated here)
# SEE ALSO SCRIPT "CUSTOM_ANALYSIS_ROOIJ"

if ('more_custom_plots' %in% desired_command) {
    
    DATASET_NAME='ROOIJonly.sp.bt_RID2l_clExtended'
    # CURRENT_RUNNAME='ALL.SP_RID2l'
    
    if (!exists('current_analysis')) {current_analysis = list()}
    current_analysis[[DATASET_NAME]] = 
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))

    # UMAP with clusters indicated (ROOIJ only)
    # ====
    mySeuratObject=current_analysis[[DATASET_NAME]]
    run_name=DATASET_NAME
    mymarkers = c('MALAT1', 'TTN', 'MYH7', 'MYH6', 'NPPA', 'NPPB', 'ACTA1','MYL2','SORBS2','CSRP3','NDUFA4','CRYAB','HSPB1', 'KCNQ1OT1')
    mypointsize=1
    
    mySeuratCommonPlots(mySeuratObject, run_name, mymarkers, mypointsize)
    # --> Show umap with annotations --> More customized/stylized version, specific size
    # --> # Distribution of patients over clusters
    
    # Some more custom UMAPs
    for (GENE in c('MYH7', 'MYH6')) {
        p=shorthand_seurat_custom_expr(seuratObject = current_analysis[[DATASET_NAME]], 
                                         gene_of_interest = GENE,
                                         textsize = 8, pointsize = 1, custom_title = GENE, mymargin = .1) 
        # p
        ggsave(filename = paste0(base_dir, 'Rplots/', DATASET_NAME, '_9_customUMAPs_',GENE,'.pdf'), plot = p,
               height=172/3-4, width=172/3-4, units='mm', device = cairo_pdf)     
    }
    # And some more
    # label:UMAPforCorrs
    for (GENE in c('NPPA', 'NPPB', 'IGFBP2', 'RTN4', 'GAPDH', 'ACTC1', 'MYL4', 'TTN', 'RYR2', 'NRAP')) {
        # GENE='RTN4'
        p=shorthand_seurat_custom_expr(seuratObject = current_analysis[[DATASET_NAME]], 
                                         gene_of_interest = GENE,
                                         textsize = 8, pointsize = 1, custom_title = GENE, mymargin = .1) 
        # p
        ggsave(filename = paste0(base_dir, 'Rplots/', DATASET_NAME, '_9_customUMAPs_',GENE,'_small.pdf'), plot = p,
               height=(PANEL_WIDTH-4)/2, width=((PANEL_WIDTH*2)-4)/3, units='mm', device = cairo_pdf)     
    }
    
    # Generate QC plots for this dataset that look pretty
    # Automatically saves plots
    newfact=gsub('R\\.','',current_analysis[[DATASET_NAME]]$annotation_patient_fct)
    current_analysis[[DATASET_NAME]]$annotation_patient_simplified_fct = factor(newfact, levels=sort(unique(newfact)))
    shorthand_plotViolinBox_custom(myseuratobjectlist = current_analysis,analysis_name = DATASET_NAME,
                                   cat_by = 'annotation_patient_simplified_fct',cat_name = 'patient',base_dir = base_dir, 
                                   manual_expr = current_analysis[[DATASET_NAME]]$percent.mt[colnames(current_analysis[[DATASET_NAME]])],
                                   manual_expr_featname='Mitochondrial percentage', type='violin', custom_ylim=c(0,100),
                                   custom_ylab = 'Mitochondrial percentage', custom_title=element_blank())
    print(paste0('mean mito pct all cells: ', mean(current_analysis[[DATASET_NAME]]$percent.mt[colnames(current_analysis[[DATASET_NAME]])])))
    print(paste0('sd mito pct all cells: ', sd(current_analysis[[DATASET_NAME]]$percent.mt[colnames(current_analysis[[DATASET_NAME]])])))
    
    # For comparison
    HU_DATASET = 'HUonly_RID2l'
    if (!(HU_DATASET %in% names(current_analysis))) {
        current_analysis[[HU_DATASET]] = 
            LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',HU_DATASET,'.h5seurat')) 
    }

    custom_levels=c('H.N1', 'H.N2', 'H.N3', 'H.N4', 'H.N5', 'H.N13', 'H.N14') 
    current_analysis[[HU_DATASET]]$annotation_patient_fct = factor(current_analysis[[HU_DATASET]]$annotation_patient_str, levels=custom_levels)
    shorthand_plotViolinBox_custom(myseuratobjectlist = current_analysis,analysis_name = HU_DATASET,
                               cat_by = 'annotation_patient_fct',cat_name = 'patient',base_dir = base_dir, 
                               manual_expr = current_analysis[[HU_DATASET]]$percent.mt[colnames(current_analysis[[HU_DATASET]])],
                               manual_expr_featname='Mitochondrial percentage', type='violin', custom_ylab = 'Mitochondrial percentage', 
                               custom_title=element_blank(), custom_ylim=c(0,100))
    
    # Now distribution of read counts
    current_analysis[[DATASET_NAME]]$nCount_rnaMitoCounts[1:10]
    current_analysis[[DATASET_NAME]]$n_counts
    current_analysis[[DATASET_NAME]]$nCount.nMT[1:10]
    ggplot(data.frame(count=log10(.1+current_analysis[[DATASET_NAME]]$nCount.nMT)), aes(x=count))+
        geom_histogram(binwidth = .1)


}    
    
################################################################################
# Some custom stuff (put somewhere else later)

if (F) {
 diff_express_clusters_save_results(DE_cluster$all_RID2l_VAR, 'all_RID2l_VAR', base_dir)
 diff_express_clusters_save_results(DE_cluster$ROOIJonly_RID2l, 'all_RID2l_VAR', base_dir)
}

################################################################################

# Probably not do this any more, way too large
if (F) {
    
    save.image(file = paste0(base_dir,'Rdata/RHL_whole_analysis.Rdata'))
    # load(file = paste0(base_dir,'Rdata/RHL_whole_analysis.Rdata'))
    
}

################################################################################
# Create overview of cell numbers

# ALL.SP_RID2l

################################################################################
# Load a dataset manually

# H5_RHL_SeuratObject_nM_sel_TEICHMANNonly.h5seurat

# DATASET_NAME='ALL.SP_RID2l'
# DATASET_NAME='TEICHMANNonly_RID2l'
# DATASET_NAME='ROOIJonly_Int1c'
# DATASET_NAME='HUonly_RID2l'
# DATASET_NAME='HUonly'
# DATASET_NAME='ROOIJonly_default'
# DATASET_NAME='ROOIJonly_RID2l'
# DATASET_NAME='ROOIJonly_RID2l_clExtended'
if (F) {
    
    # load the file
    if (!exists('current_analysis')) {current_analysis = list()}
    current_analysis[[DATASET_NAME]] =
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))
    

}
if (F) {
    # To load full pooled thing
    CURRENT_RUNNAME='all_RID2l_VAR'
    current_analysis = list()
    current_analysis[[CURRENT_RUNNAME]] = 
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_merged_noMito_sel.',CURRENT_RUNNAME,'.h5seurat'))
        # NOTE: RHL_SeuratObject_merged_noMito_sel.Rdata contains data that isn't normalized properly
    
}







