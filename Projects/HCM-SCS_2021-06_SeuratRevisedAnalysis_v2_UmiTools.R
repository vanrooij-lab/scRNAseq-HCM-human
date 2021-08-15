
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
    base_dir = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/'
    
    data_dir1 = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_UmiTools_CountTables_Wang_HCMSCS/ROOIJ/' # NULL
    data_dir2 = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_UmiTools_CountTables_Wang_HCMSCS/HU/'
        # '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_Wang_Counttables/'
} else {
    script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'
    base_dir = '/hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3/'
    
    data_dir1 = '/hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3/counttables_tsv_exclmulti_filtered/'
        # '/hpc/hub_oudenaarden/mwehrens/fastq/HCM_SCS/mapping.93.may25/'
    data_dir2 = '/hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3/counttables_tsv_exclmulti_filtered/'
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
    
    save(list = c('SCS_df_list_data'),file = paste0(base_dir,'Rdata/SCS_df_list_data.Rdata'))

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
    # load(paste0(base_dir,'Rdata/RHL_object_list.Rdata'))
    
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
    
    # actually update names    
    hgnc_names_Teich=rownames(RHL_SeuratObject_list$TEICH@assays$RNA@data)
    rownames(RHL_SeuratObject_list$TEICH@assays$RNA@data) = paste0(ensembl_IDs_check,':',hgnc_names_Teich)
    rownames(RHL_SeuratObject_list$TEICH@assays$RNA@counts) = paste0(ensembl_IDs_check,':',hgnc_names_Teich)
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
    
    save(list = c('RHL_SeuratObject_merged'),file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged.Rdata'))
    # load(paste0(base_dir,'Rdata/RHL_SeuratObject_merged.Rdata'))
    
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
    
    load(paste0(base_dir,'Rdata/RHL_SeuratObject_merged.Rdata'))
    
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
    
    # Then, subset the RNA assay (note: "noMito" is historic name, earlier only mito genes removed, now many more)
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
    RHL_SeuratObject_merged_noMito_sel <- subset(RHL_SeuratObject_merged_noMito, subset = nFeature.nMT > 50 & nFeature.nMT > 1000 & percent.mt < 80)
    # Little boilerplate, additional selection of Hu cells based on region
    # unique(RHL_SeuratObject_merged_noMito_sel[["annotation_region_str"]][RHL_SeuratObject_merged_noMito_sel[["annotation_paper_str"]]=='Hu'])
    RHL_SeuratObject_merged_noMito_sel = RHL_SeuratObject_merged_noMito_sel[,!(RHL_SeuratObject_merged_noMito_sel[["annotation_region_str"]]=='LA'&
               RHL_SeuratObject_merged_noMito_sel[["annotation_paper_str"]]=='Hu')]
    
    # Now also apply this to the stored "discarded" data
    RHL_SeuratObject_merged_noMito_sel[['rnaMitoCounts']] <- subset(RHL_SeuratObject_merged_noMito[['rnaMitoCounts']], cells = colnames(RHL_SeuratObject_merged_noMito_sel))
    if (DISCARD_GENE_OPTION) {
        RHL_SeuratObject_merged_noMito_sel[['rnaDiscardedCounts']] <- subset(RHL_SeuratObject_merged_noMito[['rnaDiscardedCounts']], cells = colnames(RHL_SeuratObject_merged_noMito_sel))
    }
    
    # table(RHL_SeuratObject_merged_sel$orig.ident)
    table(RHL_SeuratObject_merged_noMito_sel$annotation_paper_str)
    table(RHL_SeuratObject_merged_noMito_sel$annotation_patient_str)
    
    # Save it
    save(list = c('RHL_SeuratObject_merged_noMito_sel'), file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_noMito_sel.Rdata'))
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
# Now perform above analysis on the separate datasets, 
# once with RaceID2-like settings, and once with default settings

if ('split_datasets' %in% desired_command) {
    
    # Load all the data
    print('Loading')
    load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_noMito_sel.Rdata'))
    
    # Create separate subsets
    print('Splitting')
    RHL_SeuratObject_nM_sel_HUonly = subset(RHL_SeuratObject_merged_noMito_sel, subset = annotation_paper_str == 'Hu')
    RHL_SeuratObject_nM_sel_HUonly[['rnaMitoCounts']] <- subset(RHL_SeuratObject_nM_sel_HUonly[['rnaMitoCounts']], cells = colnames(RHL_SeuratObject_nM_sel_HUonly))
    if (DISCARD_GENE_OPTION) { RHL_SeuratObject_nM_sel_HUonly[['rnaDiscardedCounts']] <- subset(RHL_SeuratObject_nM_sel_HUonly[['rnaDiscardedCounts']], cells = colnames(RHL_SeuratObject_nM_sel_HUonly)) }
    
    RHL_SeuratObject_nM_sel_TEICHMANNonly = subset(RHL_SeuratObject_merged_noMito_sel, subset = annotation_paper_str == 'Teichmann')
    RHL_SeuratObject_nM_sel_TEICHMANNonly[['rnaMitoCounts']] <- subset(RHL_SeuratObject_nM_sel_TEICHMANNonly[['rnaMitoCounts']], cells = colnames(RHL_SeuratObject_nM_sel_TEICHMANNonly))
    if (DISCARD_GENE_OPTION) { RHL_SeuratObject_nM_sel_TEICHMANNonly[['rnaDiscardedCounts']] <- subset(RHL_SeuratObject_nM_sel_TEICHMANNonly[['rnaDiscardedCounts']], cells = colnames(RHL_SeuratObject_nM_sel_TEICHMANNonly)) }

    RHL_SeuratObject_nM_sel_ROOIJonly = subset(RHL_SeuratObject_merged_noMito_sel, subset = annotation_paper_str == 'vRooij')
    RHL_SeuratObject_nM_sel_ROOIJonly[['rnaMitoCounts']] <- subset(RHL_SeuratObject_nM_sel_ROOIJonly[['rnaMitoCounts']], cells = colnames(RHL_SeuratObject_nM_sel_ROOIJonly))
    if (DISCARD_GENE_OPTION) { RHL_SeuratObject_nM_sel_ROOIJonly[['rnaDiscardedCounts']] <- subset(RHL_SeuratObject_nM_sel_ROOIJonly[['rnaDiscardedCounts']], cells = colnames(RHL_SeuratObject_nM_sel_ROOIJonly)) }

    # Save those
    print('Saving')
    SaveH5Seurat(object = RHL_SeuratObject_nM_sel_HUonly, overwrite = T,
        filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_HUonly.h5seurat'))
    SaveH5Seurat(object = RHL_SeuratObject_nM_sel_TEICHMANNonly, overwrite = T,
        filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_TEICHMANNonly.h5seurat'))
        # RHL_SeuratObject_nM_sel_TEICHMANNonly = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_TEICHMANNonly.h5seurat'))
    SaveH5Seurat(object = RHL_SeuratObject_nM_sel_ROOIJonly, overwrite = T,
        filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_ROOIJonly.h5seurat'))
    
    # Perhaps a more honest comparison between Rooij and Teichmann is to 
    # select for septal cells.
    # So let's also create a Teichmann version of the data which only contains its
    # septal cells
    # 
    # Subset
    RHL_SeuratObject_nM_sel_TEICHMANN.SP.only =
        subset(RHL_SeuratObject_nM_sel_TEICHMANNonly, annotation_region_str == 'SP')
    # Save
    SaveH5Seurat(object = RHL_SeuratObject_nM_sel_TEICHMANN.SP.only, overwrite = T,
        filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_TEICHMANN.SP.only.h5seurat'))
    
}


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

if ('run_separate' %in% desired_command) {
    
    current_dataset=config$dataset
    desired_settings=config$settings
    
    seuratObject_name = paste0('RHL_SeuratObject_nM_sel_',current_dataset,'only')
    
    # load the file
    current_analysis = list()
    current_analysis[[seuratObject_name]] =
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_',seuratObject_name,'.h5seurat'))
    
    # run the seurat analysis
    if (desired_settings == 'SETTINGS_RID2l') {
        CURRENT_RUNNAME = paste0(current_dataset,'only_','RID2l')
        # run RaceID2 like settings
        current_analysis[[seuratObject_name]] = mySeuratAnalysis(current_analysis[[seuratObject_name]],
            run_name=CURRENT_RUNNAME,
            normalization.method='RC', scale.factor=median(current_analysis[[seuratObject_name]]$nCount.nMT),
            do.scale=F,do.center=F,scale.max=Inf, features_to_use_choice = 'variable') # variable because otherwise too large calculation ..
    } else {
        CURRENT_RUNNAME = paste0(current_dataset,'only_','default')
        # run default
        current_analysis[[seuratObject_name]] = mySeuratAnalysis(current_analysis[[seuratObject_name]],
            run_name=CURRENT_RUNNAME)
    } 
    
    # create respective plots 
    mySeuratCommonPlots(current_analysis[[seuratObject_name]], run_name = CURRENT_RUNNAME)
    
    # save the file
    SaveH5Seurat(object = current_analysis[[seuratObject_name]], overwrite = T,
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
# More plots

# Code in:
# source(paste0(script_dir,'Functions/MW_customized_plot_functions.R'))

if ('more_custom_plots' %in% desired_command) {
    
    CURRENT_RUNNAME='all_RID2l_VAR'
    current_analysis = list()
    current_analysis[[CURRENT_RUNNAME]] = 
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_merged_noMito_sel.',CURRENT_RUNNAME,'.h5seurat'))
    
    # Add one-letter shorthand for dataset
    name_change = c("vRooij"="R", "Hu"="H", "Teichmann"="T")
    current_analysis[[CURRENT_RUNNAME]]$annotation_paper_oneletter_fct = factor(
        name_change[current_analysis[[CURRENT_RUNNAME]]$annotation_paper_str], levels=c('R','H','T'))
    # Create custom order for annotation
    current_analysis[[CURRENT_RUNNAME]]$annotation_paper_fct = 
        factor(current_analysis[[CURRENT_RUNNAME]]$annotation_paper_str, levels=c("vRooij", "Hu", "Teichmann"))
    
    # Stress gene plots
    LIST_NAME='stress'; CAT = 'annotation_paper_oneletter_fct'; CATNAME = 'Dataset'; PLOTTYPE='boxplot'
    genes_of_interest = c('NPPA', 'NPPB', 'ACTA1', 'MYH7')#, 'MYH6')
    plot_list = lapply(genes_of_interest, 
        function(x) {shorthand_plotViolinBox_custom(current_analysis,analysis_name=CURRENT_RUNNAME,cat_by = CAT, 
                                                    gene_of_interest=x,base_dir=base_dir, cat_name = CATNAME,
                                                    type = 'violin')})
    p=wrap_plots(plot_list, nrow = 1)& theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 5, unit = "pt"))
    #p
    ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',CURRENT_RUNNAME,'_6_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'.pdf'), width = min(18.64, 4*length(genes_of_interest)), height= 4, units='cm')

    # Now create MYH6 separately, also use to add legend
    p_=shorthand_plotViolinBox_custom(current_analysis,analysis_name=CURRENT_RUNNAME,cat_by = 'annotation_paper_fct', 
                                                    gene_of_interest='MYH6',base_dir=base_dir, cat_name = CATNAME,
                                                    type = 'violin')+theme(legend.position='right')
    ggsave(plot = p_,filename = paste0(base_dir, 'Rplots/',analysis_name,'_6_',type,'_FORLEGEND.pdf'), width = 7, height= 6, units='cm')
    
    
    
    
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
# Load a dataset manually

# H5_RHL_SeuratObject_nM_sel_TEICHMANNonly.h5seurat

# DATASET_NAME='TEICHMANNonly_RID2l'
# DATASET_NAME='ROOIJonly_RID2l'
# DATASET_NAME='ROOIJonly_Int1c'
# DATASET_NAME='HUonly_RID2l'
# DATASET_NAME='HUonly'
if (F) {
    
    # load the file
    current_analysis = list()
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
