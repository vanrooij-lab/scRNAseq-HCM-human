



################################################################################

#               THIS SCRIPT WAS NOT USED, I USED CLI VERSION

################################################################################









# Trying out SCENIC

################################################################################

# Note to self for HPC use
# Run with:
# commands="prep_corr_genie-patient=R.P1-cores=10"
# commands="run_SCENIC-patient=R.P1-cores=20"

################################################################################
# Load my own main script, to be able to use my Seurat data sets
#
# Note: this also sets "base_dir", I will output to 
# base_dir/SCENIC

# set script dir (note: overwritten by loading SeuratRevisedAnalysis below)
if (exists('LOCAL')) {
    script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Projects/'
} else {
    script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'
}

desired_command='dummy'
source(paste0(script_dir, 'HCM-SCS_2021-06_SeuratRevisedAnalysis_v2_UmiTools.R'))

# SCENIC-specific libraries

library(SCENIC)

################################################################################

if (!exists('desired_command_SCENIC')) {
    cmd=command_interpreter_mw(myargs)
    desired_command_SCENIC=cmd$desired_command
    config=cmd$config
}

# Set cores and memory consistent with demands; 
# Using 36 cores (medium node), MEM usage: 
# For GENIE3 run, this was ~13079004K --> 14G    -- GENIE3 run time doesn't seem to scale with increasing cores, but nCores is an argument that's used
#           --> ask for 10 cores, 20G
# For RcisTarget, mem goes up quickly, 120G currently, keeps increasing .. 160G 5 mins in .. 164G 7 mins in .. 170 G 10 mins in
#           --> ask for 20 cores, 200G
MYMCCORES = config$cores

################################################################################
# 
# # Installing/loading SCENIC
# 
# # I'm following the tutorial at
# # http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Setup.html

# Dependencies
# ===

# ## Required
# BiocManager::install(c("AUCell", "RcisTarget")) # INSTALLED
    # FOR HPC ---
    # @ HPC, both AUCell and RcisTarget weren't installed; dependency issues
    # Note: to install RSQLite, I had to first run
    # conda install -c conda-forge boost-cpp
    # Then simply install.packages("RSQLite")
    # Afterwards, BiocManager::install("annotate") did work (automatically installed AnnotationDbi also I think)
    # 
    # -----------
# BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost # INSTALLED

## Optional (but highly recommended):
# # To score the network on cells (i.e. run AUCell):
# BiocManager::install(c("zoo", "mixtools", "rbokeh")) # INSTALLED
# # For various visualizations and perform t-SNEs:
# BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne")) # INSTALLED
# # To support paralell execution (not available in Windows):
# BiocManager::install(c("doMC", "doRNG")) # NOT INSTALLED
# # To export/visualize in http://scope.aertslab.org
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools") # NOT INSTALLED
# devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE) # NOT INSTALLED

# SCENIC itself
# ===

# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("aertslab/SCENIC") 
# packageVersion("SCENIC")

# Execute in terminal to download motif databases:
# curl https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather --output hg19-500bp-upstream-7species.mc9nr.feather &
# curl https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather --output hg19-tss-centered-10kb-7species.mc9nr.feather
#                 # mc9nr: Motif collection version 9: 24k motifs

################################################################################
# Then run SCENIC
# 
# I'm following
# http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html


# Settings for which data set to load

DATASET_NAME='ROOIJonly_RID2l'
    # all_patients = c("R.P4", "R.P5", "R.P1", "R.P2", "R.P3")

CURRENT_PATIENT = config$patient

# Load the data
# ===

#CURRENT_PATIENT = 'R.P1'
#CURRENT_PATIENT = 'R.P2'
#CURRENT_PATIENT = 'R.P3'

if ('prep_corr_genie' %in% desired_command_SCENIC) {
    
    print(paste0('Running first part, for ', CURRENT_PATIENT, ', and using ', MYMCCORES, ' cores'))
    
    # Load the data
    if (!exists('current_analysis')) {current_analysis = list()}
    current_set = paste0(CURRENT_PATIENT, 'RID2l')
    current_analysis[[current_set]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',current_set,'.h5seurat'))
    
    # Prep output dirs
    if (!(dir.exists(paste0(base_dir, 'SCENIC/',CURRENT_PATIENT,'/int/')))) {dir.create(paste0(base_dir, 'SCENIC/',CURRENT_PATIENT,'/int/'), recursive = T)}
    if (!(dir.exists(paste0(base_dir, 'SCENIC/',CURRENT_PATIENT,'/output/')))) {dir.create(paste0(base_dir, 'SCENIC/',CURRENT_PATIENT,'/output/'), recursive = T)}
    if (!(dir.exists(paste0(base_dir, 'SCENIC/',CURRENT_PATIENT,'/temp_MW/')))) {dir.create(paste0(base_dir, 'SCENIC/',CURRENT_PATIENT,'/temp_MW/'), recursive = T)}
    
    
    # Retrieve data from Seurat object 
    #
    exprMat_ <- current_analysis[[current_set]]@assays$RNA@counts
    # Simply throw out genes that have duplicate or unknown hgnc names (NA probably not interesting, duplicate symbols are very few, seem uninteresting)
    rownames_exprMat_ = shorthand_cutname(rownames(exprMat_), PART1OR2 = 2)
    rownames(exprMat_) = rownames_exprMat_
    duplicate_genes = rownames_exprMat_[(duplicated(rownames_exprMat_))]
    gene_selection_beforehand = rownames_exprMat_[!(rownames_exprMat_ %in% duplicate_genes)]
    exprMat = exprMat_[gene_selection_beforehand,]
    rm('rownames_exprMat_') # remove temporary matrix
    #
    # Acquire cell names
    # 
    cellInfo <- colnames(exprMat)
    # To use Seurat clusters as cell annotation (e.g. for visualization):
    # cellInfo <- data.frame(seuratCluster=Idents(seuratObject))

    ################################################################################
    # Preprare SCENIC object
    
    # Settings
    org <- "hgnc" # mgi or hgnc, or dmel
    if (exists('LOCAL')) { dbDir <- "/Users/m.wehrens/Data_notbacked/_Resources_BIG_SCENIC/" 
        } else { dbDir = paste0(base_dir, 'SCENIC/DATABASES/') } # RcisTarget databases location 
    myDatasetTitle <- paste0("SCS_HCM_", CURRENT_PATIENT) # choose a name for your analysis
    data(defaultDbNames)
    dbs <- defaultDbNames[[org]]
    scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=MYMCCORES) 
    # Where to store output files
    if (!(any(grepl(base_dir, scenicOptions@fileNames$output)))) {
        scenicOptions@fileNames$output[,1] = paste0(base_dir, 'SCENIC/', CURRENT_PATIENT, '/', scenicOptions@fileNames$output) # [,1] to retain names
        scenicOptions@fileNames$int[,1]    = paste0(base_dir, 'SCENIC/', CURRENT_PATIENT, '/', scenicOptions@fileNames$int)
    }
    # Additional settings
    scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
    # scenicOptions@inputDatasetInfo$colVars <- (..) #optinoal
    # Databases (I think this is covered above?)
    # scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-5kb-10species.mc8nr.feather")
    # scenicOptions@settings$db_mcVersion <- "v8"
    
    print('Starting run')
    
    # Save to use at a later time...
    saveRDS(scenicOptions, file=paste0(base_dir, 'SCENIC/',CURRENT_PATIENT,'/int/scenicOptions.Rds')) 
    
    # (Adjust minimum values according to your dataset)
    # scenicOptions=scenicOptions; minCountsPerGene=3*.01*ncol(exprMat); minSamples=ncol(exprMat)*.01
    # rownames(exprMat) %in% names(motifRankings@rankings)
    genesKept <- geneFiltering(exprMat = as.matrix(exprMat), scenicOptions=scenicOptions,
                               minCountsPerGene=3*.01*ncol(exprMat),
                               minSamples=ncol(exprMat)*.01)    
        
    # Recommended step: check whether known interesting genes are kept
    interestingGenes <- c('SORBS2', 'EPAS1', 'TTN', 'CMYA5')
    interestingGenes[which(!interestingGenes %in% genesKept)] # any missing?
    
    # Filtered expression matrix
    exprMat_filtered <- exprMat[genesKept, ]
    dim(exprMat_filtered)
    # rm('exprMat') # keep it for later use
    
    print('Finding corr mat')
    
    # Now find correlations
    runCorrelation(as.matrix(exprMat_filtered), scenicOptions)

    ################################################################################
    
    
    # Optional: add log (if it is not logged/normalized already)
    exprMat_filtered_log <- log2(exprMat_filtered+1) 
    
    print('Running GENIE3')
    
    # Run GENIE3
    runGenie3(as.matrix(exprMat_filtered_log), scenicOptions)
    
    print('Saving config and matrix..')
    
    # Save some stuff for next run
    # (Can also be done with SCENIC commands, but I find this easier)
    save(list=c('exprMat_filtered_log', 'scenicOptions'), file=paste0(base_dir, 'SCENIC/',CURRENT_PATIENT,'/temp_MW/scenicRunData.Rdata'))
    
    print('Part 1 done')
    
}


################################################################################

# Now run SCENIC using wrapper functions
# 
# Outline:
# Build the gene regulatory network: 1. Get co-expression modules 2. Get regulons (with RcisTarget): TF motif analysis)

if ('run_SCENIC' %in% desired_command_SCENIC) {
    
    print(paste0('Running actual SCENIC, for ', CURRENT_PATIENT, ', and using ', MYMCCORES, ' cores'))
    
    # Load stuff from previous part
    load(file=paste0(base_dir, 'SCENIC/',CURRENT_PATIENT,'/temp_MW/scenicRunData.Rdata'))
    
    print(paste0('Matrix loaded; size=',dim(exprMat_filtered_log)[1], 'x',dim(exprMat_filtered_log)[2]))
    
    # Some additional options
    # scenicOptions <- readRDS("int/scenicOptions.Rds")
    scenicOptions@settings$verbose <- TRUE
    scenicOptions@settings$nCores <- MYMCCORES
    scenicOptions@settings$seed <- 123
    
    # Now run
    # print('SCENIC1 ============================================================')
    # scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
        # scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) #** Only for toy run!!
    
    # print('SCENIC2 ============================================================')
    # scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
    
    print('SCENIC3 ============================================================')
    scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, as.matrix(exprMat_filtered_log))
    
    print('SCENIC done ============================================================')
    saveRDS(scenicOptions, file=paste0(base_dir, 'SCENIC/',CURRENT_PATIENT,'/int/scenicOptions.Rds')) # To save status
    
    print('Saving done ============================================================')
}

################################################################################

# Optionally, 
#   - a tSNE can be build based on the regulon activity scores
#   - regulon activity can be binarized (on/off per cell)
# 
# See http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html
# "Binarize the network activity" and "Clustering / dimensionality reduction on the regulon activity"

# # BINARIZATION
# ===
# #
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
# savedSelections <- shiny::runApp(aucellApp)
# #
# # Save the modified thresholds:
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
# saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
# #
# # scenicOptions@settings$devType="png"
# scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

# # TSNE
# nPcs <- c(5) # For toy dataset
# # nPcs <- c(5,15,50)
# scenicOptions@settings$seed <- 123 # same seed for all of them
# # Run t-SNE with different settings:
# fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
# fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# # Plot as pdf (individual files in int/):
# fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
# # Note: The toy dataset only contains ~8 regulons; using more than 8 PCs will not provide any difference…
# 
# # and to view/compare them…
# 
# par(mfrow=c(length(nPcs), 3))
# fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
# plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

# # Using only "high-confidence" regulons (normally similar)
# par(mfrow=c(3,3))
# fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
# plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)
# The chosen t-SNE can then be saved as default to use for plots (can also be “binary”, see below):
# 
# scenicOptions@settings$defaultTsne$aucType <- "AUC"
# scenicOptions@settings$defaultTsne$dims <- 5
# scenicOptions@settings$defaultTsne$perpl <- 15
# saveRDS(scenicOptions, file="int/scenicOptions.Rds")

################################################################################

# Now process the results with custom code

if (F) {
    
    SCENIC_P1_table=
        read.table('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/SCENIC/R.P1/output/Step2_regulonTargetsInfo.tsv', header=1)
    
    table(SCENIC_P1_table$coexModule)
    
    unique(SCENIC_P1_table$TF)
    
    SCENIC_P1_table[SCENIC_P1_table$TF=='GATA4',]
    
    SCENIC_P1_table[SCENIC_P1_table$highConfAnnot,]
    
    table(SCENIC_P1_table[SCENIC_P1_table$highConfAnnot,])
          
    SCENIC_P1_table_motifs=
        read.table('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/SCENIC/R.P1/output/Step2_MotifEnrichment.tsv', header = 1, sep = '\t')
    
    unique(SCENIC_P1_table_motifs$highlightedTFs)

}

################################################################################

if (F) {
    
    SCENIC_tfModules = 
        readRDS('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/SCENIC/DOWNLOAD_FULL/R.P5/int/1.6_tfModules_asDF.Rds')
    
    View(SCENIC_tfModules)
    
    SCENIC_tfModules_top50=
        SCENIC_tfModules[SCENIC_tfModules$method=='top50perTarget',]
    
    unique(SCENIC_tfModules_top50$TF)

}

################################################################################

if (F) {
    
    scenicOptions_RP1 = 
        readRDS('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/SCENIC/DOWNLOAD_FULL/R.P1/int/scenicOptions.Rds')
    
    # Update directory names to local
    scenicOptions_RP1@fileNames$output[,1] =
        gsub(pattern = '/hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3/SCENIC/', 
             replacement = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/SCENIC/DOWNLOAD_FULL/',
             x = scenicOptions_RP1@fileNames$output)
    scenicOptions_RP1@fileNames$int[,1] =
        gsub(pattern = '/hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3/SCENIC/', 
             replacement = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/SCENIC/DOWNLOAD_FULL/',
             x = scenicOptions_RP1@fileNames$int)
    
    # 
    aucell_regulonAUC <- loadInt(scenicOptions_RP1, "aucell_regulonAUC")
    
}

################################################################################

if (F) {

    regulon_table = 
        readRDS('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/SCENIC/DOWNLOAD_FULL/R.P1/int/2.6_regulons_asGeneSet.Rds')

}















