
################################################################################

# Script to (partially manually) execute at the HPC

# Note:
# HPC command to get access to node manually
# srun --nodes=1 -c 1 --mem=50G --time=1:00:00 --pty bash -i

################################################################################

# Collect commands
myargs = commandArgs(trailingOnly = T)
print(paste0('myargs=',myargs))

# Activate code from the main script
script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

# Process the commands
if (length(myargs)==0) { 
    print('No commands given.') 
} else {
    cmd=command_interpreter_mw(myargs)
}

################################################################################
# Loading HCM and Ctrl2 (Teichmann) datasets separately

if (F) { # manually executed
    
    # Now load the data
    Teichmann_all <- LoadH5Seurat('/hpc/hub_oudenaarden/mwehrens/data/Teichmann/Litvinukova_global_raw.h5seurat')
    
    # We can also load our own data easily
    HCM_DATASET_NAME='ROOIJonly_RID2l_clExtended'
    if (!exists('current_analysis')) {current_analysis = list()}
    current_analysis[[HCM_DATASET_NAME]] =
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',HCM_DATASET_NAME,'.h5seurat'))
}


################################################################################
# First update the Ctrl Teichmann dataset

if (F) { # manually executed
    
    # First, like in main script, perform some adjustments on the Teichmann data
    CalcN_out = Seurat:::CalcN(Teichmann_all)
    Teichmann_all[['nCount_RNA']]   = CalcN_out$nCount
    Teichmann_all[['nFeature_RNA']] = CalcN_out$nFeature
        
    # Update gene names to be consistent with everything
    ensembl_IDs = Teichmann_all@assays$RNA@meta.features$'gene_ids-Harvard-Nuclei'
    names(ensembl_IDs) = rownames(Teichmann_all@assays$RNA@meta.features)
        # Note that all(Teichmann_all@assays$RNA@meta.features$'gene_ids-Harvard-Nuclei' == Teichmann_all@assays$RNA@meta.features$'gene_ids-Sanger-Nuclei')
    
    ensembl_IDs_check = ensembl_IDs[rownames(Teichmann_all)]
        # this is redundant, since naming is consistent, i.e.
        # all(ensembl_IDs_check == ensembl_IDs)
        # still nice to do sanity check
    
    # actually update names    
    hgnc_names_Teich=rownames(Teichmann_all@assays$RNA@data)
    rownames(Teichmann_all@assays$RNA@data) = paste0(ensembl_IDs_check,':',hgnc_names_Teich)
    rownames(Teichmann_all@assays$RNA@counts) = paste0(ensembl_IDs_check,':',hgnc_names_Teich)
    
    # now save this updated data frame
    SaveH5Seurat(object = Teichmann_all, filename = '/hpc/hub_oudenaarden/mwehrens/data/Teichmann/Litvinukova_global_raw_MW.h5seurat')

}
    
################################################################################
# Merge the datasets

if (F) { # manually executed
    
    # Teichmann_all = LoadH5Seurat(file = '/hpc/hub_oudenaarden/mwehrens/data/Teichmann/Litvinukova_global_raw_MW.h5seurat')
    # HCM_DATASET_NAME='ROOIJonly_RID2l_clExtended'; current_analysis = list(); current_analysis[[HCM_DATASET_NAME]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',HCM_DATASET_NAME,'.h5seurat'))
    
    # Perform the merge
    dataset_TeichAll_Rooij_merged <- merge(x = Teichmann_all, y = current_analysis[[HCM_DATASET_NAME]], add.cell.ids = c('T','R'), project = "LigRec")
    object_size(RHL_SeuratObject_merged)
    
    # Save the raw merged dataset
    SaveH5Seurat(object = dataset_TeichAll_Rooij_merged, file = '/hpc/hub_oudenaarden/mwehrens/data/Teichmann/dataset_TeichAll_Rooij_merged_raw.h5seurat')
    
    # Then add annotation
    
    # Add sample information (patients in the case of Hu)
    # This just extracts the prefix added by the merge function 
    TaR_annotation_prelim=unlist(lapply(  strsplit(colnames(dataset_TeichAll_Rooij_merged), split='_'), 
        function(splitted_string) {splitted_string[1]}))
    
    # Sanity check for Teichmann
    # if (any(paste0('TEICH_',colnames(RHL_SeuratObject_list$TEICH))!=colnames(RHL_SeuratObject_merged)[TaR_annotation_prelim=='TEICH'])) {
    #     print('Annotation inconsistency.'); stop('Annotation inconsistency.')
    #     } else {print('Annotation Teichmann consistent.')}
    
    
    # Define sample annotation
    TaR_annotation_samples=rep(NA, length(TaR_annotation_prelim))
    TaR_annotation_samples[TaR_annotation_prelim=='T']=paste0('T.',dataset_TeichAll_Rooij_merged@meta.data$sample[TaR_annotation_prelim=='T'])
    TaR_annotation_samples[TaR_annotation_prelim=='R']=paste0('R.',current_analysis$ROOIJonly_RID2l_clExtended$annotation_sample_str)
    # unique(TaR_annotation_samples)
    # Determine annotation for from which paper a cell comes
    TaR_annotation_paper = rep(NA, length(TaR_annotation_prelim))
    TaR_annotation_paper[TaR_annotation_prelim=='T']='Teichmann'
    TaR_annotation_paper[TaR_annotation_prelim=='R']='vRooij'
    # Same but prettier annotation
    TaR_annotation_paper_b = rep(NA, length(TaR_annotation_prelim))
    TaR_annotation_paper_b[TaR_annotation_prelim=='T']='Ctrl'
    TaR_annotation_paper_b[TaR_annotation_prelim=='R']='HCM'
    # Get patient information for our data and Teichmann's
    TaR_annotation_patients=rep(NA, length(TaR_annotation_prelim))
    TaR_annotation_patients[TaR_annotation_prelim=='T'] = paste0('T.',dataset_TeichAll_Rooij_merged@meta.data$donor)[TaR_annotation_prelim=='T']
    TaR_annotation_patients[TaR_annotation_prelim=='R'] = current_analysis$ROOIJonly_RID2l_clExtended$annotation_patient_str
    # Now also retrieve the "region" data for the Teichmann data, which correspond for the vCMs to AX, SP, LV, RV (apex, septum, left ventricle, right v)
    # For our data, we only have septal cells
    TaR_annotation_region=TaR_annotation_prelim
    TaR_annotation_region[TaR_annotation_prelim=='R'] = 'SP'
    # For Teichmann, we need the region data
    TaR_annotation_region[TaR_annotation_prelim=='T'] = as.character(dataset_TeichAll_Rooij_merged@meta.data$region[TaR_annotation_prelim=='T'])
    
    # Add information to Seurat object
    dataset_TeichAll_Rooij_merged[["annotation_sample_fct"]]     <- as.factor(TaR_annotation_samples)
    dataset_TeichAll_Rooij_merged[["annotation_sample_str"]]     <- TaR_annotation_samples
    dataset_TeichAll_Rooij_merged[["annotation_patient_fct"]]     <- as.factor(TaR_annotation_patients)
    dataset_TeichAll_Rooij_merged[["annotation_patient_str"]]     <- TaR_annotation_patients
    dataset_TeichAll_Rooij_merged[["annotation_paper_fct"]]     <- as.factor(TaR_annotation_paper)
    dataset_TeichAll_Rooij_merged[["annotation_paper_str"]]     <- TaR_annotation_paper
    dataset_TeichAll_Rooij_merged[["annotation_region_fct"]]     <- as.factor(TaR_annotation_region)
    dataset_TeichAll_Rooij_merged[["annotation_region_str"]]     <- TaR_annotation_region
    
    # Add celltype with our cells annotated (this was added later, so not saved)
    dataset_TeichAll_Rooij_merged_noMt_sel$celltype_inclR = dataset_TeichAll_Rooij_merged_noMt_sel$cell_type
    dataset_TeichAll_Rooij_merged_noMt_sel$celltype_inclR[TaR_annotation_prelim=='R'] = 'hCM_Rooij'
    
    # Count mitochondrial reads
    # checking whether we take the right genes:
    # rownames(dataset_TeichAll_Rooij_merged)[grepl("^MT-",rownames(dataset_TeichAll_Rooij_merged))]
    dataset_TeichAll_Rooij_merged[["percent.mt"]] <- PercentageFeatureSet(dataset_TeichAll_Rooij_merged, pattern = ":MT-")
    
    # Now let's handle mito reads
    mito_genes = rownames(dataset_TeichAll_Rooij_merged)[grepl(':MT-',rownames(dataset_TeichAll_Rooij_merged))]
    # Again save the non-desired data
    dataset_TeichAll_Rooij_merged[['rnaMitoCounts']] <- CreateAssayObject(counts = dataset_TeichAll_Rooij_merged@assays$RNA@counts[mito_genes,])
    
    # Save the updated file
    SaveH5Seurat(object = dataset_TeichAll_Rooij_merged, file = '/hpc/hub_oudenaarden/mwehrens/data/Teichmann/dataset_TeichAll_Rooij_merged_annotated.h5seurat')
    
    # Then, subset the RNA assay (note: "noMito" is historic name, earlier only mito genes removed, now many more)
    dataset_TeichAll_Rooij_merged_noMt =
        subset(dataset_TeichAll_Rooij_merged, features= dataset_TeichAll_Rooij_merged@misc$desired_genes_excl_mito)
    # And add the desired stored information again
    dataset_TeichAll_Rooij_merged_noMt[['rnaMitoCounts']] = dataset_TeichAll_Rooij_merged[['rnaMitoCounts']]
        # note: to look up; don't quite understand why subsetting can't be applied to all assays at once in Seurat
    
    # New count stats
    CalcN_out_merged = Seurat:::CalcN(dataset_TeichAll_Rooij_merged_noMt)
        dataset_TeichAll_Rooij_merged_noMt[['nCount.nMT']]   = CalcN_out_merged$nCount
        dataset_TeichAll_Rooij_merged_noMt[['nFeature.nMT']] = CalcN_out_merged$nFeature
        
    # Let's select cells with 1000 reads total
    dataset_TeichAll_Rooij_merged_noMt_sel <- subset(dataset_TeichAll_Rooij_merged_noMt, subset = nCount.nMT > 1000)    
        # # Let's check some cell count stats
        # any(dataset_TeichAll_Rooij_merged_noMt$nCount.nMT <= 1000)
        # min(dataset_TeichAll_Rooij_merged_noMt$nCount.nMT[TaR_annotation_prelim=='R'])
        # min(dataset_TeichAll_Rooij_merged_noMt$nCount.nMT[TaR_annotation_prelim=='T'])
        # lowcells=colnames(dataset_TeichAll_Rooij_merged_noMt)[dataset_TeichAll_Rooij_merged_noMt$nCount.nMT==500]
        # colSums(dataset_TeichAll_Rooij_merged_noMt@assays$RNA@counts[,lowcells])
    # Also for "discarded" data
    dataset_TeichAll_Rooij_merged_noMt_sel[['rnaMitoCounts']] <- subset(dataset_TeichAll_Rooij_merged_noMt[['rnaMitoCounts']], cells = colnames(dataset_TeichAll_Rooij_merged_noMt_sel))
    
    # Now save this data frame
    SaveH5Seurat(object = dataset_TeichAll_Rooij_merged_noMt_sel, file = '/hpc/hub_oudenaarden/mwehrens/data/Teichmann/dataset_TeichAll_Rooij_merged_noMt_sel.h5seurat')
    
}


################################################################################
# Now perform first a "standard" Seurat analysis


if ('run_LigRec_stAnalysis' %in% cmd$desired_command) {
    
    CURRENT_RUNNAME='TeichAll_Rooij_dflt'
    
    dataset_TeichAll_Rooij_merged_noMt_sel = LoadH5Seurat(file = paste0(base_dir,'LR_analysis/dataset_TeichAll_Rooij_merged_noMt_sel.h5seurat'))

    # Analysis    
    dataset_TeichAll_Rooij_merged_noMt_sel_dflt = mySeuratAnalysis(dataset_TeichAll_Rooij_merged_noMt_sel, run_name = CURRENT_RUNNAME, subdir='LR_analysis/')
    
    # Save the analysis
    SaveH5Seurat(object = dataset_TeichAll_Rooij_merged_noMt_sel_dflt, overwrite = T,
        filename = paste0(base_dir,'LR_analysis/',CURRENT_RUNNAME,'.h5seurat'))
    
    # Create plots
    mySeuratCommonPlots(dataset_TeichAll_Rooij_merged_noMt_sel_dflt, run_name = CURRENT_RUNNAME, subdir='LR_analysis/')
    
}    

################################################################################
# Boiler plate fix to add cell types for us

# TaR_annotation_prelim=unlist(lapply(  strsplit(colnames(current_analysis$TeichAll_Rooij_dflt), split='_'), 
#         function(splitted_string) {splitted_string[1]}))
#     
# current_analysis$TeichAll_Rooij_dflt$celltype_inclR = current_analysis$TeichAll_Rooij_dflt$cell_type
# current_analysis$TeichAll_Rooij_dflt$celltype_inclR[TaR_annotation_prelim=='R'] = 'hCM_Rooij'
#
# 
# SaveH5Seurat(object = current_analysis[[CURRENT_RUNNAME]], overwrite = T, filename = paste0(base_dir,'LR_analysis/',CURRENT_RUNNAME,'.h5seurat'))
    

################################################################################

if ('LigRec_customMaps' %in% cmd$desired_command) {

# 
 
    # Load dataset
    CURRENT_RUNNAME='TeichAll_Rooij_dflt'
    current_analysis=list()
    current_analysis[[CURRENT_RUNNAME]] = 
        dataset_TeichAll_Rooij_merged_noMt_sel = LoadH5Seurat(file = paste0(base_dir, 'LR_analysis/', CURRENT_RUNNAME, '.h5seurat'))
    
    
    #shorthand_seurat_custom_expr(seuratObject = current_analysis[[CURRENT_RUNNAME]], 
    #                             gene_of_interest = enriched_genes_lists_clusters_ROOIJ_[[cl_name]],
    #                             textsize = 6*ZOOM_FACTOR, pointsize = .5, custom_title = cl_name, mymargin = .5*ZOOM_FACTOR, zscore = T) 
    
    custom_colors = c('#bd0020','#9d9d9c','#575756')
    
    # Another overview plot with paper
    p=DimPlot(current_analysis[[CURRENT_RUNNAME]], group.by = 'annotation_paper_fct', cols = col_vector_60, label = T, repel = T, label.size = 7/.pt, pt.size = NULL)+
        theme_void()+ggtitle(element_blank())+theme(legend.position = 'none')
    # test
    ggsave(plot = p, filename = paste0(base_dir,'LR_analysis/',CURRENT_RUNNAME,'_sourcePaper.png'), height=172/3*2-4, width=172/3*2-4, units='mm')
    
    # Another overview plot with patients
    p=DimPlot(current_analysis[[CURRENT_RUNNAME]], group.by = 'annotation_patient_fct', cols = col_vector_60, label = T, repel = T, label.size = 7/.pt, pt.size = NULL)+
        theme_void()+ggtitle(element_blank())+theme(legend.position = 'none')
    # test
    ggsave(plot = p, filename = paste0(base_dir,'LR_analysis/',CURRENT_RUNNAME,'_sourcePatient.png'), height=172/3*2-4, width=172/3*2-4, units='mm')
    
    
    # Make a plot per cell type
    p=DimPlot(current_analysis[[CURRENT_RUNNAME]], group.by = 'cell_type_inclR', cols = col_vector_60, label = T, repel = T, label.size = 7/.pt, pt.size = NULL)+
        theme_void()+ggtitle(element_blank())+theme(legend.position = 'none')
    # test
    ggsave(plot = p, filename = paste0(base_dir,'LR_analysis/',CURRENT_RUNNAME,'_celltypes.png'), height=172/3*2-4, width=172/3*2-4, units='mm')
    
    # Make a plot per cell type (with legend)
    p=DimPlot(current_analysis[[CURRENT_RUNNAME]], group.by = 'cell_type_inclR', cols = col_vector_60, label = T, repel = T, label.size = 7/.pt, pt.size = NULL)+
        theme_void()+ggtitle(element_blank())
    # test
    ggsave(plot = p, filename = paste0(base_dir,'LR_analysis/',CURRENT_RUNNAME,'_celltypes_wLegend.png'), height=172/3*2-4, width=172/3*3-4, units='mm')
       
}

################################################################################
# Now gather info about receptor and ligand expression in the large 
# datasets

if ('LigRec_Expr' %in% cmd$desired_command) {
    
    # Load dataset
    CURRENT_RUNNAME='TeichAll_Rooij_dflt'
    current_analysis=list()
    current_analysis[[CURRENT_RUNNAME]] = 
        dataset_TeichAll_Rooij_merged_noMt_sel = LoadH5Seurat(file = paste0(base_dir, CURRENT_RUNNAME, '.h5seurat'))
    
    # Load L-R of interest
    load(paste0(base_dir,'LR_analysis/LR_of_interest.Rdata'))
        # ligands_of_interest
        # receptors_of_interest
        
    # Now create convenient overview
    
    # We'll have to know cell types at some point, are they in the data already?
    
    shorthand_custom_compositeplot(seuratObject_list = current_analysis, gene_lists = ligands_of_interest, seuratObjectNameToTake = CURRENT_RUNNAME, 
                                   group.by = )
    
}










