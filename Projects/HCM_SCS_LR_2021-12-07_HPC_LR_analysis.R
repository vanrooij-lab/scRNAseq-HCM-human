
################################################################################

# Script to (partially manually) execute at the HPC

# Note:
# HPC command to get access to node manually
# srun --nodes=1 -c 1 --mem=100G --time=4:00:00 --pty bash -i

################################################################################

# Collect commands
myargs = commandArgs(trailingOnly = T)
print(paste0('myargs=',myargs))

# Activate code from the main script
script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')
base_dir_original = base_dir
base_dir = '/hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_allcells_TR/'

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
    #HCM_DATASET_NAME='ROOIJonly_RID2l_clExtended'
    HCM_DATASET_NAME='ROOIJonly.sp.bt_RID2l_clExtended'
    if (!exists('current_analysis')) {current_analysis = list()}
    current_analysis[[HCM_DATASET_NAME]] =
        LoadH5Seurat(file = paste0(base_dir_original,'Rdata/H5_RHL_SeuratObject_nM_sel_',HCM_DATASET_NAME,'.h5seurat'))
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
    SaveH5Seurat(object = Teichmann_all, filename = '/hpc/hub_oudenaarden/mwehrens/data/Teichmann/temp_Litvinukova_global_raw_MW.h5seurat')
        # Teichmann_all = LoadH5Seurat(file = '/hpc/hub_oudenaarden/mwehrens/data/Teichmann/temp_Litvinukova_global_raw_MW.h5seurat')
    
    # Note: there was a bug earlier where I used the wrong hcgn version names for gene symbols; but this is fixed now
    # I performed a quick sanity check by checking the gene "FAM129A" aka "NIBAN1", which should in both datasets be referred to as FAM129A.

}
    
################################################################################
# Merge the datasets

if (F) { # manually executed
    
    # Teichmann_all = LoadH5Seurat(file = '/hpc/hub_oudenaarden/mwehrens/data/Teichmann/Litvinukova_global_raw_MW.h5seurat')
    # HCM_DATASET_NAME='ROOIJonly_RID2l_clExtended'; current_analysis = list(); current_analysis[[HCM_DATASET_NAME]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',HCM_DATASET_NAME,'.h5seurat'))
    
    # Perform the merge
    dataset_TeichAll_Rooij_merged <- merge(x = Teichmann_all, y = current_analysis[[HCM_DATASET_NAME]], add.cell.ids = c('T','R'), project = "LigRec")
    object_size(dataset_TeichAll_Rooij_merged)
    
    # Save the raw merged dataset (since this script requires loads of memory, convenient to save intermediate state now and then)
    SaveH5Seurat(object = dataset_TeichAll_Rooij_merged, file = '/hpc/hub_oudenaarden/mwehrens/data/Teichmann/temp_dataset_TeichAll_Rooij_merged_raw.h5seurat')
        # file.exists('/hpc/hub_oudenaarden/mwehrens/data/Teichmann/temp_dataset_TeichAll_Rooij_merged_raw.h5seurat')
    
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
    TaR_annotation_samples[TaR_annotation_prelim=='R']=paste0('R.',current_analysis[[HCM_DATASET_NAME]]$annotation_sample_str)
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
    TaR_annotation_patients[TaR_annotation_prelim=='R'] = current_analysis[[HCM_DATASET_NAME]]$annotation_patient_str
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
    dataset_TeichAll_Rooij_merged$celltype_inclR = dataset_TeichAll_Rooij_merged$cell_type
    dataset_TeichAll_Rooij_merged$celltype_inclR[TaR_annotation_prelim=='R'] = 'hCM_Rooij'
    
    # Count mitochondrial reads
    # checking whether we take the right genes:
    # rownames(dataset_TeichAll_Rooij_merged)[grepl("^MT-",rownames(dataset_TeichAll_Rooij_merged))]
    dataset_TeichAll_Rooij_merged[["percent.mt"]] <- PercentageFeatureSet(dataset_TeichAll_Rooij_merged, pattern = ":MT-")
    
    # Now let's handle mito reads
    mito_genes = rownames(dataset_TeichAll_Rooij_merged)[grepl(':MT-',rownames(dataset_TeichAll_Rooij_merged))]
    # Again save the non-desired data
    dataset_TeichAll_Rooij_merged[['rnaMitoCounts']] <- CreateAssayObject(counts = dataset_TeichAll_Rooij_merged@assays$RNA@counts[mito_genes,])
    
    # Save the updated file
    SaveH5Seurat(object = dataset_TeichAll_Rooij_merged, file = '/hpc/hub_oudenaarden/mwehrens/data/Teichmann/temp_dataset_TeichAll_Rooij_merged_annotated.h5seurat')
    
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

# First move "dataset_TeichAll_Rooij_merged_noMt_sel.h5seurat" from Teichmann folder to 

if ('run_LigRec_stAnalysis' %in% cmd$desired_command) {
    
    CURRENT_RUNNAME='TeichAll_Rooij_dflt'
    
    dataset_TeichAll_Rooij_merged_noMt_sel = LoadH5Seurat(file = paste0(base_dir,'Rdata/dataset_TeichAll_Rooij_merged_noMt_sel.h5seurat'))

    # Analysis    
    dataset_TeichAll_Rooij_merged_noMt_sel_dflt = mySeuratAnalysis(dataset_TeichAll_Rooij_merged_noMt_sel, run_name = CURRENT_RUNNAME, subdir='Rplots/')
    
    # Save the analysis
    SaveH5Seurat(object = dataset_TeichAll_Rooij_merged_noMt_sel_dflt, overwrite = T,
        filename = paste0(base_dir,'Rdata/',CURRENT_RUNNAME,'.h5seurat'))
        # dataset_TeichAll_Rooij_merged_noMt_sel_dflt = LoadH5Seurat(file = paste0(base_dir,'Rdata/',CURRENT_RUNNAME,'.h5seurat'))
    
    # Create plots
    mySeuratCommonPlots(dataset_TeichAll_Rooij_merged_noMt_sel_dflt, run_name = CURRENT_RUNNAME, subdir='Rplots/')
    
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
# SaveH5Seurat(object = current_analysis[[CURRENT_RUNNAME]], overwrite = T, filename = paste0(base_dir,'Rdata/',CURRENT_RUNNAME,'.h5seurat'))
    

################################################################################

if ('LigRec_customMaps' %in% cmd$desired_command) {

# 
 
    # Load dataset
    CURRENT_RUNNAME='TeichAll_Rooij_dflt'
    current_analysis=list()
    current_analysis[[CURRENT_RUNNAME]] = LoadH5Seurat(file = paste0(base_dir, 'Rdata/', CURRENT_RUNNAME, '.h5seurat'))
    
    
    #shorthand_seurat_custom_expr(seuratObject = current_analysis[[CURRENT_RUNNAME]], 
    #                             gene_of_interest = enriched_genes_lists_clusters_ROOIJ_[[cl_name]],
    #                             textsize = 6*ZOOM_FACTOR, pointsize = .5, custom_title = cl_name, mymargin = .5*ZOOM_FACTOR, zscore = T) 
    
    custom_colors = c('#bd0020','#9d9d9c','#575756')
    
    # Another overview plot with paper
    p=DimPlot(current_analysis[[CURRENT_RUNNAME]], group.by = 'annotation_paper_fct', cols = col_vector_60, label = T, repel = T, label.size = 7/.pt, pt.size = NULL)+
        theme_void()+ggtitle(element_blank())+theme(legend.position = 'none')
    # test
    ggsave(plot = p, filename = paste0(base_dir,'Rplots/',CURRENT_RUNNAME,'_sourcePaper.png'), height=172/3*2-4, width=172/3*2-4, units='mm')
    
    # Another overview plot with patients
    p=DimPlot(current_analysis[[CURRENT_RUNNAME]], group.by = 'annotation_patient_fct', cols = col_vector_60, label = T, repel = T, label.size = 7/.pt, pt.size = NULL)+
        theme_void()+ggtitle(element_blank())+theme(legend.position = 'none')
    # test
    ggsave(plot = p, filename = paste0(base_dir,'Rplots/',CURRENT_RUNNAME,'_sourcePatient.png'), height=172/3*2-4, width=172/3*2-4, units='mm')
    
    
    # Make a plot per cell type
    p=DimPlot(current_analysis[[CURRENT_RUNNAME]], group.by = 'cell_type_inclR', cols = col_vector_60, label = T, repel = T, label.size = 7/.pt, pt.size = NULL)+
        theme_void()+ggtitle(element_blank())+theme(legend.position = 'none')
    # test
    ggsave(plot = p, filename = paste0(base_dir,'Rplots/',CURRENT_RUNNAME,'_celltypes.png'), height=172/3*2-4, width=172/3*2-4, units='mm')
    
    # Make a plot per cell type (with legend)
    p=DimPlot(current_analysis[[CURRENT_RUNNAME]], group.by = 'cell_type_inclR', cols = col_vector_60, label = T, repel = T, label.size = 7/.pt, pt.size = NULL)+
        theme_void()+ggtitle(element_blank())
    # test
    ggsave(plot = p, filename = paste0(base_dir,'Rplots/',CURRENT_RUNNAME,'_celltypes_wLegend.png'), height=172/3*2-4, width=172/3*3-4, units='mm')
       
}

################################################################################
# Now gather info about receptor and ligand expression in the large 
# datasets

if ('LigRec_Expr' %in% cmd$desired_command) {
    
    ################################################################################
    # Load & set up
    
    # Load dataset
    date()
    DATASET_NAME ='TeichAll_Rooij_dflt' # DATASET_NAME = 'ROOIJonly_RID2l_clExtended' # this is only convenient when doing some manual actions
    DATASET_NAME2='TeichAll_Rooij_dflt' 
    current_analysis=list()
    current_analysis[[DATASET_NAME]] = 
        LoadH5Seurat(file = paste0(base_dir, 'Rdata/', DATASET_NAME, '.h5seurat'))
    date()
    
    # Load L-R of interest
    load(paste0(base_dir,'LR_analysis/LR_of_interest.Rdata'))
        # ligands_of_interest
        # receptors_of_interest
    
    # Load L-R database
    source(paste0(base_dir, 'LR_analysis/human_LR_list.R'))
    # Local: source('/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Resources/human_LR_list.R')
    # Local: file.edit('/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Resources/human_LR_list.R')
    
    # Gather averages per cell type, per patient
    group.by1='annotation_patient_str'
    group.by2='cell_type_inclR'
    
    # Set up 1st list of interest
    feature_list = unlist(receptors_of_interest$HCM_up)
    feature_list_annotation = rep(names(receptors_of_interest$HCM_up), times=sapply(receptors_of_interest$HCM_up, length))
    HCM_feature_list = feature_list; HCM_feature_list_annotation = feature_list_annotation # for later use
    
    ################################################################################
    
    # Determine first the reference value, which is simply 
    # the max-value for the full expression matrix
    #
    # This is done by using the "only_return_max" option
    # Note that this approach is highly inefficient
    # TO DO: More convenient to cut it up within the function, if 
    # desired .. Makes function rather complicated though ..
    print(paste0('Calculating limits -- ',date()))
    # manual_zlims = unlist(shorthand_heatmap_feature_aggr(current_analysis=current_analysis, analysis=DATASET_NAME, 
    #                                    feature_list=feature_list, feature_list_annotation=feature_list_annotation, 
    #                                    group.by1=group.by1, group.by2=group.by2, fontsize=3, only_return_max = T))
    listp = shorthand_heatmap_feature_aggr(current_analysis=current_analysis, analysis=DATASET_NAME, 
                                        feature_list=feature_list, feature_list_annotation=feature_list_annotation, 
                                        group.by1=group.by1, group.by2=group.by2, fontsize=3, savepath = paste0(base_dir, 'Rplots/',DATASET_NAME,'__HCM_up__heatmap_feature_aggr.Rdata'))
    load(paste0(base_dir, 'Rplots/',DATASET_NAME,'__HCM_up__heatmap_feature_aggr.Rdata')) # data_LR_expression
    print(paste0('Done -- ',date()))
    
    # Functions to quickly look up related ligands
    quicky_give_L_for_R_all = function(R) {
        ligand_receptor_pairs_df[ligand_receptor_pairs_df$receptor==R,]$ligand }
    quicky_give_L_for_R_OI = function(R) {
        all_ligands_with_R = quicky_give_L_for_R_all(R);  all_ligands_with_R[all_ligands_with_R %in% ligands_of_interest$HCM_up]}
    
    ################################################################################
    # Now a second time, without redunancy with ligands
    
    feature_list_Ronly = unique(unlist(receptors_of_interest$HCM_up))
    feature_list_annotation_Ronly = rep('_', length(feature_list_Ronly))

    listp_Ronly = shorthand_heatmap_feature_aggr(current_analysis=current_analysis, analysis=DATASET_NAME, 
                                        feature_list=feature_list_Ronly, feature_list_annotation=feature_list_annotation_Ronly, 
                                        group.by1=group.by1, group.by2=group.by2, fontsize=3, 
                                        savepath = paste0(base_dir, 'Rplots/',DATASET_NAME,'__HCM_up__heatmap_feature_aggr_Ronly.Rdata'))
    load(paste0(base_dir, 'Rplots/',DATASET_NAME,'__HCM_up__heatmap_feature_aggr_Ronly.Rdata')) # data_LR_expression
    print(paste0('Done -- ',date()))
    
    
    
    all_ligands_with_R = quicky_give_L_for_R('EGFR')
    ligandsOfInterest_with_R = all_ligands_with_R[all_ligands_with_R %in% ligands_of_interest$HCM_up]
    
    ################################################################################
    # Now for ligands
    
    feature_list_Lonly = ligands_of_interest$HCM_up
    feature_list_annotation_Lonly = rep('_', length(feature_list_Lonly))

    listp_Ronly = shorthand_heatmap_feature_aggr(current_analysis=current_analysis, analysis=DATASET_NAME2, 
                                        feature_list=feature_list_Lonly, feature_list_annotation=feature_list_annotation_Lonly, 
                                        group.by1=group.by1, group.by2=group.by2, fontsize=3, 
                                        savepath = paste0(base_dir, 'Rplots/',DATASET_NAME2,'__HCM_up__heatmap_feature_aggr_Lonly.Rdata'))
    load(paste0(base_dir, 'Rplots/',DATASET_NAME2,'__HCM_up__heatmap_feature_aggr_Lonly.Rdata')) # data_LR_expression
    print(paste0('Done -- ',date()))
    
    ################################################################################
    
    # For clarity, just create manual split of this list
    if(length(receptors_of_interest$HCM_up)!=44){stop('wrong assumptions!')}
    HCM_up_sampes = list(
        HCM_up1 = receptors_of_interest$HCM_up[1:10],
        HCM_up2 = receptors_of_interest$HCM_up[11:20],
        HCM_up3 = receptors_of_interest$HCM_up[21:30],
        HCM_up4 = receptors_of_interest$HCM_up[31:40],
        HCM_up5 = receptors_of_interest$HCM_up[41:44])
    
    for (HCM_UP_smpl_name in names(HCM_up_sampes)) {
        
        print(paste0('Sample ', HCM_UP_smpl_name,' -- ',date()))
        
        current_list = HCM_up_sampes[[HCM_UP_smpl_name]]
        
        feature_list = unlist(current_list)
        feature_list_annotation = rep(names(current_list), times=sapply(current_list, length))
        
        listp=shorthand_heatmap_feature_aggr(current_analysis=current_analysis, analysis=DATASET_NAME, 
                                       feature_list=feature_list, feature_list_annotation=feature_list_annotation, 
                                       group.by1=group.by1, group.by2=group.by2, fontsize=3, 
                                       manual_zlims = manual_zlims)
        
        # Save Z-score plots
        ggsave(plot = listp$p1_Z, filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_LR_Intrx_Z_',HCM_UP_smpl_name,'_extended.png'), height=PANEL_HEIGHT*3-4, width=172-4, units='mm')
        ggsave(plot = listp$p2_Z, filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_LR_Intrx_Z_',HCM_UP_smpl_name,'_collapsed.png'), height=PANEL_HEIGHT*3-4, width=172-4, units='mm')
        
        # Save expr score plots
        ggsave(plot = listp$p1_x, filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_LR_Intrx_expr_',HCM_UP_smpl_name,'_extended.png'), height=PANEL_HEIGHT*3-4, width=172-4, units='mm')
        ggsave(plot = listp$p2_x, filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_LR_Intrx_expr_',HCM_UP_smpl_name,'_collapsed.png'), height=PANEL_HEIGHT*3-4, width=172-4, units='mm')
        
        # Save fraction expressed heatmaps
        ggsave(plot = listp$p1_f, filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_LR_Intrx_fracCells_',HCM_UP_smpl_name,'_extended.png'), height=PANEL_HEIGHT*3-4, width=172-4, units='mm')
        ggsave(plot = listp$p2_f, filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_LR_Intrx_fracCells_',HCM_UP_smpl_name,'_collapsed.png'), height=PANEL_HEIGHT*3-4, width=172-4, units='mm')
        
    }
    
    ################################################################################
    

    ################################################################################
    
    load(paste0(base_dir, 'Rplots/',DATASET_NAME2,'__HCM_up__heatmap_feature_aggr_Ronly.Rdata')) # data_LR_expression
    load(paste0(base_dir, 'Rplots/',DATASET_NAME2,'__HCM_up__heatmap_feature_aggr_Lonly.Rdata')) # data_LR_expression
    
    data_LR_expression$matrix_expr_x_g2collapse[is.na(data_LR_expression$matrix_expr_x_g2collapse)]=-1
    pheatmap(data_LR_expression$matrix_expr_x_g2collapse)
    
    # Distribution of expression
    all_expr = as.vector(data_LR_expression$matrix_expr_x_g2collapse)
    ggplot(data.frame(all_expr=all_expr), aes(x=all_expr))+
        geom_freqpoly(bins=100)+theme_bw()
    ggplot(data.frame(all_expr=all_expr[all_expr>0]), aes(x=all_expr))+
        geom_freqpoly(bins=100)+theme_bw()
    
    # Determine max. expression
    max_R_expr_forEach = apply(data_LR_expression$matrix_expr_x, 2, max)
    # Distr max. expression
    ggplot(data.frame(max_R_expr_forEach=max_R_expr_forEach), aes(x=max_R_expr_forEach))+
        geom_freqpoly(bins=100)+theme_bw()
    
    # Distribution of cell percentages
    all_f = as.vector(data_LR_expression$matrix_expr_f)
    ggplot(data.frame(all_f=all_f), aes(x=all_f))+
        geom_freqpoly(bins=100)+theme_bw()
    
    # Determine max. cell percentages 
    max_R_fraction_forEach = apply(data_LR_expression$matrix_expr_f, 2, max)
    
    # That distribution
    ggplot(data.frame(max_R_fraction_forEach=max_R_fraction_forEach), aes(x=max_R_fraction_forEach))+
        geom_freqpoly(bins=100)+theme_bw()
    
    # Outlier in terms of Z-score
    ggplot(data.frame(Z=as.vector(data_LR_expression$matrix_expr_Z_g2collapse)), aes(x=Z))+
        geom_freqpoly(bins=100)+theme_bw()
    
    max_R_Z_forEach = apply(data_LR_expression$matrix_expr_Z_g2collapse, 2, max)
    ggplot(data.frame(max_R_Z_forEach=max_R_Z_forEach), aes(x=max_R_Z_forEach))+
        geom_freqpoly(bins=100)+theme_bw()
    
    
    ################################################################################
    
    # 2nd custom heatmap, showing most relevant data (hopefully)
    
    colsNoNa = apply(data_LR_expression$matrix_expr_x, 2, function(x) {!any(is.na(x))})
    
    annot_col2 = data_LR_expression$annotation_cols[max_R_expr_forEach>2&colsNoNa,,drop=F]
    annot_col2$receptor = sapply(str_split( rownames(annot_col2), '_'), function(x){x[[1]]})
    
    pheatmap(data_LR_expression$matrix_expr_x_g2collapse[,max_R_expr_forEach>2&colsNoNa], annotation_col = annot_col2)
    pheatmap(data_LR_expression$matrix_expr_x[,max_R_expr_forEach>2&colsNoNa], annotation_row = data_LR_expression$annotation_rows)
    
    colsNoNa3 = apply(data_LR_expression$matrix_expr_Z_g2collapse, 2, function(x) {!any(is.na(x))})
    annot_col3 = data_LR_expression$annotation_cols
    annot_col3$receptor = sapply(str_split( rownames(annot_col3), '_'), function(x){x[[1]]})
    pheatmap(data_LR_expression$matrix_expr_Z_g2collapse[,max_R_Z_forEach>2&colsNoNa3], annotation_col = annot_col3[max_R_Z_forEach>2&colsNoNa3,])
    
    # using 2 params for cutoff
    colsNoNa3 = apply(data_LR_expression$matrix_expr_Z_g2collapse, 2, function(x) {!any(is.na(x))})
    annot_col3 = data_LR_expression$annotation_cols
    annot_col3$receptor = sapply(str_split( rownames(annot_col3), '_'), function(x){x[[1]]})
    p=pheatmap(data_LR_expression$matrix_expr_Z_g2collapse[,max_R_Z_forEach>1.5&max_R_fraction_forEach>.4&colsNoNa3], 
                               annotation_col = annot_col3[max_R_Z_forEach>1.5&max_R_fraction_forEach>.4&colsNoNa3,], 
                               fontsize = 7)
    ggsave(plot = p, filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_LR_Intrx_heatmap1.png'), height=150, width=300, units='mm', dpi = 600)
    # ggsave(plot = p, filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_LIGAND_heatmap1.png'), height=150, width=300, units='mm', dpi = 600)
    
    ################################################################################
    # Now also create some UMAPs
    
    # Create UMAPs of ligands
    current_ligand_list = ligands_of_interest$HCM_up
    dummy=lapply(1:length(current_ligand_list), function(idx) {
        g = current_ligand_list[[idx]]
        g_full = shorthand_seurat_fullgenename_faster(current_analysis[[DATASET_NAME]], g, return_NA = T)
        if (!is.na(g_full)) {
            
            p=FeaturePlot(current_analysis[[DATASET_NAME]], feature=g_full, cols = col_vector_60, pt.size = .2)+
                theme_void()+ggtitle(paste0('Ligand ',g))+theme(legend.position = 'none')
            # test
            ggsave(plot = p, filename = paste0(base_dir,'LR_analysis/UMAPs/',DATASET_NAME,'_UMAP_LIGAND_',g,'.png'), height=172/3*2-4, width=172/3*2-4, units='mm')
        }
    })
    
    # Create UMAPs of receptors
    current_receptor_list = unique(unlist(receptors_of_interest$HCM_up))
    dummy=lapply(1:length(current_receptor_list), function(idx) {
        g = current_receptor_list[[idx]]
        g_full = shorthand_seurat_fullgenename_faster(current_analysis[[DATASET_NAME]], g, return_NA = T)
        if (!is.na(g_full)) {
            LL = toString(  quicky_give_L_for_R_OI(g)  )
            
            p=FeaturePlot(current_analysis[[DATASET_NAME]], feature=g_full, cols = col_vector_60, pt.size = .2)+
                theme_void()+ggtitle(paste0('Expr of R: ',g,'\n(linked to L: ',LL,')'))+theme(legend.position = 'none')
            # test
            ggsave(plot = p, filename = paste0(base_dir,'LR_analysis/UMAPs/',DATASET_NAME,'_UMAP_RECEPTOR_',g,'.png'), height=172/3*2-4, width=172/3*2-4, units='mm')
        }
    })
    
    # Use:
    # feature_list
    # feature_list_annotation
    dummy=lapply(1:length(feature_list), function(idx) {
        g = feature_list[[idx]]
        g_full = shorthand_seurat_fullgenename_faster(current_analysis[[DATASET_NAME]], g, return_NA = T)
        if (!is.na(g_full)) {
            an = feature_list_annotation[[idx]]
            
            p=FeaturePlot(current_analysis[[DATASET_NAME]], feature=g_full, cols = col_vector_60, pt.size = .2)+
                theme_void()+ggtitle(paste0('Expr of R: ',g,'\n(linked to L: ',an,')'))+theme(legend.position = 'none')
            # test
            ggsave(plot = p, filename = paste0(base_dir,'LR_analysis/UMAPs/',DATASET_NAME,'_UMAP_withR-',an,'_ExprPlotFor-',g,'.png'), height=172/3*2-4, width=172/3*2-4, units='mm')
        }
    })
    
    # NOTE:
    # I should also determine percentage in which these receptors are expressed, to determine validity of the z-scores..
    
    
    
}










