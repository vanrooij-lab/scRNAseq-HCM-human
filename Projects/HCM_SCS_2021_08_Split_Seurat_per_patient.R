

################################################################################

# Seurat analyses at patient-level

# To be executed at HPC

script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

################################################################################
# Generate per patient data sets

for (DATASET_NAME in c('ROOIJonly.sp.bt_RID2l', 'HUonly.sp.bt_RID2l', 'TEICHMANNonly.sp.bt_RID2l')) {

    print(paste0('Loading data ', DATASET_NAME))
    
    # Load the data
    if (!exists('current_analysis')) {current_analysis = list()}
    current_analysis[[DATASET_NAME]] =
            LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))
    
    # Split the data
    all_patients = unique(current_analysis[[DATASET_NAME]]$annotation_patient_str)
    temp_seurat_subsets = lapply(all_patients, function(current_patient) {
                                print(paste0('Generating sep data for patient ', current_patient))
                                return(
                                    subset(current_analysis[[DATASET_NAME]], annotation_patient_str==current_patient))
                                })
    names(temp_seurat_subsets) = all_patients
    # Copy 'm into usual parameter
    for (current_seurat_subset_name in names(temp_seurat_subsets) ) {
     current_analysis[[current_seurat_subset_name]] = temp_seurat_subsets[[current_seurat_subset_name]]
    }
    # Remove temporary subset param
    rm('temp_seurat_subsets')
        
    ###
    # Perform the standard analysis
    
    # Run the analysis
    for (current_patient in all_patients) {
        
        print(paste0('Running analysis for patient ',current_patient))
        current_analysis[[paste0(current_patient,'RID2l')]] =
            mySeuratAnalysis_verybasic_part2only(mySeuratObject = current_analysis[[current_patient]], do.scale = F, do.center = F)
        
        # Save it
        SaveH5Seurat(object = current_analysis[[paste0(current_patient,'RID2l')]], overwrite = T,
            filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',paste0(current_patient,'RID2l'),'.h5seurat'))
    }
    
    # Create some plots
    print('Creating plots..')
    for (current_set in paste0(all_patients,'RID2l')) {
        mySeuratCommonPlots(mySeuratObject = current_analysis[[current_set]], run_name = current_set)
    }
    
    print('all done')

}
    
# # Save the analyses
# for (current_set in paste0(all_patients,'RID2l')) {
#     SaveH5Seurat(object = current_analysis[[current_set]], overwrite = T,
#         filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',current_set,'.h5seurat'))
# }

################################################################################
################################################################################

# Load the analyses (in case of re-running this script)
if (F) {
    DATASET_NAME='ROOIJonly.sp.bt_RID2l'
    # Load original analysis (not strictly necessary), just to get "patient annotation"
    if (!exists('current_analysis')) {current_analysis = list()}
    if (!(DATASET_NAME %in% names(current_analysis))) { 
        current_analysis[[DATASET_NAME]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat')) 
    }
    # Patient names
    all_patients = unique(current_analysis[[DATASET_NAME]]$annotation_patient_str)
    
    for (current_set in paste0(all_patients,'RID2l')) {
        current_analysis[[current_set]] =
            LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',current_set,'.h5seurat'))
    }
}

################################################################################
################################################################################

# Some additional custom plots

if (F) {
    for (GENE in c('TTN','XIRP2','CMYA5')) {
        
        #GENE='TTN'
        #GENE='XIRP2'
        #GENE='CMYA5'
    
        # Make plot
        plot_list = lapply(sort(all_patients), function(current_patient) {
            # current_patient = 'R.P1'
            shorthand_seurat_custom_expr(seuratObject = current_analysis[[paste0(current_patient,'RID2l')]], gene_of_interest = GENE, custom_title = paste0(GENE,' (',gsub('^R\\.','',current_patient),')'), mymargin = 1, add_box = T)
            })
         
        p=wrap_plots(plot_list, nrow=1)    
        p
        # Save
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/R.PALL.RooijPerPatient_Custom_',GENE,'.pdf'), 
               width = (PANEL_WIDTH*3-4), height= (PANEL_WIDTH*3-4)/5, units='mm', device = cairo_pdf)
    
    }

}
















