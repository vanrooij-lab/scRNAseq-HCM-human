

################################################################################

# Seurat analyses at patient-level

################################################################################

# Load the data
DATASET_NAME='ROOIJonly_RID2l'
# DATASET_NAME='TEICHMANNonly_RID2l'
# DATASET_NAME='HUonly_RID2l'


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

################################################################################

# Run the analysis
for (current_patient in all_patients) {
    print(paste0('Running analysis for patient ',current_patient))
    current_analysis[[paste0(current_patient,'RID2l')]] =
        mySeuratAnalysis_verybasic_part2only(mySeuratObject = current_analysis[[current_patient]], do.scale = F, do.center = F)
}
print('all done')

# Save the analyses
for (current_set in paste0(all_patients,'RID2l')) {
    SaveH5Seurat(object = current_analysis[[current_set]], overwrite = T,
        filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',current_set,'.h5seurat'))
}

# Create some plots
for (current_set in paste0(all_patients,'RID2l')) {
    mySeuratCommonPlots(mySeuratObject = current_analysis[[current_set]], run_name = current_set)
}


################################################################################
################################################################################

# Some additional custom plots


for (GENE in c('TTN','XIRP2','CMYA5')) {
    
    #GENE='TTN'
    #GENE='XIRP2'
    #GENE='CMYA5'

    # Make plot
    plot_list = lapply(sort(all_patients), function(current_patient) {
        shorthand_seurat_custom_expr(seuratObject = current_analysis[[paste0(current_patient,'RID2l')]], gene_of_interest = GENE, custom_title = paste0(GENE,' (',current_patient,')'))})
    p=wrap_plots(plot_list, nrow=1)    
    p
    # Save
    ggsave(plot = p,filename = paste0(base_dir, 'Rplots/R.PALL.RooijPerPatient_Custom_',GENE,'.pdf'), 
           width = min(184.6, 40*5), height= 40, units='mm')

}


















