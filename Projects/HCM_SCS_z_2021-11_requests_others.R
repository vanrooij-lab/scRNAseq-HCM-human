

# Requests by Anne
gene_list_custom = list()
gene_list_custom$request_anne = c('LYPLAL1', 'DES', 'XIRP1', 'MYL12A', 'PSMA5', 'ANKRD2', 'NDUFS7', 'GLRX2', 'DPH5', 'LTBP1', 'PPP1R12B', 'FHL2', 'TECRL', 'IDH2', 'HSPB1', 'LDB3', 'CRYAB', 'ACTC1', 'MYL3')


################################################################################
################################################################################

# TO BE EXECUTED AT HPC

################################################################################
# Loading data

# Load the main script
script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

# --- load / set up data
    
# CURRENT_RUNNAME='all_RID2l_VAR'
CURRENT_RUNNAME='ALL.SP_RID2l_clExtended'

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
# Create custom order for annotation
current_analysis[[CURRENT_RUNNAME]]$annotation_paper_fct = 
    factor(current_analysis[[CURRENT_RUNNAME]]$annotation_paper_str, levels=c("vRooij", "Hu", "Teichmann"))

# --- end of load / set up data

################################################################################

# UMAPs per gene
for (GENE in gene_list_custom$request_anne) {

    # GENE='NPPA'
    p=shorthand_seurat_custom_expr(seuratObject = current_analysis[[CURRENT_RUNNAME]], 
                                         gene_of_interest = GENE,
                                         textsize = 8, pointsize = .1, custom_title = GENE, mymargin = .1) 
    # p
    ggsave(filename = paste0(base_dir, 'Rplots/anne/', CURRENT_RUNNAME, '_9_customUMAPs_',GENE,'_small.pdf'), plot = p,
           height=(172/3-4), width=(172/3-4), units='mm', device = cairo_pdf)

}

# combined
for (listname in names(gene_list_custom)) {
        
    p=shorthand_seurat_custom_expr(seuratObject = current_analysis[[CURRENT_RUNNAME]], 
                                         gene_of_interest = gene_list_custom[[listname]],
                                         textsize = 6, pointsize = .5, custom_title = listname, mymargin = .1, zscore = T) 
    ggsave(filename = paste0(base_dir, 'Rplots/anne/', CURRENT_RUNNAME, '_9_customUMAPs_COMPOSITE-',listname,'_small.pdf'), plot = p,
               height=(172/3-4), width=(172/3-4), units='mm', device = cairo_pdf)
}

# Summary plots
# Note: output here also goes to standard "Rplots" dir

# Box plots
shorthand_custom_boxplot(seuratObject_list=current_analysis, 
                         gene_lists=gene_list_custom, 
                         seuratObjectNameToTake=CURRENT_RUNNAME, 
                         group.by='annotation_paper_beatified', 
                         topX=Inf, mylimits=.01) 

# Summary plots
p_lists=shorthand_custom_compositeplot(seuratObject_list=current_analysis, 
                         gene_lists=gene_list_custom, 
                         seuratObjectNameToTake=CURRENT_RUNNAME, 
                         group.by='annotation_paper_beatified', 
                         group.by2='annotation_patient_fct',
                         zscore=T) 
p=wrap_plots(p_lists$p_violin_list, nrow=1)
ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customCOMPOSITE_ANNE.pdf'), plot = p,
       height=(PANEL_WIDTH-4), width=(PANEL_WIDTH-4), units='mm', device=cairo_pdf)

p=wrap_plots(p_lists$p_bar_list_g2, nrow=1)
ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customCOMPOSITE_ANNE_g2.pdf'), plot = p,
       height=(PANEL_WIDTH-4), width=(PANEL_WIDTH-4),  units='mm', device=cairo_pdf)

################################################################################



################################################################################
################################################################################

# TO BE EXECUTED LOCALLY
# (On Rooij maps)

################################################################################
# Load data locally

LOCAL=1; script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Projects/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

################################################################################

DATASET_NAME='ROOIJonly_RID2l_clExtended'
if (!exists('current_analysis')) {current_analysis = list()}
current_analysis[[DATASET_NAME]] =
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))

# UMAPs per gene
for (GENE in gene_list_custom$request_anne) {

    # GENE='NPPA'
    p=shorthand_seurat_custom_expr(seuratObject = current_analysis[[DATASET_NAME]], 
                                         gene_of_interest = GENE,
                                         textsize = 8, pointsize = 1, custom_title = GENE, mymargin = .1) 
    # p
    ggsave(filename = paste0(base_dir, 'Rplots/anne/', DATASET_NAME, '_9_customUMAPs_',GENE,'_small.pdf'), plot = p,
           height=(172/3-4), width=(172/3-4), units='mm', device = cairo_pdf)

}












