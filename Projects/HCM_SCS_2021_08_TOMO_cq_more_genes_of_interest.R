
script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'; desired_command='dummy'; source(paste0(script_dir, 'HCM-SCS_2021-06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

DATASET_NAME='ALL.SP_RID2l'

if (!exists('current_analysis')) {current_analysis = list()}
    current_analysis[[DATASET_NAME]] =
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))

# Super redundant with previous, but here again in case you're doing manual run for regulon plots only
if (!('annotation_paper_oneletter_fct' %in% names(current_analysis[[DATASET_NAME]]))) {
    name_change = c("vRooij"="R", "Hu"="H", "Teichmann"="T")
    current_analysis[[DATASET_NAME]]$annotation_paper_oneletter_fct = 
        factor(name_change[current_analysis[[DATASET_NAME]]$annotation_paper_str], levels=c('R','H','T'))
}    
    
# Additional plots for Tomo
# Box plots
tomo_gene_list = list(TOMO=c('ACTC1', 'MYL12A', 'LTBP1', 'PSMA5', 'LYPLAL1', 'GLRX2', 'IDH2', 'FHL2', 'DPH5', 'PPP1R12B',  'ANKRD2', 'TECRL', 'NDUFS7', 'HSPB1', 'LDB3', 'MYL3', 'XIRP1', 'CRYAB', 'DES'))
shorthand_custom_boxplot(seuratObject_list=current_analysis, 
                                 gene_lists=tomo_gene_list, 
                                 seuratObjectNameToTake=DATASET_NAME, 
                                 group.by='annotation_paper_oneletter_fct', 
                                 topX=Inf, mylimits=.01) 

# UMAPs
p_tomo_umap = lapply(tomo_gene_list$TOMO, function(gene_name) {
                    shorthand_seurat_custom_expr(seuratObject = current_analysis[[DATASET_NAME]], 
                                     gene_of_interest = gene_name,
                                     textsize = 6, pointsize = .5, custom_title = gene_name, mymargin = .1, zscore = T) 
            })
p=wrap_plots(p_tomo_umap, ncol=4)
ggsave(filename = paste0(base_dir, 'Rplots/', DATASET_NAME, '_9_customUMAPs_TOMO.pdf'), plot = p,
               height=ceiling(length(tomo_gene_list$TOMO)/4)*30, width=min(184.6-4, 30*4), units='mm')     
ggsave(filename = paste0(base_dir, 'Rplots/', DATASET_NAME, '_9_customUMAPs_TOMO.png'), plot = p,
               height=ceiling(length(tomo_gene_list$TOMO)/4)*30, width=min(184.6-4, 30*4), units='mm')     

 
