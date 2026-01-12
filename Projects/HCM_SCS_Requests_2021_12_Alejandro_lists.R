
# Creating some UMAPs for Alejandro

# Script to convert mouse genes to human
source('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Functions/Conversion_mouse_human_genes.R')



big_list_Alejandro = list(
    Ligands_MP = c('C1qb', 'C1qa', 'Apoe', 'Igfbp4', 'Ly86', 'F13a1', 'Cp', 'Gas6', 'Timp2'), 
    Receptors = c('Cd93', 'Cr1l', 'Cspg4', 'Itga4', 'Itga9', 'Itgb1', 'Cd180', 'Fzd8', 'Lrp6', 'Chrna4', 'Ldlr', 'Lrp1', 'Lrp2', 'Lrp5', 'Lrp8', 'Scarb1', 'Sorl1', 'Vldlr'), 
    Exercise = c('Sox17', 'Foxp1'),
    TAB = c('Filip1', 'Chd2', and 'Sorbs2', 'Nrap', 'Xirp1', 'Xirp2', 'Ankrd1', 'Lmod2', 'Myh6', 'Myh7', 'Nppa', 'Nppb', 'Des', 'Hspb7', 'Dusp27', 'Csrp3', 'Thbs4', 'Synpo2l', 'Gapdh', 'Hprt'))
big_list_Alejandro_sorted = lapply(big_list_Alejandro, sort)

# Humanize the lists
#big_list_Alejandro_humanized = 
#    lapply(big_list_Alejandro, convertMouseGeneList)
#big_list_Alejandro_humanized_sorted = lapply(big_list_Alejandro_humanized, sort)
# lapply(big_list_Alejandro_humanized_sorted, toString)

# convertMouseGeneList(c('Filip1', 'Chd2', 'Sorbs2'))

big_list_Alejandro_humanized_sorted = list(
    Ligands_MP = c('APOE', 'C1QA', 'C1QB', 'CP', 'F13A1', 'GAS6', 'IGFBP4', 'LY86', 'TIMP2'),
    Receptors = c('CD180', 'CD93', 'CHRNA4', 'CR1', 'CSPG4', 'FZD8', 'ITGA4', 'ITGA9', 'ITGB1', 'LDLR', 'LRP1', 'LRP2', 'LRP5', 'LRP6', 'LRP8', 'SCARB1', 'SORL1', 'VLDLR'),
    Exercise = c('FOXP1', 'SOX17'),
    TAB = c("FILIP1", "CHD2",   "SORBS2", 'ANKRD1', 'CSRP3', 'DES', 'GAPDH', 'HPRT1', 'HSPB7', 'LMOD2', 'MYH6', 'MYH7', 'NPPA', 'NPPB', 'NRAP', 'STYXL2', 'SYNPO2L', 'THBS4', 'XIRP1', 'XIRP2'))



##########

# Make the plots @ HPC

script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

DATASET_NAME='TeichAll_Rooij_dflt' 
current_analysis=list()
current_analysis[[DATASET_NAME]] = 
    LoadH5Seurat(file = paste0(base_dir, 'LR_analysis/', DATASET_NAME, '.h5seurat'))
    
# Now run the plots
for (current_ligand_list in big_list_Alejandro_humanized_sorted) {
    
    dummy=lapply(1:length(current_ligand_list), function(idx) {
        g = current_ligand_list[[idx]]
        g_full = shorthand_seurat_fullgenename_faster(current_analysis[[DATASET_NAME]], g, return_NA = T)
        
        if (!is.na(g_full)) {
            
            p=FeaturePlot(current_analysis[[DATASET_NAME]], feature=g_full, cols = col_vector_60, pt.size = .2)+
                theme_void()+ggtitle(paste0('Gene: ',g))+theme(legend.position = 'none')
            # test
            ggsave(plot = p, filename = paste0(base_dir,'LR_analysis/UMAPs/Alejandro/',DATASET_NAME,'_UMAP_gene_',g,'.png'), height=172/3*2-4, width=172/3*2-4, units='mm')
        }
        
    })
    
}

# Now also create violins
for (current_ligand_list in big_list_Alejandro_humanized_sorted) {
    
    dummy=lapply(1:length(current_ligand_list), function(idx) {
        g = current_ligand_list[[idx]]
        g_full = shorthand_seurat_fullgenename_faster(current_analysis[[DATASET_NAME]], g, return_NA = T)
        
        if (!is.na(g_full)) {
            
            # Violin plot
            p=VlnPlot(object = current_analysis[[DATASET_NAME]],features = g_full, cols = col_vector_60, 
                      group.by = 'cell_type_inclR', pt.size = 0)+
                ggtitle(paste0('Gene: ',g))+theme(legend.position = 'none')
                        
            # save
            ggsave(plot = p, filename = paste0(base_dir,'LR_analysis/UMAPs/Alejandro/',DATASET_NAME,'_VLN_gene_',g,'.pdf'), 
                   height=172/3*2-4, width=172/3*3-4, units='mm', device = cairo_pdf)
            
        }
        
    })
    
}






