
########################################################################

# For convenience, on HPC, this script can be sourced by:
# script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')
#
# Local:
# LOCAL=1; script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Projects/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

########################################################################

# Code copied from main script
if (F) {
    
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
    
}

########################################################################

if (F) {
    
    ROOIJ_DATASET = 'ROOIJonly.sp.bt_RID2l'
    
    dir.create(paste0(base_dir, 'Rplots/2022_10-extra/'))
    
    # Get TF names
    load(file=paste0(base_dir,'Rdata/',ROOIJ_DATASET,'__SCENIC_reg_top_genes_sorted_full.Rdata')) # SCENIC_reg_top_genes_sorted_full
    SCENIC_TF_list = names(SCENIC_reg_top_genes_sorted_full)
    print(paste0('TF names from file: ',toString(SCENIC_TF_list))) # ZEB1, FOXN3, NFE2L1, MAFK, CLOCK, REST, CEBPZ, MEF2D, MITF, MEF2A, ETV1, CREB1, MXI1, JUND, YY1, SREBF2, CEBPB, SRF, NR3C1, ESRRA, USF2, ZNF91
    
    # Expression of the SCENIC TFs
    LIST_NAME='SCENIC_TFs'; CAT = 'annotation_paper_beatified'; CATNAME = 'Dataset'; 
    genes_of_interest = SCENIC_TF_list
    plot_list_violins = lapply(genes_of_interest, 
        function(x) {shorthand_plotViolinBox_custom(current_analysis,analysis_name=CURRENT_RUNNAME,cat_by = CAT, 
                                                    gene_of_interest=x,base_dir=base_dir, cat_name = CATNAME,
                                                    type = 'violin')})
    plot_list_boxplots = lapply(genes_of_interest, 
        function(x) {shorthand_plotViolinBox_custom(current_analysis,analysis_name=CURRENT_RUNNAME,cat_by = CAT, 
                                                    gene_of_interest=x,base_dir=base_dir, cat_name = CATNAME,
                                                    type = 'boxplot')})
    
    plot_list_prchigh = lapply(genes_of_interest, 
        function(x) {shorthand_plotPercentageCellsHighForGene(myseuratobjectlist=current_analysis, analysis_name=CURRENT_RUNNAME,cat_by = CAT, 
                                                              gene_of_interest=x, cat_name = CATNAME)})
    plot_list_prchigh90 = lapply(genes_of_interest, 
        function(x) {shorthand_plotPercentageCellsHighForGene(myseuratobjectlist=current_analysis, analysis_name=CURRENT_RUNNAME,cat_by = CAT, 
                                                              gene_of_interest=x, cat_name = CATNAME, PERCENTILE = .9)})
    
    plot_list_freq = lapply(genes_of_interest, 
        function(x) {shorthand_plotFreqPoly_custom(current_analysis,analysis_name=CURRENT_RUNNAME,cat_by = CAT, 
                                                    gene_of_interest=x,base_dir=base_dir, cat_name = CATNAME)})
    
    for (p_idx in 1:length(genes_of_interest)) {
        
        
        currentgene=genes_of_interest[p_idx]
        
        # Save the violin
        p=plot_list_violins[[p_idx]]
        PLOTTYPE='violin'
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'_',currentgene,'_separated.pdf'), width = 172/3-4, height= 172/3-4, units='mm', device=cairo_pdf)    
        p=p+give_better_textsize_plot(7)
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'__small_',currentgene,'_separated.pdf'), width = 172/4-4, height= 172/4-4, units='mm', device=cairo_pdf)
        
        # Save the boxplot
        p=plot_list_boxplots[[p_idx]]
        PLOTTYPE='boxplot'
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'_',currentgene,'_separated.pdf'), width = 172/3-4, height= 172/3-4, units='mm', device=cairo_pdf)    
        p=p+give_better_textsize_plot(7)
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'__small_',currentgene,'_separated.pdf'), width = 172/4-4, height= 172/4-4, units='mm', device=cairo_pdf)
        
        # Save the percentage high plots
        p=plot_list_prchigh[[p_idx]]
        PLOTTYPE='prchigh'
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'_',currentgene,'_separated.pdf'), width = 172/3-4, height= 172/3-4, units='mm', device=cairo_pdf)    
        p=p+give_better_textsize_plot(7)
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'__small_',currentgene,'_separated.pdf'), width = 172/4-4, height= 172/4-4, units='mm', device=cairo_pdf)
        # Save the percentage high plots for 90%
        p=plot_list_prchigh90[[p_idx]]
        PLOTTYPE='prchigh90'
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'_',currentgene,'_separated.pdf'), width = 172/3-4, height= 172/3-4, units='mm', device=cairo_pdf)    
        p=p+give_better_textsize_plot(7)
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'__small_',currentgene,'_separated.pdf'), width = 172/4-4, height= 172/4-4, units='mm', device=cairo_pdf)
        
        # Save the freq plots
        p=plot_list_freq[[p_idx]]
        PLOTTYPE='histogram'
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'_',currentgene,'_separated.pdf'), width = 172/3*2-4, height= 172/3-4, units='mm', device=cairo_pdf)    
        p=p+give_better_textsize_plot(7)
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'__small_',currentgene,'_separated.pdf'), width = 172/4*2-4, height= 172/4-4, units='mm', device=cairo_pdf)
        
        
        print(paste0(currentgene, ' done ..'))
    }
    
    # From other code:
    # 
    # # Small 2x2 version
    # p=wrap_plots(plot_list, nrow = 2)& theme(plot.margin = margin(t = .5, r = .5, b = .5, l = 2, unit = "mm"))
    # ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',CURRENT_RUNNAME,'_6_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'.pdf'), width = 172/3-4, height= 172/3-4, units='mm', device=cairo_pdf)
    # #p
    # p=wrap_plots(plot_list, nrow = 1)& theme(plot.margin = margin(t = .5, r = .5, b = .5, l = 2, unit = "mm"))
    # # ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',CURRENT_RUNNAME,'_6_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'-style2.pdf'), width = 172-4, height= (172-4)/4, units='mm', device=cairo_pdf)
    # #ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',CURRENT_RUNNAME,'_6_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'-style2.pdf'), width = (172/3*2)-4, height= ((172/3*2)-4)/4*1.5, units='mm', device=cairo_pdf)
    # ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',CURRENT_RUNNAME,'_6_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'-style2.pdf'), width = (172/3*2)-4, height= PANEL_HEIGHT, units='mm', device=cairo_pdf)
    # 
    
}

########################################################################
# More expression comparison between the sets for hand-picked genes
    
if (F) {

    genes_of_interest_batch2_ =    
        # ±Correlated with size
        c('CALM3', 'HSP90AA1', 'PKM',
        # Additions from correlated with NPPA
        'RTN4', 'PSAP', 'MFGE8', 'PROS1', 'YY1', 'PIN1', 'TSC22D1', 'TGM2',
        # additions from DE in mice n=1
        "CCN2", "CALR", 
        "NFE2L1", "DBP", "ANKRD1",   "SFPQ", 
        # Additions from HCM vs. Ctrl1(±2)
        "GNAS", "LPL", "NFIC",
        # pos con
        "SRF")
            # conversion_mouse_to_human_symbols[c('Ccn2', 'Nppb', 'Rtn4', 'Calr')] # "CCN2" "NPPB" "RTN4" "CALR"
            #conversion_mouse_to_human_symbols[c('Nfe2l1', 'Dbp', 'Ankrd1', 'Sfpq')] # c("NFE2L1",    "DBP", "ANKRD1",   "SFPQ" )
        # genes_of_interest_batch2_[genes_of_interest_batch2_ %in% shorthand_cutname( rownames(current_analysis[[CURRENT_RUNNAME]]) )]
        # CCN2 is not detected in human data ..
        # conversion_ens_sym_107manual_df[conversion_ens_sym_107manual_df$Gene.name=='CCN2',]
        'ENSG00000118523' %in% shorthand_cutname( rownames(current_analysis[[CURRENT_RUNNAME]]), PART1OR2 = 2)
    print(paste0('Couldnt find: ', 
                 toString(   genes_of_interest_batch2_[!(genes_of_interest_batch2_ %in% shorthand_cutname( rownames(current_analysis[[CURRENT_RUNNAME]]))) ]  ) 
                          ))
    genes_of_interest_batch2 = genes_of_interest_batch2_[genes_of_interest_batch2_ %in% shorthand_cutname( rownames(current_analysis[[CURRENT_RUNNAME]]) )]
    
    CAT = 'annotation_paper_beatified'; CATNAME = 'Dataset'; 
    # genes_of_interest_batch2
    genes_of_interest=genes_of_interest_batch2
    plot_list_freq = lapply(genes_of_interest, 
            function(x) {shorthand_plotFreqPoly_custom(current_analysis,analysis_name=CURRENT_RUNNAME,cat_by = CAT, 
                                                        gene_of_interest=x,base_dir=base_dir, cat_name = CATNAME)})
        
    # save the plots
    LIST_NAME='GOI_all_r2'
    for (p_idx in 1:length(genes_of_interest)) {
        
        currentgene=genes_of_interest[p_idx]
        
        # Save the freq plots
        p=plot_list_freq[[p_idx]]
        PLOTTYPE='histogram'
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'_',currentgene,'_separated.pdf'), width = 172/3*2-4, height= 172/3-4, units='mm', device=cairo_pdf)    
        p=p+give_better_textsize_plot(7)
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'__small_',currentgene,'_separated.pdf'), width = 172/4*2-4, height= 172/4-4, units='mm', device=cairo_pdf)
        
        
        print(paste0(currentgene, ' done ..'))
    }

}

################################################################################
# Genes of interest batch 3

  
if (F) {

    dir.create(paste0(base_dir, 'Rplots/2022_10-extra/r3/'))
    
    # These genes are either in top 10 module 4 or top 10 cluster
    genes_of_interest_batch3_ =    
        c('CMYA5', 'XIRP2', 'ZNF106', 'MAP4', 'NRAP', 'ANKRD1', 'TTN', 'HIPK2', 'DES', 'TECRL', 'TCAP', 'MYOM1',
          # c(
          # module 4 TFs homer
           'SCL', 'TEAD4', 'TCF7L2', 'GATA1', 'MEF2C',
          # cl 2 TFs 
          # c(
          'MEF2C','TEAD4','MEF2B','TCF3','ARID3B','GATA4','TBX5','CEBPB','AR','DSP','MYO18B','NRAP'
          )
            

    # check names        
    print(paste0('Couldnt find: ', 
                 toString(   genes_of_interest_batch3_[!(genes_of_interest_batch3_ %in% shorthand_cutname( rownames(current_analysis[[CURRENT_RUNNAME]]))) ]  ) 
                          ))
    genes_of_interest_batch3 = genes_of_interest_batch3_[genes_of_interest_batch3_ %in% shorthand_cutname( rownames(current_analysis[[CURRENT_RUNNAME]]) )]
    
    
    
    CAT = 'annotation_paper_beatified'; CATNAME = 'Dataset'; 
    # genes_of_interest_batch2
    genes_of_interest=genes_of_interest_batch3
    plot_list_freq = lapply(genes_of_interest, 
            function(x) {shorthand_plotFreqPoly_custom(current_analysis,analysis_name=CURRENT_RUNNAME,cat_by = CAT, 
                                                        gene_of_interest=x,base_dir=base_dir, cat_name = CATNAME)})
        
    # save the plots
    LIST_NAME='GOI_all_r3'
    for (p_idx in 1:length(genes_of_interest)) {
        
        currentgene=genes_of_interest[p_idx]
        
        # Save the freq plots
        p=plot_list_freq[[p_idx]]
        PLOTTYPE='histogram'
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/r3/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'_',currentgene,'_separated.pdf'), width = 172/3*2-4, height= 172/3-4, units='mm', device=cairo_pdf)    
        p=p+give_better_textsize_plot(7)
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/r3/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'__small_',currentgene,'_separated.pdf'), width = 172/4*2-4, height= 172/4-4, units='mm', device=cairo_pdf)
        
        
        print(paste0(currentgene, ' done ..'))
    }

}

################################################################################

if (F) {

    
    dir.create(paste0(base_dir, 'Rplots/2022_10-extra/r4_Ale/'))
    
    # Alejandro's genes from proteomics on TAB mice
    genes_of_interest_batch4_ =    
        c('XIRP2', 'XIRP1', 'NRAP', 'NPPB', 'NPPA', 'MYH7', 'LMOD2', 'HSPB7', 'DUSP27', 'DES', 'CSRP3', 'ANKRD1', 'SYNPO2L', 'SORBS2', 'FILIP1', 'CDH2')
            

    # check names        
    print(paste0('Couldnt find: ', 
                 toString(   genes_of_interest_batch4_[!(genes_of_interest_batch4_ %in% shorthand_cutname( rownames(current_analysis[[CURRENT_RUNNAME]]))) ]  ) 
                          ))
    genes_of_interest_batch4 = genes_of_interest_batch4_[genes_of_interest_batch4_ %in% shorthand_cutname( rownames(current_analysis[[CURRENT_RUNNAME]]) )]
    
    
    
    CAT = 'annotation_paper_beatified'; CATNAME = 'Dataset'; 
    # genes_of_interest_batch2
    genes_of_interest=genes_of_interest_batch4
    plot_list_freq = lapply(genes_of_interest, 
            function(x) {shorthand_plotFreqPoly_custom(current_analysis,analysis_name=CURRENT_RUNNAME,cat_by = CAT, 
                                                        gene_of_interest=x,base_dir=base_dir, cat_name = CATNAME)})
        
    # save the plots
    LIST_NAME='GOI_all_r4'
    for (p_idx in 1:length(genes_of_interest)) {
        
        currentgene=genes_of_interest[p_idx]
        
        # Save the freq plots
        p=plot_list_freq[[p_idx]]
        PLOTTYPE='histogram'
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/r4_Ale/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'_',currentgene,'_separated.pdf'), width = 172/3*2-4, height= 172/3-4, units='mm', device=cairo_pdf)    
        p=p+give_better_textsize_plot(7)
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/2022_10-extra/r4_Ale/',CURRENT_RUNNAME,'_9_',PLOTTYPE,'_by-', CAT, '_',LIST_NAME,'__small_',currentgene,'_separated.pdf'), width = 172/4*2-4, height= 172/4-4, units='mm', device=cairo_pdf)
        
        
        print(paste0(currentgene, ' done ..'))
    }
    

}
    
################################################################################
# General enrichments of TFs

ANALYSIS_NAME_ALL = 'ALL.SP_btypSel_RID2l_clExtended'
load(paste0(base_dir,'Rdata/DE_cluster__',ANALYSIS_NAME_ALL,'.Rdata')) # DE_cluster

View(DE_cluster$ALL.SP_btypSel_RID2l_clExtended$`4`)

genes_DE_allData_cluster4 = DE_cluster$ALL.SP_btypSel_RID2l_clExtended$`4`

genes_DE_allData_cluster4$gene = rownames(genes_DE_allData_cluster4)
genes_DE_allData_cluster4$ens = shorthand_cutname( rownames(genes_DE_allData_cluster4), PART1OR2 = 1)
genes_DE_allData_cluster4$sym = shorthand_cutname( rownames(genes_DE_allData_cluster4), PART1OR2 = 2)

table_to_export_TF_allpooled = genes_DE_allData_cluster4[genes_DE_allData_cluster4$sym %in% human_TFs_combined_temp,]
table_to_export_TF_allpooled=table_to_export_TF_allpooled[order(table_to_export_TF_allpooled$avg_log2FC, decreasing = T),c('sym','avg_log2FC','p_val','p_val_adj')]
View(table_to_export_TF_allpooled[table_to_export_TF_allpooled$avg_log2FC>0,])

# For highlighting purposes:
SCENIC_TF_list[SCENIC_TF_list %in% table_to_export_TF_allpooled[table_to_export_TF_allpooled$avg_log2FC>0,]$sym]

# Now also show ligands (note that I did this earlier too)
table_to_export_lig_allpooled = genes_DE_allData_cluster4[genes_DE_allData_cluster4$ens %in% ligand_list_unique_ENS,]
table_to_export_lig_allpooled=table_to_export_lig_allpooled[order(table_to_export_lig_allpooled$avg_log2FC, decreasing = T),c('sym','avg_log2FC','p_val','p_val_adj')]
View(table_to_export_lig_allpooled[table_to_export_lig_allpooled$avg_log2FC>0,])


########################################################################
# Some customized stuff
# (Note: this was abandoned)

if (F) {
    
    # From main file
    TF_SCENIC_favorites_target_expression = R_values_ordered$Group.2[R_values_ordered$x_compmax>.2]
    
    # By eye based on above Violin plots    
    TF_SCENIC_manual_favorites_expression = 
        c('YY1', 'SREBF2', 'REST', 'NFE2L1','MEF2D','MAFK','ETV1','CREB1','CEBPZ','CEBPB')
        # typo sanity check
        # all(manual_favorite_TF_list %in% R_values_ordered$Group.2)
    toString(TF_SCENIC_manual_favorites_expression)
    
    toString(TF_SCENIC_favorites_target_expression[TF_SCENIC_favorites_target_expression %in% TF_SCENIC_manual_favorites_expression])
    
    plot_Venn_MW_2lists(TF_SCENIC_manual_favorites_expression,TF_SCENIC_favorites_target_expression, 
                        'TF_expression', 'target_expression')
    
}


########################################################################


load(file = paste0(base_dir,'Rdata/FSCA__correlations_FSCA_per_patient_combined.Rdata')) # correlations_FSCA_per_patient_combined


# correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$gene_symbol=='YY1',]

# Most of the SCENIC TFs are not in the FSCA data
TF_SCENIC_manual_favorites_expression
TF_SCENIC_manual_favorites_expression[TF_SCENIC_manual_favorites_expression %in% correlations_FSCA_per_patient_combined$gene_symbol]
toString(SCENIC_TF_list[SCENIC_TF_list %in% correlations_FSCA_per_patient_combined$gene_symbol])

# And those that are don't seem to be significant
View(correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$gene_symbol %in% TF_SCENIC_manual_favorites_expression,])
View(correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$gene_symbol %in% SCENIC_TF_list,])

max(correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$gene_symbol %in% SCENIC_TF_list,]$P4_cor,
    correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$gene_symbol %in% SCENIC_TF_list,]$P5_cor)
min(correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$gene_symbol %in% SCENIC_TF_list,]$P4_p.adjust,
    correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$gene_symbol %in% SCENIC_TF_list,]$P5_p.adjust)










