
################################################################################

# Ligand-Receptor analysis of the HCM SCS data
# =====
# Approach:
# 1. Identify ligand genes that are enriched in 
#       - HCM data vs others
#       - Our HCM clusters
#       --> both can be done using the existing FC tables, 
#           by filtering for ligands (XXX et al.).

################################################################################
# Libs etc

# Load libraries etc using the main HCM SCS script
LOCAL=1; script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Projects/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

# First load a list of receptor-ligand interactions
source('/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Resources/human_LR_list.R')
    # file.edit('/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Resources/human_LR_list.R')
    # ligand_list
    # ligand_receptor_pairs_df

# Convenient later
source('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Functions/Conversion_mouse_human_genes.R')
    # file.edit('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Functions/Conversion_mouse_human_genes.R')

################################################################################

# Load the HCM-specific dataset
DATASET_NAME='ROOIJonly_RID2l_clExtended'
if (!exists('current_analysis')) {current_analysis = list()}
current_analysis[[DATASET_NAME]] =
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))

# Load the FC-tables for pooled data
load(paste0(base_dir, 'Rdata/DE_cluster__ALL.SP_btypSel_RID2l_clExtended.Rdata')) # DE_cluster
# load(paste0(base_dir, 'Rdata/DE_cluster__ALL.SP_RID2l_clExtended.Rdata')) # DE_cluster
DE_cluster_ALL.SP_btypSel_RID2l_clExtended = DE_cluster$ALL.SP_btypSel_RID2l_clExtended
    # Note that cluster 4 is the one which contains all our cells
# load(paste0(base_dir, 'Rdata/enriched_genes_lists_clusters_ALL.SP__ALL.SP_RID2l_clExtended.Rdata'))

# And now also for the HCM-specific clusters
load(paste0(base_dir, 'Rdata/DE_cluster__ROOIJonly.sp.bt_RID2l_clExtended.Rdata')) # DE_cluster
DE_cluster_ROOIJonly_RID2l_clExtended = DE_cluster$ROOIJonly_RID2l_clExtended

# First gather the appropriate list of FC genes
# View(DE_cluster$ALL.SP_RID2l_clExtended$`4`)
FC_pooled_cluster_4 = DE_cluster_ALL.SP_btypSel_RID2l_clExtended$`4` 
FC_pooled_cluster_4$genename_short = shorthand_cutname( rownames(DE_cluster_ALL.SP_btypSel_RID2l_clExtended$`4`) )
FC_pooled_cluster_4_sel = FC_pooled_cluster_4[FC_pooled_cluster_4$p_val_adj<.05&FC_pooled_cluster_4$avg_log2FC>0,]

# Now make an overview of ligands only
FC_pooled_cluster_4_sel_ligands = FC_pooled_cluster_4_sel[FC_pooled_cluster_4_sel$genename_short %in% ligand_list,]
FC_pooled_cluster_4_sel_ligands_ordered = FC_pooled_cluster_4_sel_ligands[order(FC_pooled_cluster_4_sel_ligands$avg_log2FC,decreasing = T),]
FC_pooled_cluster_4_sel_ligands_ordered$avg_FC = 2^FC_pooled_cluster_4_sel_ligands_ordered$avg_log2FC

# Save the ligand list
if (!dir.exists( paste0(base_dir,'LR_analysis/'))) {dir.create( paste0(base_dir,'LR_analysis/'))}
openxlsx::write.xlsx(x=FC_pooled_cluster_4_sel_ligands_ordered, file = paste0(base_dir,'LR_analysis/','DE_cluster_ALL.SP_btypSel_RID2l_clExtended__FC_pooled_cluster_4_sel_ligands_ordered.xlsx'), overwrite = T) 


# Make list of top-10 ligands
ligands_HCMenriched_top10 = FC_pooled_cluster_4_sel_ligands_ordered$genename_short[1:10]
ligands_HCMenriched_all   = FC_pooled_cluster_4_sel_ligands_ordered$genename_short

# Now get the list that Louk & Bas made
louk_bas_ischemia_list = sort(c('Nppb', 'Calr', 'Mfge8', 'B2m', 'Adam10', 'Hmgb1', 'Hspg2', 'Sema7a', 'Vim', 'Gnai2', 'Gnas', 'Hsp90aa1', 'Hbegf', 'Tgm2', 'Pkm'))
# louk_bas_ischemia_list_humanized = convertMouseGeneList(louk_bas_ischemia_list)
louk_bas_ischemia_list_humanized = sort(c('HSP90AA1','HMGB1', 'TGM2', 'VIM', 'GNAI2', 'ADAM10', 'SEMA7A', 'GNAS', 'MFGE8', 'NPPB', 'PKM', 'B2M', 'HBEGF', 'CALR', 'HSPG2'))
    # Only Hsp90aa1 wasn't found, through "convertMouseGeneList", but I noticed HSP90AA1 is present in the dataset
ligands_both_hcm_ischemia = intersect(louk_bas_ischemia_list_humanized, ligands_HCMenriched_all)
ligands_both_only_hcm = ligands_HCMenriched_all[!(ligands_HCMenriched_all %in% louk_bas_ischemia_list_humanized)]
venn_simple_plot_mw(list(ischemia=louk_bas_ischemia_list_humanized,hcm=ligands_HCMenriched_all))
toString(ligands_both_hcm_ischemia)
toString(ligands_both_only_hcm)

# Create overview plot
FC_pooled_cluster_4_sel_ligands_orderedAlph = FC_pooled_cluster_4_sel_ligands_ordered[order(FC_pooled_cluster_4_sel_ligands_ordered$genename_short),]
p_list = lapply(1:(dim(FC_pooled_cluster_4_sel_ligands_orderedAlph)[1]), function(idx) {
    g =  FC_pooled_cluster_4_sel_ligands_orderedAlph$genename_short[idx]
    fc = 2^FC_pooled_cluster_4_sel_ligands_orderedAlph$avg_log2FC[idx]
    if (FC_pooled_cluster_4_sel_ligands_orderedAlph %in% ligands_both_hcm_ischemia) {mark=''} else {mark='*'}
    shorthand_seurat_custom_expr(seuratObject = current_analysis$ROOIJonly_RID2l_clExtended, gene_of_interest = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended, g), 
                                 pointsize = .5, custom_title = paste0(mark, g, '\n(FC=',round(fc,2),')'))})
wrap_plots(p_list)

VlnPlot(current_analysis$ROOIJonly_RID2l_clExtended, features = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended, 'GNAS'))
# FeaturePlot(current_analysis$ROOIJonly_RID2l_clExtended, features = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended, 'GNAS'))


##########

# Now also check for ligands that are differentially expressed between our clusters

DE_cluster_ROOIJonly_RID2l_clExtended

# Perform selection of both ligands and significant 
DE_cluster_ROOIJonly_RID2l_clExtended_sel_L = 
    lapply(1:length(DE_cluster_ROOIJonly_RID2l_clExtended), function(idx) {
        
        tbl=DE_cluster_ROOIJonly_RID2l_clExtended[[idx]]
    
        # First select sign. enriched genes
        tbl$genename_short = shorthand_cutname( rownames(tbl) )
        tbl_sel = tbl[tbl$p_val_adj<.05&tbl$avg_log2FC>0,]
        
        # Then sort ligands out of these
        tbl_sel_L = tbl_sel[tbl_sel$genename_short %in% ligand_list,]
        tbl_sel_L$fc = 2^tbl_sel_L$avg_log2FC
        tbl_sel_L$cluster = names(DE_cluster_ROOIJonly_RID2l_clExtended)[[idx]]
        
        return(tbl_sel_L)
        
    })

# Create overview table
DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined = 
    Reduce(rbind, DE_cluster_ROOIJonly_RID2l_clExtended_sel_L)
    # View(DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined[DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined$cluster!=6,]) 
    # View(DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined[!(DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined$cluster %in% c(5,6)),]) 

sapply(DE_cluster_ROOIJonly_RID2l_clExtended_sel_L, function(tbl) {
  paste0(tbl$genename_short, ' (FC=', round(2^tbl$avg_log2FC,1) , ')')
})
sapply(DE_cluster_ROOIJonly_RID2l_clExtended_sel_L, function(tbl) {
  tbl$genename_short  
})
 
VlnPlot(current_analysis$ROOIJonly_RID2l_clExtended, features = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended, 'VIM'))
VlnPlot(current_analysis$ROOIJonly_RID2l_clExtended, features = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended, 'SORBS1'))
VlnPlot(current_analysis$ROOIJonly_RID2l_clExtended, features = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended, 'LPL'))

shorthand_plotViolinBox_custom(current_analysis, analysis_name = 'ROOIJonly_RID2l_clExtended', cat_by = 'clusters_custom', cat_name = 'cluster', gene_of_interest = 'LPL', base_dir = base_dir)
VlnPlot(current_analysis$ROOIJonly_RID2l_clExtended, features = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended, 'LPL'))

# Also check ligands under-represented in cluster 5 (because those might be stress-related)
DE_cluster_ROOIJonly_RID2l_clExtended_5low = DE_cluster_ROOIJonly_RID2l_clExtended$'5'[DE_cluster_ROOIJonly_RID2l_clExtended$`5`$avg_log2FC<0&DE_cluster_ROOIJonly_RID2l_clExtended$`5`$p_val_adj<0.05,]
DE_cluster_ROOIJonly_RID2l_clExtended_5low$genename_short = shorthand_cutname(rownames(DE_cluster_ROOIJonly_RID2l_clExtended_5low))
DE_cluster_ROOIJonly_RID2l_clExtended_5low_L = DE_cluster_ROOIJonly_RID2l_clExtended_5low[DE_cluster_ROOIJonly_RID2l_clExtended_5low$genename_short %in% ligand_list,]
DE_cluster_ROOIJonly_RID2l_clExtended_5low_L$fc = 2^DE_cluster_ROOIJonly_RID2l_clExtended_5low_L$avg_log2FC
DE_cluster_ROOIJonly_RID2l_clExtended_5low_L$genename_short

# Venn cl.5 low and hcm.enriched
venn_simple_plot_mw( list(cl5.low = DE_cluster_ROOIJonly_RID2l_clExtended_5low_L$genename_short, hcm.enriched = ligands_HCMenriched_all) )

# Create overview plot of cl.5 reduced ligands
DE_cluster_ROOIJonly_RID2l_clExtended_5low_L_alph = DE_cluster_ROOIJonly_RID2l_clExtended_5low_L[order(DE_cluster_ROOIJonly_RID2l_clExtended_5low_L$genename_short),]
p_list_cl5.low = lapply(1:(dim(DE_cluster_ROOIJonly_RID2l_clExtended_5low_L_alph)[1]), function(idx) {
    g =  DE_cluster_ROOIJonly_RID2l_clExtended_5low_L_alph$genename_short[idx]
    fc = 2^DE_cluster_ROOIJonly_RID2l_clExtended_5low_L_alph$avg_log2FC[idx]
    if (DE_cluster_ROOIJonly_RID2l_clExtended_5low_L_alph %in% ligands_both_hcm_ischemia) {mark=''} else {mark='*'}
    shorthand_seurat_custom_expr(seuratObject = current_analysis$ROOIJonly_RID2l_clExtended, gene_of_interest = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended, g), 
                                 pointsize = .5, custom_title = paste0(mark, g, '\n(FC=',round(fc,2),')'))})
p=wrap_plots(p_list_cl5.low, nrow=2)
ggsave(plot=p, filename = paste0(base_dir,'LR_analysis/','ROOIJonly_RID2l_clExtended','ligands_cl5low.pdf'), 
                width=172/3-4, height=172/3/1.5-4, units='mm', device = cairo_pdf)

# Create overview of cluster-specific ligands
DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined_sortClus = DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined[order(DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined$cluster),]
DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined_sortClus_noCl6 = DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined_sortClus[DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined_sortClus$cluster!=6,]
p_list_cl5.low = lapply(1:(dim(DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined_sortClus_noCl6)[1]), function(idx) {
    g =  DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined_sortClus_noCl6$genename_short[idx]
    fc = 2^DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined_sortClus_noCl6$avg_log2FC[idx]
    cl = DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined_sortClus_noCl6$cluster[idx]
    if (DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined_sortClus_noCl6 %in% ligands_both_hcm_ischemia) {mark=''} else {mark='*'}
    shorthand_seurat_custom_expr(seuratObject = current_analysis$ROOIJonly_RID2l_clExtended, gene_of_interest = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended, g), 
                                 pointsize = .5, custom_title = paste0('Cl. ',cl, '; ', mark, g, '\n(FC=',round(fc,2),')'))})
p=wrap_plots(p_list_cl5.low)
p
ggsave(plot=p, filename = paste0(base_dir,'LR_analysis/','ROOIJonly_RID2l_clExtended','ligands_high_perCl.pdf'), 
                width=172-4, height=172-4, units='mm', device = cairo_pdf)

for (cl in unique(DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined_sortClus_noCl6$cluster)) {
   unique(DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined_sortClus_noCl6$cluster), function(cl) {
    glist = list(cl.X.up = DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined_sortClus_noCl6$genename_short[DE_cluster_ROOIJonly_RID2l_clExtended_sel_L_joined_sortClus_noCl6$cluster==cl], hcm.enriched = ligands_HCMenriched_all)
    names(glist) = c(paste0('cl.',cl,'.up'), 'HCM')
    venn_simple_plot_mw( glist )
}

################################################################################
# OK, so now the most important question, are receptors for these ligands also
# expressed in other cell types of the heart?
    
# First define ligand lists of interest
ligands_of_interest = lapply(DE_cluster_ROOIJonly_RID2l_clExtended_sel_L[1:5], function(X){sort(X$genename_short)})
names(ligands_of_interest) = paste0('cl.',1:5)
ligands_of_interest$cl.5.down = DE_cluster_ROOIJonly_RID2l_clExtended_5low_L_alph$genename_short
ligands_of_interest$HCM_up = sort(ligands_HCMenriched_all)

# Then look up receptors for those clusters
receptors_of_interest =
    lapply(ligands_of_interest, function(lst)
              sapply(lst, function(L) {ligand_receptor_pairs_df$receptor[ligand_receptor_pairs_df$ligand==L]})
    )
    

# Save the data of interest
save(list=c('ligands_of_interest', 'receptors_of_interest'), file = paste0(base_dir,'LR_analysis/LR_of_interest.Rdata'))


# Then go to HPC to load all data, and show bar-plots to gauge receptor activity    
# \/  \/  \/    
    
################################################################################    
    
    
    
    
    
    
    
    


