
venn_simple_plot_mw(list(
    TTN=gene_lists_customcorrelated_reorganized$`posCorrWith_ENSG00000155657:TTN`,
    XIRP2=gene_lists_customcorrelated_reorganized$`posCorrWith_ENSG00000163092:XIRP2`))

sum(gene_lists_customcorrelated_reorganized$`posCorrWith_ENSG00000155657:TTN` %in% gene_lists_customcorrelated_reorganized$`posCorrWith_ENSG00000163092:XIRP2`)
# sanity: sum(gene_lists_customcorrelated_reorganized$`posCorrWith_ENSG00000163092:XIRP2` %in% gene_lists_customcorrelated_reorganized$`posCorrWith_ENSG00000155657:TTN`)
length(gene_lists_customcorrelated_reorganized$`posCorrWith_ENSG00000163092:XIRP2`)


################################################################################
# Comparison between the modules and regulons

# Now load the SCENIC regulons and custom modules
# SCENIC regulons
load(paste0(base_dir, 'Rdata/SCENIC_reg_top_genes_sorted_full.Rdata')) # SCENIC_reg_top_genes_sorted_full
# Custom modules
ANALYSIS_NAME = "ROOIJonly_RID2l"
load(paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_sorted_shortname.Rdata')) # core_regulons_sorted_shortname

###

# 
# Note, for HPC custom analyses, regulon script can be sourced as
# desired_command_regulon='dummy'; source('/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/HCM_SCS_2021_06_SeuratRevised_Regulon_v3.R')
#
# LOCAL=1; desired_command_regulon='dummy'; source('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Projects/HCM_SCS_2021_06_SeuratRevised_Regulon_v3.R')

###

core_regulons_sorted_shortname_=core_regulons_sorted_shortname
names(core_regulons_sorted_shortname_) = paste0(gsub('s\\.R','_Module',names(core_regulons_sorted_shortname_)),'_')
currently_pooled_regulons = c(SCENIC_reg_top_genes_sorted_full, core_regulons_sorted_shortname_)
regulon_overlap_heatmap(currently_pooled_regulons, base_dir, 'custom_pool_SCENIC_customModules', MYTREECUTTINGHEIGHT=2, myfontsize=8, makeallheatmaps=T, cutree_k=NULL)
    
    
    
    
    