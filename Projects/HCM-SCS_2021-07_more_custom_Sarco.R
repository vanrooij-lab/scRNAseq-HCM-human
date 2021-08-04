
library(scales)
library(patchwork)





sarcomere_GO_table_ = read.table('/Users/m.wehrens/Data/_2019_02_HCM_SCS/misc_analyses/GO-sarcomere_associated_genes.txt', sep = '\t')
names(sarcomere_GO_table_)=c('bioentity_label','GO','organism')
    # source: http://amigo.geneontology.org/amigo/term/GO:0030017

#sarcomere_GO_table = sarcomere_GO_table_[sarcomere_GO_table_$organism=='NCBITaxon:9606'&sarcomere_GO_table_$GO=='GO:0030017',]
sarcomere_GO_table = sarcomere_GO_table_[sarcomere_GO_table_$organism=='NCBITaxon:9606',]

sarcomere_genes = unique(sarcomere_GO_table$bioentity_label)
sarcomere_genes
toString(sarcomere_genes)

sarcomere_genes_maya = sort(c('MYBPC3', 'MYH6', 'TTN', 'RYR2', 'CMYA5', 'XIRP2', 'ACTN1', 'TNNT2', 'TPM1', 'MYL2', 'MYL3', 'TNNI3', 'ACTC1')) # CMYA3-->XIRP2

sarcomere_genes_maya[!(sarcomere_genes_maya %in% sarcomere_genes)]
serca or atp2a
# c('SERCA', 'ATP2A') %in% sarcomere_genes

venn_simple_plot_mw(venn_list = list(sarco1=sarcomere_genes, sarco2=sarcomere_genes_maya))

shorthand_seurat_fullgenename(sarcomere_genes_maya)

plot_list_sarco=lapply(sarcomere_genes_maya, function(x) {shorthand_seurat_custom_expr(x)})

wrap_plots(plot_list_sarco)*theme(legend.position='none')

shorthand_seurat_custom_expr('MALAT1')

#####

VlnPlot(current_analysis$ROOIJonly_RID2l, features = shorthand_seurat_fullgenename('MALAT1'), )

#####


# Plot some of anne's genes
wrap_plots(lapply(c('ATP2A2'), function(x) {shorthand_seurat_custom_expr(x)}))*theme(legend.position='none')





