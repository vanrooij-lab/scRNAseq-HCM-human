
library(scales)
library(patchwork)


################################################################################


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



################################################################################
# Checking cell cycle genes

# Gene Set: KEGG_CELL_CYCLE
# https://www.gsea-msigdb.org/gsea/msigdb/cards/KEGG_CELL_CYCLE
cellcycle_genes = c('ABL1', 'ANAPC1', 'ANAPC10', 'ANAPC11', 'ANAPC13', 'ANAPC2', 'ANAPC4', 'ANAPC5', 'ANAPC7', 'ATM', 'ATR', 'BUB1', 'BUB1B', 'BUB3', 'CCNA1', 'CCNA2', 'CCNB1', 'CCNB2', 'CCNB3', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CCNE2', 'CCNH', 'CDC14A', 'CDC14B', 'CDC16', 'CDC20', 'CDC23', 'CDC25A', 'CDC25B', 'CDC25C', 'CDC26', 'CDC27', 'CDC45', 'CDC6', 'CDC7', 'CDK1', 'CDK2', 'CDK4', 'CDK6', 'CDK7', 'CDKN1A', 'CDKN1B', 'CDKN1C', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CDKN2D', 'CHEK1', 'CHEK2', 'CREBBP', 'CUL1', 'DBF4', 'E2F1', 'E2F2', 'E2F3', 'E2F4', 'E2F5', 'EP300', 'ESPL1', 'FZR1', 'GADD45A', 'GADD45B', 'GADD45G', 'GSK3B', 'HDAC1', 'HDAC2', 'MAD1L1', 'MAD2L1', 'MAD2L2', 'MCM2', 'MCM3', 'MCM4', 'MCM5', 'MCM6', 'MCM7', 'MDM2', 'MYC', 'ORC1', 'ORC2', 'ORC3', 'ORC4', 'ORC5', 'ORC6', 'PCNA', 'PKMYT1', 'PLK1', 'PRKDC', 'PTTG1', 'PTTG2', 'RAD21', 'RB1', 'RBL1', 'RBL2', 'RBX1', 'SFN', 'SKP1', 'SKP1P2', 'SKP2', 'SMAD2', 'SMAD3', 'SMAD4', 'SMC1A', 'SMC1B', 'SMC3', 'STAG1', 'STAG2', 'TFDP1', 'TFDP2', 'TGFB1', 'TGFB2', 'TGFB3', 'TP53', 'TTK', 'WEE1', 'WEE2', 'YWHAB', 'YWHAE', 'YWHAG', 'YWHAH', 'YWHAQ', 'YWHAZ', 'ZBTB17')

p=shorthand_seurat_custom_expr(seuratObject = current_analysis[[ANALYSIS_NAME]], 
                                                    gene_of_interest = cellcycle_genes, textsize=8, pointsize=1, 
                                                    custom_title = 'Cell cycle', mymargin = .5, zscore = T)

ggsave(filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_9_custom_CellCycle.pdf'), 
        plot = p, width=30, height=30, units='mm') # 184.6/3*2-4
    
