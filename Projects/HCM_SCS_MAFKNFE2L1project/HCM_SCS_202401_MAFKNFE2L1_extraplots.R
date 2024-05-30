
################################################################################

# Load our data locally
#
# See also the file 
# ~/Documents/git_repos/SCS_More_analyses/Projects/HCM_SCS_2021_07_custom_analysis_ROOIJ.R
# for loading the object
#
if (F) {
    LOCAL=1; script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Projects/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')
    
    library(clustree)
    library(pheatmap)
    library(wesanderson) # only for some custom plots

        
    # DATASET_NAME='ROOIJonly.sp.bt_RID2l'
    DATASET_NAME='ROOIJonly.sp.bt_RID2l_clExtended' # includes some extra analysis
    DATASET_NAME_prev = 'ROOIJonly.sp.bt_RID2l'
    
    # load the file
    #if (!exists('current_analysis')) {current_analysis = list()}
    #current_analysis[[DATASET_NAME]] =
    #    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))
    # load the clExtended
    if (!exists('current_analysis')) {current_analysis=list()}
    current_analysis[[DATASET_NAME]] =
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))

    # OR ALTERNATIVELY
    # current_analysis[[DATASET_NAME]] = readRDS(paste0(base_dir,'Rdata/2024/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.Rds'))
    
    # Get TF names
    load(file=paste0(base_dir,'Rdata/',DATASET_NAME_prev,'__SCENIC_reg_top_genes_sorted_full.Rdata')) # SCENIC_reg_top_genes_sorted_full
    SCENIC_TF_list = names(SCENIC_reg_top_genes_sorted_full)
    print(paste0('TF names from file: ',toString(SCENIC_TF_list))) # ZEB1, FOXN3, NFE2L1, MAFK, CLOCK, REST, CEBPZ, MEF2D, MITF, MEF2A, ETV1, CREB1, MXI1, JUND, YY1, SREBF2, CEBPB, SRF, NR3C1, ESRRA, USF2, ZNF91

}

# Create dirs if applicable
if(!dir.exists(paste0(base_dir,'Rplots/2024/MAFKNFE2L1/'))){dir.create(paste0(base_dir,'Rplots/2024/MAFKNFE2L1/'), recursive=T)}

# Additional directories I'm using
HFTF_datadir = '/Users/m.wehrens/Data/2023_08_HFTFs/analysis_202308/'

################################################################################
# Some administrative things that only needed doing once (you can ignore this)

if (F) {
    # Let's now also save things we use to Rds format
    if(!dir.exists(paste0(base_dir,'Rdata/2024/'))){dir.create(paste0(base_dir,'Rdata/2024/'))}
    saveRDS(current_analysis[[DATASET_NAME]], paste0(base_dir,'Rdata/2024/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.Rds'))
}


################################################################################
# Some custom functions

shorthand24_histogram = function(SO, G_full)  {
    
    df_expr = data.frame(expression_count = SO@assays$RNA@data[G_full,])
    
    ggplot(df_expr, aes(x=expression_count))+
        geom_freqpoly()+theme_bw()+give_better_textsize_plot(8)+
        xlab('Gene read count (normalized)')+ylab('Times observed')+ggtitle(shorthand_cutname( G ) )
    
}

shorthand24_expression = function(SO, G_full)  {
    
    expression_count = SO@assays$RNA@data[G_full,]
    return(expression_count)
    
}


################################################################################
# Now make some plots that we were interested in..

# Show MAFK/NFE2L1 expression (again)
G='MAFK'
G_full = shorthand_seurat_fullgenename(current_analysis[[DATASET_NAME]], G)
p=FeaturePlot(current_analysis[[DATASET_NAME]], G_full, 
                    cols = pals::kovesi.diverging_rainbow_bgymr_45_85_c67(100), pt.size = .25)+ # pals::viridis(100)
                theme_void()+theme(legend.position = 'none', panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
                ggtitle(shorthand_cutname(G))
p
ggsave(plot=p, filename = paste0(base_dir,'Rplots/2024/MAFKNFE2L1/',ANALYSIS_NAME,'_expression_',G,'.pdf'), 
                height=50, width=50, units='mm', device=cairo_pdf)
# Show MAFK/NFE2L1 expression (again)
G='NFE2L1'
G_full = shorthand_seurat_fullgenename(current_analysis[[DATASET_NAME]], G)
p=FeaturePlot(current_analysis[[DATASET_NAME]], G_full, 
                    cols = pals::kovesi.diverging_rainbow_bgymr_45_85_c67(100), pt.size = .25)+theme_void()+
                theme(legend.position = 'none', panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
                ggtitle(shorthand_cutname(G))
p
ggsave(plot=p, filename = paste0(base_dir,'Rplots/2024/MAFKNFE2L1/',ANALYSIS_NAME,'_expression_',G,'.pdf'), 
                height=50, width=50, units='mm', device=cairo_pdf)
print(paste0(base_dir,'Rplots/2024/MAFKNFE2L1/'))


### UMAPs | # Now show this on somewhat customized UMAP
# First load data
UMAP1=current_analysis[[DATASET_NAME]]@reductions$umap@cell.embeddings[,1]
UMAP2=current_analysis[[DATASET_NAME]]@reductions$umap@cell.embeddings[,2]
G='MAFK'
G_full = shorthand_seurat_fullgenename(current_analysis[[DATASET_NAME]], G)
expr_MAFK = shorthand24_expression(current_analysis[[DATASET_NAME]], G_full)
threshold_MAFKhigh = quantile(expr_MAFK, .9)
G='NFE2L1'
G_full = shorthand_seurat_fullgenename(current_analysis[[DATASET_NAME]], G)
expr_NFE2L1 = shorthand24_expression(current_analysis[[DATASET_NAME]], G_full)
threshold_NFE2L1high = quantile(expr_NFE2L1, .9)

### Histograms | Check MAFK expression distribution
G='MAFK'
G_full = shorthand_seurat_fullgenename(current_analysis[[DATASET_NAME]], G)
p=shorthand24_histogram(current_analysis[[DATASET_NAME]], G_full)+
    geom_vline(xintercept = threshold_MAFKhigh, linetype='dashed')
p
ggsave(plot=p, filename = paste0(base_dir,'Rplots/2024/MAFKNFE2L1/',ANALYSIS_NAME,'_histexpr_',G,'.pdf'), 
                height=50, width=50, units='mm', device=cairo_pdf)
# And NFE2L1
G='NFE2L1'
G_full = shorthand_seurat_fullgenename(current_analysis[[DATASET_NAME]], G)
p=shorthand24_histogram(current_analysis[[DATASET_NAME]], G_full)+
    geom_vline(xintercept = threshold_NFE2L1high, linetype='dashed')
p
ggsave(plot=p, filename = paste0(base_dir,'Rplots/2024/MAFKNFE2L1/',ANALYSIS_NAME,'_histexpr_',G,'.pdf'), 
                height=50, width=50, units='mm', device=cairo_pdf)

# Show umap
expr_overview = rep(NA, length(expr_MAFK))
expr_overview[(expr_MAFK<=threshold_MAFKhigh ) & (expr_NFE2L1<=threshold_NFE2L1high)]='none'
expr_overview[(expr_MAFK>threshold_MAFKhigh  ) & (expr_NFE2L1<=threshold_NFE2L1high)]='MAFK'
expr_overview[(expr_NFE2L1>threshold_NFE2L1high) & (expr_MAFK<=threshold_MAFKhigh)]='NFE2L1'
expr_overview[(expr_NFE2L1>threshold_NFE2L1high) & (expr_MAFK>threshold_MAFKhigh)]='both'
df_toplot = data.frame(UMAP1=UMAP1, UMAP2=UMAP2, expr_overview=as.factor(expr_overview))
p=ggplot(df_toplot, aes(x=UMAP1, y=UMAP2, color=expr_overview))+
    geom_point(size=.25)+theme_bw()+
    scale_colour_manual(# name = "XXX", 
                      #labels = c('none','MAFK','NFE2L1','both'),
                      values = c(none="#cccccc",MAFK="mediumvioletred",NFE2L1="mediumturquoise",both="black"))+
    give_better_textsize_plot(8)+
    theme_void()+theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
    
p
ggsave(plot=p, filename = paste0(base_dir,'Rplots/2024/MAFKNFE2L1/',ANALYSIS_NAME,'_expression_','custom_NFE2L1_n_MAFK','.pdf'), 
                height=50, width=80, units='mm', device=cairo_pdf)

# How many cells have both?
sum(expr_overview=='both')
# Expected based on random:
length(expr_overview)*.1^2
# and for .8 treshold; length(expr_overview)*.2^2

# What is the distribution of cells? (scatter plot MAFK vs. NFE2L1)
# note that normalization of data didn't involve logfold transformatino (see methods my paper, RC normaliztaion was used)
ggplot(data.frame(expr_MAFK=expr_MAFK, expr_NFE2L1=expr_NFE2L1), aes(x=log10(.1+expr_NFE2L1), y=log10(.1+expr_MAFK)))+
    geom_point()+theme_bw()

table(expr_overview)
round(table(expr_overview)/length(expr_overview)*100,1)
.077^2

################################################################################
# Now 

# Identify cells of interest
cellpopulation_MAFK_n_NFE2L1 = rep(NA, length(expr_overview))
cellpopulation_MAFK_n_NFE2L1[expr_overview=='both']='positive'
cellpopulation_MAFK_n_NFE2L1[expr_overview!='both']='negative'
current_analysis[[DATASET_NAME]]$cellpopulation_MAFK_n_NFE2L1 = cellpopulation_MAFK_n_NFE2L1

# Identify markers from NFE2L1+ MAFK+ cells
FindMarkers_out_NFE2L1_n_MAFK_positive = 
              FindMarkers(current_analysis[[DATASET_NAME]],
                          ident.1 = 'positive',
                          ident.2 = 'negative',
                          group.by = 'cellpopulation_MAFK_n_NFE2L1')

# What about any NFE2L1+?
cellpopulation_NFE2L1 = rep(NA, length(expr_overview))
cellpopulation_NFE2L1[expr_NFE2L1 >  threshold_NFE2L1high]='positive'
cellpopulation_NFE2L1[expr_NFE2L1 <= threshold_NFE2L1high]='negative'
current_analysis[[DATASET_NAME]]$cellpopulation_NFE2L1 = cellpopulation_NFE2L1
FindMarkers_out_NFE2L1_positive = 
              FindMarkers(current_analysis[[DATASET_NAME]],
                          ident.1 = 'positive',
                          ident.2 = 'negative',
                          group.by = 'cellpopulation_NFE2L1')

# What about any MAFK+?
cellpopulation_MAFK = rep(NA, length(expr_overview))
cellpopulation_MAFK[expr_NFE2L1 >  threshold_MAFKhigh]='positive'
cellpopulation_MAFK[expr_NFE2L1 <= threshold_MAFKhigh]='negative'
current_analysis[[DATASET_NAME]]$cellpopulation_MAFK = cellpopulation_MAFK
FindMarkers_out_MAFK_positive = 
              FindMarkers(current_analysis[[DATASET_NAME]],
                          ident.1 = 'positive',
                          ident.2 = 'negative',
                          group.by = 'cellpopulation_MAFK')

# Now the two alone
current_analysis[[DATASET_NAME]]$overview_MAFK_NFE2L1_highlow = expr_overview
FindMarkers_out_MAFKonly_positive = 
              FindMarkers(current_analysis[[DATASET_NAME]],
                          ident.1 = 'MAFK',
                          ident.2 = 'none',
                          group.by = 'overview_MAFK_NFE2L1_highlow')
FindMarkers_out_NFE2L1only_positive = 
              FindMarkers(current_analysis[[DATASET_NAME]],
                          ident.1 = 'NFE2L1',
                          ident.2 = 'none',
                          group.by = 'overview_MAFK_NFE2L1_highlow')



# Now overlap those with the hits from the RNA-seq OE of NFE2L1/MAFK
genes_upin_NFE2L1_n_MAFK = shorthand_cutname( rownames(
    FindMarkers_out_NFE2L1_n_MAFK_positive[FindMarkers_out_NFE2L1_n_MAFK_positive$p_val_adj<0.05,]) )
genes_upin_NFE2L1 = shorthand_cutname( rownames(
    FindMarkers_out_NFE2L1_positive[FindMarkers_out_NFE2L1_positive$p_val_adj<0.05,]) )
genes_upin_MAFK = shorthand_cutname( rownames(
    FindMarkers_out_MAFK_positive[FindMarkers_out_MAFK_positive$p_val_adj<0.05,]) )
genes_upin_NFE2L1only = shorthand_cutname( rownames(
    FindMarkers_out_NFE2L1only_positive[FindMarkers_out_NFE2L1only_positive$p_val_adj<0.05,]) )
genes_upin_MAFKonly = shorthand_cutname( rownames(
    FindMarkers_out_MAFKonly_positive[FindMarkers_out_MAFKonly_positive$p_val_adj<0.05,]) )

FindMarkers_out_NFE2L1_n_MAFK_positive$sym = shorthand_cutname(rownames(FindMarkers_out_NFE2L1_n_MAFK_positive))
FindMarkers_out_NFE2L1_positive$sym = shorthand_cutname(rownames(FindMarkers_out_NFE2L1_positive))
FindMarkers_out_MAFK_positive$sym = shorthand_cutname(rownames(FindMarkers_out_MAFK_positive))
FindMarkers_out_NFE2L1only_positive$sym = shorthand_cutname(rownames(FindMarkers_out_NFE2L1only_positive))
FindMarkers_out_MAFKonly_positive$sym = shorthand_cutname(rownames(FindMarkers_out_MAFKonly_positive))
openxlsx::write.xlsx(FindMarkers_out_NFE2L1_n_MAFK_positive[order(FindMarkers_out_NFE2L1_n_MAFK_positive$avg_log2FC,decreasing = T),], 
                     paste0(base_dir,'Rplots/2024/MAFKNFE2L1/', 'xlsx_FindMarkers_out_NFE2L1_n_MAFK_positive.xlsx'), row.names=T)
openxlsx::write.xlsx(FindMarkers_out_NFE2L1_positive[order(FindMarkers_out_NFE2L1_positive$avg_log2FC,decreasing = T),], 
                     paste0(base_dir,'Rplots/2024/MAFKNFE2L1/', 'xlsx_FindMarkers_out_NFE2L1_positive.xlsx'), row.names=T)
openxlsx::write.xlsx(FindMarkers_out_MAFK_positive[order(FindMarkers_out_MAFK_positive$avg_log2FC,decreasing = T),], 
                     paste0(base_dir,'Rplots/2024/MAFKNFE2L1/', 'xlsx_FindMarkers_out_MAFK_positive.xlsx'), row.names=T)
openxlsx::write.xlsx(FindMarkers_out_NFE2L1only_positive[order(FindMarkers_out_NFE2L1only_positive$avg_log2FC,decreasing = T),], 
                     paste0(base_dir,'Rplots/2024/MAFKNFE2L1/', 'xlsx_FindMarkers_out_NFE2L1only_positive.xlsx'), row.names=T)
openxlsx::write.xlsx(FindMarkers_out_MAFKonly_positive[order(FindMarkers_out_MAFKonly_positive$avg_log2FC,decreasing = T),], 
                     paste0(base_dir,'Rplots/2024/MAFKNFE2L1/', 'xlsx_FindMarkers_out_MAFKonly_positive.xlsx'), row.names=T)

# Comparing the two conditions from the scRNA-seq data
venn_simple_plot_mw(list(both= genes_upin_NFE2L1_n_MAFK, NFE2L1=genes_upin_NFE2L1))
venn_simple_plot_mw(list(both= genes_upin_NFE2L1_n_MAFK, MAFK=genes_upin_MAFK))
venn_simple_plot_mw(list(NFE2L1=genes_upin_NFE2L1, MAFK=genes_upin_MAFK))
venn_simple_plot_mw(list(NFE2L1only=genes_upin_NFE2L1only, MAFKonly=genes_upin_MAFKonly))

# Looking at how MAFK+ and NFE2L1+ alone compare to MAFK+NFE2L1+ selected cells
genes_upin_either = unique(c(genes_upin_MAFK, genes_upin_NFE2L1))
venn_simple_plot_mw(list(both= genes_upin_NFE2L1_n_MAFK, either=genes_upin_either))
    
genes_upin_eitheralone = unique(c(genes_upin_MAFKonly, genes_upin_NFE2L1only))
venn_simple_plot_mw(list(both= genes_upin_NFE2L1_n_MAFK, either_alone=genes_upin_eitheralone))
genes_unique_NFE2L1_n_MAFK = genes_upin_NFE2L1_n_MAFK[!(genes_upin_NFE2L1_n_MAFK %in% genes_upin_eitheralone)]

# Compare with MAFK+NFE2L1+ unique genes from OE experiment
genes_interxMAFKNFE2L1_OEexp =
    readRDS(paste0(HFTF_datadir,'lists/', 'genes_interxMAFKNFE2L1.Rds'))

# comparing genes "unique" to MAFK+NFE2L1+ cells (by measure of enrichcment, might depend cutoff)
venn_simple_plot_mw( list(HCM=genes_unique_NFE2L1_n_MAFK, OE=genes_interxMAFKNFE2L1_OEexp) )

# Comparing genes that are simply up in the MAFK+NFE2L1+ cells in scRNA-seq, compared OE experiment
venn_simple_plot_mw( list(HCM=genes_upin_NFE2L1_n_MAFK, OE=genes_interxMAFKNFE2L1_OEexp) )

# More extensive lists
df_collected_deseq2_results_OEexp = 
    readRDS(paste0(HFTF_datadir,'Rdata/HFTFs-202308_DESEQ2_list.Rds')) # df_collected_deseq2_results
df_DESEQ_MAFK_n_NFE2L1 = df_collected_deseq2_results_OEexp$`MAFK+NFE2L1.vs.GFP`
df_DESEQ_MAFK_n_NFE2L1_sel = df_DESEQ_MAFK_n_NFE2L1[df_DESEQ_MAFK_n_NFE2L1$padj<0.05&df_DESEQ_MAFK_n_NFE2L1$log2FoldChange>0,]
df_DESEQ_MAFK_n_NFE2L1_sel_order = df_DESEQ_MAFK_n_NFE2L1_sel[order(df_DESEQ_MAFK_n_NFE2L1_sel$log2FoldChange, decreasing = T), ]

genes_MAFK_n_NFE2L1_OEexp = 
    shorthand_cutname(rownames(df_DESEQ_MAFK_n_NFE2L1_sel_order))#[1:100])

# Genes of interest detected from both analyses
# Comparing MAFK+NFE2L1+ cells from RNA-seq analysis to genes up in the OE experiment
venn_simple_plot_mw( list(HCM=genes_upin_NFE2L1_n_MAFK, OE=genes_MAFK_n_NFE2L1_OEexp) )
genes_upin_NFE2L1_n_MAFK[genes_upin_NFE2L1_n_MAFK %in% genes_MAFK_n_NFE2L1_OEexp]
toString(genes_upin_NFE2L1_n_MAFK[genes_upin_NFE2L1_n_MAFK %in% genes_MAFK_n_NFE2L1_OEexp])

# 15k genes (findMarkers performs the test on the "data" slot),
# however, only genes are tested that have are found in >1% of cells
# (per default)
dim(current_analysis$ROOIJonly.sp.bt_RID2l_clExtended@assays$RNA@data)
sum(rowMeans(current_analysis$ROOIJonly.sp.bt_RID2l_clExtended@assays$RNA@data>0)>0.01)
# so 11k genes remain..
#
# From OE experiment analysis, data_container$groups$PCR.10$dt$ndata, 
# 2242 genes detected after QC

# WRONG phrasing:
# when drawing 175 genes from 2292 genes, what's the probability
# distribution of expected hits, given that there's (115/15240=)
# 0.0075 fraction of positive genes.

# Perhaps a better assessemnt:
# The bulk RNA-seq data only contained 2242 genes after QC, and only 1039 
# after DESeq2 stringent cutoff, of those 175 (17%) were deemed to show 
# significant differential expression. Of the ±11K genes considered in the 
# single cell findmarkers output, only 115 were considered significant (±1%). 
# What is the probability the overlap between those two is 12 if random genes 
#were drawn? Given that all 175 OE genes are present in the complete HCM data, 
# the probability to draw ≥12 of those from the HCM set in 115 trials is  3.6e-7. 
# However, this might not hold true, as 

# This shows what happens if we perform 100 experiments, what the probability
# is that there are x successes.
ggplot(data.frame(x=seq(0,20), p=dbinom(seq(0,20), 100, .1)), aes(x=x,y=p))+
    geom_point()
# Now we can do something similar
ggplot(data.frame(x=seq(0,50), p=dbinom(seq(0,50), 115, 175/11e3)), aes(x=x,y=p))+
    geom_point()
qbinom(.95, 175, 115/11e3) # meaning that there's a <5% chance we have 3+ hits

# what's the chance on 12 or more successes,
# that equals the chance on NOT getting between 0 to 11 successes.
1-sum(dbinom(0:11, 115, 175/11e3))
1-pbinom(11, 115, 175/11e3)
# so that chance is pretty small
# allthough assuming higher numbers, ie from the bulk seq chance
1-pbinom(11, 115, 175/2242)

# chance on no hits in 12 trials:
# (but this is something else than asking how
# many hits you'll get after 175 tries.
(1-0.0075)^12

# Sanity check, which indeed agrees with the internet calculator
# https://stattrek.com/online-calculator/binomial using the same numbers
# "Cumulative probability: P(X≥12)"
1-pbinom(11, 175, .1)

################################################################################
# Relative abundance of MAF and NFE2L genes

# Let's check the abundances of the MAF genes in my data;
myMAFgenes = rownames(current_analysis$ROOIJonly.sp.bt_RID2l_clExtended)[grepl(':MAF',rownames(current_analysis$ROOIJonly.sp.bt_RID2l_clExtended))]
df_MAF_expression = as.data.frame(apply(current_analysis$ROOIJonly.sp.bt_RID2l_clExtended@assays$RNA@data[myMAFgenes,], 1, mean))
colnames(df_MAF_expression) = 'mean_expression'
df_MAF_expression$gene_sym = shorthand_cutname(rownames(df_MAF_expression))
p=ggplot(df_MAF_expression, aes(x=gene_sym, y=mean_expression))+
    geom_bar(stat="identity")+theme_bw()+give_better_textsize_plot(8)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                                    strip.text.x = element_text(angle = 90))
p
ggsave(plot=p, filename = paste0(base_dir,'Rplots/2024/MAFKNFE2L1/',ANALYSIS_NAME,'_meanexpression_','MAFgenes','.pdf'), 
                height=50, width=50, units='mm', device=cairo_pdf)

# Let's check the abundances of the MAF genes in my data;
myNFELgenes = rownames(current_analysis$ROOIJonly.sp.bt_RID2l_clExtended)[grepl(':NFE2L',rownames(current_analysis$ROOIJonly.sp.bt_RID2l_clExtended))]
df_NFE2L_expression = as.data.frame(apply(current_analysis$ROOIJonly.sp.bt_RID2l_clExtended@assays$RNA@data[myNFELgenes,], 1, mean))
colnames(df_NFE2L_expression) = 'mean_expression'
df_NFE2L_expression$gene_sym = shorthand_cutname(rownames(df_NFE2L_expression))
p=ggplot(df_NFE2L_expression, aes(x=gene_sym, y=mean_expression))+
    geom_bar(stat="identity")+theme_bw()+give_better_textsize_plot(8)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                                    strip.text.x = element_text(angle = 90))
p
ggsave(plot=p, filename = paste0(base_dir,'Rplots/2024/MAFKNFE2L1/',ANALYSIS_NAME,'_meanexpression_','NFE2genes','.pdf'), 
                height=50, width=30, units='mm', device=cairo_pdf)

################################################################################
# Let's check how stress genes relate to MAFK/NFE2L1 expressoin

G='NPPA'
G_full = shorthand_seurat_fullgenename(current_analysis[[DATASET_NAME]], G)
expr_NPPA = shorthand24_expression(current_analysis[[DATASET_NAME]], G_full)
G='MYH7'
G_full = shorthand_seurat_fullgenename(current_analysis[[DATASET_NAME]], G)
expr_MYH7 = shorthand24_expression(current_analysis[[DATASET_NAME]], G_full)
G='ACTA1'
G_full = shorthand_seurat_fullgenename(current_analysis[[DATASET_NAME]], G)
expr_ACTA1 = shorthand24_expression(current_analysis[[DATASET_NAME]], G_full)

# Create df
df_expr_stress=data.frame(expr_MAFK=expr_MAFK, expr_NFE2L1=expr_NFE2L1, 
                          expr_NPP=expr_NPPA, expr_MYH7=expr_MYH7, expr_ACTA1=expr_ACTA1)

# NPPA plot
p=ggplot(df_expr_stress, 
       aes(x=log10(.1+expr_NFE2L1), y=log10(.1+expr_MAFK), color=log10(.1+expr_NPPA)))+
    geom_point(size=.5)+theme_bw()+
    scale_colour_viridis_c()+give_better_textsize_plot(8)+
    theme(legend.key.size = unit(2, 'mm'))
ggsave(plot=p, filename = paste0(base_dir,'Rplots/2024/MAFKNFE2L1/',ANALYSIS_NAME,'_custom_MAFK_NFE2L1_scatter_NPPA.pdf'), 
                height=50, width=100, units='mm', device=cairo_pdf)
    #scale_colour_gradientn(colours = c("blue", "white", "red"))
# MYH7 plot
p=ggplot(df_expr_stress, 
       aes(x=log10(.1+expr_NFE2L1), y=log10(.1+expr_MAFK), color=log10(.1+expr_MYH7)))+
    geom_point(size=.5)+theme_bw()+
    scale_colour_viridis_c()+give_better_textsize_plot(8)+
    theme(legend.key.size = unit(2, 'mm'))
ggsave(plot=p, filename = paste0(base_dir,'Rplots/2024/MAFKNFE2L1/',ANALYSIS_NAME,'_custom_MAFK_NFE2L1_scatter_MYH7.pdf'), 
                height=50, width=100, units='mm', device=cairo_pdf)
    #scale_colour_gradientn(colours = c("blue", "white", "red"))
# ACTA1 plot
p=ggplot(df_expr_stress, 
       aes(x=log10(.1+expr_NFE2L1), y=log10(.1+expr_MAFK), color=log10(.1+expr_ACTA1)))+
    geom_point(size=.5)+theme_bw()+
    scale_colour_viridis_c()+give_better_textsize_plot(8)+
    theme(legend.key.size = unit(2, 'mm'))
ggsave(plot=p, filename = paste0(base_dir,'Rplots/2024/MAFKNFE2L1/',ANALYSIS_NAME,'_custom_MAFK_NFE2L1_scatter_ACTA1.pdf'), 
                height=50, width=100, units='mm', device=cairo_pdf)
    #scale_colour_gradientn(colours = c("blue", "white", "red"))

# custom levels
current_analysis$ROOIJonly.sp.bt_RID2l_clExtended$overview_MAFK_NFE2L1_highlow_fct=
    factor(current_analysis$ROOIJonly.sp.bt_RID2l_clExtended$overview_MAFK_NFE2L1_highlow, 
           levels=c('none','NFE2L1','MAFK','both'))
MYPLOTSIZEFACTOR=5
G='NPPA'
G_full = shorthand_seurat_fullgenename(current_analysis[[DATASET_NAME]], G)
p=VlnPlot(object = current_analysis$ROOIJonly.sp.bt_RID2l_clExtended, features = G_full, group.by='overview_MAFK_NFE2L1_highlow_fct')+
    give_better_textsize_plot(8*MYPLOTSIZEFACTOR)
ggsave(plot=p, filename = paste0(base_dir,'Rplots/2024/MAFKNFE2L1/',ANALYSIS_NAME,'_custom_MAFK_NFE2L1_Vln_NPPA.pdf'), 
                height=50*MYPLOTSIZEFACTOR, width=100*MYPLOTSIZEFACTOR, units='mm', device=cairo_pdf)
G='ACTA1'
G_full = shorthand_seurat_fullgenename(current_analysis[[DATASET_NAME]], G)
p=VlnPlot(object = current_analysis$ROOIJonly.sp.bt_RID2l_clExtended, features = G_full, group.by='overview_MAFK_NFE2L1_highlow_fct')+
    give_better_textsize_plot(8*MYPLOTSIZEFACTOR)
ggsave(plot=p, filename = paste0(base_dir,'Rplots/2024/MAFKNFE2L1/',ANALYSIS_NAME,'_custom_MAFK_NFE2L1_Vln_ACTA1.pdf'), 
                height=50*MYPLOTSIZEFACTOR, width=100*MYPLOTSIZEFACTOR, units='mm', device=cairo_pdf)
G='MYH7'
G_full = shorthand_seurat_fullgenename(current_analysis[[DATASET_NAME]], G)
p=VlnPlot(object = current_analysis$ROOIJonly.sp.bt_RID2l_clExtended, features = G_full, group.by='overview_MAFK_NFE2L1_highlow_fct')+
    give_better_textsize_plot(8*MYPLOTSIZEFACTOR)
ggsave(plot=p, filename = paste0(base_dir,'Rplots/2024/MAFKNFE2L1/',ANALYSIS_NAME,'_custom_MAFK_NFE2L1_Vln_MYH7.pdf'), 
                height=50*MYPLOTSIZEFACTOR, width=100*MYPLOTSIZEFACTOR, units='mm', device=cairo_pdf)



#library(rgl)
#plot3d( expr_MAFK, expr_NFE2L1, expr_NPPA, type = "s", radius = .2 )



