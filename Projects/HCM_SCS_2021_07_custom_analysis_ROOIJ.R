

# Further analysis of clustering.

# So with the first iteration of the Seurat clustering tool,
# we found multiple clusters that have similar DEs.
#
# We can lower the resolution, but then one outlier clsuter
# gets merged into bigger ones, which introduces artifacts.
# 
# So we use clustree to identify hierarchy, which allows us
# to decrease "redundancy" in the clustering, but keep the
# outlier cluster separate.

####################################################################################################

# Execute the Seurat lib loading part of main script first.

# install.packages('clustree')

library(clustree)
library(pheatmap)
#library(scales)

####################################################################################################

# For future reference: using SC3
#library(SC3)
# SCE <- as.SingleCellExperiment(DietSeurat(Seurat_Obj))

####################################################################################################

DATASET_NAME='ROOIJonly_RID2l'

# load the file
if (!exists('current_analysis')) {current_analysis = list()}
current_analysis[[DATASET_NAME]] =
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))

####################################################################################################
# Create different number of clusters

# clustree is a simple visualization tool, which requires different levels of clustering to be present
current_analysis$ROOIJonly_RID2l_clExtended=FindClusters(current_analysis$ROOIJonly_RID2l, resolution = seq(.1,1,.1))

# Illustration of different clusterings
p=    DimPlot(current_analysis$ROOIJonly_RID2l_clExtended, group.by = 'RNA_snn_res.0.1')+
      DimPlot(current_analysis$ROOIJonly_RID2l_clExtended, group.by = 'RNA_snn_res.0.4')+
      DimPlot(current_analysis$ROOIJonly_RID2l_clExtended, group.by = 'RNA_snn_res.1')
ggsave(plot=p, filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_Clustering_multiLevels.pdf'), 
                height=50, width=200, units='mm')

# Visualize hierarchy using clustree
clustree_out = clustree(current_analysis$ROOIJonly_RID2l_clExtended)
ggsave(plot=clustree_out, filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_Clustree_analysis.pdf'), 
                height=250, width=200, units='mm')

####################################################################################################

# A simpler approach is to look at overlapping DE genes

# First do for max. resolution
DimPlot(current_analysis$ROOIJonly_RID2l_clExtended, group.by = 'RNA_snn_res.1')
DE_cluster$ROOIJonly_RID2l_clExtended =
    diff_express_clusters(mySeuratObject = current_analysis$ROOIJonly_RID2l_clExtended, mc.cores = MYMCCORES, custom_ident = 'RNA_snn_res.1')

# Then look at the overlap of clusters

# Create pairs to compare
df_compare = tidyr::expand_grid(x=names(DE_cluster$ROOIJonly_RID2l_clExtended), y=names(DE_cluster$ROOIJonly_RID2l_clExtended))

# Calculate overlaps
df_compare$overlap = sapply(1:dim(df_compare)[1], function(X) {
    df1=DE_cluster$ROOIJonly_RID2l_clExtended[[     df_compare$x[[X]]     ]]
    df2=DE_cluster$ROOIJonly_RID2l_clExtended[[     df_compare$y[[X]]     ]]
    topX_1=rownames(df1[order(df1$avg_log2FC, decreasing = T),][1:20,])
    topX_2=rownames(df2[order(df2$avg_log2FC, decreasing = T),][1:20,])
    sum(topX_1 %in% topX_2) /
        min(length(topX_1), length(topX_1))
    })
df_compare$overlap_N = sapply(1:dim(df_compare)[1], function(X) {
    df1=DE_cluster$ROOIJonly_RID2l_clExtended[[     df_compare$x[[X]]     ]]
    df2=DE_cluster$ROOIJonly_RID2l_clExtended[[     df_compare$y[[X]]     ]]
    topX_1=rownames(df1[order(df1$avg_log2FC, decreasing = T),][1:20,])
    topX_2=rownames(df2[order(df2$avg_log2FC, decreasing = T),][1:20,])
    sum(topX_1 %in% topX_2)
    })

# Create matrix
matrix_compare <- reshape2::acast(df_compare, x~y, value.var="overlap")
hclust_out=hclust(d = as.dist(1-matrix_compare), method = 'ward.D2')
plot(hclust_out, hang=-1)

# Pretty matrix
ggplot(df_compare, aes(x=factor(x, levels=hclust_out$labels[hclust_out$order]), y=factor(y, levels=hclust_out$labels[hclust_out$order]), fill=overlap))+
    geom_tile()+
    geom_text(aes(label=overlap_N), color='white')+theme_minimal()+
    xlab('Cluster')+ylab('Cluster')+
    scale_fill_gradientn(colors=rainbow_colors)+give_better_textsize_plot(8)+theme(legend.position = 'none')

# Based on this analysis, we can merge certain clusters
cutree_out=cutree(hclust_out, k = 5)
df_annotation=data.frame(row.names = names(cutree_out), group=factor(cutree_out))
pheatmap(matrix_compare, cluster_cols = hclust_out, cluster_rows = hclust_out, annotation_col = df_annotation)

# re-assign clusters
reassignment_lookup = cutree_out
current_analysis$ROOIJonly_RID2l_clExtended[['clusters_custom']]=factor(cutree_out[current_analysis$ROOIJonly_RID2l_clExtended$RNA_snn_res.1], 
    levels=1:length(reassignment_lookup))

# Show it
DimPlot(current_analysis$ROOIJonly_RID2l_clExtended, group.by = 'clusters_custom')

####################################################################################################

# Now run diff. expr.
# First set identities to the custom ones
Idents(current_analysis$ROOIJonly_RID2l_clExtended) <- current_analysis$ROOIJonly_RID2l_clExtended$clusters_custom
DimPlot(current_analysis$ROOIJonly_RID2l_clExtended)
DE_cluster$ROOIJonly_RID2l_clExtended =
    diff_express_clusters(mySeuratObject = current_analysis$ROOIJonly_RID2l_clExtended, mc.cores = MYMCCORES)
table_topDE = diff_express_clusters_save_results(all_markers = DE_cluster$ROOIJonly_RID2l_clExtended, run_name = 'ROOIJonly_RID2l_clExtended', base_dir = base_dir, topX = 30)


####################################################################################################

################################################################################
# Run cluster plotting script again
ANALYSIS_NAME_clExtended='ROOIJonly_RID2l_clExtended'
NR_GENES_PER_CLUSTER=4
# Retrieve 4 top genes
gene_symbol_table_clExt=diff_express_clusters_save_results(DE_cluster$ROOIJonly_RID2l_clExtended, run_name=ANALYSIS_NAME_cl0.4, base_dir, 
    topX=NR_GENES_PER_CLUSTER, easy_names=T, save=F) 
list_p=list()
for (cluster_idx in 1:dim(gene_symbol_table_clExt)[2]) {

    # cluster_idx=1
    # cluster_idx=2
    
    current_genes = gene_symbol_table_clExt[,cluster_idx][1:NR_GENES_PER_CLUSTER]
    
    list_p[[cluster_idx]]=wrap_plots(lapply(current_genes, 
                        function(x) {shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME_cl0.4]], x, textsize=6, pointsize=.3)}
                        ), 
                 nrow=NR_GENES_PER_CLUSTER)*theme(legend.position='none')#+
        #plot_annotation(title = paste0('Cluster ', cluster_idx-1))
    print(list_p[[cluster_idx]])
    
    ggsave(plot=list_p[[cluster_idx]], filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME_cl0.4,'_7_clExt_ClusterMarkers_cl_',cluster_idx,'.pdf'), 
                height=40, width=16, units='mm')
    #VlnPlot(current_analysis$ROOIJonly_RID2l, features = shorthand_seurat_fullgenename(current_genes))

}
# Save one additional plot with a legend to use for figures
p=list_p[[cluster_idx]][[1]]+theme(legend.position = 'right', )
ggsave(plot=p, filename = paste0(base_dir,'Rplots/_0_LEGENDBAR_rainbow.pdf'), height=40, width=16, units='mm')

################################################################################
# Print new umap

p=DimPlot(current_analysis$ROOIJonly_RID2l_clExtended, cols = rep(col_vector_60,2), label = T, repel = T, label.size = 7)+
    theme_void()+ggtitle(element_blank())+theme(legend.position = 'none')
# Save files
ggsave(plot = p, filename = paste0(base_dir,'Rplots/',run_name,'_2_umapLabeled_by_',current_annotation,'.pdf'), height=5.5, width=5.5, units='cm')

################################################################################
# New custom plot with patient - cluster distr.

# Custom plot that shows distribution of patients over clusters
# mymaxy=1.5*max(table(Idents(mySeuratObject)))
p=ggplot(data.frame( cluster = Idents(mySeuratObject),
                Donor = mySeuratObject$annotation_patient_fct))+
    geom_bar(aes(x=cluster, fill=Donor))+theme_bw()+
    xlab('Cluster')+ylab('Number of cells')+
    give_better_textsize_plot(8)+
    theme(legend.position = 'right', legend.key.size = unit(3, "mm"),  
            plot.margin = margin(t = 0, r = 0, b = 0, l = 5, unit = "pt"))
    #theme(legend.position=c(.99,.99), legend.justification = c(1,1), legend.key.size = unit(3, "mm"))
    # ylim(c(0,mymaxy)); p
# p
ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_5_Barplot_PatientCluster_ClExt_distr.pdf'), 
    plot = p, height=40, width=42, units='mm')


####################################################################################################
# Also some other gene lists

ANALYSIS_NAME='ROOIJonly_RID2l'

SETNAME='sarcomere'
sarcomere_genes_maya = sort(c('MYBPC3', 'MYH6', 'TTN', 'RYR2', 'CMYA5', 'XIRP2', 'ACTN1', 'TNNT2', 'TPM1', 'MYL2', 'MYL3', 'TNNI3', 'ACTC1')) # CMYA3-->XIRP2
sarcomere_genes_maya = sarcomere_genes_maya[!(sarcomere_genes_maya %in% gene_symbol_table_clExt)]
length(sarcomere_genes_maya)
p=wrap_plots(lapply(sarcomere_genes_maya, 
                        function(x) {shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], x, textsize=6, pointsize=.3)}
                        ), 
                 ncol=3)*theme(legend.position='none')
ggsave(plot=p, filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_GenesOfInterest_',SETNAME,'.pdf'), 
                width=3*16, height=4*10, units='mm')
    # 16*10mm per plot

##

SETNAME='NPPA-related'
# selected_genes
# "DKK3"   "HSPB1"  "MYL7"   "RTN4"   "IGFBP2" "NPPB"   "CD63"   "TPM3"   "ACTA1"  "ACTC1"  "ATP5MG" "XIRP1"  "MYH6"
NPPA_related_genes = c('IGFBP2', 'RTN4', 'CRYAB', 'MYH6', 'MYL7','GAPDH','XIRP1')
p=wrap_plots(lapply(NPPA_related_genes, 
                        function(x) {shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], x, textsize=6, pointsize=.3)}
                        ), 
                 nrow=2)*theme(legend.position='none')
# p
ggsave(plot=p, filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_GenesOfInterest_',SETNAME,'.pdf'), 
                height=50, width=184/3*2-(4*1), units='mm') 

##

SETNAME='TTN-related'
# Note: hand picked; I omitted top genes that were already printed
TTN_related_genes = c('DST', 'NEBL', 'MYOM1', 'COX7A1')
p=wrap_plots(lapply(TTN_related_genes, 
                        function(x) {shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], x, textsize=6, pointsize=.3)}
                        ), 
                 nrow=2)*theme(legend.position='none')
# p
ggsave(plot=p, filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_GenesOfInterest_',SETNAME,'.pdf'), 
                height=184/3-4, width=184/3-4, units='mm') 


##

SETNAME='TFs'
TF_genes = sort(c('MEF2A', 'MEF2C', 'MEF2D', 'NR2F2', 'CAMTA1', 'CAMTA2', 'FLII', 'DR1', 'ARNT', # 'MEF2B'
                                'IRX5', 'DR1', 'SOX4', 'SNAI2', 'RARG', 'NR3C1', 'CENPA',
                                'GATA2'))#'GATA1', 'GFI1B'
p=wrap_plots(lapply(TF_genes, 
                        function(x) {shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], x, textsize=6, pointsize=.3)}
                        ), 
                 nrow=1)*theme(legend.position='none')
ggsave(plot=p, filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_GenesOfInterest_',SETNAME,'.pdf'), 
                height=20, width=184.6, units='mm')

####################################################################################################
# Some more UMAPs of special interest

CURRENT_GENE='NPPA'
p=shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], CURRENT_GENE, textsize=8, pointsize=1)
ggsave(plot=p, filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_GenesOfInterest_',CURRENT_GENE,'.pdf'), 
                height=60, width=60, units='mm')

CURRENT_GENE='TTN'
p=shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], CURRENT_GENE, textsize=8, pointsize=1)
ggsave(plot=p, filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_GenesOfInterest_',CURRENT_GENE,'.pdf'), 
                height=60, width=60, units='mm')

# Old code below
####################################################################################################
####################################################################################################
####################################################################################################

# old code
# Extra marker plots
#genes_of_interest=c('NPPA', 'FHL1', 'XIRP2', 'NEAT1', 'FHL2', 'CMYA5', 'NPPB', 'PCDH7', 'MYOM1', 'MALAT1', 'CKM', 'XIRP2')
#genes_of_interest_fullname=rownames(current_analysis$ROOIJonly_RID2l)[grepl(paste0(paste0(':',genes_of_interest,'$'),collapse = '|'),rownames(current_analysis$ROOIJonly_RID2l))]
#FeaturePlot(current_analysis$ROOIJonly_RID2l, genes_of_interest_fullname, cols=rainbow_colors)



# run_name='ROOIJonly_RID2l'
#current_analysis[[ANALYSIS_NAME]]
#gene_symbol_table=DE_cluster$ROOIJonly_RID2l
ANALYSIS_NAME='ROOIJonly_RID2l'
load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/Rdata/DE_cluster_ROOIJonly_RID2l.Rdata')
NR_GENES_PER_CLUSTER=4
gene_symbol_table=diff_express_clusters_save_results(DE_cluster$ROOIJonly_RID2l, run_name, base_dir, topX=NR_GENES_PER_CLUSTER, easy_names=T, save=F) 
list_p=list()
for (cluster_idx in 1:dim(gene_symbol_table)[2]) {

    # cluster_idx=1
    # cluster_idx=2
    
    current_genes = gene_symbol_table[,cluster_idx][1:NR_GENES_PER_CLUSTER]
    
    list_p[[cluster_idx]]=wrap_plots(lapply(current_genes, 
                        function(x) {shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], x, textsize=6, pointsize=.3)}
                        ), 
                 nrow=NR_GENES_PER_CLUSTER)*theme(legend.position='none')#+
        #plot_annotation(title = paste0('Cluster ', cluster_idx-1))
    print(list_p[[cluster_idx]])
    
    ggsave(plot=list_p[[cluster_idx]], filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_ClusterMarkers_cl_',cluster_idx,'.pdf'), 
                height=40, width=16, units='mm')
    #VlnPlot(current_analysis$ROOIJonly_RID2l, features = shorthand_seurat_fullgenename(current_genes))

}
# Save one additional plot with a legend to use for figures
p=list_p[[cluster_idx]][[1]]+theme(legend.position = 'right', )
ggsave(plot=p, filename = paste0(base_dir,'Rplots/_0_LEGENDBAR_rainbow.pdf'), height=40, width=16, units='mm')




VlnPlot(current_analysis[[ANALYSIS_NAME]], features = shorthand_seurat_fullgenename(current_analysis[[ANALYSIS_NAME]], TF_genes))

shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], 'MEF2A', textsize=12, pointsize=2)
shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], 'MEF2C', textsize=12, pointsize=2)
shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], 'NR2F2', textsize=12, pointsize=2)
shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], 'CAMTA1', textsize=12, pointsize=2)
shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], 'FLII', textsize=12, pointsize=2)
shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], 'DR1', textsize=12, pointsize=2)
shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], 'ARNT', textsize=12, pointsize=2)

shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], 'ARNT', textsize=12, pointsize=2)


################################################################################
# What's the nature of those outlier cells

FeaturePlot(current_analysis$ROOIJonly_RID2l, features = 'percent.mt', cols=rainbow_colors)
FeaturePlot(current_analysis$ROOIJonly_RID2l, features = 'nFeature.nMT', cols=rainbow_colors)
FeaturePlot(current_analysis$ROOIJonly_RID2l, features = 'nCount.nMT', cols=rainbow_colors)

################################################################################

current_analysis$ROOIJonly_RID2l = RunTSNE(current_analysis$ROOIJonly_RID2l)
DimPlot(current_analysis$ROOIJonly_RID2l, reduction = 'tsne')

FeaturePlot(current_analysis$ROOIJonly_RID2l, reduction = 'tsne', features = shorthand_seurat_fullgenename('TTN'), cols=rainbow_colors)

################################################################################
# What is we focus on the main cluster?

DimPlot(current_analysis$ROOIJonly_RID2l)
current_analysis$ROOIJonly_RID2l_sub = 
    subset(current_analysis$ROOIJonly_RID2l, cells=colnames(current_analysis$ROOIJonly_RID2l)[Idents(current_analysis$ROOIJonly_RID2l)!=3])

current_analysis$ROOIJonly_RID2l_sub <- ScaleData(current_analysis$ROOIJonly_RID2l_sub , verbose = FALSE, do.scale = F, do.center = F)
current_analysis$ROOIJonly_RID2l_sub <- RunPCA(current_analysis$ROOIJonly_RID2l_sub , npcs = 30, verbose = FALSE)
current_analysis$ROOIJonly_RID2l_sub <- RunUMAP(current_analysis$ROOIJonly_RID2l_sub, reduction = "pca", dims = 1:20)
current_analysis$ROOIJonly_RID2l_sub <- FindNeighbors(current_analysis$ROOIJonly_RID2l_sub, reduction = "pca", dims = 1:30)
current_analysis$ROOIJonly_RID2l_sub = FindClusters(current_analysis$ROOIJonly_RID2l_sub)
#current_analysis$ROOIJonly_RID2l_sub = FindClusters(current_analysis$ROOIJonly_RID2l_sub, resolution = .2)
current_analysis$ROOIJonly_RID2l_sub = FindClusters(current_analysis$ROOIJonly_RID2l_sub)

DimPlot(current_analysis$ROOIJonly_RID2l_sub)

FeaturePlot(current_analysis$ROOIJonly_RID2l_sub, features = shorthand_seurat_fullgenename('TTN'), cols=rainbow_colors)

df_clustersDE_subclust =
    diff_express_clusters(mySeuratObject = current_analysis$ROOIJonly_RID2l_sub, mc.cores = MYMCCORES)
diff_express_clusters_save_results(all_markers = df_clustersDE_subclust, run_name = 'ROOIJsubclust', base_dir = base_dir, topX = 30)


################################################################################
# Clustering shows somewhat weird patterns, so perhaps we should adjust some parameters
# (Weird, e.g.: cluster 0 shows some outliers across the map)

# library("leiden")
# library(leidenAlg)
# library(reticulate)
# use_python('/Users/m.wehrens/anaconda3/bin/python')

# Leidenalg doesn't seem to be working
current_analysis$ROOIJonly_RID2l_clusterPlaying = 
    FindClusters(current_analysis$ROOIJonly_RID2l, algorithm = 4, resolution = .6)
DimPlot(current_analysis$ROOIJonly_RID2l_clusterPlaying)
df_clustersDE_leiden =
    diff_express_clusters(mySeuratObject = current_analysis$ROOIJonly_RID2l_clusterPlaying, mc.cores = MYMCCORES)
diff_express_clusters_save_results(all_markers = df_clustersDE_leiden, run_name = 'ROOIJalgLeiden', base_dir = base_dir, topX = 30)

# Algorithm 3
current_analysis$ROOIJonly_RID2l_clusterPlaying = 
    FindClusters(current_analysis$ROOIJonly_RID2l, algorithm = 3)
DimPlot(current_analysis$ROOIJonly_RID2l_clusterPlaying)
df_clustersDE_alg3 =
    diff_express_clusters(mySeuratObject = current_analysis$ROOIJonly_RID2l_clusterPlaying, mc.cores = MYMCCORES)
diff_express_clusters_save_results(all_markers = df_clustersDE_alg3, run_name = 'ROOIJalg3', base_dir = base_dir, topX = 30)
# View(df_clustersDE_alg3)

# Algorithm 2
current_analysis$ROOIJonly_RID2l_clusterPlaying = 
    FindClusters(current_analysis$ROOIJonly_RID2l, algorithm = 2)
DimPlot(current_analysis$ROOIJonly_RID2l_clusterPlaying)
df_clustersDE_alg2 =
    diff_express_clusters(mySeuratObject = current_analysis$ROOIJonly_RID2l_clusterPlaying, mc.cores = MYMCCORES)
diff_express_clusters_save_results(all_markers = df_clustersDE_alg2, run_name = 'ROOIJalg2', base_dir = base_dir, topX = 30)

# Algorithm 1, resolution tweaked
current_analysis$ROOIJonly_RID2l_clusterPlaying = 
    FindClusters(current_analysis$ROOIJonly_RID2l, resolution = 1)
DimPlot(current_analysis$ROOIJonly_RID2l_clusterPlaying)
df_clustersDE_alg1_res1 =
    diff_express_clusters(mySeuratObject = current_analysis$ROOIJonly_RID2l_clusterPlaying, mc.cores = MYMCCORES)
diff_express_clusters_save_results(all_markers = df_clustersDE_alg1_res1, run_name = 'ROOIJalg1_res1', base_dir = base_dir, topX = 30)

# Comparison previous
DimPlot(current_analysis$ROOIJonly_RID2l)
    # Conclusion --> resolution 1 is better, since .5 assigns far-away cells to "left" cluster

# TTN expression for reference (in earlier analysis, TTN turned up as important feature)
FeaturePlot(current_analysis$ROOIJonly_RID2l_clusterPlaying, features = shorthand_seurat_fullgenename(gene_names = 'TTN', seuratObject = current_analysis$ROOIJonly_RID2l_clusterPlaying), cols=rainbow_colors)

# Algorithm 1, resolution tweaked
current_analysis$ROOIJonly_RID2l_cluster_0.4 = 
    FindClusters(current_analysis$ROOIJonly_RID2l, resolution = .4)
DimPlot(current_analysis$ROOIJonly_RID2l_cluster_0.4)
DE_cluster$ROOIJonly_RID2l_cluster_0.4 =
    diff_express_clusters(mySeuratObject = current_analysis$ROOIJonly_RID2l_cluster_0.4, mc.cores = MYMCCORES)
diff_express_clusters_save_results(all_markers = DE_cluster$df_clustersDE_alg1_res04, run_name = 'ROOIJalg1_res04', base_dir = base_dir, topX = 30)


################################################################################
# What happens if we focus on the large main cluster?

current_analysis$ROOIJonly_RID2l_sub_0.4 =         
    FindClusters(current_analysis$ROOIJonly_RID2l_sub, resolution = .4)
DimPlot(current_analysis$ROOIJonly_RID2l_sub_0.4)
df_clustersDE_alg1_sub_res04 =
    diff_express_clusters(mySeuratObject = current_analysis$ROOIJonly_RID2l_sub_0.4, mc.cores = MYMCCORES)
diff_express_clusters_save_results(all_markers = df_clustersDE_alg1_sub_res04, run_name = 'ROOIJalg1_sub_res04', base_dir = base_dir, topX = 30)
FeaturePlot(current_analysis$ROOIJonly_RID2l_sub_0.4, features = shorthand_seurat_fullgenename(gene_names = 'TTN', seuratObject = current_analysis$ROOIJonly_RID2l_clusterPlaying), cols=rainbow_colors)



