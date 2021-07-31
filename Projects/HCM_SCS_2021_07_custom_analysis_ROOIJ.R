library(scales)

# old code
# Extra marker plots
#genes_of_interest=c('NPPA', 'FHL1', 'XIRP2', 'NEAT1', 'FHL2', 'CMYA5', 'NPPB', 'PCDH7', 'MYOM1', 'MALAT1', 'CKM', 'XIRP2')
#genes_of_interest_fullname=rownames(current_analysis$ROOIJonly_RID2l)[grepl(paste0(paste0(':',genes_of_interest,'$'),collapse = '|'),rownames(current_analysis$ROOIJonly_RID2l))]
#FeaturePlot(current_analysis$ROOIJonly_RID2l, genes_of_interest_fullname, cols=rainbow_colors)

# function to plot expression of a gene
shorthand_seurat_custom_expr = function(gene_of_interest) {
    
    # gene_of_interest = 'KRT6A'
    
    gene_of_interest_fullname =rownames(current_analysis$ROOIJonly_RID2l)[grepl(paste0(':',gene_of_interest,'$'),rownames(current_analysis$ROOIJonly_RID2l))]
    
    # retrieve current expression
    current_expr = current_analysis$ROOIJonly_RID2l@assays$RNA@data[gene_of_interest_fullname,]
    expr_limits=calc_limits(current_expr, percentile = .03)
    
    ggplot(data.frame(UMAP_1=current_analysis$ROOIJonly_RID2l@reductions$umap@cell.embeddings[,1],
                      UMAP_2=current_analysis$ROOIJonly_RID2l@reductions$umap@cell.embeddings[,2],
                      expr=current_expr),
            mapping = aes(x=UMAP_1, y=UMAP_2, color=expr))+
        geom_point()+scale_color_gradientn(colours=rainbow_colors, limits=c(0,expr_limits[2]), oob=squish)+
        give_better_textsize_plot(12)+ggtitle(gene_of_interest)+
        theme_void()#+theme(legend.position = 'none')
       
}
# shorthand to get full names
shorthand_seurat_fullgenename = function(gene_names) {
    gene_of_interest_fullname =rownames(current_analysis$ROOIJonly_RID2l)[grepl(paste0(paste0(':',gene_names,'$'),collapse = '|'),rownames(current_analysis$ROOIJonly_RID2l))]
    return(gene_of_interest_fullname)
}


NR_GENES_PER_CLUSTER=6
for (cluster_idx in 1:dim(gene_symbol_table)[2]) {

    # cluster_idx=1
    # cluster_idx=2
    
    current_genes = gene_symbol_table[,cluster_idx][1:NR_GENES_PER_CLUSTER]
    
    p=wrap_plots(lapply(current_genes, function(x) {shorthand_seurat_custom_expr(x)}), nrow=NR_GENES_PER_CLUSTER)+
        plot_annotation(title = paste0('Cluster ', cluster_idx))
    print(p)
    
    #VlnPlot(current_analysis$ROOIJonly_RID2l, features = shorthand_seurat_fullgenename(current_genes))
       
}

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








