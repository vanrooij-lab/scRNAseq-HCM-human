
DimPlot(current_analysis$ROOIJonly_Int1c)

DE_cluster$ROOIJonly_Int1c =
    diff_express_clusters(mySeuratObject = current_analysis$ROOIJonly_Int1c, mc.cores = MYMCCORES)
diff_express_clusters_save_results(all_markers = DE_cluster$ROOIJonly_Int1c, run_name = 'ROOIJonly_Int1c', base_dir = base_dir, topX = 30)

# TTN
FeaturePlot(current_analysis$ROOIJonly_Int1c, features = shorthand_seurat_fullgenename(gene_names = 'TTN', seuratObject = current_analysis$ROOIJonly_RID2l_clusterPlaying), cols=rainbow_colors)
# NPPA
FeaturePlot(current_analysis$ROOIJonly_Int1c, features = shorthand_seurat_fullgenename(gene_names = 'NPPA', seuratObject = current_analysis$ROOIJonly_Int1c), cols=rainbow_colors)
