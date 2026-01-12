VlnPlot(current_analysis$ROOIJonly_RID2l_clExtended, features = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended,'MYH7'))

VlnPlot(current_analysis$ROOIJonly_RID2l_clExtended, features = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended,'TTN'))

VlnPlot(current_analysis$ROOIJonly_RID2l_clExtended, features = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended,'MALAT1'))

VlnPlot(current_analysis$ROOIJonly_RID2l_clExtended, features = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended,'RYR2'))

VlnPlot(current_analysis$ROOIJonly_RID2l_clExtended, features = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended,c('MYL2','FHL2')))


FeaturePlot(current_analysis$ROOIJonly_RID2l_clExtended, features = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended,c('RPL32')), cols = rainbow_colors)+
    theme_void()
VlnPlot(current_analysis$ROOIJonly_RID2l_clExtended, features = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended,c('RPL32')), cols = rainbow_colors)+
    theme_void()

FeaturePlot(current_analysis$ROOIJonly_RID2l_clExtended, features = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended,c('GAPDH')), cols = rainbow_colors)+
    theme_void()
VlnPlot(current_analysis$ROOIJonly_RID2l_clExtended, features = shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l_clExtended,c('GAPDH')), cols = rainbow_colors)+
    theme_void()
