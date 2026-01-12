

ggplot(data.frame(DAPI=current_analysis$ROOIJonly_RID2l_clExtended_FSCA$DAPI), aes(x=DAPI))+
    geom_freqpoly()

FeaturePlot(current_analysis$ROOIJonly_RID2l_clExtended_FSCA, 'DAPI')
VlnPlot(current_analysis$ROOIJonly_RID2l_clExtended_FSCA, 'DAPI')
