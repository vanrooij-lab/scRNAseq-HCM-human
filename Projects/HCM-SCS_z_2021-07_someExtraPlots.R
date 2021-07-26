


CURRENT_RUNNAME="ROOIJonly_RID2l"

# Compare enrichment
current_markers_a <- FindMarkers(current_analysis[[CURRENT_RUNNAME]], ident.1 = 0, min.pct = .05, assay = 'RNAlog10')
current_markers_b <- FindMarkers(current_analysis[[CURRENT_RUNNAME]], ident.1 = 0, min.pct = .05, assay = 'RNA')

# Some features of interest
FeaturePlot(current_analysis[[CURRENT_RUNNAME]], features = 'percent.mt', cols = rainbow_colors)
FeaturePlot(current_analysis[[CURRENT_RUNNAME]], features = 'nCount_RNA', cols = rainbow_colors) # full count
FeaturePlot(current_analysis[[CURRENT_RUNNAME]], features = 'nCount_nMT_RNA', cols = rainbow_colors) # non-mito count

# Regulons
# core 1
FeaturePlot(current_analysis[[CURRENT_RUNNAME]], features = c('CFLAR', 'MYH7B', 'HOOK2', 'DDX17'), cols = rainbow_colors)
# core 2
FeaturePlot(current_analysis[[CURRENT_RUNNAME]], features = c('RPS12' ,'RPL13A' ,'PTMA' ,'RPL14'), cols = rainbow_colors)
# core 3
FeaturePlot(current_analysis[[CURRENT_RUNNAME]], features = c('VIM', 'ACTB', 'B2M', 'HLA-C'), cols = rainbow_colors)
# core 4
FeaturePlot(current_analysis[[CURRENT_RUNNAME]], features = c('CRYAB', 'MYL2', 'TNNI3', 'TPM1'), cols = rainbow_colors)
# core 5
FeaturePlot(current_analysis[[CURRENT_RUNNAME]], features = c('MAP4', 'MYH7', 'ZNF106', 'XIRP2'), cols = rainbow_colors)




