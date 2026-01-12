################################################################################

# List of convenient data files
# These are available locally at:
# srv-lnx-varo1:/opt/backup_wehrens/data/Wehrens2022/Rdata-important
# (to do: make them available online)
# (note to self: these data are on my local laptop at: /Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3b/Rdata-important)

# LIST OF DATA FILES
# 
# Data from myectomy samples as Seurat object, both in Rds and h5seurat format.
# H5_RHL_SeuratObject_nM_sel_ROOIJonly.sp.bt_RID2l_clExtended.Rds
# H5_RHL_SeuratObject_nM_sel_ROOIJonly.sp.bt_RID2l_clExtended.h5seurat
# 
# Data from all cells as Seurat object, both in Rds and h5seurat format.
# H5_RHL_SeuratObject_nM_sel_ALL.SP_btypSel_RID2l_clExtended.Rds
# H5_RHL_SeuratObject_nM_sel_ALL.SP_btypSel_RID2l_clExtended.h5seurat
# 
# List of marker genes per cluster for myectomy cells (figure 1 of paper); access object by, e.g. for cluster 2: DE_cluster$ROOIJonly.sp.bt_RID2l$`2`.
# DE_cluster__ROOIJonly.sp.bt_RID2l_clExtended.Rdata
# 
# List of regulons (figure 4 of the paper); object name: SCENIC_reg_top_genes_sorted_full.
# ROOIJonly.sp.bt_RID2l__SCENIC_reg_top_genes_sorted_full.Rdata
# 
# List of gene modules (figure 5 of the paper); object name: core_regulons_sorted.
# ROOIJonly.sp.bt_RID2l_core_regulons_sorted.Rdata

################################################################################
# Example of how to load and use the data

data_dir = '/Users/m.wehrens/Data_Hubrecht/_2019_02_HCM_SCS/2021_HPC_analysis.3b/Rdata-important/'

# Load the data
Seurat_RooijOnly = readRDS(paste0(data_dir, 'H5_RHL_SeuratObject_nM_sel_ROOIJonly.sp.bt_RID2l_clExtended.Rds'))
Seurat_RooijAll  = readRDS(paste0(data_dir, 'H5_RHL_SeuratObject_nM_sel_ALL.SP_btypSel_RID2l_clExtended.Rds'))

# Show the data with default clusters determined in the paper
DimPlot(Seurat_RooijOnly)
DimPlot(Seurat_RooijAll)

# Plot expression of a gene:
# Genes have dual names, look up one like:
rownames(Seurat_RooijOnly)[    grepl('NFE2L1',rownames(Seurat_RooijOnly))    ]
FeaturePlot(object = Seurat_RooijOnly, features = 'ENSG00000082641:NFE2L1')

################################################################################

seurat_alldata = readRDS('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2024_datasharing/H5_RHL_SeuratObject_nM_sel_ALL.SP_btypSel_RID2l_clExtended.Rds')



