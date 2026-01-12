
########################################################################

# Extended analysis; 
# - looking at analysis with default Seurat settings


########################################################################

# For convenience, source on HPC first the main script;
# script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')
#
# or locally:
# LOCAL=1; script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Projects/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')


####################################################################################################

DATASET_NAME_NORM= 'ROOIJonly_default'
DATASET_NAME_Rcl = 'ROOIJonly_RID2l_clExtended'
DATASET_NAME_Int1 = 'ROOIJonly_Int1c'

# load the file with the default analysis
if (!exists('current_analysis')) {current_analysis = list()}
current_analysis[[DATASET_NAME_NORM]] =
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME_NORM,'.h5seurat'))

# Load the integrated data set
current_analysis[[DATASET_NAME_NORM]] =
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME_NORM,'.h5seurat'))


# Also load the RID2l analysis for reference
current_analysis[[DATASET_NAME_Int1]] =
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME_Int1,'.h5seurat'))


####################################################################################################
# Compare the UMAP with the old data

# Sanity check
max(rowMeans(current_analysis$ROOIJonly_default@assays$RNA@scale.data))
min(rowMeans(current_analysis$ROOIJonly_default@assays$RNA@scale.data))
apply(current_analysis$ROOIJonly_default@assays$RNA@scale.data,1,sd)
    # seems to be normalized

# Project RID2l clusters on top
current_analysis[[paste0(DATASET_NAME_NORM, '_clExtended')]]$clusters_RID2l = 
    Idents(current_analysis[[DATASET_NAME_Rcl]])[colnames(current_analysis[[paste0(DATASET_NAME_NORM, '_clExtended')]])]
DimPlot(current_analysis[[paste0(DATASET_NAME_NORM, '_clExtended')]], group.by = 'clusters_RID2l')

# Look at some of our favorite genes
GENES_OF_INTEREST = c('XIRP2','NPPA', 'NEAT1','KRT6A')
currentfullgenenames=shorthand_seurat_fullgenename_faster(current_analysis[[paste0(DATASET_NAME_NORM, '_clExtended')]], GENES_OF_INTEREST)
# FeaturePlot(current_analysis[[paste0(DATASET_NAME_NORM, '_clExtended')]], features = currentfullgenenames, cols = rainbow_colors)&theme_void()&theme(legend.position = 'none')
# custom fn
plist=lapply(GENES_OF_INTEREST, 
        function(gene) {shorthand_seurat_custom_expr(seuratObject = current_analysis[[paste0(DATASET_NAME_NORM, '_clExtended')]], 
                                                     gene_of_interest = gene, textsize = 10, pointsize = .1)})
p=wrap_plots(plist, nrow=2)
p
ggsave(plot=p, filename = paste0(base_dir,'Rplots/','ROOIJonly_default_clExtended_','customplots_genesOfInterest.pdf'), width=PANEL_WIDTH-4, height=PANEL_HEIGHT, units='mm', device = cairo_pdf)


# Compare clustering with alluvial
# Create alluvial plot
comparing_cluster_df = data.frame(clustering_default = Idents(current_analysis$ROOIJonly_default_clExtended),
        clustering_RID2l = Idents(current_analysis$ROOIJonly_RID2l_clExtended), freq=1)

alluvial_comparing_cluster_df = aggregate(comparing_cluster_df$freq, 
    by=list(clustering_default=comparing_cluster_df$clustering_default,
            clustering_RID2l=comparing_cluster_df$clustering_RID2l),
    FUN=sum)

library(ggalluvial)
p=ggplot(alluvial_comparing_cluster_df,
       aes(y = x, axis1 = clustering_RID2l, axis2 = clustering_default)) +
  geom_alluvium(aes(fill = as.factor(clustering_RID2l)), width = 1/12)+
  geom_stratum(width = 1/12, fill = "grey", color = "black") +
  geom_label(stat = "stratum", infer.label = TRUE)+give_better_textsize_plot(10)+theme_void()+
    ggtitle('No FN vs. FN clusters')+theme(legend.position = 'none')
ggsave(plot=p, filename = paste0(base_dir,'Rplots/customplots_alluvial_ClDef-vs-ClRidl2.pdf'), width=PANEL_WIDTH-4, height=PANEL_HEIGHT*2, units='mm', device = cairo_pdf)

####################################################################################################
# Now adjust the clustering resolution of the feature-normalized dataset

# first perform multiple resolutions
current_analysis[[paste0(DATASET_NAME_NORM, '_clExtended')]]=FindClusters(current_analysis[[DATASET_NAME_NORM]], resolution = seq(.1,1,.1))

# Inspect results
p=    DimPlot(current_analysis[[paste0(DATASET_NAME_NORM, '_clExtended')]], group.by = 'RNA_snn_res.0.1')+
      DimPlot(current_analysis[[paste0(DATASET_NAME_NORM, '_clExtended')]], group.by = 'RNA_snn_res.0.3')+
      DimPlot(current_analysis[[paste0(DATASET_NAME_NORM, '_clExtended')]], group.by = 'RNA_snn_res.0.4')+
      DimPlot(current_analysis[[paste0(DATASET_NAME_NORM, '_clExtended')]], group.by = 'RNA_snn_res.1') & theme_void()
p

# Update the clusters with the chosen resolution
new_cls_def=as.numeric(current_analysis[[paste0(DATASET_NAME_NORM, '_clExtended')]]$RNA_snn_res.0.3)
Idents(current_analysis[[paste0(DATASET_NAME_NORM, '_clExtended')]]) <- factor(new_cls_def, levels=min(new_cls_def):max(new_cls_def))
DimPlot(current_analysis[[paste0(DATASET_NAME_NORM, '_clExtended')]])

# Now run enrichment analysis
DE_cluster=list()
DE_cluster$ROOIJonly_default_clExtended =
    diff_express_clusters(mySeuratObject = current_analysis$ROOIJonly_default_clExtended, mc.cores = MYMCCORES)

# Export results
DE_out_ROOIJonly_RID2l_clExtended = diff_express_clusters_save_results(
  all_markers = DE_cluster$ROOIJonly_default_clExtended, run_name = 'ROOIJonly_default_clExtended', base_dir = base_dir, topX = 30, extendedOutput = T, FC_cutoff = 1.1, pval_cutoff = .05)

####################################################################################################
# Compare RID2l with feature-norm data (aka default)

p=  (DimPlot(current_analysis$ROOIJonly_default_clExtended, pt.size = .1)+#, label = T, label.size = 10)+
    (DimPlot(current_analysis$ROOIJonly_default_clExtended, pt.size = .1, group.by = 'clusters_RID2l', cols = colors_distinguishable)+#, label = T, label.size = 10)+
         ggtitle(element_blank()))&
    theme_void()&theme(legend.position = 'none'))+
      plot_annotation(title = 'FS and no FS clustering',
                  theme = theme(plot.title = element_text(size = 10*2)))
p
ggsave(plot=p, filename = paste0(base_dir,'Rplots/customplots_UMAP_ClDef-vs-ClRidl2.pdf'), width=2*PANEL_WIDTH-4, height=PANEL_WIDTH, units='mm', device = cairo_pdf)

####################################################################################################
####################################################################################################
# Now let's take a look at the integrated dataset

DimPlot(current_analysis$ROOIJonly_Int1c)
DimPlot(current_analysis$ROOIJonly_Int1c)

# Project RID2l clusters on top
current_analysis$ROOIJonly_Int1c$clusters_RID2l = 
    Idents(current_analysis[[DATASET_NAME_Rcl]])[colnames(current_analysis$ROOIJonly_Int1c)]
p=DimPlot(current_analysis$ROOIJonly_Int1c, group.by = 'clusters_RID2l', label = T, label.size = 10/.pt, pt.size = .1, repel = T, label.box = T)+ggtitle('Integrated; original clusters')+give_better_textsize_plot(8)+theme_void()+theme(legend.position = 'none', plot.title = element_text(size=10))
p
ggsave(plot=p, filename = paste0(base_dir,'Rplots/','ROOIJonly_Int1c','_customplots_UMAP_ClRID2l.pdf'), width=PANEL_WIDTH-4, height=PANEL_WIDTH, units='mm', device = cairo_pdf)

# Look at some of our favorite genes
GENES_OF_INTEREST = c('XIRP2','NPPA', 'NEAT1','KRT6A')
currentfullgenenames=shorthand_seurat_fullgenename_faster(current_analysis[[paste0(DATASET_NAME_NORM, '_clExtended')]], GENES_OF_INTEREST)
# FeaturePlot(current_analysis$ROOIJonly_Int1c, features = currentfullgenenames, cols = rainbow_colors)&theme_void()&theme(legend.position = 'none')
# custom fn
plist=lapply(GENES_OF_INTEREST, 
        function(gene) {shorthand_seurat_custom_expr(seuratObject = current_analysis$ROOIJonly_Int1c, 
                                                     gene_of_interest = gene, textsize = 10, pointsize = .1)})
p=wrap_plots(plist[!is.na(plist)], nrow=2)
p
#ggsave(plot=p, filename = paste0(base_dir,'Rplots/','ROOIJonly_default_clExtended_','customplots_genesOfInterest.pdf'), width=PANEL_WIDTH-4, height=PANEL_HEIGHT, units='mm', device = cairo_pdf)


####################################################################################################
# So clusters don't seem to hold after integrated data
# However, do they hold in single patients?

all_patients = paste0('R.P',1:5,'RID2l')

for (DATASET_PT in all_patients) {

    # DATASET_PT = 'R.P1RID2l'
    
    current_analysis[[DATASET_PT]] =
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_PT,'.h5seurat'))
}

for (DATASET_PT in all_patients) {
    # Project RID2l clusters on top
    current_analysis[[DATASET_PT]]$clusters_RID2l = 
        Idents(current_analysis[[DATASET_NAME_Rcl]])[colnames(current_analysis[[DATASET_PT]])]
}

# Now create plot
FONTSIZE=6
plist = lapply(1:5, function(PT_NR) {
    DATASET_PT = paste0('R.P',PT_NR,'RID2l')
    p=DimPlot(current_analysis[[DATASET_PT]], group.by = 'clusters_RID2l', label = T, label.size = FONTSIZE/.pt, pt.size = .1, repel = T, label.box = T)+
        ggtitle(paste0('P',PT_NR,'; original clusters'))+give_better_textsize_plot(FONTSIZE)+theme_void()+theme(legend.position = 'none', plot.title = element_text(size=FONTSIZE))
    if (PT_NR==3) {p=p+xlim(c(-4,4))} # ugly boiler plate to avoid 7 outlier cells to screw up visibility
    return(p)})
    # DimPlot(current_analysis[[paste0('R.P',3,'RID2l')]])+xlim(c(-4,4))

p=wrap_plots(plist[!is.na(plist)], nrow=2)
p
ggsave(plot=p, filename = paste0(base_dir,'Rplots/','ROOIJonly_default_clExtended_','customplots_genesOfInterest.pdf'), 
       width=(3/2)*PANEL_WIDTH-4, height=(2/2)*PANEL_WIDTH, units='mm', device = cairo_pdf)




