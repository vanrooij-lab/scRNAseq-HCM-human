


# TRIAGE analysis

RTS_table = 
    read.table('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Resources/TRIAGE_human_rts.txt', head=0)
colnames(RTS_table) = c('gene','RTS','TF')
rownames(RTS_table) = RTS_table$gene
# View(RTS_table)

# Load previous analysis
DATASET_NAME='ROOIJonly_RID2l'
if (!exists('current_analysis')) {current_analysis = list()}
current_analysis[[DATASET_NAME]] =
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))

# Copy original analysis as basis
current_analysis$ROOIJonly_TRIAGE = current_analysis$ROOIJonly_RID2l

# Now re-run parts of the analysis, but applying the TRIAGE     
current_analysis$ROOIJonly_TRIAGE

# Normalize data (log(exp/cellTotal*median_cellSum+1)
medianCellCount = median(colSums(current_analysis$ROOIJonly_TRIAGE@assays$RNA@counts))
current_analysis$ROOIJonly_TRIAGE <- NormalizeData(current_analysis$ROOIJonly_TRIAGE, normalization.method = 'LogNormalize', scale.factor = medianCellCount) # values are defaults

# Perform triage transformation manually

# First select genes that are present in the TRIAGE RTS list
# current_analysis$ROOIJonly_TRIAGE@assays$RNA@scale.data
gene_names_RooijRID2l = data.frame(gene_short=shorthand_cutname(rownames(current_analysis$ROOIJonly_TRIAGE@assays$RNA@data)),
                                    gene=rownames(current_analysis$ROOIJonly_TRIAGE@assays$RNA@data))
# First selection of genes
gene_names_RooijRID2l$sel1 = gene_names_RooijRID2l$gene_short %in% RTS_table$gene
# Removing duplicate names (only minor fraction where ENS couldn't be linked symobl and TBCE and POLR2J3)
g1 = gene_names_RooijRID2l$gene_short[gene_names_RooijRID2l$sel1] # gene selection 1
gene_names_RooijRID2l$sel2 = F
gene_names_RooijRID2l$sel2[gene_names_RooijRID2l$sel1] = !(g1 %in% g1[duplicated(g1)])
# Final selection, both long and short names
g2      = gene_names_RooijRID2l$gene_short[gene_names_RooijRID2l$sel2]
g2_long = gene_names_RooijRID2l$gene[gene_names_RooijRID2l$sel2]
exp_matrix_triage =
    current_analysis$ROOIJonly_TRIAGE@assays$RNA@data[g2_long,]

# Now actually perform the transormation
RTS = RTS_table[g2,]$RTS # bit redundant, but sanity check
exp_matrix_triage = apply(exp_matrix_triage, 2, function(X) {X*RTS})

# Just some checking
if (F) {
    View(exp_matrix_triage)
    exp_matrix_triage_RS=data.frame(total=rowSums(exp_matrix_triage))
    View(exp_matrix_triage_RS)
}

# Now put back in Seurat again, and run rest of analysis
current_analysis$ROOIJonly_TRIAGE@assays$RNA@scale.data = exp_matrix_triage

# Also export original matrix, because I'd like to -- as sanity check -- run it through the
# original TRIAGE script 
# Note: this is log10(1+exp) transformed already, so will tell TRIAGE not to log-transform
matrix_for_export = current_analysis$ROOIJonly_TRIAGE@assays$RNA@data[g2_long,]
rownames(matrix_for_export) = g2
write.table(file = paste0(base_dir, 'TRIAGE/exp_matrix_counts.tsv'), x = matrix_for_export, quote = F, sep='\t')
    # See below, I'll also import the result again to compare

# 
# Usually you'd do scaling by Seurat, but that is now not necessary since i did it manually
current_analysis$ROOIJonly_TRIAGE =
    mySeuratAnalysis_verybasic_part2only(mySeuratObject = current_analysis$ROOIJonly_TRIAGE, 
                                         skip_scaling = T, features_to_use_choice = 'all', cluster_resolution = 1)

# all_features = rownames(current_analysis$ROOIJonly_TRIAGE@assays$RNA@scale.data)
# current_analysis$ROOIJonly_TRIAGE <- RunPCA(object=current_analysis$ROOIJonly_TRIAGE, npcs = 30, verbose = FALSE, features=all_features)
# current_analysis$ROOIJonly_TRIAGE <- RunUMAP(current_analysis$ROOIJonly_TRIAGE, reduction = "pca", dims = 1:30)
# current_analysis$ROOIJonly_TRIAGE <- FindNeighbors(current_analysis$ROOIJonly_TRIAGE, reduction = "pca", dims = 1:30)
# current_analysis$ROOIJonly_TRIAGE <- FindClusters(current_analysis$ROOIJonly_TRIAGE, resolution = 1)

DimPlot(current_analysis$ROOIJonly_TRIAGE)
DimPlot(current_analysis$ROOIJonly_TRIAGE, group.by='annotation_patient_fct')

all_markers_TRIAGE = 
    diff_express_clusters(mySeuratObject = current_analysis$ROOIJonly_TRIAGE, mc.cores = 1)
top10_markers = diff_express_clusters_save_results(all_markers = all_markers_TRIAGE, run_name = 'ROOIJonly_TRIAGE', base_dir = base_dir, topX = 20, easy_names = T)

####################################################################################################
# Which clusters overlap?

top10_markers_list = lapply(1:ncol(top10_markers), function(i) {top10_markers[,i][!is.na(top10_markers[,i])]})
names(top10_markers_list) = paste0('cl.',names(all_markers_TRIAGE))
regulon_overlap_heatmap(pooled_regulons = top10_markers_list, base_dir = base_dir, run_name = 'all_markers_TRIAGE', makeallheatmaps = F)

####################################################################################################
# Do it the same way as for "main" analysis

DE_cluster=list()

# Run different clustering resolutions
current_analysis$ROOIJonly_TRIAGE_clExt=FindClusters(current_analysis$ROOIJonly_TRIAGE, resolution = seq(.1,1,.1))

# First do for max. resolution
DimPlot(current_analysis$ROOIJonly_TRIAGE_clExt, group.by = 'RNA_snn_res.1')
DE_cluster$ROOIJonly_TRIAGE_clExt =
    diff_express_clusters(mySeuratObject = current_analysis$ROOIJonly_TRIAGE_clExt, mc.cores = MYMCCORES, custom_ident = 'RNA_snn_res.1')

# Then look at the overlap of clusters

# Create pairs to compare
df_compare = tidyr::expand_grid(x=names(DE_cluster$ROOIJonly_TRIAGE_clExt), y=names(DE_cluster$ROOIJonly_TRIAGE_clExt))

# Calculate overlaps
df_compare$overlap = sapply(1:dim(df_compare)[1], function(X) {
    df1=DE_cluster$ROOIJonly_TRIAGE_clExt[[     df_compare$x[[X]]     ]]
    df2=DE_cluster$ROOIJonly_TRIAGE_clExt[[     df_compare$y[[X]]     ]]
    topX_1=rownames(df1[order(df1$avg_log2FC, decreasing = T),][1:min(dim(df1)[1],20),])
    topX_2=rownames(df2[order(df2$avg_log2FC, decreasing = T),][1:min(dim(df2)[1],20),])
    sum(topX_1 %in% topX_2) /
        min(length(topX_1), length(topX_2))
    })
df_compare$overlap_N = sapply(1:dim(df_compare)[1], function(X) {
    df1=DE_cluster$ROOIJonly_TRIAGE_clExt[[     df_compare$x[[X]]     ]]
    df2=DE_cluster$ROOIJonly_TRIAGE_clExt[[     df_compare$y[[X]]     ]]
    topX_1=rownames(df1[order(df1$avg_log2FC, decreasing = T),][1:min(dim(df1)[1],20),])
    topX_2=rownames(df2[order(df2$avg_log2FC, decreasing = T),][1:min(dim(df2)[1],20),])
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
current_analysis$ROOIJonly_TRIAGE_clExt[['clusters_custom']]=factor(cutree_out[current_analysis$ROOIJonly_TRIAGE_clExt$RNA_snn_res.1], 
    levels=1:length(reassignment_lookup))

# Show it
DimPlot(current_analysis$ROOIJonly_TRIAGE_clExt, group.by = 'clusters_custom')


# Coarse grain cluster duplicates out
# Merge 0 and 2, 1 and 4.

################################################################################

# Now run diff. expr.
# First set identities to the custom ones
Idents(current_analysis$ROOIJonly_TRIAGE_clExt) <- current_analysis$ROOIJonly_TRIAGE_clExt$clusters_custom
DimPlot(current_analysis$ROOIJonly_TRIAGE_clExt)
DE_cluster$ROOIJonly_TRIAGE_clExt =
    diff_express_clusters(mySeuratObject = current_analysis$ROOIJonly_TRIAGE_clExt, mc.cores = MYMCCORES)
table_topDE = diff_express_clusters_save_results(all_markers = DE_cluster$ROOIJonly_TRIAGE_clExt, run_name = 'ROOIJonly_TRIAGE_clExt', base_dir = base_dir, topX = 30)

####################################################################################################
# Show patients

p=DimPlot(object = current_analysis$ROOIJonly_TRIAGE_clExt, group.by = 'annotation_patient_fct', 
          label = F, repel = T, label.size = 7, pt.size = .1)+
                theme_void()+
                give_better_textsize_plot(6)+ggtitle(element_blank())+
                theme(legend.position = 'right',axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                        legend.key.height = unit(1, 'mm'))+
                guides(colour = guide_legend(override.aes = list(size=.5)))
        
p
ggsave(plot=p, filename = paste0(base_dir,'Rplots/',TRIAGE_DATASET,'_9_Patients.pdf'), 
                height=40, width=40, units='mm')

################################################################################
# Now compare with previous analysis

current_analysis$ROOIJonly_TRIAGE_clExt$clusters_custom

current_analysis$ROOIJonly_RID2l$Cl_TRIAGE = current_analysis$ROOIJonly_TRIAGE_clExt$clusters_custom[colnames(current_analysis$ROOIJonly_RID2l)]
p=DimPlot(object = current_analysis$ROOIJonly_RID2l, group.by = 'Cl_TRIAGE', label = F, repel = T, label.size = 7, pt.size = .25)+
                theme_void()+
                give_better_textsize_plot(6)+ggtitle(element_blank())+
                theme(legend.position = 'right',axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                        legend.key.height = unit(1, 'mm'))+
                guides(colour = guide_legend(override.aes = list(size=.5)))
        
p
ggsave(plot=p, filename = paste0(base_dir,'Rplots/',TRIAGE_DATASET,'_9_ClustersProjectedOnRID2l.pdf'), 
                height=40, width=40, units='mm')
    


# And the converse

current_analysis$ROOIJonly_TRIAGE_clExt$Cl_RID2l = Idents(current_analysis$ROOIJonly_RID2l)[colnames(current_analysis$ROOIJonly_TRIAGE_clExt)]
DimPlot(current_analysis$ROOIJonly_TRIAGE_clExt, group.by='Cl_RID2l')

################################################################################
# Save this analysis

TRIAGE_DATASET='ROOIJonly_TRIAGE_clExt'

TRIAGE_DATASET='ROOIJonly_TRIAGE_clExt'
save(list='DE_cluster', file=paste0(base_dir, 'Rdata/DE_cluster__ROOIJonly_TRIAGE_clExt.Rdata'))
SaveH5Seurat(object = current_analysis[[TRIAGE_DATASET]], filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',TRIAGE_DATASET,'.h5seurat'), overwrite = T)

################################################################################
# Now make the standard plots

TRIAGE_DATASET='ROOIJonly_TRIAGE_clExt'
TRIAGE_topXMarkers = as.vector(table_topDE[1:2,])[!is.na(as.vector(table_topDE[1:2,]))]
mySeuratCommonPlots(mySeuratObject = current_analysis[[TRIAGE_DATASET]], run_name = TRIAGE_DATASET, mymarkers = TRIAGE_topXMarkers[1], mypointsize = .5)

plot_list = lapply(TRIAGE_topXMarkers, function(g) {
    shorthand_seurat_custom_expr(current_analysis[[TRIAGE_DATASET]], gene_of_interest = g, 
                                 pointsize = .5, textsize = 8, mymargin = 1, zscore = F)})
wrap_plots(plot_list, nrow=2)            
                 
extra_to_plot = c('HAND1','GATA4','TBX3','IER3')

plot_list2 = lapply(TRIAGE_topXMarkers, function(g) {
    shorthand_seurat_custom_expr(current_analysis[[TRIAGE_DATASET]], gene_of_interest = g, 
                                 pointsize = .5, textsize = 8, mymargin = 1, zscore = F)})


list_p=list(); list_p2=list()
NR_GENES_PER_CLUSTER=2
for (cluster_idx in 1:dim(table_topDE)[2]) {

    # cluster_idx=1
    # cluster_idx=2
    
    current_genes = table_topDE[,cluster_idx][1:NR_GENES_PER_CLUSTER]
    
    # On triage dataset
    # ====
    
    list_p[[cluster_idx]]=wrap_plots(lapply(current_genes, 
                        function(x) {shorthand_seurat_custom_expr(current_analysis[[TRIAGE_DATASET]], x, textsize=6, pointsize=.3)}
                        ), 
                 nrow=NR_GENES_PER_CLUSTER)*theme(legend.position='none')#+
        #plot_annotation(title = paste0('Cluster ', cluster_idx-1))
    print(list_p[[cluster_idx]])
    
    ggsave(plot=list_p[[cluster_idx]], filename = paste0(base_dir,'Rplots/',TRIAGE_DATASET,'_7_clExt_ClusterMarkers_cl_',cluster_idx,'.pdf'), 
                height=20, width=16, units='mm')
    #VlnPlot(current_analysis$ROOIJonly_RID2l, features = shorthand_seurat_fullgenename(current_genes))
    
    # On original dataset (plots little bit bigger)
    list_p2[[cluster_idx]]=wrap_plots(lapply(current_genes, 
                        function(x) {shorthand_seurat_custom_expr(current_analysis[[DATASET_NAME]], x, textsize=8, pointsize=.6)}
                        ), 
                 nrow=NR_GENES_PER_CLUSTER)*theme(legend.position='none')#+
        #plot_annotation(title = paste0('Cluster ', cluster_idx-1))
    print(list_p2[[cluster_idx]])
    
    ggsave(plot=list_p2[[cluster_idx]], filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_PROJECTIONoriginal_clExt_ClusterMarkers_cl_',cluster_idx,'.pdf'), 
                height=40, width=32, units='mm')

}
# Save one additional plot with a legend to use for figures
p=list_p[[cluster_idx]][[1]]+theme(legend.position = 'right', )
ggsave(plot=p, filename = paste0(base_dir,'Rplots/_0_LEGENDBAR_rainbow.pdf'), height=40, width=16, units='mm')

########################################################

# 
# SETNAME='TTN-related'
# # Note: hand picked; I omitted top genes that were already printed
# TTN_related_genes = c('DST', 'NEBL', 'MYOM1', 'COX7A1')
# p=wrap_plots(lapply(TTN_related_genes, 
#                         function(x) {shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], x, textsize=6, pointsize=.3)}
#                         ), 
#                  nrow=2)*theme(legend.position='none')
# # p
# ggsave(plot=p, filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_GenesOfInterest_',SETNAME,'.pdf'), 
#                 height=184/3-4, width=184/3-4, units='mm') 



########################################################
# Sanity check, is my transformation same as theirs?

triage_python_matrix = read.table(paste0(base_dir, '/TRIAGE/exp_matrix_triage_out.tsv'))
rownames(triage_python_matrix) = shorthand_seurat_fullgenename_faster(seuratObject = current_analysis[[TRIAGE_DATASET]], 
                                                               gene_names = rownames(triage_python_matrix))
current_genes = rownames(triage_python_matrix)
current_cells = colnames(triage_python_matrix)

all((triage_python_matrix[current_genes,current_cells] - exp_matrix_triage[current_genes,current_cells])<.01)
all((triage_python_matrix[current_genes,current_cells] - exp_matrix_triage[current_genes,current_cells])<.000001)

shorthand_seurat_fullgenename(seuratObject = current_analysis[[TRIAGE_DATASET]], 
                                                               gene_names = 'TBX5')
shorthand_seurat_fullgenename(seuratObject = current_analysis[[TRIAGE_DATASET]], 
                                                               gene_names = 'TTN')

GENE='ENSG00000089225:TBX5'
GENE="ENSG00000155657:TTN"
triage_python_matrix[current_genes,current_cells][GENE,1:20]
exp_matrix_triage[current_genes,current_cells][GENE,1:20]




######################################################################

# export mini
matrix_for_export_mini = current_analysis$ROOIJonly_TRIAGE@assays$RNA@data[g2_long,][1:30,1:30]
rownames(matrix_for_export_mini) = g2[1:30]
write.table(file = paste0(base_dir, 'TRIAGE/exp_matrix_counts_mini.tsv'), x = matrix_for_export_mini, quote = F, sep='\t')

View(as.matrix(matrix_for_export_mini))











