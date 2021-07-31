
# Note that PCA requires feature normalization
# However, feature normalization also emphasizes differences between patients;
# so variation might come from patient-to-patient variation.

# Standard PCA is done on variable features, but those don't include TTN
# (Which has a very high expression, and is a bit of an outlier in the 
# expression vs. var plot)
#

# Create copy of analysis
customnewname=paste0(DATASET_NAME, '_addCustom')
current_analysis[[customnewname]] = current_analysis[[DATASET_NAME]]

# Identify some highly expressed genes that we'd like to take along
average_gene_expressions = rowSums(current_analysis[[customnewname]]@assays$RNA@data)
View(rowSums(current_analysis[[customnewname]]@assays$RNA@data))
top30_genes = names(sort(average_gene_expressions, decreasing = T)[1:30])

plot1 <- VariableFeaturePlot(current_analysis[[customnewname]])
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# First scale features
features_of_interest = unique(c(top30_genes, rownames(current_analysis[[DATASET_NAME]])))
current_analysis[[customnewname]] = 
    ScaleData(current_analysis[[DATASET_NAME]], features = features_of_interest)
current_analysis[[customnewname]] = 
    RunPCA(current_analysis[[customnewname]], features = features_of_interest)

# Now look at the components
# View(current_analysis[[customnewname]]@reductions$pca)
# View(current_analysis[[customnewname]]@reductions$pca@feature.loadings)

# Some of the Seurat standard plots
ElbowPlot(current_analysis[[customnewname]])
VizDimLoadings(current_analysis[[customnewname]], dims = 1:4, reduction = "pca")+give_better_textsize_plot(12)

DimPlot(current_analysis[[customnewname]], reduction = 'pca')
DimPlot(current_analysis[[customnewname]], reduction = 'pca', group.by = 'annotation_patient_str')
FeaturePlot(current_analysis[[customnewname]], reduction = 'pca', features = 'percent.mt')

print(current_analysis[[customnewname]][["pca"]], dims = 1:5, nfeatures = 10)

#View(current_analysis[[customnewname]][["pca"]])
#View(current_analysis[[customnewname]][["pca"]]@feature.loadings)

# Some of 
current_pc_loadings=current_analysis[[customnewname]][["pca"]]@feature.loadings
rownames(current_pc_loadings) = sapply(str_split(rownames(current_pc_loadings), pattern=':'), function(x) {if (x[[2]]=='') {x[[1]]} else {x[[2]]} })

# Create little overview
df_pc_overview = data.frame(pc=numeric(), loading=numeric(), rank=numeric(), sign=character(), gene=character())
NR_GENES=30
for (PC in 1:5) {
    current_genes_neg =  sort(current_pc_loadings[,PC], decreasing = F)[1:NR_GENES]
    current_genes_pos =  sort(current_pc_loadings[,PC], decreasing = T)[1:NR_GENES]
    df_pc_overview = rbind(df_pc_overview, 
        data.frame(pc=PC, rank=c(1:NR_GENES,1:NR_GENES), loading=c(current_genes_neg, current_genes_pos), 
                    sign=rep(c('-','+'),each=NR_GENES), 
                    gene=c(names(current_genes_neg), names(current_genes_pos))))
    
    #lapply(1:5, function(X) {sort(current_pc_loadings[,X], decreasing = T)[1:30]})
    #lapply(1:5, function(X) {sort(current_pc_loadings[,X])[1:30]})
}
# positive loadings
ggplot(df_pc_overview[df_pc_overview$sign=='+',], aes(x=pc, y=rank, fill=loading))+
    scale_y_reverse()+
    #scale_y_continuous(breaks=1:30, labels=1:30, limits = rev )+
    #scale_y_discrete(breaks= 1:30)+
    geom_tile()+theme_minimal()+
    geom_text(aes(label=gene),color='white')+
    ggtitle('Top PC loadings (positive)')
    #+scale_fill_gradientn(colors=rainbow_colors)+
# negative loadings
ggplot(df_pc_overview[df_pc_overview$sign=='-',], aes(x=pc, y=rank, fill=loading))+
    geom_tile()+theme_minimal()+
    geom_text(aes(label=gene),color='white')+
    scale_y_reverse()+ggtitle('Top PC loadings (negative)')#+scale_fill_gradientn(colors=rainbow_colors)
    
# Try to make a neat little table like for the clusters
table_pc_pos=
            sapply(1:5, function(x) {
                sapply(1:30, function(y) {
                    df_pc_overview[df_pc_overview$pc==x&df_pc_overview$rank==y&df_pc_overview$sign=='+',]$gene
                })
            })
colnames(table_pc_pos) = paste0('PC_',1:5)
table_pc_neg=
            sapply(1:5, function(x) {
                sapply(1:30, function(y) {
                    df_pc_overview[df_pc_overview$pc==x&df_pc_overview$rank==y&df_pc_overview$sign=='-',]$gene
                })
            })
colnames(table_pc_neg) = paste0('PC_',1:5)
openxlsx::write.xlsx(x=list('PC_loadings_positive'=table_pc_pos,'PC_loadings_negative'=table_pc_neg),
                     file=paste0(base_dir,'Rplots/custom_PCA_toploadings.xlsx'))

# Project PCs on umap
current_pc_loadings['TTN',]

View(current_analysis[[customnewname]][["pca"]]@cell.embeddings)


current_analysis[[customnewname]]$PC_1=
    current_analysis[[customnewname]][["pca"]]@cell.embeddings[,'PC_1']
current_analysis[[customnewname]]$PC_2=
    current_analysis[[customnewname]][["pca"]]@cell.embeddings[,'PC_2']
current_analysis[[customnewname]]$PC_3=
    current_analysis[[customnewname]][["pca"]]@cell.embeddings[,'PC_3']
current_analysis[[customnewname]]$PC_4=
    current_analysis[[customnewname]][["pca"]]@cell.embeddings[,'PC_4']
current_analysis[[customnewname]]$PC_5=
    current_analysis[[customnewname]][["pca"]]@cell.embeddings[,'PC_5']

FeaturePlot(current_analysis[[customnewname]], features = 'PC_1', cols=rainbow_colors)
FeaturePlot(current_analysis[[customnewname]], features = 'PC_2', cols=rainbow_colors)
FeaturePlot(current_analysis[[customnewname]], features = 'PC_3', cols=rainbow_colors)
FeaturePlot(current_analysis[[customnewname]], features = 'PC_4', cols=rainbow_colors)
FeaturePlot(current_analysis[[customnewname]], features = 'PC_5', cols=rainbow_colors)

library(patchwork)
(FeaturePlot(current_analysis[[customnewname]], features = 'PC_1', cols=rainbow_colors)+theme(legend.position = 'none'))+
(FeaturePlot(current_analysis[[customnewname]], features = 'PC_2', cols=rainbow_colors)+theme(legend.position = 'none'))+
(FeaturePlot(current_analysis[[customnewname]], features = 'PC_3', cols=rainbow_colors)+theme(legend.position = 'none'))+
(FeaturePlot(current_analysis[[customnewname]], features = 'PC_4', cols=rainbow_colors)+theme(legend.position = 'none'))+
plot_layout(ncol = 2)


