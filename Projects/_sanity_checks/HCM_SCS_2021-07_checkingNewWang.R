
################################################################################

library(pheatmap)
library(ggrepel)
library(patchwork)

################################################################################

# This processes the original Wang data into a Seurat object

if (F) {
        
    # OK now we can just compare this with the count file that was provided by Hu lab
    
    dataset_list_original_counttablesWANG = c(GSE109816='/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/Original_counttables/GSE109816_normal_heart_umi_matrix.csv',
                                              GSE121893='/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/Original_counttables/GSE121893_human_heart_sc_umi.csv')
    
    SCS_df_list_data_raw_WANGORIGINAL = loadData_MW_parallel(dataset_list_original_counttablesWANG, 
                                                         mc.cores=MYMCCORES, prefix=F, sep=',')
    names(SCS_df_list_data_raw_WANGORIGINAL) = names(dataset_list_original_counttablesWANG)
    
    # Create Seurat objects
    RHL_SeuratObject_list_WANGORIGINAL = mclapply(1:length(SCS_df_list_data_raw_WANGORIGINAL), 
                function(idx) {
                    object=CreateSeuratObject(counts = SCS_df_list_data_raw_WANGORIGINAL[[idx]], project = names(SCS_df_list_data_raw_WANGORIGINAL)[idx])
                    print(paste0(names(SCS_df_list_data_raw_WANGORIGINAL)[idx],' done .'))
                    return(object)
                    }, mc.cores = 1)
    names(RHL_SeuratObject_list_WANGORIGINAL) = names(SCS_df_list_data_raw_WANGORIGINAL)
    object_size(RHL_SeuratObject_list_WANGORIGINAL) 
    SaveH5Seurat(object = RHL_SeuratObject_list_WANGORIGINAL[[1]], overwrite = T,
            filename = paste0(base_dir,'Rdata/RHL_SeuratObject_list_WANGORIGINAL_1.h5seurat'))
    SaveH5Seurat(object = RHL_SeuratObject_list_WANGORIGINAL[[2]], overwrite = T,
            filename = paste0(base_dir,'Rdata/RHL_SeuratObject_list_WANGORIGINAL_2.h5seurat'))
    
    ###    
    
    # merge objects
    RHL_SeuratObject_merged_WANGORIGINAL <- merge(RHL_SeuratObject_list_WANGORIGINAL[[1]], y = unlist(RHL_SeuratObject_list_WANGORIGINAL)[2:length(RHL_SeuratObject_list_WANGORIGINAL)], 
                                            add.cell.ids = names(RHL_SeuratObject_list_WANGORIGINAL), project = "H")
    dim(RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@counts)    
    
    ###
    
    # Now add some annotation
    currentcellnames=gsub(pattern = 'GSE[0-9]+_',replacement = '',colnames(RHL_SeuratObject_merged_WANGORIGINAL))
    RHL_SeuratObject_merged_WANGORIGINAL[['annotation_sample_str']] = metadata_Wang_full_table[currentcellnames,]$plate_nr
    RHL_SeuratObject_merged_WANGORIGINAL[['annotation_sample_fct']] = as.factor(RHL_SeuratObject_merged_WANGORIGINAL$annotation_sample_str)
    RHL_SeuratObject_merged_WANGORIGINAL[['annotation_patient_str']] = metadata_Wang_full_table[currentcellnames,]$sample
    RHL_SeuratObject_merged_WANGORIGINAL[['annotation_patient_fct']] = as.factor(RHL_SeuratObject_merged_WANGORIGINAL$annotation_patient_str)
    RHL_SeuratObject_merged_WANGORIGINAL[['annotation_paper_str']] = rep('Hu', dim(RHL_SeuratObject_merged_WANGORIGINAL)[2])
    RHL_SeuratObject_merged_WANGORIGINAL[['annotation_paper_fct']] = factor(RHL_SeuratObject_merged_WANGORIGINAL$annotation_paper_str)
    
    # Save final one
    SaveH5Seurat(object = RHL_SeuratObject_merged_WANGORIGINAL, overwrite = T,
            filename = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_WANGORIGINAL.h5seurat'))
    # RHL_SeuratObject_merged_WANGORIGINAL = LoadH5Seurat(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_WANGORIGINAL.h5seurat'))
    
    # Remove working file
    rm('RHL_SeuratObject_merged_WANGORIGINAL')
    
    #test_table = read.table('/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/Original_counttables/GSE121893_human_heart_sc_umi.csv', header = 1, row.names = 1, sep = ',')
    #rm('test_table')
    
    # Load the dataset
    # RHL_SeuratObject_merged_WANGORIGINAL= LoadH5Seurat(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_WANGORIGINAL.h5seurat'))
    #

}



################################################################################

# Now compare with latest Seurat Wang/Hu object
# If necessary, load
# current_analysis=list()
# current_analysis$HUonly_RID2l = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))
# current_analysis$RHL_SeuratObject_merged_WANGORIGINAL = LoadH5Seurat(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_WANGORIGINAL.h5seurat'))

if (F) {
    # To load this files, use 
    if (!exists('current_analysis')) {current_analysis = list()}
    current_analysis$RHL_SeuratObject_merged_WANGORIGINAL = LoadH5Seurat(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_WANGORIGINAL.h5seurat'))
    current_analysis$HUonly_RID2l = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_','HUonly_RID2l','.h5seurat'))
}    

##
# First load name conversion
load(file = paste0(base_dir,'Rdata/metadata_Wang_full_table_selection.Rdata'))
load(file = paste0(base_dir,'Rdata/metadata_Wang_full_table.Rdata')) # also load full table for later making overview
    # loads metadata_Wang_full_table_selection
    # metadata_Wang_full_table_selection$ID_MW; metadata_Wang_full_table_selection$ID
# Make convenient conversion table
conversion_table_old_to_new = metadata_Wang_full_table_selection$ID_MW
names(conversion_table_old_to_new) = metadata_Wang_full_table_selection$ID

names_newWANG = colnames(current_analysis$HUonly_RID2l@assays$RNA@counts)
names_oldWANG_ = colnames(current_analysis$RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@counts)
names_oldWANG = sapply(str_split(names_oldWANG_, pattern = '_'), function(x) {paste0(x[2:length(x)],collapse='_')})
# sanity check
length(names_oldWANG) == length(unique(names_oldWANG)) # TRUE

current_analysis$RHL_SeuratObject_merged_WANGORIGINAL[['ID_new']] = conversion_table_old_to_new[colnames(current_analysis$RHL_SeuratObject_merged_WANGORIGINAL)]

# Now perform two tests: 
# 1. see whether correlation can identify old-new pairs (use seurat to get compare count tables)
# 2. plot together on tsne map, showing connections (use seurat)

# First merge
# Create new object to create consistent gene names
temp_counttable_new=current_analysis$HUonly_RID2l@assays$RNA@counts
rownames(temp_counttable_new) = shorthand_cutname(gene_names = rownames(temp_counttable_new))
genes_new = rownames(temp_counttable_new)
colnames(temp_counttable_new) = sapply(str_split(colnames(temp_counttable_new),pattern = '\\.'), function(x) {paste0(c(x[[1]],x[[2]],x[[5]]),collapse='-')})
cells_new=colnames(temp_counttable_new)
    # "N1.97474.a_R0.C60.AGAACGCCATC" --> "N6-92563-AACCTTATCGG"
temp_Seurat_new = CreateSeuratObject(temp_counttable_new)
    # Warning: Non-unique features (rownames) present in the input matrix, making unique
# For old data, do same to get cellnames correct
temp_counttable_old=current_analysis$RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@counts
temp_counttable_old = temp_counttable_old[,names_oldWANG %in% metadata_Wang_full_table_selection$ID]
colnames(temp_counttable_old) = sapply(str_split(colnames(temp_counttable_old), pattern = '_'), function(x) {paste0(x[2:length(x)],collapse='_')})
colnames(temp_counttable_old) = conversion_table_old_to_new[colnames(temp_counttable_old)]
temp_counttable_old = temp_counttable_old[,colnames(temp_counttable_old) %in% cells_new]
cells_old = colnames(conversion_table_old_to_new)
genes_old=rownames(temp_counttable_old)
temp_Seurat_old = CreateSeuratObject(temp_counttable_old)
# 
current_analysis$WANG_old_new=
    merge(temp_Seurat_old, temp_Seurat_new,
                                            add.cell.ids = c('old','new'), project = "check")
rm('temp_counttable_new','temp_Seurat_new','temp_counttable_old','temp_Seurat_old')

current_analysis$WANG_old_new <- NormalizeData(current_analysis$WANG_old_new) 
current_analysis$WANG_old_new <- FindVariableFeatures(current_analysis$WANG_old_new)
current_analysis$WANG_old_new <- ScaleData(current_analysis$WANG_old_new) # features = all_variable_features is default
current_analysis$WANG_old_new <- RunPCA(object=current_analysis$WANG_old_new, npcs = 30)
current_analysis$WANG_old_new <- RunUMAP(current_analysis$WANG_old_new, reduction = "pca", dims = 1:30)
current_analysis$WANG_old_new <- FindNeighbors(current_analysis$WANG_old_new, reduction = "pca", dims = 1:30)
current_analysis$WANG_old_new <- FindClusters(current_analysis$WANG_old_new)#, resolution = .4)

current_analysis$WANG_old_new[['version']]= sapply(str_split(colnames(current_analysis$WANG_old_new),pattern = '_'), function(x){x[[1]]})

DimPlot(current_analysis$WANG_old_new, group.by='version')

venn_simple_plot_mw(list(old=genes_old,new=genes_new))

gene_intersect = intersect(genes_old, genes_new)
current_analysis$WANG_old_new_geneIntersect=
    subset(current_analysis$WANG_old_new, features=gene_intersect)
current_analysis$WANG_old_new_geneIntersect <- NormalizeData(current_analysis$WANG_old_new_geneIntersect)
current_analysis$WANG_old_new_geneIntersect <- FindVariableFeatures(current_analysis$WANG_old_new_geneIntersect)
current_analysis$WANG_old_new_geneIntersect <- ScaleData(current_analysis$WANG_old_new_geneIntersect, features =gene_intersect) # features = all_variable_features is default
current_analysis$WANG_old_new_geneIntersect <- RunPCA(object=current_analysis$WANG_old_new_geneIntersect, npcs = 30, features =gene_intersect)
current_analysis$WANG_old_new_geneIntersect <- RunUMAP(current_analysis$WANG_old_new_geneIntersect, reduction = "pca", features =gene_intersect)
current_analysis$WANG_old_new_geneIntersect <- FindNeighbors(current_analysis$WANG_old_new_geneIntersect, reduction = "pca", dims = 1:30)
current_analysis$WANG_old_new_geneIntersect <- FindClusters(current_analysis$WANG_old_new_geneIntersect)#, resolution = .4)

################################################################################
# Save the resulting Seurat analyses

SaveH5Seurat(object = current_analysis$WANG_old_new, overwrite = T,
            filename = paste0(base_dir,'Rdata/RHL_WANG_old_new.h5seurat'))
            # current_analysis$WANG_old_new = LoadH5Seurat(file = paste0(base_dir,'Rdata/RHL_WANG_old_new.h5seurat'))

SaveH5Seurat(object = current_analysis$WANG_old_new_geneIntersect, overwrite = T,
            filename = paste0(base_dir,'Rdata/RHL_WANG_old_new_geneIntersect.h5seurat'))
            # current_analysis$WANG_old_new_geneIntersect = LoadH5Seurat(file = paste0(base_dir,'Rdata/RHL_WANG_old_new_geneIntersect.h5seurat'))

################################################################################

p= DimPlot(current_analysis$WANG_old_new_geneIntersect, group.by='version')+give_better_textsize_plot(8)+theme_void()
p
ggsave(filename = paste0(base_dir, 'Rplots/_0checks_HU_newoldcell_mixing.pdf'), plot=p, width=50, height=50, units='mm')

# PCA analysis
p=DimPlot(current_analysis$WANG_old_new_geneIntersect, group.by='version', reduction = 'pca')+give_better_textsize_plot(8)+theme_void()+ggtitle('Re-mapping Wang et al. data')
p
ggsave(filename = paste0(base_dir, 'Rplots/_0checks_HU_newoldcell_mixing_PCA.pdf'), plot=p, width=50, height=50, units='mm')

#DimPlot(current_analysis$WANG_old_new_geneIntersect, reduction = 'pca')+give_better_textsize_plot(8)+theme_void()

current_analysis$WANG_old_new_geneIntersect[['cellname_noON']] = 
    as.vector(sapply(colnames(current_analysis$WANG_old_new_geneIntersect), function(x){gsub(pattern = 'old_|new_', replacement = '', x=x)}))

df_plot=data.frame(
    pc1=current_analysis$WANG_old_new_geneIntersect@reductions$pca@cell.embeddings[,1],
    pc2=current_analysis$WANG_old_new_geneIntersect@reductions$pca@cell.embeddings[,2],
    cellname=as.vector(current_analysis$WANG_old_new_geneIntersect[['cellname_noON']]),
    version=current_analysis$WANG_old_new_geneIntersect$version)
p1=ggplot(df_plot, aes(x=pc1, y=pc2, color=version))+
    theme_minimal()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())+
    geom_point()+theme(legend.position='none')+xlab('PC1')+ylab('PC2')
p1
p2=ggplot(df_plot, aes(x=pc1, y=pc2, color=cellname_noON))+
    theme_minimal()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())+
    geom_line()+theme(legend.position='none')+xlab('PC1')+ylab('PC2')
p2
p=(p1+p2+plot_layout(nrow = 1))
p
ggsave(filename = paste0(base_dir, 'Rplots/_0checks_HU_newoldcell_mixing_PCA-lines.pdf'), plot=p, width=100, height=50, units='mm')

# Which genes are in PC1?
View(current_analysis$WANG_old_new_geneIntersect@reductions$pca@feature.loadings)


# Note that this plot has an extreme legend, so don't print that ..
p=ggplot(df_plot)+   
  theme_minimal()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())+
  geom_line(aes(x=pc1, y=pc2, color=cellname_noON), size=.1, alpha=.5)+
  geom_point(aes(x=pc1, y=pc2, shape=version), color='black', size=.5, alpha=.5)+
  theme(legend.position='none')+
  give_better_textsize_plot(8)
p
ggsave(filename = paste0(base_dir, 'Rplots/_0checks_HU_newoldcell_mixing_PCA-lines-Combined.pdf'), plot=p, width=PANEL_WIDTH-4, height=PANEL_WIDTH-4, units='mm', device = cairo_pdf)    


####
# What about correlations?
WITHSET='WANG_old_new_geneIntersect'
WITHSET='WANG_old_new'
oldcellnames=colnames(current_analysis[[WITHSET]]@assays$RNA@counts[,as.vector(current_analysis[[WITHSET]][['version']]=='old')])
newcellnames=colnames(current_analysis[[WITHSET]]@assays$RNA@counts[,as.vector(current_analysis[[WITHSET]][['version']]=='new')])
joined_datamatrix=as.matrix(current_analysis[[WITHSET]]@assays$RNA@data)
oldname=oldcellnames[1]

# Try for 1 cell
cor_out = t(cor(x = joined_datamatrix[,oldname], y=joined_datamatrix[,current_analysis[[WITHSET]][['version']]=='new']))
oldname
# View(cor_out)

# Now calculate full matrix
cor_matrix_out = cor(x=joined_datamatrix[,sort(oldcellnames)],y=joined_datamatrix[,sort(newcellnames)])
dim(cor_matrix_out)

# Show part of cor matrix
pheatmap(cor_matrix_out[1:10,1:10], cluster_rows = F, cluster_cols = F)

# Now look for each cell, which other cell (in other set) matches best
matches_x_to_y = apply(cor_matrix_out, 1, function(x) {order(x,decreasing = T)[1]})
# Reverse (redundant because t(M)=M)
matches_y_to_x = apply(cor_matrix_out, 2, function(x) {order(x,decreasing = T)[1]})
# Strongest test, do they all pair (old-new) correctly?
all(matches_x_to_y == 1:length(matches_x_to_y))
all(matches_y_to_x == 1:length(matches_x_to_y))

# Partial heatmap (full = too big to display)
p=pheatmap(cor_matrix_out[1:100,1:100], cluster_rows = F, cluster_cols = F, fontsize = 2)
ggsave(paste0(base_dir,'Rplots/_0checks_HU_newcelloldcellpairing_corrmat100.pdf'), plot = p, width=200, height=200, units='mm')
# Illustrating how well it matches for all
p=ggplot(data.frame(cell_nr=1:dim(cor_matrix_out)[1], matched_nr = matches_x_to_y))+
    geom_point(aes(x=cell_nr,y=matched_nr))+theme_bw()
# p
ggsave(paste0(base_dir,'Rplots/_0checks_HU_newcelloldcellpairing_cellNrMatch.pdf'), plot = p, width=50, height=50, units='mm')

################################################################################
# Overview tables of samples based metadata_Wang_full_table & metadata_Wang_full_table_sel

# Let's also create a table that's a bit more clear
metadata_Wang_full_table$group
metadata_Wang_full_table$sample
metadata_Wang_full_table$Type
metadata_Wang_full_table$ident
metadata_Wang_full_table$plate_nr

metadata_Wang_full_table_selFields = metadata_Wang_full_table[,c('group','sample','Type','plate_nr')] # ,'ident'
    # Leaving out ident, since it creates much more complex table

metadata_Wang_full_table_selFields_summaryTable = data.frame(table(metadata_Wang_full_table_selFields))
View(metadata_Wang_full_table_selFields_summaryTable)

sum(metadata_Wang_full_table_selFields_summaryTable$Freq>0)

###

metadata_Wang_full_table_selection_selFields = metadata_Wang_full_table_selection[,c('group','sample','Type','plate_nr', 'ident')]
metadata_Wang_full_table_selection_selFields_summaryTable = data.frame(table(metadata_Wang_full_table_selection_selFields))

View(metadata_Wang_full_table_selection_selFields_summaryTable[metadata_Wang_full_table_selection_selFields_summaryTable$Freq>0,])
dim(metadata_Wang_full_table_selection_selFields_summaryTable[metadata_Wang_full_table_selection_selFields_summaryTable$Freq>0,])

metadata_Wang_full_table_selection_selFields_summaryTable_final = metadata_Wang_full_table_selection_selFields_summaryTable[metadata_Wang_full_table_selection_selFields_summaryTable$Freq>0,]
openxlsx::write.xlsx(x = as.data.frame(metadata_Wang_full_table_selection_selFields_summaryTable_final), 
                     file = paste0(base_dir,'Rplots/mapped_Wang_cells_overview_stats.xlsx'), overwrite=T)

##MARKERNOTEREMOVE-XXXX

################################################################################

# A better cor matrix might perhaps be on normalized expression, as it is not
# biased by highly expressed genes (that might also be more susceptible to
# mitochondrial mapping artifacts etc)

# Previously, I used scale.data, but this is biased towards genes that are different between
# the two sets, since it selects for max var genes; as such this decreases the correlations
# 
# joined_datamatrix_featnorm = as.matrix(current_analysis[[WITHSET]]@assays$RNA@scale.data)

# Now use custom normalized count matrix
joined_datamatrix_featnorm2 = t(scale(t(as.matrix(current_analysis[[WITHSET]]@assays$RNA@data[rowSums(current_analysis[[WITHSET]]@assays$RNA@counts>0)>10,]))))
cor_matrix_featnorm_out = cor(x=joined_datamatrix_featnorm2[,sort(oldcellnames)],y=joined_datamatrix_featnorm2[,sort(newcellnames)])
dim(cor_matrix_featnorm_out)

pheatmap(cor_matrix_featnorm_out[1:10,1:10], cluster_rows = F, cluster_cols = F)

# Partial heatmap (full = too big to display)
p=pheatmap(cor_matrix_featnorm_out[1:10,1:10], cluster_rows = F, cluster_cols = F, fontsize = 2)
ggsave(paste0(base_dir,'Rplots/_0checks_HU_newcelloldcellpairing_corrmat100_featnorm_10x10.pdf'), plot = p, width=200, height=200, units='mm')
p=pheatmap(cor_matrix_featnorm_out[1:100,1:100], cluster_rows = F, cluster_cols = F, fontsize = 2)
ggsave(paste0(base_dir,'Rplots/_0checks_HU_newcelloldcellpairing_corrmat100_featnorm.pdf'), plot = p, width=200, height=200, units='mm')
# Illustrating how well it matches for all
matches_x_to_y_featnorm = apply(cor_matrix_featnorm_out, 1, function(x) {order(x,decreasing = T)[1]})
p=ggplot(data.frame(cell_nr=1:dim(cor_matrix_featnorm_out)[1], matched_nr = matches_x_to_y_featnorm))+
    geom_point(aes(x=cell_nr,y=matched_nr))+theme_bw()
# p
ggsave(paste0(base_dir,'Rplots/_0checks_HU_newcelloldcellpairing_cellNrMatch.pdf'), plot = p, width=50, height=50, units='mm')

# Sanity check showing read counts
hist(log10(.1+colSums(current_analysis[[WITHSET]]@assays$RNA@counts)))

# Distribution of 1st and 2nd best R-values
Rmatch1_featnorm = apply(cor_matrix_featnorm_out, 1, function(x) {x[order(x,decreasing = T)][1]})
Rmatch2_featnorm = apply(cor_matrix_featnorm_out, 1, function(x) {x[order(x,decreasing = T)][2]})
p=ggplot(data.frame(Rmatch=c(Rmatch1_featnorm, Rmatch2_featnorm), matchtype=rep(c('1st_match','2nd_best'), each=length(Rmatch1_featnorm))))+
    geom_freqpoly(aes(x=Rmatch, color=matchtype), binwidth=.05)+theme_bw()+give_better_textsize_plot(8)+xlab('Correlation coefficient')
p+give_better_textsize_plot(12)
ggsave(paste0(base_dir,'Rplots/_0checks_HU_newcelloldcellpairing_Rvals_matches.pdf'), plot = p, width=50, height=50, units='mm', device = cairo_pdf)

################################################################################

# Another simple check would be whether cells cluster in similar ways

class_old_or_new = rep(NA, length(Cluster_assignments))
class_old_or_new[grepl('^old_',names(Cluster_assignments))] = 'old'
class_old_or_new[grepl('^new_',names(Cluster_assignments))] = 'new'

current_analysis$WANG_old_new_geneIntersect$source_old_or_new = class_old_or_new
p=DimPlot(current_analysis$WANG_old_new_geneIntersect, group.by = 'source_old_or_new')+theme_void()+ggtitle(element_blank())
p
ggsave(paste0(base_dir,'Rplots/_0checks_HU_newcelloldcellpairing_clusterAssignments_UMAP_source.pdf'), plot = p, width=75, height=75, units='mm', device = cairo_pdf)


p=DimPlot(current_analysis$WANG_old_new_geneIntersect)+theme_void()
ggsave(paste0(base_dir,'Rplots/_0checks_HU_newcelloldcellpairing_clusterAssignments_UMAP.pdf'), plot = p, width=75, height=75, units='mm', device = cairo_pdf)

Cluster_assignments = Idents(current_analysis$WANG_old_new_geneIntersect)

Cluster_assignments_old = Cluster_assignments[grepl('^old_',names(Cluster_assignments))]
names(Cluster_assignments_old) = gsub('^old_','',names(Cluster_assignments_old))
Cluster_assignments_new = Cluster_assignments[grepl('^new_',names(Cluster_assignments))]
names(Cluster_assignments_new) = gsub('^new_','',names(Cluster_assignments_new))

compare_clusters_assignments_old_new = 
    data.frame(cell_name=names(Cluster_assignments_old), cluster_old=Cluster_assignments_old[names(Cluster_assignments_old)], cluster_new=Cluster_assignments_new[names(Cluster_assignments_old)])

# compare_clusters_assignments_old_new_melted = 
#     data.frame(cell_name=c(names(Cluster_assignments_old),names(Cluster_assignments_new)), 
#                cluster=c(Cluster_assignments_old[names(Cluster_assignments_old)], Cluster_assignments_new[names(Cluster_assignments_old)]),
#                old_or_new=rep(c('old','new'), each=length(Cluster_assignments_old)))

compare_clusters_assignments_old_new_table = table(compare_clusters_assignments_old_new[c('cluster_old','cluster_new')])

library(ggalluvial) # install.packages('ggalluvial')

p=ggplot(as.data.frame(compare_clusters_assignments_old_new_table),
       aes(y = Freq, axis1 = cluster_old, axis2 = cluster_new)) +
  #geom_alluvium(aes(fill = cluster_old), width = 1/12)
  geom_alluvium(aes(fill = cluster_old)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  #scale_x_discrete(limits = c("Survey", "Response"),
  #                 expand = c(0.15, 0.05)) +
  theme_void()+theme(legend.position = 'none', 
                     axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())+
    give_better_textsize_plot(7)
p
ggsave(paste0(base_dir,'Rplots/_0checks_HU_newcelloldcellpairing_clusterAssignments.pdf'), plot = p, width=50, height=100, units='mm', device = cairo_pdf)

# Perform cluster enrichment
DE_Wang_New_old = diff_express_clusters(current_analysis$WANG_old_new_geneIntersect, mc.cores = 1)
# Export results
DE_Wang_New_old_Tables = diff_express_clusters_save_results(
  all_markers = DE_Wang_New_old, run_name = 'DE_Wang_New_old', base_dir = base_dir, topX = 60, extendedOutput = T, FC_cutoff = 1.1, pval_cutoff = .05, easy_names = F)
    
DE_Wang_New_old_OldVsNew = diff_express_clusters(current_analysis$WANG_old_new_geneIntersect, mc.cores = 1, custom_ident = factor(grepl('^new_',names(Cluster_assignments))*1, levels=c(0,1)))
DE_Wang_New_old_OldVsNew_Tables = diff_express_clusters_save_results(
  all_markers = DE_Wang_New_old_OldVsNew, run_name = 'DE_Wang_New_old_OldVsNew', base_dir = base_dir, topX = 60, extendedOutput = T, FC_cutoff = 1.1, pval_cutoff = .05, easy_names = F)

################################################################################

# Project original Wang cluster assignments on the data
cellnames_WANG_old_new_geneIntersect_ = gsub(pattern = '^old_|^new_', replacement = '', colnames(current_analysis$WANG_old_new_geneIntersect))
metadata_Wang_full_table_selection_ = metadata_Wang_full_table_selection
rownames(metadata_Wang_full_table_selection_) = metadata_Wang_full_table_selection_$ID_MW
current_analysis$WANG_old_new_geneIntersect[['ident_wang']] = metadata_Wang_full_table_selection_[cellnames_WANG_old_new_geneIntersect_,]$ident
p=DimPlot(current_analysis$WANG_old_new_geneIntersect, group.by = 'ident_wang')+theme_void()+theme()+give_better_textsize_plot(10)+theme_void_extramw_removeTickText()+ggtitle(element_blank())
p
ggsave(paste0(base_dir,'Rplots/_0checks_HU_newcelloldcellpairing_UMAPoriginalWANGIdents.pdf'), plot = p, width=75, height=75, units='mm', device = cairo_pdf)

# Now also project Wang sample annotation on the data
current_analysis$WANG_old_new_geneIntersect[['sample_annot_wang']] = metadata_Wang_full_table_selection_[cellnames_WANG_old_new_geneIntersect_,]$Type
p=DimPlot(current_analysis$WANG_old_new_geneIntersect, group.by = 'sample_annot_wang')+theme_void()+theme()+give_better_textsize_plot(10)+theme_void_extramw_removeTickText()+ggtitle(element_blank())
p
ggsave(paste0(base_dir,'Rplots/_0checks_HU_newcelloldcellpairing_UMAPoriginalWANG_SampleType.pdf'), plot = p, width=75, height=75, units='mm', device = cairo_pdf)


################################################################################

# Just some more sanity checking

# just sanity checking naming consistency
rownames(metadata_Wang_full_table_selection)[1]
colnames(current_analysis$HUonly_RID2l)[1]

# Retrieving some extra annotation from metadata table
cellnames_HUonly_RID2l_ = gsub(pattern = '\\.a|\\.b', replacement = '', colnames(current_analysis$HUonly_RID2l))
current_analysis$HUonly_RID2l[['ident_wang']] = metadata_Wang_full_table_selection[cellnames_HUonly_RID2l_,]$ident
current_analysis$HUonly_RID2l[['enrichment_group']] = metadata_Wang_full_table_selection[cellnames_HUonly_RID2l_,]$group
current_analysis$HUonly_RID2l[['Individual']] = metadata_Wang_full_table_selection[cellnames_HUonly_RID2l_,]$Individual    

# Some plots
DimPlot(current_analysis$HUonly_RID2l)
p=DimPlot(current_analysis$HUonly_RID2l, group.by = 'ident_wang')+theme_void()+theme()+give_better_textsize_plot(10)+theme_void_extramw_removeTickText()+ggtitle(element_blank())
p
ggsave(paste0(base_dir,'Rplots/_0checks_HU_cellIdentitiesWangOnMyUMAP.pdf'), plot = p, width=50, height=100, units='mm', device = cairo_pdf)
DimPlot(current_analysis$HUonly_RID2l, group.by = 'enrichment_group')
DimPlot(current_analysis$HUonly_RID2l, group.by = 'annotation_patient_fct')
DimPlot(current_analysis$HUonly_RID2l, group.by = 'annotation_region_fct')
DimPlot(current_analysis$HUonly_RID2l, group.by = 'Individual')

FeaturePlot(current_analysis$HUonly_RID2l, features = 'nCount.nMT')


# checking patient source consistency
all(
sapply(current_analysis$HUonly_RID2l[['annotation_patient_str']], function(x) {gsub(pattern = 'H\\.', replacement = '', x = x )}) == 
    current_analysis$HUonly_RID2l[['Individual']])

# Checking expression of some (marker) genes
some_markers = shorthand_seurat_fullgenename(seuratObject = current_analysis$HUonly_RID2l, c('TTN','NPPA'))
FeaturePlot(current_analysis$HUonly_RID2l, features = some_markers)







################################################################################
# Older stuff




# Compare some arbitrary cells
names_new = colnames(current_analysis$HUonly_RID2l@assays$RNA@counts)
names_old = colnames(current_analysis$RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@counts)
shared_gene_names = unique(c(names(cellcounts_new),names(cellcounts_old)))



venn_simple_plot_mw(list(genes_new=names(cellcounts_new),genes_old=names(cellcounts_old)))
cells_to_sample = sample(1:length(names_new), 16)

listP=list(); listP_NAlab=list(); counter=0
for (cell_idx in cells_to_sample) {
    
    counter=counter+1

    cellcounts_new = current_analysis$HUonly_RID2l@assays$RNA@counts[,names_new[cell_idx]]
    cellcounts_old = current_analysis$RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@counts[,wang_convert_name_new_to_old(names_new[cell_idx], add_GSE = T)]
    
    # for display
    old_cell_name = wang_convert_name_new_to_old(names_new[cell_idx], add_GSE = T)
    
    compare_cell_df = data.frame(new_UMI=cellcounts_new[shared_gene_names], 
                                 old_UMI=cellcounts_old[shared_gene_names],
                                 gene=shared_gene_names, row.names = shared_gene_names)
    compare_cell_df[is.na(compare_cell_df)]=-1
    
    compare_cell_df_sel = compare_cell_df[compare_cell_df$new_UMI>0|compare_cell_df$old_UMI>0,]
    
    # determine top 10 genes that are 0 in other dataset
    compare_cell_df_sel_sel0new=compare_cell_df_sel[compare_cell_df_sel$new_UMI==0,]
    compare_cell_df_sel_sel0old=compare_cell_df_sel[compare_cell_df_sel$old_UMI==0,]
    genes0high_1 = 
        compare_cell_df_sel_sel0new[order(compare_cell_df_sel_sel0new$old_UMI, decreasing = T),]$gene[1:10]
    genes0high_2 = 
        compare_cell_df_sel_sel0old[order(compare_cell_df_sel_sel0old$new_UMI, decreasing = T),]$gene[1:10]
    # determine top 5 genes that are absent in other dataset
    compare_cell_df_sel_selNAnew=compare_cell_df_sel[compare_cell_df_sel$new_UMI==-1,]
    compare_cell_df_sel_selNAold=compare_cell_df_sel[compare_cell_df_sel$old_UMI==-1,]
    genesNAhigh_1 = 
        compare_cell_df_sel_selNAnew[order(compare_cell_df_sel_selNAnew$old_UMI, decreasing = T),]$gene[1:5]
    genesNAhigh_2 = 
        compare_cell_df_sel_selNAold[order(compare_cell_df_sel_selNAold$new_UMI, decreasing = T),]$gene[1:5]
    
    
    listP[[counter]] = 
        ggplot(compare_cell_df_sel, aes(x=log10(1.1+old_UMI), y=log10(1.1+new_UMI)))+
            geom_point()+theme_bw()+
            geom_abline(slope = 1, intercept = 0)+
            geom_text_repel(data=compare_cell_df_sel[c(genes0high_1, genes0high_2),],
                    aes(x=log10(1.1+old_UMI), y=log10(1.1+new_UMI), label=gene), color='red', size=3, max.overlaps = 100)+
        ggtitle(paste0(old_cell_name,'\n',names_new[cell_idx]))+give_better_textsize_plot(8)
    
    listP_NAlab[[counter]] = 
        ggplot(compare_cell_df_sel, aes(x=log10(1.1+old_UMI), y=log10(1.1+new_UMI)))+
                geom_point()+theme_bw()+
                geom_abline(slope = 1, intercept = 0)+
                geom_text_repel(data=compare_cell_df_sel[c(genesNAhigh_1, genesNAhigh_2),],
                        aes(x=log10(1.1+old_UMI), y=log10(1.1+new_UMI), label=gene), color='purple', size=3, max.overlaps = 100)+#, angle=45)+
            ggtitle(paste0(old_cell_name,'\n',names_new[cell_idx]))+give_better_textsize_plot(8)

    # Example of gene where synonyms seem to be the problem
    # ATP5J (old name) / ATP5PF (new name)

}
p_combined = wrap_plots(listP)
p_combined_NAlab = wrap_plots(listP_NAlab)
print(p_combined)
NR=2
ggsave(filename = paste0(base_dir,'Rplots/','MAPPING_CHECK','_0_comparing_newWangOldWang_Expression16Cells_',NR,'.png'), plot = p_combined, height=35, width=35, units = 'cm')
ggsave(filename = paste0(base_dir,'Rplots/','MAPPING_CHECK','_0_comparing_newWangOldWang_Expression16Cells_NAlab_',NR,'.png'), plot = p_combined_NAlab, height=35, width=35, units = 'cm')

################################################################################

# OLD MAPPING
# Let's now compare the vCM's we wanted to select between 
# their mapping and our mapping

old_vCM_cellnames = sapply(colnames(PREL_SeuratObject_merged_WANGonly), wang_convert_name_new_to_old, add_GSE=T)

RHL_SeuratObject_merged_WANGORIGINAL_vCM = subset(RHL_SeuratObject_merged_WANGORIGINAL, cells = old_vCM_cellnames)

RHL_SeuratObject_merged_WANGORIGINAL_vCM = 
    mySeuratAnalysis(mySeuratObject = RHL_SeuratObject_merged_WANGORIGINAL_vCM, run_name = 'WANGORIGINAL_vCM')
mySeuratCommonPlots(RHL_SeuratObject_merged_WANGORIGINAL_vCM, run_name = 'WANGORIGINAL_vCM')

# NEW MAPPING
# Let's now compare the vCM's we wanted to select between 
# their mapping and our mapping

PREL_SeuratObject_merged_WANGonly = 
    mySeuratAnalysis(mySeuratObject = PREL_SeuratObject_merged_WANGonly, run_name = 'PRELWANGMAPPEDAGAIN_vCM')
mySeuratCommonPlots(PREL_SeuratObject_merged_WANGonly, run_name = 'PRELWANGMAPPEDAGAIN_vCM')


################################################################################

# What do they do about pseudogenes?

mito_names = 
    rownames(RHL_SeuratObject_merged_WANGORIGINAL)[grepl('^MT', rownames(RHL_SeuratObject_merged_WANGORIGINAL))]
        # They also do seem to have pseudogenes in there; let's check their expression

mito_related_totals = 
    rowSums(RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@data[mito_names,])
    # View(mito_related_totals)

# e.g.
compare_mito_df = data.frame(MTND6P3 = RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@data['MTND6P3',], 
                             MT.ND6 = RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@data['MT-ND6',])
#library(zoo)
#y1=rollmean(compare_mito_df$expr_real, 100)
ggplot(compare_mito_df,aes(x=MT.ND6,y=MTND6P3))+
    geom_point(alpha=.25)+theme_bw()+geom_abline(intercept = 0, slope = 1)+
    ggtitle('MT-ND6 real vs. pseudo expr')
    #geom_line(aes(y=rollmean(expr_pseudo, 30, na.pad=T)), color='red')

# e.g.
compare_mito_df2 = data.frame(MT.ATP6 = RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@data['MT-ATP6',], 
                             MTATP6P1 = RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@data['MTATP6P1',])
#library(zoo)
#y1=rollmean(compare_mito_df$expr_real, 100)
ggplot(compare_mito_df2,aes(x=MT.ATP6,y=MTATP6P1))+
    geom_point(alpha=.25)+theme_bw()+geom_abline(intercept = 0, slope = 1)+
    ggtitle('MT-ATP6 real vs. pseudo expr')


################################################################################









