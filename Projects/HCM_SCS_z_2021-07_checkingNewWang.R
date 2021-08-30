
library(ggrepel)
library(patchwork)

################################################################################

# This loads (an older version) of the new wang data, don't use this.

if (F) {
    # Let's just do some checks on the Wang data beforehand
    SCS_df_list_data_raw_WANGonly = loadData_MW_parallel(dataset_list_paths_WANG, mc.cores=MYMCCORES, prefix=F)
    SCS_df_list_data_WANGonly = lapply(SCS_df_list_data_raw_WANGonly, preprocess_convertAAnames_toSymbol, revert_to_hgnc=T, script_dir=script_dir)
    PREL_SeuratObject_list_WANGOnly = mclapply(1:length(SCS_df_list_data_WANGonly), 
                function(idx) {
                    object=CreateSeuratObject(counts = SCS_df_list_data_WANGonly[[idx]], project = names(SCS_df_list_data_WANGonly)[idx])
                    print(paste0(names(SCS_df_list_data_WANGonly)[idx],' done .'))
                    return(object)
                    }, mc.cores = 1)
    names(PREL_SeuratObject_list_WANGOnly) = names(SCS_df_list_data_WANGonly)
    object_size(PREL_SeuratObject_list_WANGOnly) 
    
    # merged object
    PREL_SeuratObject_merged_WANGonly <- merge(PREL_SeuratObject_list_WANGOnly[[1]], y = unlist(PREL_SeuratObject_list_WANGOnly)[2:length(PREL_SeuratObject_list_WANGOnly)], 
                                            add.cell.ids = names(PREL_SeuratObject_list_WANGOnly), project = "H")
    dim(PREL_SeuratObject_merged_WANGonly@assays$RNA@counts)        
    
    # Add some annotation
    currentcellnames=colnames(PREL_SeuratObject_merged_WANGonly)
    PREL_SeuratObject_merged_WANGonly[['annotation_sample_str']] = sapply(str_split(string = currentcellnames, pattern = '_'),function(x){x[[2]]})
    PREL_SeuratObject_merged_WANGonly[['annotation_sample_fct']] = as.factor(PREL_SeuratObject_merged_WANGonly$annotation_sample_str)
    PREL_SeuratObject_merged_WANGonly[['annotation_patient_str']] = sapply(str_split(string = currentcellnames, pattern = '_'),function(x){x[[1]]})
    PREL_SeuratObject_merged_WANGonly[['annotation_patient_fct']] = as.factor(PREL_SeuratObject_merged_WANGonly$annotation_patient_str)
    PREL_SeuratObject_merged_WANGonly[['annotation_paper_str']] = rep('Hu', dim(PREL_SeuratObject_merged_WANGonly)[2])
    PREL_SeuratObject_merged_WANGonly[['annotation_paper_fct']] = factor(PREL_SeuratObject_merged_WANGonly$annotation_paper_str)
    
    SaveH5Seurat(object = PREL_SeuratObject_merged_WANGonly, overwrite = T,
            filename = paste0(base_dir,'Rdata/PREL_SeuratObject_merged_WANGonly.h5seurat'))

}

################################################################################

# This processes the original Wang data into a Seurat object

if (F) {
        
    # OK now we can just compare this with the count file that was provided by Hu lab
    
    dataset_list_original_counttablesWANG = c(GSE109816='/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/Original_counttables/GSE109816_normal_heart_umi_matrix.csv',
                                              GSE121893='/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/Original_counttables/GSE121893_human_heart_sc_umi.csv')
    
    SCS_df_list_data_raw_WANGORIGINAL = loadData_MW_parallel(dataset_list_original_counttablesWANG, 
                                                         mc.cores=MYMCCORES, prefix=F, sep=',')
    names(SCS_df_list_data_raw_WANGORIGINAL) = names(dataset_list_original_counttablesWANG)
    
    # repeat above
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
    
    
    test_table = read.table('/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/Original_counttables/GSE121893_human_heart_sc_umi.csv', header = 1, row.names = 1, sep = ',')
    rm('test_table')
    
    # Load the dataset
    # RHL_SeuratObject_merged_WANGORIGINAL= LoadH5Seurat(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_WANGORIGINAL.h5seurat'))
    #

}

if (F) {
    # To load this file, use 
    current_analysis$RHL_SeuratObject_merged_WANGORIGINAL = LoadH5Seurat(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_WANGORIGINAL.h5seurat'))
}    

##########

# Now compare with latest Seurat Wang/Hu object
# If necessary, load
# current_analysis=list()
# current_analysis$HUonly_RID2l = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))
# current_analysis$RHL_SeuratObject_merged_WANGORIGINAL = LoadH5Seurat(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_WANGORIGINAL.h5seurat'))

##
# First load name conversion
load(file = paste0(base_dir,'Rdata/metadata_Wang_full_table_selection.Rdata'))
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
# current_analysis$WANG_old_new <- FindNeighbors(current_analysis$WANG_old_new, reduction = "pca", dims = 1:30)
# current_analysis$WANG_old_new <- FindClusters(current_analysis$WANG_old_new, resolution = cluster_resolution)

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

p= DimPlot(current_analysis$WANG_old_new_geneIntersect, group.by='version')+give_better_textsize_plot(8)+theme_void()
ggsave(filename = paste0(base_dir, 'Rplots/_0checks_HU_newoldcell_mixing.pdf'), plot=p, width=50, height=50, units='mm')

# PCA analysis
DimPlot(current_analysis$WANG_old_new_geneIntersect, group.by='version', reduction = 'pca')+give_better_textsize_plot(8)+theme_void()
DimPlot(current_analysis$WANG_old_new_geneIntersect, reduction = 'pca')+give_better_textsize_plot(8)+theme_void()

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
p2=ggplot(df_plot, aes(x=pc1, y=pc2, color=cellname_noON))+
    theme_minimal()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())+
    geom_line()+theme(legend.position='none')+xlab('PC1')+ylab('PC2')
p=(p1+p2+plot_layout(nrow = 1))
ggsave(filename = paste0(base_dir, 'Rplots/_0checks_HU_newoldcell_mixing_PCA-lines.pdf'), plot=p, width=100, height=50, units='mm')

# Which genes are in PC1?
View(current_analysis$WANG_old_new_geneIntersect@reductions$pca@feature.loadings)

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
DimPlot(current_analysis$HUonly_RID2l, group.by = 'ident_wang')
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









