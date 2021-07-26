
library(ggrepel)
library(patchwork)

################################################################################

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


################################################################################

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

# merged object
RHL_SeuratObject_merged_WANGORIGINAL <- merge(RHL_SeuratObject_list_WANGORIGINAL[[1]], y = unlist(RHL_SeuratObject_list_WANGORIGINAL)[2:length(RHL_SeuratObject_list_WANGORIGINAL)], 
                                        add.cell.ids = names(RHL_SeuratObject_list_WANGORIGINAL), project = "H")
dim(RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@counts)    

################################################################################

# Now add some annotation
currentcellnames=gsub(pattern = 'GSE[0-9]+_',replacement = '',colnames(RHL_SeuratObject_merged_WANGORIGINAL))
RHL_SeuratObject_merged_WANGORIGINAL[['annotation_sample_str']] = metadata_Wang_full_table[currentcellnames,]$plate_nr
RHL_SeuratObject_merged_WANGORIGINAL[['annotation_sample_fct']] = as.factor(RHL_SeuratObject_merged_WANGORIGINAL$annotation_sample_str)
RHL_SeuratObject_merged_WANGORIGINAL[['annotation_patient_str']] = metadata_Wang_full_table[currentcellnames,]$sample
RHL_SeuratObject_merged_WANGORIGINAL[['annotation_patient_fct']] = as.factor(RHL_SeuratObject_merged_WANGORIGINAL$annotation_patient_str)
RHL_SeuratObject_merged_WANGORIGINAL[['annotation_paper_str']] = rep('Hu', dim(RHL_SeuratObject_merged_WANGORIGINAL)[2])
RHL_SeuratObject_merged_WANGORIGINAL[['annotation_paper_fct']] = factor(RHL_SeuratObject_merged_WANGORIGINAL$annotation_paper_str)


SaveH5Seurat(object = RHL_SeuratObject_merged_WANGORIGINAL, overwrite = T,
        filename = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_WANGORIGINAL.h5seurat'))

test_table = read.table('/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/Original_counttables/GSE121893_human_heart_sc_umi.csv', header = 1, row.names = 1, sep = ',')
rm('test_table')

# Load the dataset
# RHL_SeuratObject_merged_WANGORIGINAL= LoadH5Seurat(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged_WANGORIGINAL.h5seurat'))
#


##########

# Compare some arbitrary cells
names_new = colnames(PREL_SeuratObject_merged_WANGonly@assays$RNA@counts)
names_old = colnames(RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@counts)
shared_gene_names = unique(c(names(cellcounts_new),names(cellcounts_old)))

venn_simple_plot_mw(list(genes_new=names(cellcounts_new),genes_old=names(cellcounts_old)))
cells_to_sample = sample(1:length(names_new), 16)

listP=list(); listP_NAlab=list(); counter=0
for (cell_idx in cells_to_sample) {
    
    counter=counter+1

    cellcounts_new = PREL_SeuratObject_merged_WANGonly@assays$RNA@counts[,names_new[cell_idx]]
    cellcounts_old = RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@counts[,wang_convert_name_new_to_old(names_new[cell_idx], add_GSE = T)]
    
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








