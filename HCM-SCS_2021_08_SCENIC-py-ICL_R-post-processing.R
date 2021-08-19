
################################################################################

# Scenic script version 2
# ===
# Using Python/command line on HPC
# This script just prepares loom files with expression matrices of the respective patients
# Note that matrices will be filtered, but not log-transformed 
# I'm following the SCENIC example/walkthrough PBMC10k_SCENIC-protocol-CLI here.

################################################################################
# Load my own main script, to be able to use my Seurat data sets
#
# Note: this also sets "base_dir", I will output to 
# base_dir/SCENIC

# set script dir (note: overwritten by loading SeuratRevisedAnalysis below)
if (exists('LOCAL')) {
    script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Projects/'
} else {
    script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'
}

desired_command='dummy'
source(paste0(script_dir, 'HCM-SCS_2021-06_SeuratRevisedAnalysis_v2_UmiTools.R'))

# Also get some fns from my own regulons script to do post-processing
desired_command_regulon='dummy'; source(paste0(script_dir, if (exists('LOCAL')) { 'Projects/' } else { '' }, 'HCM_SCS_2021_06_SeuratRevised_Regulon_v3.R'))


################################################################################
# Libs

library(pheatmap)

################################################################################

if (!dir.exists(paste0(base_dir, 'Ldata/'))) { dir.create(paste0(base_dir, 'Ldata/')) }

DATASET_NAME='ROOIJonly_RID2l'
ALL_PATIENTS = paste0("R.P", 1:5)

for (CURRENT_PATIENT in ALL_PATIENTS) {
# for (CURRENT_PATIENT in ALL_PATIENTS[2:5]) {
    
    # Load respective dataset
    if (!exists('current_analysis')) {current_analysis = list()}
    current_set = paste0(CURRENT_PATIENT, 'RID2l')
    current_analysis[[current_set]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',current_set,'.h5seurat'))
    
    # Retrieve data from Seurat object 
    # ===
    exprMat_ <- current_analysis[[current_set]]@assays$RNA@counts
    # Simply throw out genes that have duplicate or unknown hgnc names (NA probably not interesting, duplicate symbols are very few, seem uninteresting)
    rownames_exprMat_ = shorthand_cutname(rownames(exprMat_), PART1OR2 = 2)
    rownames(exprMat_) = rownames_exprMat_
    duplicate_genes = rownames_exprMat_[(duplicated(rownames_exprMat_))]
    gene_selection_beforehand = rownames_exprMat_[!(rownames_exprMat_ %in% duplicate_genes)]
    exprMat = exprMat_[gene_selection_beforehand,]
        dim(exprMat)
    rm('rownames_exprMat_') # remove temporary matrix
        
    # Quality control is in principle already done, but many genes might be present because
    # they're expressed in one of the other datasets/patients.
    # So let's here again create a gene filter
    genesStatPositiveCells = rowMeans(exprMat>0)
    exprMat_filtered = exprMat[genesStatPositiveCells>.1,]
        dim(exprMat_filtered) # 5833 x 306
        
    # Creating a loom file
    print(paste0('Writing loom for ', current_set))
    lfile <- loomR::create(filename = paste0(base_dir, 'Ldata/', current_set, '.loom'), 
                           data = exprMat_filtered, overwrite = T) 
                           #feature.attrs = list(gene_info=), 
                           #cell.attrs = list(cell=colnames(current_analysis$R.P1RID2l@assays$RNA@data)))
    
}


################################################################################
# Now collect all the result matrices per patient
# Also construct regulon lists using this info

ALL_PATIENTS = paste0("R.P", 1:5)

loom_file_list = list()
scenic_regulons_collected_all_patients_list=list()
scenic_regulons_collected_all_patients_list_regnames=list()
for (CURRENT_PATIENT in ALL_PATIENTS) {

    # collect in list
    loom_file_list[[CURRENT_PATIENT]] <- loomR::connect(filename=paste0(base_dir, '/SCENIC/',CURRENT_PATIENT,'/output/SCENIC_out_',CURRENT_PATIENT,'RID2l.loom'), 
                        mode = 'r')

    df_regulons_compare = as.data.frame(loom_file_list[[CURRENT_PATIENT]]$row.attrs$Regulons[])
    row.names(df_regulons_compare)=loom_file_list[[CURRENT_PATIENT]]$row.attrs$Gene[]
    
    scenic_regulons_collected_all_patients_list[[CURRENT_PATIENT]] = 
        lapply(1:ncol(df_regulons_compare), function(reg_idx) {rownames(df_regulons_compare)[df_regulons_compare[,reg_idx]==1]})
    names(scenic_regulons_collected_all_patients_list[[CURRENT_PATIENT]]) = colnames(df_regulons_compare) # paste0(CURRENT_PATIENT,'.',colnames(df_regulons_compare))
    scenic_regulons_collected_all_patients_list_regnames[[CURRENT_PATIENT]] = colnames(df_regulons_compare)
}

# scenic_regulons_collected_all_patients_list -- holds a list per patient
# scenic_regulons_collected_all_patients_list_regnames -- holds the regulon names per patient

# now create also flattened list
scenic_regulons_collected_all_patients = unlist(scenic_regulons_collected_all_patients_list, recursive = F)
scenic_regulons_collected_all_patients_regnames = 
    unlist(sapply(1:length(scenic_regulons_collected_all_patients_list), function(i) {names(scenic_regulons_collected_all_patients_list[[i]])}), recursive = F)
names(scenic_regulons_collected_all_patients_regnames) = names(scenic_regulons_collected_all_patients)

################################################################################
# Now some analyis of this data

length(scenic_regulons_collected_all_patients)
    
regulon_overlap_heatmap(pooled_regulons = scenic_regulons_collected_all_patients, base_dir = base_dir, run_name = 'SCENIC_regulons_collected_all_patients', MYTREECUTTINGHEIGHT = 5.14, myfontsize = 2)

regulon_overlap_heatmap(pooled_regulons = scenic_regulons_collected_all_patients_list_regnames, base_dir = base_dir, run_name = 'SCENIC_regulons_overlap_regnames', MYTREECUTTINGHEIGHT = .1, myfontsize = 2)

all_regulons = Reduce(x = scenic_regulons_collected_all_patients_list_regnames, f = unique)

library(UpSetR)
# Manual upset df to keep rownames
# upsetR_df = UpSetR::fromList(scenic_regulons_collected_all_patients_list_regnames)
upset_df = data.frame(lapply(scenic_regulons_collected_all_patients_list_regnames, function(tf_list) {all_regulons %in% tf_list}), row.names=all_regulons)
UpSetR::upset(upsetR_list)

pheatmap(1*upset_df, fontsize = 4)


current_selection_regulons_SCENIC = rownames(upset_df[apply(upset_df, 1, sum)>=4,])
# potential other selections:
rownames(upset_df[apply(upset_df, 1, sum)>=5,])
rownames(upset_df[apply(upset_df, 1, sum)>=3,])

# Find out gene overlap between those ones
scenic_regulons_collected_all_patients_overlapping_ones = 
    scenic_regulons_collected_all_patients[scenic_regulons_collected_all_patients_regnames %in% current_selection_regulons_SCENIC]

# regulon_overlap_heatmap(pooled_regulons = scenic_regulons_collected_all_patients_overlapping_ones, base_dir = base_dir, run_name = 'SCENIC_regulons_collected_all_patients_overlapping', MYTREECUTTINGHEIGHT = .1, myfontsize = 2)

# Better calculate the overlap for each respective regulon between the patietns
overlapping_scores_regulons = 
    sapply(current_selection_regulons_SCENIC, function(tf) {
        # tf = current_selection_regulons_SCENIC[1]
        current_reg_multi_pt = scenic_regulons_collected_all_patients[scenic_regulons_collected_all_patients_regnames == tf]
        all_targets = Reduce(x = current_reg_multi_pt, f = unique)
        upset_df_tf = data.frame(lapply(current_reg_multi_pt, function(target_list) {all_targets %in% target_list}), row.names=all_targets)
        current_overlap = mean(apply(upset_df_tf, 1, sum)>=3)
    })
# rownames(upset_df_tf[apply(upset_df_tf, 1, sum)>=2,])
overlapping_scores_regulons_df = data.frame(overlap=overlapping_scores_regulons, gene=names(overlapping_scores_regulons))
overlapping_scores_regulons_df$gene = factor(overlapping_scores_regulons_df$gene, levels=overlapping_scores_regulons_df$gene[order(overlapping_scores_regulons_df$overlap)])

ggplot(overlapping_scores_regulons_df, aes(x=gene, y=overlap))+
    geom_bar(stat='identity')+coord_flip()+theme_bw()+give_better_textsize_plot(8)+ylim(c(0,1))

# Now a secondary thing, let's calculate potential overlap with each other
# First determine core regulons
SCENIC_regulons_core_genes = 
    sapply(current_selection_regulons_SCENIC[overlapping_scores_regulons>.1], function(tf) { # [overlapping_scores_regulons>.1]
        # tf = current_selection_regulons_SCENIC[1]
        current_reg_multi_pt = scenic_regulons_collected_all_patients[scenic_regulons_collected_all_patients_regnames == tf]
        all_targets = Reduce(x = current_reg_multi_pt, f = unique)
        upset_df_tf = data.frame(lapply(current_reg_multi_pt, function(target_list) {all_targets %in% target_list}), row.names=all_targets)
        current_overlap = rownames(upset_df_tf)[apply(upset_df_tf, 1, sum)>=3]
    })
# Take any with >=10 genes from these
SCENIC_regulons_core_genes_sel = SCENIC_regulons_core_genes[sapply(SCENIC_regulons_core_genes, length)>=10]

regulon_overlap_heatmap(pooled_regulons = SCENIC_regulons_core_genes_sel, base_dir = base_dir, run_name = 'SCENIC_regulons_overlap_core_sel', MYTREECUTTINGHEIGHT = 1.8, myfontsize = 2)

################################################################################
# Now determine the relateness of their regulons to my regulons

# Load my regulons
load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/Rplots/ROOIJonly_RID2l_core_regulons_sorted_shortname.Rdata')
    # core_regulons_sorted_shortname

# Now compare matrix
SCENIC_MW_regulon_comparison = 
    sapply(core_regulons_sorted_shortname, function(core_mw) { # [overlapping_scores_regulons>.1]
        sapply(SCENIC_regulons_core_genes, function(core_scenic) { # [overlapping_scores_regulons>.1]
            overlap = mean(core_scenic %in% core_mw)
        })
    })

pheatmap(SCENIC_MW_regulon_comparison)

################################################################################
# And additionally we can project their regulons on our UMAP

ANALYSIS_NAME = "ROOIJonly_RID2l"

if (F) {
    
    # Load analysis if necessary
    if (!exists('current_analysis')) {current_analysis = list()}
    if (!(ANALYSIS_NAME %in% names(current_analysis))) {
        current_analysis[[ANALYSIS_NAME]] =
            LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',ANALYSIS_NAME,'.h5seurat'))
    }
    
    
    # X=shorthand_cutname(core_regulons$s.R.4)
    
    plotlist=lapply(1:length(SCENIC_regulons_core_genes),
        function(X) {shorthand_seurat_custom_expr(seuratObject = current_analysis[[ANALYSIS_NAME]], 
                                                    gene_of_interest = SCENIC_regulons_core_genes[[X]], textsize=6, pointsize=1, 
                                                    custom_title = names(SCENIC_regulons_core_genes)[X], mymargin = .5, zscore = T)})
    plots_per_row=round(184.6/30)
    n_rows=length(plotlist)/plots_per_row
    p=wrap_plots(plotlist, nrow = ceiling(n_rows))
    
    #p
    ggsave(filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_RegulonsSCENIC_UMAP_compositeExpr.pdf'), 
        plot = p, width=30*plots_per_row, height=30*n_rows, units='mm') # 184.6/3*2-4
    
    # 2nd version of plots
    p=wrap_plots(plotlist, nrow = 1)
    ggsave(filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_RegulonsSCENIC_UMAP_compositeExpr.pdf'), 
        plot = p, width=184/3*2-4, height=20*1, units='mm') # 184.6/3*2-4
    
    
}

################################################################################

# Now look at results
lfile <- loomR::connect(filename='/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/SCENIC/R.P1/output/SCENIC_out_R.P1RID2l.loom', 
                        mode = 'r')

lfile[["matrix"]][1:5, 1:5]

full.matrix <- lfile$matrix[, ]
dim(x = full.matrix) 
    # 306 x 5833
    # note: these are the cells x genes again


# So, e.g. 
lfile$row.attrs$Gene[1:3]
    # This shows genes 1-3 of the original table
dim(lfile$row.attrs$Regulons[1:3])
    # This shows for each of the 127 regulons, whether the gene is a member or not (I suppose)

# So let's create a heatmap of this
df_regulons_compare = as.data.frame(lfile$row.attrs$Regulons[])
row.names(df_regulons_compare)=lfile$row.attrs$Gene[]

# Let's also gather the regulon members
scenic_regulons_collected = lapply(1:ncol(df_regulons_compare), function(reg_idx) {rownames(df_regulons_compare)[df_regulons_compare[,reg_idx]==1]})
names(scenic_regulons_collected) = colnames(df_regulons_compare)

# p=pheatmap(as.matrix(df_regulons_compare))


regulon_overlap_heatmap(pooled_regulons = scenic_regulons_collected, base_dir = base_dir, run_name = 'SCENIC_test', MYTREECUTTINGHEIGHT = 1, myfontsize = 2)

    
    

# Conversely, we also have regulon "activity" over the cells
# (±composite expression)
# lfile$col.attrs$



# And 
lfile[["row_attrs/Regulons"]][1:3]

rownames(full.matrix) = lfile$row.attrs$Regulons[[]]
lfile[["row_attrs/Regulons"]][1]
lfile[["row_attrs/Regulons"]][2]
lfile[["row_attrs/Regulons"]][2]

lfile[["row_attrs/Genes"]]



lfile$col.attrs

df_1 = 
    data.frame(auc=lfile$col.attrs$RegulonsAUC, 
            cellid=lfile$col.attrs$CellID)

regulons = lfile$row.attrs$Regulons

lfile$row.attrs$Gene


#####

attrs <- c("nUMI", "nGene", "orig.ident")
attr.df <- lfile$get.attribute.df(MARGIN = 2, attribute.names = attrs)
head(x = attr.df)



