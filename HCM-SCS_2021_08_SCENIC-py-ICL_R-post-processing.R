
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
# Libs/more config

library(pheatmap)

DATASET_NAME='ROOIJonly_RID2l'

################################################################################
# Export patients for analysis

if (!dir.exists(paste0(base_dir, 'Ldata/'))) { dir.create(paste0(base_dir, 'Ldata/')) }

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

################################################################################
# Upset like heatmap
# ===

upset_df = data.frame(lapply(scenic_regulons_collected_all_patients_list_regnames, function(tf_list) {all_regulons %in% tf_list}), row.names=all_regulons)

# Upset plot not so clear
# library(UpSetR)
# upsetR_df = UpSetR::fromList(scenic_regulons_collected_all_patients_list_regnames)
# UpSetR::upset(upsetR_list)

# Heatmap of upset df
cellsize=200/nrow(upset_df)
p=pheatmap(1*upset_df, fontsize = cellsize*.pt, cluster_cols = F, cluster_rows = T, color = c('white','black'), legend = F, 
           treeheight_row = 0, cellwidth = cellsize*.pt, cellheight = cellsize*.pt, border_color = NA) # .pt = points/mm
p
ggsave(filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_RegulonsSCENIC_overlapHeatmap.pdf'), 
        plot = p, width=15+ncol(upset_df)*cellsize, height=15+nrow(upset_df)*cellsize, units='mm') # 184.6/3*2-4
    
################################################################################
# Now retreive the NES scores

myRegScores=list()
for (CURRENT_PATIENT in ALL_PATIENTS) {
 
    # Gather some info about the reliability of the link between gene sets and their TF
    add.info.meta = read.csv(paste0(base_dir,'/SCENIC/',CURRENT_PATIENT,'/reg.csv'), head=F, nrows=3)
    add.info = read.csv(paste0(base_dir,'/SCENIC/',CURRENT_PATIENT,'/reg.csv'), head=F, skip = 3)
    colnames(add.info)=c(add.info.meta[3, 1:2], add.info.meta[2, 3:length(add.info)])
        # View(add.info)
    
    # Now summarize this to have NES scores available per regulon
    myRegScores[[CURRENT_PATIENT]] = aggregate(list(NES=add.info$NES), by = list(TF=add.info$TF), FUN=median)
    rownames(myRegScores[[CURRENT_PATIENT]]) = myRegScores[[CURRENT_PATIENT]]$TF
       
}

# Slightly more sophisticated "upset_df"
upset_df_NES = data.frame(lapply(names(scenic_regulons_collected_all_patients_list_regnames), 
     function(current_patient)  {
        tf_list = scenic_regulons_collected_all_patients_list_regnames[[current_patient]]
        regulons_of_interest = all_regulons[all_regulons %in% tf_list]
        regulons_of_interest_  = gsub(pattern = '\\([+-]\\)',replacement = '', regulons_of_interest)
        values = myRegScores[[current_patient]][regulons_of_interest_,]$NES
        column_w_score = rep(NA,length(all_regulons))
        column_w_score[all_regulons %in% tf_list] = values
        return(column_w_score)
    }), row.names=all_regulons)
colnames(upset_df_NES) = names(scenic_regulons_collected_all_patients_list_regnames)

# Make heatmap again
cellsize=200/nrow(upset_df_NES)
upset_df_NES_ = upset_df_NES; upset_df_NES_[is.na(upset_df_NES)]=0
upset_df_NES_ = Reduce(rbind, rev(lapply(split(upset_df_NES_, apply(upset_df[rownames(upset_df_NES_),], 1, sum)), function(X) {
        # order method of the separate parts 
        # X[hclust(dist(X))$order,] # clustering
        X[order(apply(X, 1, mean), decreasing = T),]
    })))
p=pheatmap(upset_df_NES_, fontsize = cellsize*.pt, cluster_cols = F, cluster_rows = F, legend = F, 
           treeheight_row = 0, cellwidth = cellsize*.pt, cellheight = cellsize*.pt, border_color = NA) # .pt = points/mm
p
ggsave(filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_RegulonsSCENIC_overlapHeatmap_NES.pdf'), 
        plot = p, width=15+ncol(upset_df)*cellsize, height=15+nrow(upset_df)*cellsize, units='mm') # 184.6/3*2-4


################################################################################
# Delve further in selected regulons, select first by requiring minimally present in 4 patients
NR_PATIENTS_SCEN_FILTER = 4

current_selection_regulons_SCENIC = rownames(upset_df[apply(upset_df, 1, sum)>=NR_PATIENTS_SCEN_FILTER,])
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
        # tf = 'STAT6(+)' 
        # tf = 'MEF2A(+)' # current_selection_regulons_SCENIC[1]
        # tf = 'MLX(+)'
        current_reg_multi_pt = scenic_regulons_collected_all_patients[scenic_regulons_collected_all_patients_regnames == tf]
        all_targets = unique(Reduce(x = current_reg_multi_pt, f = c)) # sort(Reduce(x = current_reg_multi_pt, f = unique))
        upset_df_tf = data.frame(lapply(current_reg_multi_pt, function(target_list) {all_targets %in% target_list}), row.names=all_targets)
        # current_overlap = mean(apply(upset_df_tf, 1, sum)>=3)
        rsize_per_patient=apply(upset_df_tf, 2, sum)
        effective_max_hits = min(sort(rsize_per_patient[rsize_per_patient>0])[1:3]) 
            # rationale: if we were to have a regulon present in 4 patients, and one of those patient-regulons is very small, but the others are big,
            # i'd like to only compare overlap for the three biggest ones
            # for now, something that happens in 3 patients is considered above the magic treshold
        current_overlap = min(1, sum(apply(upset_df_tf, 1, sum)>=3)/effective_max_hits)
        
    })
# rownames(upset_df_tf[apply(upset_df_tf, 1, sum)>=2,])
overlapping_scores_regulons_df = data.frame(overlap=overlapping_scores_regulons, gene=names(overlapping_scores_regulons))
overlapping_scores_regulons_df$gene = factor(overlapping_scores_regulons_df$gene, levels=overlapping_scores_regulons_df$gene[order(overlapping_scores_regulons_df$overlap)])
# overlapping_scores_regulons_df$score = overlapping_scores_regulons_df$gene
# Determine NES scores per TF
overlapping_scores_regulons_df$NES = apply(upset_df_NES, 1, function(x) {median(x[x>0])})[overlapping_scores_regulons_df$gene]

# Now show gene overlap between those regulons that were identified in the separate patients 
p=ggplot(overlapping_scores_regulons_df, aes(x=gene, y=overlap))+
    geom_bar(stat='identity')+coord_flip()+theme_bw()+give_better_textsize_plot(8)+ylim(c(0,1))+ylab('Fraction identical genes\nbetween patients')+xlab('Regulon')
# p
ggsave(filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_RegulonsSCENIC_MemberGeneConsistencyPatients.pdf'), 
        plot = p, width=50, height=8/.pt*nrow(overlapping_scores_regulons_df), units='mm') # 184.6/3*2-4

# Now including NES color
p=ggplot(overlapping_scores_regulons_df, aes(x=gene, y=overlap, fill=NES))+
    geom_bar(stat='identity')+coord_flip()+theme_bw()+give_better_textsize_plot(8)+ylim(c(0,1))+ylab('Fraction identical genes\nbetween patients')+xlab('Regulon')+
    # scale_fill_gradientn(colours = rainbow_colors)+
    scale_fill_gradientn(colours = c('orange','red'))+theme(legend.position='none')
# p
ggsave(filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_RegulonsSCENIC_MemberGeneConsistencyPatients_NES.pdf'), 
        plot = p, width=50, height=8/.pt*nrow(overlapping_scores_regulons_df), units='mm') # 184.6/3*2-4


################################################################################
# Now a secondary thing, let's calculate potential overlap with each other
# First determine core regulons, genes present in three regulons, only consider regulons that show >50% gene overlap, finally only take regulons >10 genes

PERCENTAGE_TRESHOLD=.5
SCENIC_regulons_core_genes = 
    sapply(current_selection_regulons_SCENIC[overlapping_scores_regulons>PERCENTAGE_TRESHOLD], function(tf) { # [overlapping_scores_regulons>.1]
        # tf = current_selection_regulons_SCENIC[1]
        current_reg_multi_pt = scenic_regulons_collected_all_patients[scenic_regulons_collected_all_patients_regnames == tf]
        all_targets = unique(Reduce(x = current_reg_multi_pt, f = c))
        upset_df_tf = data.frame(lapply(current_reg_multi_pt, function(target_list) {all_targets %in% target_list}), row.names=all_targets)
        current_overlap = rownames(upset_df_tf)[apply(upset_df_tf, 1, sum)>=3]
    })
# Take any with >=10 genes from these
SCENIC_regulons_core_genes_sel = SCENIC_regulons_core_genes[sapply(SCENIC_regulons_core_genes, length)>=10]

regulon_overlap_heatmap(pooled_regulons = SCENIC_regulons_core_genes_sel, base_dir = base_dir, run_name = 'SCENIC_regulons_overlap_core_sel', MYTREECUTTINGHEIGHT = 1.8, myfontsize = 8, makeallheatmaps=T)

length(SCENIC_regulons_core_genes_sel)

save(list='SCENIC_regulons_core_genes_sel', file=paste0(base_dir, 'Rdata/SCENIC_regulons_core_genes_sel.Rdata'))

################################################################################
# Create lists of top genes for the selected core regulons

# First load and organize weights of the genes

# Now also get gene weights
RELEVANT_REGULONS = gsub(pattern = '\\([+-]\\)', replacement = '', x = names(SCENIC_regulons_core_genes_sel))
regulon_gene_importance = list()
for (CURRENT_PATIENT in ALL_PATIENTS) {
    print(paste0('Looking up selected genes for ',CURRENT_PATIENT))
    
    adj.info = read.csv(paste0(base_dir, '/SCENIC/',CURRENT_PATIENT,'/adj.csv'), head=1)
    adj.info = adj.info[adj.info$TF %in% RELEVANT_REGULONS,]
    
    # Importance: higher = more important
    regulon_gene_importance[[CURRENT_PATIENT]] = 
        lapply(unique(adj.info$TF), function(tf) {adj.info[adj.info$TF==tf,c('target', 'importance')]})
    names(regulon_gene_importance[[CURRENT_PATIENT]]) = unique(adj.info$TF)
    
}
rm('adj.info')


# Now first create a summary parameter for each of the lists 
SCENIC_reg_top_genes = 
    lapply(names(SCENIC_regulons_core_genes_sel), function(tf) {
    
        # tf = 'ZEB1(+)'
        # tf = 'MEF2A(+)'
        # tf = 'USF2(+)'
        current_member_list = SCENIC_regulons_core_genes_sel[[tf]]
    
        tf_ = gsub('\\([+-]\\)', '', tf)
        
        relevant_weight_values_per_pt = lapply(regulon_gene_importance, function(X) {X[[tf_]]})
        
        importance_df = data.frame(lapply(relevant_weight_values_per_pt, function(X) {rownames(X) = X$target; X[current_member_list,]$importance}), row.names=current_member_list)
        importance_df$median_nona = apply(importance_df, 1, function(x) {median(x[!is.na(x)])})
        top_genes = rownames(importance_df[order(importance_df$median_nona, decreasing = T),])[1:10]
    })
names(SCENIC_reg_top_genes) = gsub('\\([+-]\\)', '', names(SCENIC_regulons_core_genes_sel))
SCENIC_reg_top_genes_df = data.frame(SCENIC_reg_top_genes)

openxlsx::write.xlsx(x= SCENIC_reg_top_genes_df[,sort(colnames(SCENIC_reg_top_genes_df))], file = paste0(base_dir,'Rplots/',DATASET_NAME,'_SCENIC_top10_regulons.xlsx'), overwrite = T)        
        
# NOTE: ALSO GIVE LENGTHS!!
SCENIC_regulon_lengths = 
    data.frame(lapply(SCENIC_regulons_core_genes_sel[sort(names(SCENIC_regulons_core_genes_sel))], length))
colnames(SCENIC_regulon_lengths)=gsub('\\([+-]\\)','',sort(names(SCENIC_regulons_core_genes_sel)))
SCENIC_regulon_lengths[2,] = paste0(colnames(SCENIC_regulon_lengths), ' (', SCENIC_regulon_lengths[1,],')')
openxlsx::write.xlsx(x= SCENIC_regulon_lengths, file = paste0(base_dir,'Rplots/',DATASET_NAME,'_SCENIC_regulon_lengths.xlsx'), overwrite = T)        

regulon_gene_importance$R.P1$tf

SCENIC_regulons_core_genes_sel


# Better calculate the overlap for each respective regulon between the patietns
lists_of_regulon_genes = 
    sapply(current_selection_regulons_SCENIC, function(tf) {
        
        # tf="MEF2A(+)"
        
        current_reg_multi_pt = scenic_regulons_collected_all_patients[scenic_regulons_collected_all_patients_regnames == tf]
        
        
        
    })


################################################################################
# Create new regulon overview where all genes are sorted

# First whole df
SCENIC_reg_top_genes_sorted_full_dfs = 
    lapply(names(SCENIC_regulons_core_genes_sel), function(tf) {
    
        # tf = 'ZEB1(+)'
        # tf = 'MEF2A(+)'
        # tf = 'USF2(+)'
        current_member_list = SCENIC_regulons_core_genes_sel[[tf]]
    
        tf_ = gsub('\\([+-]\\)', '', tf)
        
        relevant_weight_values_per_pt = lapply(regulon_gene_importance, function(X) {X[[tf_]]})
        
        importance_df = data.frame(lapply(relevant_weight_values_per_pt, function(X) {rownames(X) = X$target; X[current_member_list,]$importance}), row.names=current_member_list)
        importance_df$median_nona = apply(importance_df, 1, function(x) {median(x[!is.na(x)])})
        
        # top_genes = rownames(importance_df[order(importance_df$median_nona, decreasing = T),])[1:100]
        importance_df[order(importance_df$median_nona, decreasing = T),]
    })
names(SCENIC_reg_top_genes_sorted_full_dfs) = gsub('\\([+-]\\)', '', names(SCENIC_regulons_core_genes_sel))

# Now just lists of names
SCENIC_reg_top_genes_sorted_full =
    lapply(names(SCENIC_regulons_core_genes_sel), function(tf) {
            tf_ = gsub('\\([+-]\\)', '', tf)
            rownames(SCENIC_reg_top_genes_sorted_full_dfs[[tf_]])
        })
names(SCENIC_reg_top_genes_sorted_full) = gsub('\\([+-]\\)', '', names(SCENIC_regulons_core_genes_sel))


# Save the list of SCENIC regulon genes
save(x='SCENIC_reg_top_genes_sorted_full', file=paste0(base_dir,'Rdata/',DATASET_NAME,'__SCENIC_reg_top_genes_sorted_full.Rdata'))

#####

lapply(SCENIC_reg_top_genes_sorted_full, function(X) {which(X %in% 'MYL2')})
lapply(SCENIC_reg_top_genes_sorted_full, length)


lapply(SCENIC_reg_top_genes_top100, function(X) {which(X %in% 'MYL2')})

# Sanity check:
GENE='SRF'
target_genes= shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l, rownames(SCENIC_reg_top_genes_sorted_full_dfs[[GENE]]))
gene_of_interest=shorthand_seurat_fullgenename_faster(current_analysis$ROOIJonly_RID2l, GENE)
correlations=
    sapply(target_genes[target_genes!=gene_of_interest], function(t) {
        #out=cor.test(current_analysis$ROOIJonly_RID2l@assays$RNA@data[gene_of_interest,], current_analysis$ROOIJonly_RID2l@assays$RNA@data[t,])
        #out$estimate
        cor(current_analysis$ROOIJonly_RID2l@assays$RNA@data[gene_of_interest,], current_analysis$ROOIJonly_RID2l@assays$RNA@data[t,])
        })

df_toplot=data.frame(corr=correlations, importance=SCENIC_reg_top_genes_sorted_full_dfs[[GENE]][names(target_genes)[target_genes!=gene_of_interest],]$median_nona, name=names(correlations))
high_corr=rownames(df_toplot[order(df_toplot$corr, decreasing = T),][1:15,])
high_importance=rownames(df_toplot[1:15,])
ggplot(df_toplot,aes(x=corr, y=importance))+
    geom_point()+
    geom_text_repel(data=df_toplot[unique(c(high_corr,high_importance)),], aes(label=name), color='red')+
    geom_smooth(method='lm')+theme_bw()+ggtitle(GENE)

################################################################################
# Now determine the relateness of their regulons to my regulons

# Load my regulons
load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/Rplots/ROOIJonly_RID2l_core_regulons_sorted_shortname.Rdata')
    # core_regulons_sorted_shortname

# Now compare matrix
SCENIC_MW_regulon_comparison = 
    sapply(core_regulons_sorted_shortname, function(core_mw) { # [overlapping_scores_regulons>.1]
        sapply(SCENIC_regulons_core_genes, function(core_scenic) { # [overlapping_scores_regulons>.1]
            overlap = sum(core_scenic %in% core_mw)/min(length(core_scenic), length(core_mw))
        })
    })

pheatmap(SCENIC_MW_regulon_comparison)

# Heatmap again
myfontsize=8
p=pheatmap(SCENIC_MW_regulon_comparison, cluster_rows = T,cluster_cols = F, 
            fontsize = myfontsize, fontsize_col = myfontsize, fontsize_row = myfontsize, treeheight_row = 0, treeheight_col = 0, legend = F, limits=c(0,1), border_color = NA)
            #annotation_colors = list(colors=col_Dark2[1:max(cutree_out)])))
print(p)
ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_7_regulon_SCENIC_MW_gene_overlap-v1b-L.pdf'), 
    plot = p, width=ncol(SCENIC_MW_regulon_comparison)*myfontsize/.pt*1.1+20, height=nrow(SCENIC_MW_regulon_comparison)*myfontsize/.pt*1.1+15, units='mm', limitsize = F)

################################################################################
# And additionally we can project their regulons on our UMAP

DATASET_NAME = "ROOIJonly_RID2l"

if (F) {
    
    # Load analysis if necessary
    if (!exists('current_analysis')) {current_analysis = list()}
    if (!(DATASET_NAME %in% names(current_analysis))) {
        current_analysis[[DATASET_NAME]] =
            LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))
    }
    
    
    # X=shorthand_cutname(core_regulons$s.R.4)
    
    plotlist=lapply(1:length(SCENIC_regulons_core_genes),
        function(X) {shorthand_seurat_custom_expr(seuratObject = current_analysis[[DATASET_NAME]], 
                                                    gene_of_interest = SCENIC_regulons_core_genes[[X]], textsize=6, pointsize=.5, 
                                                    custom_title = names(SCENIC_regulons_core_genes)[X], mymargin = .5, zscore = T)})
    plots_per_row=round((184.6/2)/20)
    n_rows=length(plotlist)/plots_per_row
    p=wrap_plots(plotlist, nrow = ceiling(n_rows))
    
    #p
    ggsave(filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_RegulonsSCENIC_UMAP_compositeExpr.pdf'), 
        plot = p, width=20*plots_per_row, height=20*n_rows, units='mm') # 184.6/3*2-4
    ggsave(filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_RegulonsSCENIC_UMAP_compositeExpr-v2.pdf'), 
        plot = p, width=100, height=75, units='mm') # 184.6/3*2-4
    
}

################################################################################
# And also the underlying TFs


regulon_controllers = sort(c('CEBPB', 'ESRRA', 'ETV1', 'FOXN3', 'JUND', 'MAFK', 'MEF2A', 'MEF2D', 'MXI1', 'NR3C1', 'SREBF2', 'SRF', 'USF2', 'YY1', 'ZEB1', 'ZNF91'))


p_list=lapply(regulon_controllers, function(gene) {
    shorthand_seurat_custom_expr(seuratObject = current_analysis[[ANALYSIS_NAME]], 
                                                    gene_of_interest = gene, textsize=4, pointsize=.25, 
                                                    custom_title = paste0(gene), mymargin = .5, zscore = T)})
p=wrap_plots(p_list, nrow=1)#ceiling(sqrt(length(p_list))))
ggsave(filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_9_custom_SCENIC_tfs.pdf'), 
        plot = p, width=9.5*16, height=10, units='mm') # 184.6/3*2-4

################################################################################
# Now also sort their regulons' genes, just using the same method as for my
# own regulon analysis
# 
# df_core_regulon_overview_list=list()
# 
# function(core_regulon_list, current_analysis, DATASET_NAME) {
#     
#     for (core_reg_name in names(core_regulon_list)) {
#         
#         if (length(core_regulon_list[[core_reg_name]])<1) { break }
#         
#         # core_reg_name='s.R.4'
#         
#         all_patients = unique(current_analysis[[DATASET_NAME]]$annotation_patient_str)
#         #all_patients = names(collected_regulon_objects$TEICHMANNonly_RID2l)
#         df_reg_genescores = data.frame(gene=character(), corr=numeric())
#         for (current_patient in all_patients) {
#         
#             print(paste0('Calculating ',core_reg_name, ' - ',current_patient))
#             
#             patient_selection = current_analysis[[DATASET_NAME]]$annotation_patient_str==current_patient
#             
#             expr_mat = as.matrix(current_analysis[[DATASET_NAME]]@assays$RNA@data)
#             rownames(expr_mat) = shorthand_cutname(rownames(expr_mat))
#             
#             # current_analysis[[DATASET_NAME]]@assays$RNA@counts[core_regulon_list[[core_reg_name]], ]
#             # Let's do this per patient again
#             cor_out = cor(t(expr_mat[core_regulon_list[[core_reg_name]], patient_selection]))
#             cor_reg_avg = colMeans(cor_out)
#             #View(cor_reg_avg)
#             
#             df_reg_genescores=rbind(df_reg_genescores, data.frame(gene=names(cor_reg_avg), corr=cor_reg_avg))
#             
#         }
#         df_reg_genescores_sum       = aggregate(x = list(median_corr=df_reg_genescores$corr), by = list(gene=df_reg_genescores$gene), median)
#         
#         # !!!! CONTINUE EDITING HERE!!!! 
#         XXXXX
#         df_reg_genescores_sum$nr_patients = XXXX core_regulons_gene_info_table[[core_reg_name]][df_reg_genescores_sum$gene,]$total_occurence XXXX
#         XXXXX
#         # !!!! CONTINUE EDITING HERE!!!!
#         rownames(df_reg_genescores_sum) = df_reg_genescores_sum$gene
#         
#         # Sort first by patients, then by correlation (median)
#         df_reg_genescores_sum$nr_patients <- ordered(df_reg_genescores_sum$nr_patients, levels = sort(unique(df_reg_genescores_sum$nr_patients), decreasing = T))
#         df_core_regulon_overview_list[[core_reg_name]]=
#             do.call(rbind, lapply(split(df_reg_genescores_sum, df_reg_genescores_sum$nr_patients), function(x) x[order(x$median_corr,decreasing = T),]))
#         
#     }
#     
#     # Create output
#     core_regulons_sorted = lapply(df_core_regulon_overview_list, function(x) {x$gene})
#     matrix_core_regulons_padded = t(plyr::ldply(core_regulons_sorted, rbind))
#     df_core_regulons_padded = data.frame(matrix_core_regulons_padded[-1,])
#     colnames(df_core_regulons_padded)=matrix_core_regulons_padded[1,]
#     openxlsx::write.xlsx(x= df_core_regulons_padded, file = paste0(base_dir,'Rplots/',DATASET_NAME,'_core_regulons.xlsx'), overwrite = T)
#     
#     df_core_regulons_padded_shortname=as.data.frame(shorthand_cutname_table(df_core_regulons_padded))
#     openxlsx::write.xlsx(x= df_core_regulons_padded_shortname, file = paste0(base_dir,'Rplots/',DATASET_NAME,'_core_regulons_shortname.xlsx'), overwrite = T)
#     
# }
# 
# if (F) {
# 
#     core_regulon_list=SCENIC_regulons_core_genes
#     
#     
#     
# 
# }    
    
################################################################################

    
    
    
    
################################################################################
# Code graveyard
    
# # Now look at results
# lfile <- loomR::connect(filename='/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/SCENIC/R.P1/output/SCENIC_out_R.P1RID2l.loom', 
#                         mode = 'r')
# 
# lfile[["matrix"]][1:5, 1:5]
# 
# full.matrix <- lfile$matrix[, ]
# dim(x = full.matrix) 
#     # 306 x 5833
#     # note: these are the cells x genes again
# 
# 
# # So, e.g. 
# lfile$row.attrs$Gene[1:3]
#     # This shows genes 1-3 of the original table
# dim(lfile$row.attrs$Regulons[1:3])
#     # This shows for each of the 127 regulons, whether the gene is a member or not (I suppose)
# 
# # So let's create a heatmap of this
# df_regulons_compare = as.data.frame(lfile$row.attrs$Regulons[])
# row.names(df_regulons_compare)=lfile$row.attrs$Gene[]
# 
# # Let's also gather the regulon members
# scenic_regulons_collected = lapply(1:ncol(df_regulons_compare), function(reg_idx) {rownames(df_regulons_compare)[df_regulons_compare[,reg_idx]==1]})
# names(scenic_regulons_collected) = colnames(df_regulons_compare)
# 
# # p=pheatmap(as.matrix(df_regulons_compare))
# 
# 
# regulon_overlap_heatmap(pooled_regulons = scenic_regulons_collected, base_dir = base_dir, run_name = 'SCENIC_test', MYTREECUTTINGHEIGHT = 1, myfontsize = 2)
# 
#     
#     
# 
# # Conversely, we also have regulon "activity" over the cells
# # (±composite expression)
# # lfile$col.attrs$
# 
# 
# 
# # And 
# lfile[["row_attrs/Regulons"]][1:3]
# 
# rownames(full.matrix) = lfile$row.attrs$Regulons[[]]
# lfile[["row_attrs/Regulons"]][1]
# lfile[["row_attrs/Regulons"]][2]
# lfile[["row_attrs/Regulons"]][2]
# 
# lfile[["row_attrs/Genes"]]
# 
# 
# 
# lfile$col.attrs
# 
# df_1 = 
#     data.frame(auc=lfile$col.attrs$RegulonsAUC, 
#             cellid=lfile$col.attrs$CellID)
# 
# regulons = lfile$row.attrs$Regulons
# 
# lfile$row.attrs$Gene
# 
# 
# #####
# 
# attrs <- c("nUMI", "nGene", "orig.ident")
# attr.df <- lfile$get.attribute.df(MARGIN = 2, attribute.names = attrs)
# head(x = attr.df)



