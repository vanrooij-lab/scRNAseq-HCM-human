

# Basic function to convert mouse to human gene names
source(paste0(script_dir, '/Functions/Conversion_mouse_human_genes.R'))


################################################################################

# Let's obtain the nomura data (see below for loading humanized data)
Nomura_data = openxlsx::read.xlsx('/Users/m.wehrens/Data/_2019_02_HCM_SCS/Nomura_data/41467_2018_6639_MOESM7_ESM.xlsx')
genes_in_Nomura_hypertr_modules = Nomura_data[Nomura_data$assigned.module %in% c('M1', 'M2', 'M5', 'M11', 'M16'),]$gene.name

human_eq_genes_in_Nomura_hypertr_modules = convertMouseGeneList(genes_in_Nomura_hypertr_modules)

save(list='human_eq_genes_in_Nomura_hypertr_modules', file=paste0(base_dir, 'Rdata/human_eq_genes_in_Nomura_hypertr_modules.Rdata'))
# load(file=paste0(base_dir, 'Rdata/human_eq_genes_in_Nomura_hypertr_modules.Rdata')) # human_eq_genes_in_Nomura_hypertr_modules

# Now per Nomura module

hyplist = c('M1', 'M2', 'M5', 'M11', 'M16')
names(hyplist) = hyplist
Nomura_hyp_lists = lapply( hyplist , function(X) {Nomura_data[Nomura_data$assigned.module == X,]$gene.name})

Nomura_hyp_lists_humanized = lapply(Nomura_hyp_lists, convertMouseGeneList)

lapply(Nomura_hyp_lists_humanized, length)

save(list='Nomura_hyp_lists_humanized', file=paste0(base_dir, 'Rdata/Nomura_hyp_lists_humanized.Rdata'))
# load(file=paste0(base_dir, 'Rdata/Nomura_hyp_lists_humanized.Rdata')) # Nomura_hyp_lists_humanized

################################################################################

# Now load the SCENIC regulons and custom modules
# SCENIC regulons
# load(paste0(base_dir, 'Rdata/SCENIC_regulons_core_genes_sel.Rdata')) # SCENIC_regulons_core_genes_sel
load(paste0(base_dir, 'Rdata/SCENIC_reg_top_genes_sorted_full.Rdata')) # SCENIC_reg_top_genes_sorted_full
# Custom modules
ANALYSIS_NAME = "ROOIJonly.sp.bt_RID2l"
ANALYSIS_NAME_clExtended = paste0(ANALYSIS_NAME, '_clExtended')
# load(paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_sorted_shortname.Rdata')) # core_regulons_sorted_shortname
# load(file=paste0(base_dir,'Rdata/',DATASET_NAME,'__SCENIC_reg_top_genes_sorted_full.Rdata')) # SCENIC_reg_top_genes_sorted_full
load(file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_core_regulons_sorted.Rdata')) # core_regulons_sorted
core_regulons_sorted_shortname = lapply(core_regulons_sorted, shorthand_cutname)


################################################################################
# Calculate overlaps with specific modules/regulons

quick_overlap = function(A,B) {
    sum(A %in% B)/
        min(  length(A), length(B) )
}

# For SCENIC

overlap_df_SCENIC = expand.grid(names(Nomura_hyp_lists_humanized), names(SCENIC_reg_top_genes_sorted_full))

overlap_df_SCENIC$overlap = sapply(1:(dim(overlap_df_SCENIC)[1]), function(idx) {quick_overlap(Nomura_hyp_lists_humanized[[   overlap_df_SCENIC$Var1[idx]   ]], 
                                                                                   SCENIC_reg_top_genes_sorted_full[[  overlap_df_SCENIC$Var2[idx]  ]])})

ggplot(overlap_df_SCENIC, aes(x=Var1, y=Var2, fill=overlap))+
    geom_tile()

# Matrix

overlap_matrix_SCENIC = matrix(rep(NA, dim(overlap_df_SCENIC)[1]), ncol = length(Nomura_hyp_lists_humanized))
rownames(overlap_matrix_SCENIC)=names(SCENIC_reg_top_genes_sorted_full)
colnames(overlap_matrix_SCENIC)=names(Nomura_hyp_lists_humanized)
for (idx in 1:nrow(overlap_df_SCENIC)) {
    overlap_matrix_SCENIC[overlap_df_SCENIC$Var2[idx],  overlap_df_SCENIC$Var1[idx]] =  overlap_df_SCENIC$overlap[idx]
}

pheatmap(overlap_matrix_SCENIC, cluster_cols = F)
sums_overlap_matrix_SCENIC = apply(overlap_matrix_SCENIC, 1, sum)
    
View(data.frame(sums_overlap_matrix_SCENIC))

# View(overlap_df_SCENIC)

# For the modules
# ====

overlap_df_MODULES = expand.grid(names(Nomura_hyp_lists_humanized), names(core_regulons_sorted_shortname))

overlap_df_MODULES$overlap = sapply(1:(dim(overlap_df_MODULES)[1]), function(idx) {quick_overlap(Nomura_hyp_lists_humanized[[   overlap_df_MODULES$Var1[idx]   ]], 
                                                                                           core_regulons_sorted_shortname[[  overlap_df_MODULES$Var2[idx]  ]])})

ggplot(overlap_df_MODULES, aes(x=Var1, y=Var2, fill=overlap))+
    geom_tile()

# View(overlap_df_MODULES)

# Matrix

overlap_matrix_MODULES = matrix(rep(NA, dim(overlap_df_MODULES)[1]), ncol = length(Nomura_hyp_lists_humanized))
rownames(overlap_matrix_MODULES)=names(core_regulons_sorted_shortname)
colnames(overlap_matrix_MODULES)=names(Nomura_hyp_lists_humanized)
for (idx in 1:nrow(overlap_df_MODULES)) {
    overlap_matrix_MODULES[overlap_df_MODULES$Var2[idx],  overlap_df_MODULES$Var1[idx]] =  overlap_df_MODULES$overlap[idx]
}

pheatmap(overlap_matrix_MODULES, cluster_cols = F, cluster_rows = F)

################################################################################




# now look at some stats
sum(core_regulons_sorted_shortname[[1]] %in% human_eq_genes_in_Nomura_hypertr_modules)
length(shared_regulon_genes_list[[1]])
module_nomura_overlap = sapply(1:5,function(X){
    sum(core_regulons_sorted_shortname[[X]] %in% human_eq_genes_in_Nomura_hypertr_modules)/
        min(  length(core_regulons_sorted_shortname[[X]]), length(human_eq_genes_in_Nomura_hypertr_modules) )  })
names(module_nomura_overlap) = names(core_regulons_sorted_shortname)
module_nomura_overlap
round(module_nomura_overlap*100,0)

SCENIC_nomura_overlap = sapply(1:length(SCENIC_reg_top_genes_sorted_full),function(X){
    sum(SCENIC_reg_top_genes_sorted_full[[X]] %in% human_eq_genes_in_Nomura_hypertr_modules)/
        min(   length(human_eq_genes_in_Nomura_hypertr_modules), length(SCENIC_reg_top_genes_sorted_full[[X]])   ) })
names(SCENIC_nomura_overlap)=names(SCENIC_reg_top_genes_sorted_full)
SCENIC_nomura_overlap
View(data.frame(names(SCENIC_reg_top_genes_sorted_full), SCENIC_nomura_overlap))
names(SCENIC_nomura_overlap) = names(SCENIC_reg_top_genes_sorted_full)
names(SCENIC_reg_top_genes_sorted_full)[max(SCENIC_nomura_overlap)==SCENIC_nomura_overlap]
sapply(SCENIC_reg_top_genes_sorted_full, length)
round(sort(SCENIC_nomura_overlap, decreasing = T)*100,0)

ggplot(data.frame(regulon=names(module_nomura_overlap), overlap=module_nomura_overlap*100), 
       aes(x=regulon, y=overlap))+
    geom_bar(stat='identity')+theme_bw()+coord_flip()

ggplot(data.frame(regulon=names(SCENIC_nomura_overlap), overlap=SCENIC_nomura_overlap*100), 
       aes(x=regulon, y=overlap))+
    geom_bar(stat='identity')+theme_bw()+coord_flip()

################################################################################

# Create Nomura overlap sets of interest

# Set 1: module 2 - nomura overlapping genes
genelist_module2_nomura_overlap = 
    core_regulons_sorted_shortname[[2]][core_regulons_sorted_shortname[[2]] %in% human_eq_genes_in_Nomura_hypertr_modules]
genelist_module2_nomura_overlap

save(list='genelist_module2_nomura_overlap', file = paste0(base_dir,'Rdata/zcustom__genelist_module2_nomura_overlap.Rdata'))
# load(paste0(base_dir,'Rdata/zcustom__genelist_module2_nomura_overlap.Rdata')) # genelist_module2_nomura_overlap

# Set 2: SCENIC top overlapping - nomura overlapping genes
# (Note: in earlier analysis, FOXN3 was top-overlapping gene w/ Nomura, now it's CEBPZ (FOXN3 is now 2nd))
genelist_SCENIC.FOXN3_nomura_overlap = 
    SCENIC_reg_top_genes_sorted_full[['FOXN3']][SCENIC_reg_top_genes_sorted_full[['FOXN3']] %in% human_eq_genes_in_Nomura_hypertr_modules]
genelist_SCENIC.FOXN3_nomura_overlap

save(list='genelist_SCENIC.FOXN3_nomura_overlap', file = paste0(base_dir,'Rdata/zcustom__genelist_SCENIC.FOXN3_nomura_overlap.Rdata'))
# load(paste0(base_dir,'Rdata/zcustom__genelist_SCENIC.FOXN3_nomura_overlap.Rdata')) # genelist_SCENIC.FOXN3_nomura_overlap


# Set 2b: SCENIC top overlapping (CEBPZ) members with Nomura genes
genelist_SCENIC.CEBPZ_nomura_overlap = 
    SCENIC_reg_top_genes_sorted_full[['CEBPZ']][SCENIC_reg_top_genes_sorted_full[['CEBPZ']] %in% human_eq_genes_in_Nomura_hypertr_modules]
genelist_SCENIC.CEBPZ_nomura_overlap

save(list='genelist_SCENIC.CEBPZ_nomura_overlap', file = paste0(base_dir,'Rdata/zcustom__genelist_SCENIC.CEBPZ_nomura_overlap.Rdata'))






