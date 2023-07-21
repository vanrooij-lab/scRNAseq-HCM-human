


# Follow up code 2023-06-05





if (!dir.exists(paste0(base_dir, 'Rtxt/'))) {dir.create(paste0(base_dir, 'Rtxt/'))}


# Note that there's already excel files with this too
current_TF = "NFE2L1"
write.table(x = data.frame(targets=SCENIC_reg_top_genes_sorted_full[[current_TF]]), 
            file = paste0(base_dir, 'Rtxt/',current_TF,'_targets.csv'), quote = F, col.names = F, row.names = F)


# Check some things
length(SCENIC_reg_top_genes_sorted_full$NFE2L1)
length(SCENIC_reg_top_genes_sorted_full$MAFK)

core_regulons_sorted_s.R.4_short = shorthand_cutname( core_regulons_sorted$s.R.4 )

# Intersect between module 4 and NFE2L1
venn_simple_plot_mw(list(module=core_regulons_sorted_s.R.4_short, reg_NFE2L1=SCENIC_reg_top_genes_sorted_full$NFE2L1))
toString(intersect(core_regulons_sorted_s.R.4_short, SCENIC_reg_top_genes_sorted_full$NFE2L1))
# Intersect between module 4 and MAFK
venn_simple_plot_mw(list(module=core_regulons_sorted_s.R.4_short, reg_MAFK=SCENIC_reg_top_genes_sorted_full$MAFK))
toString(intersect(core_regulons_sorted_s.R.4_short, SCENIC_reg_top_genes_sorted_full$MAFK))
toString(intersect(core_regulons_sorted_s.R.4_short, SCENIC_reg_top_genes_sorted_full$MAFK)[1:10])

venn_simple_plot_mw(list(reg_NFE2L1=SCENIC_reg_top_genes_sorted_full$NFE2L1, reg_MAFK=SCENIC_reg_top_genes_sorted_full$MAFK))

toString(intersect(SCENIC_reg_top_genes_sorted_full$NFE2L1, SCENIC_reg_top_genes_sorted_full$MAFK))
length(intersect(SCENIC_reg_top_genes_sorted_full$NFE2L1, SCENIC_reg_top_genes_sorted_full$MAFK))

 
intersect_MAFK_NFE2L1 =
    intersect(SCENIC_reg_top_genes_sorted_full$NFE2L1, SCENIC_reg_top_genes_sorted_full$MAFK)
toString(intersect(intersect_MAFK_NFE2L1, core_regulons_sorted_s.R.4_short))
length(intersect(intersect_MAFK_NFE2L1, core_regulons_sorted_s.R.4_short))

# Very little of MAFK+NFE2L1 are related to top 20 module 4
intersect(SCENIC_reg_top_genes_sorted_full$MAFK, core_regulons_sorted_s.R.4_short[1:20])
intersect(SCENIC_reg_top_genes_sorted_full$NFE2L1, core_regulons_sorted_s.R.4_short[1:20])
    # Note that NFE2L1 has the most interesting targets..

# write table with intersect of all 10 genes
write.table(x = 
    data.frame(intersect=intersect(intersect_MAFK_NFE2L1, core_regulons_sorted_s.R.4_short)),
    file = paste0(base_dir, 'Rtxt/','overlap_MAFK_NFE2L1_mod4.csv'), quote = F, col.names = F, row.names = F
    )
    
################################################################################
# Small more convenient plot..

genes_allfrom3sets = unique(c(core_regulons_sorted_s.R.4_short, SCENIC_reg_top_genes_sorted_full$NFE2L1, SCENIC_reg_top_genes_sorted_full$MAFK))
length(genes_allfrom3sets)

ranks.mod4   = 1:length(core_regulons_sorted_s.R.4_short)
names(ranks.mod4) = core_regulons_sorted_s.R.4_short
ranks.NFE2L1 = 1:length(SCENIC_reg_top_genes_sorted_full$NFE2L1)
names(ranks.NFE2L1) = SCENIC_reg_top_genes_sorted_full$NFE2L1
ranks.MAFK   = 1:length(SCENIC_reg_top_genes_sorted_full$MAFK)
names(ranks.MAFK) = SCENIC_reg_top_genes_sorted_full$MAFK

df_membergenes_3sets =
    data.frame(row.names = genes_allfrom3sets, 
               mod.4  = genes_allfrom3sets %in% core_regulons_sorted_s.R.4_short,
               NFE2L1 = genes_allfrom3sets %in% SCENIC_reg_top_genes_sorted_full$NFE2L1,
               MAFK   = genes_allfrom3sets %in% SCENIC_reg_top_genes_sorted_full$MAFK)

df_membergenes_3sets_ranked =
    data.frame(row.names = genes_allfrom3sets, 
               mod.4  = ranks.mod4[genes_allfrom3sets],
               NFE2L1 = ranks.NFE2L1[genes_allfrom3sets],
               MAFK   = ranks.MAFK[genes_allfrom3sets])
df_membergenes_3sets_ranked[is.na(df_membergenes_3sets_ranked)] = 500

pheatmap(1*df_membergenes_3sets, cluster_cols = F, cluster_rows = F)
pheatmap(1*df_membergenes_3sets, cluster_cols = F, cluster_rows = T, fontsize = 5, treeheight_row = 0, color = c('white','black'))
p=pheatmap(df_membergenes_3sets_ranked, cluster_cols = F, cluster_rows = T, fontsize = 3, treeheight_row = 0, color = c('white','grey'), legend = F)
p
if (!dir.exists(paste0(base_dir,'Rplots/NFE2L1MAFK/'))){dir.create(paste0(base_dir,'Rplots/NFE2L1MAFK/'))}
ggsave(plot=p, filename = paste0(base_dir,'Rplots/NFE2L1MAFK/heatmap_members_comparision.pdf'), 
                height=172, width=172/5, units='mm')

breaksList=1:160
p=pheatmap(df_membergenes_3sets_ranked, cluster_cols = F, cluster_rows = T, fontsize = 3, treeheight_row = 0, legend = T, color=rev(viridis(100)), breaks = breaksList)
ggsave(plot=p, filename = paste0(base_dir,'Rplots/NFE2L1MAFK/heatmap_members_comparision_ranked.pdf'), 
                height=172, width=172/5, units='mm')

df_membergenes_3sets_ranked_forexport = df_membergenes_3sets_ranked
df_membergenes_3sets_ranked_forexport[df_membergenes_3sets_ranked_forexport==500] = "-"

write.table(x = 
    df_membergenes_3sets_ranked_forexport,
    file = paste0(base_dir, 'Rtxt/','overlap_MAFK_NFE2L1_mod4_overviewAll.tsv'), quote = F, col.names = T, row.names = T, sep = "\t")









