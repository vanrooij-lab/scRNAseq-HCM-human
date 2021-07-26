
################################################################################
# Calculate mito totals in original Wang/Hu data

gene_names = rownames(RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@data)
mito_gene_names = gene_names[grepl('^MT-',rownames(RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@data))]

mito_totals=
    rowSums(RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@data[mito_gene_names,])

ggplot(data.frame(gene=names(mito_totals), expression=mito_totals))+
    geom_bar(aes(x=gene, y=expression), stat='identity')+
    coord_flip()+theme_bw()

################################################################################
# Now compare with my current (WANG4) mapping

load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/Rdata/RHL_SeuratObject_nM_sel_HUonly.Rdata')

mito_table_HU_WANG4 = RHL_SeuratObject_nM_sel_HUonly@misc$mito.counts[,colnames(RHL_SeuratObject_nM_sel_HUonly)]

mito_totals_HU_WANG4 =
    rowSums(mito_table_HU_WANG4)

ggplot(data.frame(gene=names(mito_totals_HU_WANG4), expression=mito_totals_HU_WANG4))+
    geom_bar(aes(x=gene, y=expression), stat='identity')+
    coord_flip()+theme_bw()

################################################################################
# Actual scatter comparison

ggplot(data.frame(gene=names(mito_totals), 
                    expression_new=mito_totals_HU_WANG4[names(mito_totals)],
                    expression_old=mito_totals),
        aes(x=expression_old, y=expression_new))+
    geom_point()+
    geom_text_repel(aes(label=gene))+
    theme_bw()+geom_abline(slope = 1,intercept = 0)


################################################################################
# Compare some of the pseudogene counts

pseudo_gene_selection=c('MTND4P12', 'MTATP6P1', 'MTND1P23', 'MTND2P28', 'MTCO1P12', 'MTRNR2L12', 'MTND4P35', 'MTND6P4', 'MTCO1P40')

psuedo_expr_ori = rowSums(RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@data[pseudo_gene_selection[pseudo_gene_selection%in% rownames(RHL_SeuratObject_merged_WANGORIGINAL@assays$RNA@data)],])
psuedo_expr_new = rowSums(RHL_SeuratObject_nM_sel_HUonly@assays$RNA@data[pseudo_gene_selection,])

df_psuedo = data.frame(expression_new=psuedo_expr_new[names(psuedo_expr_new)],
                    expression_old=psuedo_expr_ori[names(psuedo_expr_new)],
                    gene=names(psuedo_expr_new))
df_psuedo[is.na(df_psuedo)]=0


mylim=c(0,max(c(psuedo_expr_ori,psuedo_expr_new)))
ggplot(df_psuedo,
        aes(x=expression_old, y=expression_new))+
    geom_point()+
    geom_text_repel(aes(label=gene))+
    theme_bw()+geom_abline(slope = 1,intercept = 0)+coord_fixed(ratio = 1, xlim = mylim, ylim = mylim)



