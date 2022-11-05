
################################################################################

# Looking more into NPPA-correlated genes
# Load this volcano by execyting ..geneCorrelationsOfInterest.. file, or 
# directly load this volcano
CURRENT_GENE='ENSG00000175206:NPPA'
CURRENT_DATASET="ROOIJonly.sp.bt_RID2l"
load(file = paste0(base_dir,'Rdata/Volcano_df_gene_pt__',CURRENT_GENE,'_',CURRENT_DATASET,'.Rdata'))
Volcano_df_NPPA_perpatient=Volcano_df_gene_pt[[CURRENT_GENE]][[CURRENT_DATASET]]
rm('Volcano_df_gene_pt')
# However, some convenient aggregate dataframes where already created; so we can just use those
# Execute earlier code in ..geneCorrelationsOfInterest.. to obtain 
# df_corr_collection$melt$`ENSG00000175206:NPPA`$ROOIJonly.sp.bt_RID2l
Volcano_df_NPPA_melted = df_corr_collection$melt$`ENSG00000175206:NPPA`$ROOIJonly.sp.bt_RID2l
        
################################################################################

# determine limits for color scale

p=ggplot(Volcano_df_NPPA_melted)+
    geom_freqpoly(aes(x=corr))+
    theme_bw()+give_better_textsize_plot(8)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
ggsave(filename = paste0(base_dir,'Rplots/2022_10-extra/','ROOIJonly.sp.bt_RID2l','_9_DistributionCorrsWNppa.pdf'), 
                             plot = p, height=30, width=30, units = 'mm', device=cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5

calc_limits(Volcano_df_NPPA_melted$corr,.01)[1]
# -0.113439
calc_limits(Volcano_df_NPPA_melted$corr,.01)[2]
# 0.194202

LIMITS_COLORSCALE = c(-.25, .25)

# Now select genes that are significantly correlated, and are a ligand
# Volcano_df_NPPA_melted_selLigSig 
genes_selLigSig = unique(Volcano_df_NPPA_melted[Volcano_df_NPPA_melted$pval.adj.sign & (shorthand_cutname(Volcano_df_NPPA_melted$gene_name, PART1OR2 = 1) %in% ligand_list_unique_ENS),]$gene_name)
genes_selLigSig_corrsum = aggregate(list(sum_corr=Volcano_df_NPPA_melted_selLigSig$corr), by=list(gene=Volcano_df_NPPA_melted_selLigSig$gene_name_short), FUN=sum)
genes_selLigSig_ranked = genes_selLigSig_corrsum$gene[order(genes_selLigSig_corrsum$sum_corr, decreasing = F)]
Volcano_df_NPPA_melted_selLigSig = Volcano_df_NPPA_melted[Volcano_df_NPPA_melted$gene_name %in% genes_selLigSig,]
Volcano_df_NPPA_melted_selLigSig$gene_name_short_sortedFct = factor(Volcano_df_NPPA_melted_selLigSig$gene_name_short, levels=genes_selLigSig_ranked)
p=ggplot(Volcano_df_NPPA_melted_selLigSig, aes(x=donor, y=gene_name_short_sortedFct, fill=corr))+
    geom_tile()+
    geom_point(data=Volcano_df_NPPA_melted_selLigSig[Volcano_df_NPPA_melted_selLigSig$pval.adj.sign,], size=2, color='black')+
    #scale_fill_gradientn(colors=rev(pals::viridis(n=100)))+
    scale_fill_gradientn(colours = rainbow_colors, breaks=seq(-1,1,.05), limits=LIMITS_COLORSCALE, oob=squish)+
    theme_bw()+give_better_textsize_plot(8)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab('Patient')+ylab('Gene')+ggtitle('Corr w/ NPPA, ligs')
p
ncol_effective=length(unique(Volcano_df_NPPA_melted_selLigSig$donor))
nrow_effective=length(unique(Volcano_df_NPPA_melted_selLigSig$gene_name_short))
ggsave(filename = paste0(base_dir,'Rplots/2022_10-extra/','ROOIJonly.sp.bt_RID2l','_9_LigandsCorrelatedWithNPPA_withLegend.pdf'), 
                             plot = p, height=min(172,3*nrow_effective+20), width=1/3*172-4, units = 'mm', device=cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5
p2=p+theme(legend.position = 'none')
ggsave(filename = paste0(base_dir,'Rplots/2022_10-extra/','ROOIJonly.sp.bt_RID2l','_9_LigandsCorrelatedWithNPPA.pdf'), 
                             plot = p2, height=min(172,3*nrow_effective+20), width=1/3*172-4, units = 'mm', device=cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5

# Repeat for TFs           
# Now select genes that are significantly correlated, and are a *transcription factor*
# Volcano_df_NPPA_melted_selTFSig 
genes_selTFSig = unique(Volcano_df_NPPA_melted[Volcano_df_NPPA_melted$pval.adj.sign & (shorthand_cutname(Volcano_df_NPPA_melted$gene_name, PART1OR2 = 2) %in% human_TFs_combined_temp),]$gene_name)
genes_selTFSig_corrsum = aggregate(list(sum_corr=Volcano_df_NPPA_melted[Volcano_df_NPPA_melted$gene_name %in% genes_selTFSig,]$corr), by=list(gene=Volcano_df_NPPA_melted[Volcano_df_NPPA_melted$gene_name %in% genes_selTFSig,]$gene_name_short), FUN=sum)
genes_selTFSig_ranked = genes_selTFSig_corrsum$gene[order(genes_selTFSig_corrsum$sum_corr, decreasing = F)]
Volcano_df_NPPA_melted_selTFSig = Volcano_df_NPPA_melted[Volcano_df_NPPA_melted$gene_name %in% genes_selTFSig,]
Volcano_df_NPPA_melted_selTFSig$gene_name_short_sortedFct = factor(Volcano_df_NPPA_melted_selTFSig$gene_name_short, levels=genes_selTFSig_ranked)
p=ggplot(Volcano_df_NPPA_melted_selTFSig, aes(x=donor, y=gene_name_short_sortedFct, fill=corr))+
    geom_tile()+
    geom_point(data=Volcano_df_NPPA_melted_selTFSig[Volcano_df_NPPA_melted_selTFSig$pval.adj.sign,], size=2, color='black')+
    #scale_fill_gradientn(colors=rev(pals::viridis(n=100)))+
    scale_fill_gradientn(colours = rainbow_colors, breaks=seq(-1,1,.05), limits=LIMITS_COLORSCALE, oob=squish)+
    theme_bw()+give_better_textsize_plot(8)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab('Patient')+ylab('Gene')+ggtitle('Corr w/ NPPA, TFs')
p
ncol_effective=length(unique(Volcano_df_NPPA_melted_selTFSig$donor))
nrow_effective=length(unique(Volcano_df_NPPA_melted_selTFSig$gene_name_short))
ggsave(filename = paste0(base_dir,'Rplots/2022_10-extra/','ROOIJonly.sp.bt_RID2l','_9_TFsCorrelatedWithNPPA_withLegend.pdf'), 
                             plot = p, height=min(172,3*nrow_effective+20), width=1/3*172-4, units = 'mm', device=cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5
p2=p+theme(legend.position = 'none')
ggsave(filename = paste0(base_dir,'Rplots/2022_10-extra/','ROOIJonly.sp.bt_RID2l','_9_TFsCorrelatedWithNPPA.pdf'), 
                             plot = p2, height=min(172,3*nrow_effective+20), width=1/3*172-4, units = 'mm', device=cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5
 
################################################################################
# Now look at genes enriched in a cluster

load(paste0(base_dir,'Rdata/DE_cluster__',ANALYSIS_NAME_clExtended,'.Rdata')) # DE_cluster

View(DE_cluster$ROOIJonly.sp.bt_RID2l_clExtended$`3`)
genes_DE_cluster3_df = DE_cluster$ROOIJonly.sp.bt_RID2l_clExtended$`3`

genes_DE_cluster3_df$gene = rownames(genes_DE_cluster3_df)
genes_DE_cluster3_df$ens = shorthand_cutname( rownames(genes_DE_cluster3_df), PART1OR2 = 1)
genes_DE_cluster3_df$sym = shorthand_cutname( rownames(genes_DE_cluster3_df), PART1OR2 = 2)

View(genes_DE_cluster3_df[genes_DE_cluster3_df$ens %in% ligand_list_unique_ENS,])

table_to_export = genes_DE_cluster3_df[genes_DE_cluster3_df$ens %in% ligand_list_unique_ENS,]
table_to_export=table_to_export[order(table_to_export$avg_log2FC, decreasing = T),c('sym','avg_log2FC','p_val','p_val_adj')]
View(table_to_export)

# Let's do the same for transription factors

table_to_export_TF = genes_DE_cluster3_df[genes_DE_cluster3_df$sym %in% human_TFs_combined_temp,]
table_to_export_TF=table_to_export_TF[order(table_to_export_TF$avg_log2FC, decreasing = T),c('sym','avg_log2FC','p_val','p_val_adj')]
View(table_to_export_TF)

human_TFs_combined_temp








