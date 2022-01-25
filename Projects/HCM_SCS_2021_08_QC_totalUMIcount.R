


    

script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

current_analysis$Rooij_raw = 
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_Rooij_raw.h5seurat'))

current_analysis$Hu_raw = 
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_Hu_raw.h5seurat'))
        
current_analysis$Rooij_raw

mito_genes = rownames(current_analysis$Rooij_raw)[grepl(':MT-',rownames(current_analysis$Rooij_raw))]
current_analysis$Rooij_raw@misc$desired_genes_excl_mito = 
            rownames(current_analysis$Rooij_raw)[!(rownames(current_analysis$Rooij_raw) %in% mito_genes)]
current_analysis$Rooij_raw_nMt =
        subset(current_analysis$Rooij_raw, features= current_analysis$Rooij_raw@misc$desired_genes_excl_mito)
cell_counts = colSums(current_analysis$Rooij_raw_nMt@assays$RNA@counts)

# spit out some stats
median_UMI_count = median(cell_counts[cell_counts>1e3])
median_UMI_count # 2201.5
cells_included = sum(cell_counts>1e3)
cells_included # 2292

# Make histogram with color-coding
hist_out = hist(log10(.1+cell_counts), breaks = seq(-1.1,5,.1), plot=F)
plot_df = data.frame(mids=hist_out$mids, counts=hist_out$counts)
plot_df$Status = 'Excluded'
plot_df$Status[plot_df$mids>log10(.1+1000)] = 'Included'
plot_df$Status=factor(plot_df$Status,levels = c('Included','Excluded'))

p=ggplot(plot_df)+
    geom_bar(aes(x=mids, y=counts, fill=Status), stat='Identity')+
    #geom_histogram(aes(x=log10(.1+cell_counts)), bins=50)+
    geom_vline(xintercept = log10(.1+1000))+
    theme_bw()+ylab('Number of cells')+xlab('log10(UMI count + 0.1)')+ggtitle('Total UMI count distribution')+
    give_better_textsize_plot(8)+theme(legend.position='none')
p
ggsave(filename = paste0(base_dir,'Rplots/','QC_Rooij_raw','_0_totalUMI_perCell.pdf'), 
    plot = p, height=172/3-4, width=172/3-4, units='mm', device = cairo_pdf)
p=p+theme(legend.key.size = unit(2, 'mm'), legend.position=c(0.1, 0.9))
ggsave(filename = paste0(base_dir,'Rplots/','QC_Rooij_raw','_0_totalUMI_perCell_Legend.pdf'), 
    plot = p, height=172/3-4, width=172/3-4, units='mm', device = cairo_pdf)



ggplot(data.frame(UMI=cell_counts))+
    geom_histogram(aes(x=UMI), fill='lightblue')+
    geom_freqpoly(data=data.frame(UMI=cell_counts[cell_counts>1000]), aes(x=UMI), color='red', alpha=.5)+
    theme_bw()+
    xlim(0,10000)+geom_vline(xintercept = median(cell_counts[cell_counts>1000]))
