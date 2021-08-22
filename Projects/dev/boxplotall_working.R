
library(reshape2)


# Some extra comparative plots between bulk level data
load(file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_sorted.Rdata')) # core_regulons_sorted

current_genes = core_regulons_sorted$s.R.2[1:10]

VlnPlot(current_analysis$ROOIJonly_RID2l, group.by = 'annotation_paper_fct', features = core_regulons_sorted$s.R.2[1:10])

current_analysis$ROOIJonly_RID2l$dummy = round(runif(length(current_analysis$ROOIJonly_RID2l$annotation_paper_fct))*2)

VlnPlot(current_analysis$ROOIJonly_RID2l, group.by = 'dummy', features = core_regulons_sorted$s.R.2[1:10])

shorthand_custom_boxplot = function(seuratObject_list, gene_lists, seuratObjectNameToTake) {
    
    for (current_list_name in gene_lists) {
        
        current_genes = gene_lists[[current_list_name]]
        
        current_exprs = as.matrix(current_analysis$ROOIJonly_RID2l@assays$RNA@data[current_genes,]) # removes dgMatrix class if there
        current_exprs = t(scale(t(current_exprs), center=F))
        rownames(current_exprs) = shorthand_cutname(rownames(current_exprs))
         
        df_ = data.frame(t(current_exprs))
        df_$dataset = current_analysis$ROOIJonly_RID2l$annotation_paper_fct
        df_$dataset = as.factor(current_analysis$ROOIJonly_RID2l$dummy) # REMOVE, ONLY FOR TESTING
        df = reshape2::melt(data = df_, id.vars = 'dataset', value.name = 'expression', variable.name='gene')
        
        p=ggplot(df, mapping=aes(x=gene, y=expression, fill=dataset))+
            #geom_violin(position='dodge')+
            geom_boxplot(position='dodge')+
            xlab('Gene')+ylab('UMI count (normalized)')+
            give_better_textsize_plot(8)+theme_bw()+
            theme(legend.position = 'right', legend.key.size = unit(3, "mm"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        ggsave(plot = p, paste0(base_dir, 'Rplots/', seuratObjectNameToTake, '_9_customBoxplotGenes_', current_list_name))
        
        ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_5_Barplot_PatientCluster_distr.pdf'), 
            plot = p, height=33.28, width=33.28, units='mm')
        
    }
}


    