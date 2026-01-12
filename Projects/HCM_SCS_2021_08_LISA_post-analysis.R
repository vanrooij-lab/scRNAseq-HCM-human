library(ggrepel)


lisa_results_collected = list(); lisa_results_collected_extrasel=list(); p_list = list()
lisa_results_collected_extrasel_top5=list(); lisa_results_collected_top5=list()
lisa_results_collected_extrasel_top5_new=list()
lisa_joined_list=list()
merged_list2_with_ann=c(); merged_list2_listCl=list()
for (cl_idx in 1:6) {
    
    lisa_header = 
        read.table(paste0(base_dir, 'LISA/oneshot_cl.',cl_idx,'_out.lisa.tsv'), sep='\t', nrows=1)
    lisa_result = 
        read.table(paste0(base_dir, 'LISA/oneshot_cl.',cl_idx,'_out.lisa.tsv'), sep='\t', skip=1)
    colnames(lisa_result) = lisa_header
    
    lisa_header_down = 
        read.table(paste0(base_dir, 'LISA/oneshot_Downcl.',cl_idx,'_out.lisa.tsv'), sep='\t', nrows=1)
    lisa_result_down = 
        read.table(paste0(base_dir, 'LISA/oneshot_Downcl.',cl_idx,'_out.lisa.tsv'), sep='\t', skip=1)
    colnames(lisa_result_down) = lisa_header_down
    
    lisa_joined = 
        dplyr::left_join(lisa_result, lisa_result_down, by = 'sample_id', suffix = c('_up','_down'))
    lisa_joined = lisa_joined[order(lisa_joined$summary_p_value_up, decreasing = F),]
    # View(lisa_joined)
    
    lisa_joined_list[[cl_idx]] = lisa_joined
    
    # lisa_joined$summary_p_value_up_adjusted = p.adjust(lisa_joined$summary_p_value_up)
    #ggplot(data.frame(p=lisa_joined$summary_p_value_up), aes(x=p))+
    #           geom_histogram(bins=50)
    # Note: I can't find p-val adjustment in the paper.
    # I don't think they apply it (also looking at code).
    # but from the histogram, it doesn't appear necessary. 
    # Usually you correct when there are multiple samples, 
    # and there's high noise that leads to false positives.
    # The latter might not be the case here;
    # Not sure about the role of the p-value in their strategy.
    # Perhaps could be viewed more as a measure of goodness
    # of fit of a certain TF.
    
    # First select identified Lisa factors
    lisa_joined_selection = lisa_joined[lisa_joined$summary_p_value_up<.01,]
    lisa_joined_sel_perFactor = split(lisa_joined_selection, f = lisa_joined_selection$factor_up)
    lisa_perfactor_celltypes = lapply(lisa_joined_sel_perFactor, function(X) {
        length(unique(X$cell_type_up))
    })
    lisa_joined_selection$celltype_count = lisa_perfactor_celltypes[lisa_joined_selection$factor_up]
    
    # Then select Lisa factors that are detected specifically in CMs or up in multiple cell types
    lisa_joined_selection_sel2 = lisa_joined_selection[lisa_joined_selection$celltype_count>1|lisa_joined_selection$cell_type_up=='Cardiomyocyte',]
    
    # Now resp. collect top hits & special hits, which are CM-specific or multi-up 
    lisa_results_collected[[cl_idx]] = unique(lisa_joined_selection$factor_up)
    lisa_results_collected_extrasel[[cl_idx]] = unique(lisa_joined_selection_sel2$factor_up)
    
    # Then save these respective two
    lisa_results_collected_top5[[cl_idx]] = lisa_results_collected[[cl_idx]][1:(min(5, length(lisa_results_collected[[cl_idx]])))]
    lisa_results_collected_extrasel_top5[[cl_idx]] = lisa_results_collected_extrasel[[cl_idx]][1:(min(5, length(lisa_results_collected_extrasel[[cl_idx]])))]
    
    # And also create a little extra list, showing which factors the special hits identified in addition (if any)
    lisa_results_collected_extrasel_top5_new[[cl_idx]] = lisa_results_collected_extrasel_top5[[cl_idx]][    !(lisa_results_collected_extrasel_top5[[cl_idx]] %in% lisa_results_collected_top5[[cl_idx]])       ]
    
    # Now create a second extra list, which is the merge of the two respective lists, but with indication which one is member of 2nd list
    merged_list2 = unique(c(lisa_results_collected_top5[[cl_idx]], lisa_results_collected_extrasel_top5[[cl_idx]]))
    merged_list2_annotation = sapply(merged_list2 %in% lisa_results_collected_extrasel_top5[[cl_idx]], function(X) {if (X) {'*'} else {''}})
    merged_list2_with_ann[cl_idx] = toString(paste0(merged_list2,merged_list2_annotation))
    merged_list2_listCl[[cl_idx]] = merged_list2
    
    # current_lisa_results_top20 = lisa_results_collected[[cl_idx]][1:(min(20, length(lisa_results_collected[[cl_idx]])))]
    
    # Some complicated filtering that wasn't useful
    
    # print(paste0('Cl ', cl_idx, ': ', 
    #     toString(unique(lisa_joined[lisa_joined$cell_type_up=='Cardiomyocyte'&lisa_joined$"summary_p_value_up"<.05,]$factor_up))))
    # 
    # interesting_results_up = unique(lisa_joined[lisa_joined$summary_p_value_up<1e-4,]$factor_up)
    # down_stats_per_factor = 
    #     sapply(interesting_results_up, function(f) {min(lisa_joined[lisa_joined$factor_down %in% f,]$summary_p_value_down)})
    # interesting_results_up_filtered = interesting_results_up[down_stats_per_factor>.05]
    # 
    # lisa_results_collected[[cl_idx]] = interesting_results_up_filtered
    # current_lisa_results_top20 = lisa_results_collected[[cl_idx]][1:(min(20, length(lisa_results_collected[[cl_idx]])))]
    
}

lisa_export_simple =
    data.frame(cluster_TFs=sapply(lisa_results_collected_top5, toString), 
               cluster_TFs_additional_stricter=sapply(lisa_results_collected_extrasel_top5_new, toString),
               merge=paste0(sapply(lisa_results_collected_top5, toString), ' (', sapply(lisa_results_collected_extrasel_top5_new, toString),')'),
               cluster_TFs_stricter=sapply(lisa_results_collected_extrasel_top5, toString),
               merge_ann = merged_list2_with_ann)


openxlsx::write.xlsx(x=lisa_export_simple, file=paste0(base_dir,'Rplots/Lisa_clusters_short_summary.xlsx'), overwrite = T)

for (cl_idx in 1:5) {    
    
    lisa_joined = lisa_joined_list[[cl_idx]]
    
    p_list[[cl_idx]]=ggplot(lisa_joined, aes(x=-log10(summary_p_value_up), y=-log10(summary_p_value_down)))+
        geom_point()+
        
        geom_point(data=lisa_joined[lisa_joined$cell_type_up == 'Cardiomyocyte'|lisa_joined$cell_type_down == 'Cardiomyocyte',], size=4, color='gray')+
        geom_point(data=lisa_joined[lisa_joined$factor_up %in% merged_list2_listCl[[cl_idx]],], aes(color=factor_up))+
        #geom_text_repel(data=lisa_joined[lisa_joined$factor_up %in% current_lisa_results_top20,], 
        #                aes(x=-log10(summary_p_value_up), y=-log10(summary_p_value_down), label=factor_up), max.overlaps = Inf, force = 30, color='red')+
        theme_bw()+ggtitle(paste0('Cluster ',cl_idx))+
        give_better_textsize_plot(8)

    
}

wrap_plots(p_list)




