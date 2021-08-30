
lisa_results_collected = list(); p_list = list()
for (cl_idx in 1:5) {
    
    lisa_header = 
        read.table(paste0('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/LISA/strictbg_oneshot_cl.',cl_idx,'_out.lisa.tsv'), sep='\t', nrows=1)
    lisa_result = 
        read.table(paste0('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/LISA/strictbg_oneshot_cl.',cl_idx,'_out.lisa.tsv'), sep='\t', skip=1)
    colnames(lisa_result) = lisa_header
    
    lisa_header_down = 
        read.table(paste0('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/LISA/strictbg_oneshot_Downcl.',cl_idx,'_out.lisa.tsv'), sep='\t', nrows=1)
    lisa_result_down = 
        read.table(paste0('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/LISA/strictbg_oneshot_Downcl.',cl_idx,'_out.lisa.tsv'), sep='\t', skip=1)
    colnames(lisa_result_down) = lisa_header_down
    
    lisa_joined = 
        dplyr::left_join(lisa_result, lisa_result_down, by = 'sample_id', suffix = c('_up','_down'))
    lisa_joined = lisa_joined[order(lisa_joined$summary_p_value_up, decreasing = F),]
    
    print(paste0('Cl ', cl_idx, ': ', 
        toString(unique(lisa_joined[lisa_joined$cell_type_up=='Cardiomyocyte'&lisa_joined$"summary_p_value_up"<.05,]$factor_up))))
    
    interesting_results_up = unique(lisa_joined[lisa_joined$summary_p_value_up<1e-4,]$factor_up)
    down_stats_per_factor = 
        sapply(interesting_results_up, function(f) {min(lisa_joined[lisa_joined$factor_down %in% f,]$summary_p_value_down)})
    interesting_results_up_filtered = interesting_results_up[down_stats_per_factor>.05]
    
    lisa_results_collected[[cl_idx]] = interesting_results_up_filtered
    current_lisa_results_top20 = lisa_results_collected[[cl_idx]][1:(min(20, length(lisa_results_collected[[cl_idx]])))]
    
    library(ggrepel)
    p_list[[cl_idx]]=ggplot(lisa_joined, aes(x=-log10(summary_p_value_up), y=-log10(summary_p_value_down)))+
        geom_point()+
        geom_point(data=lisa_joined[lisa_joined$factor_up %in% current_lisa_results_top20,], aes(color=factor_up))+
        #geom_text_repel(data=lisa_joined[lisa_joined$factor_up %in% current_lisa_results_top20,], 
        #                aes(x=-log10(summary_p_value_up), y=-log10(summary_p_value_down), label=factor_up), max.overlaps = Inf, force = 30, color='red')+
        theme_bw()+
        give_better_textsize_plot(8)
    
    
    
}

wrap_plots(p_list)




