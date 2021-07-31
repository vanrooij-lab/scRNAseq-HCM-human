

# Violin plot
# Use lapply and wrap_plots to create a list of plots
shorthand_plotViolinBox_custom = function(myseuratobjectlist, analysis_name, cat_by='seurat_clusters', cat_name = 'category', gene_of_interest,base_dir, percentile=.02, type='boxplot') {

    # get full gene names
    gene_of_interest_fullname = rownames(myseuratobjectlist[[analysis_name]])[grepl(paste0(paste0(':',gene_of_interest,'$'),collapse='|'),rownames(myseuratobjectlist[[analysis_name]]))]
    
    # create df for plotting
    #plot_df = data.frame(count=numeric(), cat=character(), gene=character())
    #for (gene in gene_of_interest_fullname) {
        # retrieve current expression
        #current_expr = myseuratobjectlist[[analysis_name]]@assays$RNA@data[gene,]
        #expr_limits=calc_limits(current_expr, percentile = .03)

        #plot_df = rbind(plot_df, data.frame(count=current_expr, cat=myseuratobjectlist[[analysis_name]][[cat_by]], gene=gene))
    #}
    
    # retrieve expression
    current_expr = myseuratobjectlist[[analysis_name]]@assays$RNA@data[gene_of_interest_fullname,]
    print(paste0('Retrieved gene ',gene_of_interest_fullname))
    
    # create dataframe to plot
    plot_df = data.frame(count=current_expr, cat=myseuratobjectlist[[analysis_name]][[cat_by]], gene=gene_of_interest_fullname)
    colnames(plot_df) = c('count',cat_name)
    
    # decide on limits
    # per group, determine where 2% percentile lies to set that as maximum to the plot
    my_max_limit=max(sapply(unique(plot_df[[cat_name]]), function(x) { calc_limits(current_expr[plot_df[[cat_name]]==x], percentile = .02)[2] } ))

    # start plotting w/ desired options
    p=ggplot(plot_df, aes_string(x=cat_name, y='count',fill=cat_name))
    
    if (type=='violin') {p=p+geom_violin()}
    if (type=='boxplot') {p=p+geom_boxplot()}
    
    p=p+
        #geom_bar(aes_string(y='count', fill=as.factor(cat_name)),stat='summary',fun='mean',fill='grey',color='black')+
        #geom_boxplot(alpha=0)+
        ggtitle(gene_of_interest)+theme_bw()+give_better_textsize_plot(10)+
        theme(legend.position='none')+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    # save 1
    ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',analysis_name,'_6_violin_by-', cat_name, '_',gene_of_interest,'.pdf'), width = 4, height= 6, units='cm')
    
    # add ylim and save 2
    p=p+ylim(c(0,my_max_limit))
    ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',analysis_name,'_6_violin_ylim98_by-', cat_name, '_',gene_of_interest,'.pdf'), width = 4, height= 6, units='cm')
    
    # return plot
    return(p)
    
}
