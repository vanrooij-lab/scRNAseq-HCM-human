

# For convenience, little fn to change text size ggplot in one go
give_better_textsize_plot <- function(TEXTSIZE){
  theme(#legend.position="none",
        text = element_text(size=TEXTSIZE),
        axis.text = element_text(size=TEXTSIZE),
        plot.title = element_text(size=TEXTSIZE))
}



# Violin plot
# Use lapply and wrap_plots to create a list of plots
shorthand_plotViolinBox_custom = function(myseuratobjectlist, analysis_name, cat_by='seurat_clusters', 
    cat_name = 'category', gene_of_interest,base_dir, percentile=.02, type='boxplot') {

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
    
    if (type=='violin') {p=p+geom_violin(scale='width')}
    if (type=='boxplot') {p=p+geom_boxplot()}
    
    p=p+
        #geom_bar(aes_string(y='count', fill=as.factor(cat_name)),stat='summary',fun='mean',fill='grey',color='black')+
        #geom_boxplot(alpha=0)+
        ggtitle(gene_of_interest)+theme_bw()+
        theme(legend.position='none')+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        xlab(element_blank())+give_better_textsize_plot(8)
    
    # save 1
    ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',analysis_name,'_6_',type,'_by-', cat_name, '_',gene_of_interest,'.pdf'), width = 4, height= 6, units='cm')
    
    # add ylim and save 2
    p=p+ylim(c(0,my_max_limit))
    ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',analysis_name,'_6_',type,'_ylim98_by-', cat_name, '_',gene_of_interest,'.pdf'), width = 4, height= 6, units='cm')
    
    # return plot
    return(p)
    
}

shorthand_cutname = function(gene_names, PART1OR2=2) {
    
  # if something to be split
  if (any(grepl(gene_names, pattern = ':'))) {
    # split, also for genes w/o 2nd annotation, just return 1st part instead
    return(sapply(str_split(gene_names, pattern = ':'), function(Y) { if (Y[[2]]=="") {Y[[1]]} else {Y[[PART1OR2]]} }))
  } else {return(gene_names)}
    
}

shorthand_cutname_table = function(gene_table, PART1OR2=2) {
    
 sapply(gene_table, function(X) {
            sapply(str_split(X, pattern = ':'), function(Y) { if (is.na(Y)) {return(NA)} else if (Y[[2]]=="") {Y[[1]]} else {Y[[PART1OR2]]} })
            })   
}

################################################################################

# function to plot expression of a gene
shorthand_seurat_custom_expr = function(seuratObject, gene_of_interest, textsize=8, pointsize=1, mypercentile=0.03, custom_title=NULL,mymargin=0.1,zscore=F) {
    
    if (length(gene_of_interest)>1){
      print('You gave >1 genes, assuming you want composite expression.')  
      print('(Otherwise, use lapply & wrap_plots.)')
    }
    if (is.null(custom_title)) {
      mytitle = gene_of_interest
    } else { mytitle = custom_title }
  
    # gene_of_interest = 'KRT6A'
    print(paste0('Retrieving ',gene_of_interest))
    gene_of_interest_fullname = shorthand_seurat_fullgenename(seuratObject, gene_of_interest)
    
    # retrieve current expression
    current_expr = as.matrix(seuratObject@assays$RNA@data[gene_of_interest_fullname,])
    
    # Create composite expression
    if (length(gene_of_interest)>1){
      print('Calculating composite..')
      if (zscore) { scaled_expr = t(scale(t(current_expr),center = T, scale = T))
      } else { scaled_expr = t(scale(t(current_expr),center = F, scale = T)) }
      # current_expr=colSums(scaled_expr)
      current_expr=colMeans(scaled_expr)
    }
    
    expr_limits=calc_limits(current_expr, percentile = mypercentile)
    
    ggplot(data.frame(UMAP_1=seuratObject@reductions$umap@cell.embeddings[,1],
                      UMAP_2=seuratObject@reductions$umap@cell.embeddings[,2],
                      expr=current_expr),
            mapping = aes(x=UMAP_1, y=UMAP_2, color=expr))+
        geom_point(size = pointsize, stroke = 0, shape = 16)+
        scale_color_gradientn(colours=rainbow_colors, limits=c(0,expr_limits[2]), oob=squish)+
        #+ggtitle(gene_of_interest)+
        annotate("text", -Inf, Inf, label = mytitle, hjust = 0, vjust = 1, size=textsize / .pt)+
        theme_void()+
        #give_better_textsize_plot(textsize)+
        theme(legend.position = 'none', plot.margin = margin(mymargin,mymargin,0,0,'mm'))
        
        
}
# shorthand to get full names
shorthand_seurat_fullgenename = function(seuratObject, gene_names) {
    
    gene_of_interest_fullname =rownames(seuratObject)[grepl(paste0(paste0(':',gene_names,'$'),collapse = '|'),rownames(seuratObject))]
    return(gene_of_interest_fullname)
    
}






