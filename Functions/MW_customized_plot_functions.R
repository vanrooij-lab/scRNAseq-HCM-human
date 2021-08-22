

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
    # Retrieving full ensembl names if necessary
    print(paste0('Retrieving ',gene_of_interest))
    if ((substr(x = gene_of_interest[1], 1, 3)=='ENS')) {
      gene_of_interest_fullname = gene_of_interest
    } else {
      gene_of_interest_fullname = shorthand_seurat_fullgenename(seuratObject, gene_of_interest)
    }
    
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



################################################################################



# Some extra comparative plots between bulk level data
# 
# load(file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_sorted.Rdata')) # core_regulons_sorted
# current_analysis$ROOIJonly_RID2l$dummyXX = round(runif(length(current_analysis$ROOIJonly_RID2l$annotation_paper_fct))*2)
#
# shorthand_custom_boxplot(seuratObject_list, gene_lists, seuratObjectNameToTake, group.by='annotation_paper_fct', topX=10, mylimits=.01)
shorthand_custom_boxplot = function(seuratObject_list, gene_lists, seuratObjectNameToTake, group.by='annotation_paper_fct', topX=10, mylimits=.01) {
    
    for (current_list_name in names(gene_lists)) {
        
        current_genes = gene_lists[[current_list_name]]
        current_genes = current_genes[1:(min(length(current_genes), topX))]
        
        # Retrieving full ensembl names if necessary
        print(paste0('Retrieving ',current_genes))
        if ((substr(x = current_genes[1], 1, 3)=='ENS')) {
          gene_of_interest_fullname = current_genes
        } else {
          gene_of_interest_fullname = shorthand_seurat_fullgenename(seuratObject_list[[seuratObjectNameToTake]], current_genes)
        }
        
        current_exprs = as.matrix(seuratObject_list[[seuratObjectNameToTake]]@assays$RNA@data[gene_of_interest_fullname,]) # removes dgMatrix class if there
        current_exprs = t(scale(t(current_exprs), center=T))
        rownames(current_exprs) = shorthand_cutname(rownames(current_exprs))
         
        df_ = data.frame(t(current_exprs))
        df_[[group.by]] = as.factor(current_analysis[[seuratObjectNameToTake]][[group.by, drop=T]])
        df = reshape2::melt(data = df_, id.vars = group.by, value.name = 'expression', variable.name='gene')
        
        agr_list=list(df[[group.by]], df$gene); names(agr_list) = c(group.by, 'gene')
        df_agr = aggregate(x = list(expression=df$expression), by = agr_list, FUN = mean)
        df_agr$sd = aggregate(x = list(sd=df$expression), by = agr_list, FUN = sd)$sd
        df_agr$ymin=df_agr$expression-df_agr$sd
        df_agr$ymax=df_agr$expression+df_agr$sd
        
        mylims=calc_limits(values = df$expression, percentile = mylimits)
        
        p=ggplot(df, mapping=aes_string(x='gene', y='expression', fill=group.by))+
            #geom_bar(data=df_agr, stat='identity',position='dodge')+
            geom_boxplot(position=position_dodge(width=1), outlier.size = .1, lwd=.25)+
            #geom_violin(alpha=.3, scale = 'width', width=.5, position = position_dodge(width = 1))+ # alpha=0.2, 
            geom_point(data=df_agr, stat='identity', position=position_dodge(width=1), shape=23, size=.5)+
            xlab('Gene')+ylab('Expression')+
            theme_bw()+give_better_textsize_plot(6)+
            theme(legend.position = 'none', legend.key.size = unit(3, "mm"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            #ylim(c(min(df$expression),mylims[2]))
            ylim(c(max(-3,mylims[1]),min(3,mylims[2])))
        p 
        
        # p=ggplot(df_agr, mapping=aes_string(x='gene', ymin='ymin', ymax='ymax', y='expression', fill=group.by))+
        #     geom_bar(stat='identity',position='dodge')+
        #     geom_errorbar(stat='identity',position='dodge')+
        #     #geom_violin(position='dodge', draw_quantiles = T)+ # alpha=0.2, 
        #     xlab('Gene')+ylab('Expression')+
        #     theme_bw()+give_better_textsize_plot(6)+
        #     theme(legend.position = 'none', legend.key.size = unit(3, "mm"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        #     #ylim(c(min(df$expression),mylims[2]))
        #     ylim(c(-1,min(3,mylims[2])))
        # p 
        
        ggsave(filename = paste0(base_dir, 'Rplots/', seuratObjectNameToTake, '_9_customBoxplotGenes_', current_list_name,'.pdf'), plot = p,
               height=42, width=42, units='mm')
        
        # If 1st plot, also print legend for reference
        if (current_list_name == names(gene_lists)[1]) {
        p=p+theme(legend.position='right')
        ggsave(filename = paste0(base_dir, 'Rplots/', seuratObjectNameToTake, '_9_customBoxplotGenes_LEGEND.pdf'), plot = p,
               height=42, width=60, units='mm')      
        }
    }
}

shorthand_custom_compositeplot = function(seuratObject_list, gene_lists, seuratObjectNameToTake, group.by='annotation_paper_fct', group.by2=NULL, zscore=T, myfontsize=6, custom_legend=NULL, mymargin=2, mypointsize=.25) {
    
  print('Note: this function assumes input names are long format.')
  
    p_violin_list=list(); p_bar_list = list(); p_bar_list_g2 = list()
    for (current_list_name in names(gene_lists)) {
        
        # Overall mean
        current_genes = gene_lists[[current_list_name]]
        
        # Retrieving full ensembl names if necessary
        print(paste0('Retrieving ',current_genes))
        if ((substr(x = current_genes[1], 1, 3)=='ENS')) {
          gene_of_interest_fullname = current_genes
        } else {
          gene_of_interest_fullname = shorthand_seurat_fullgenename(seuratObject_list[[seuratObjectNameToTake]], current_genes)
        }
        
        current_exprs = as.matrix(seuratObject_list[[seuratObjectNameToTake]]@assays$RNA@data[gene_of_interest_fullname,]) # removes dgMatrix class if there
        
        # Create composite expression
        if (zscore) { 
          scaled_exprs = t(scale(t(current_exprs),center = T, scale = T))
        } else { 
          scaled_exprs = t(scale(t(current_exprs),center = F, scale = T)) 
        }
        
        # current_expr=colSums(scaled_expr)
        current_expr=colMeans(scaled_exprs)
        
        # Plotting data frames
        df = data.frame(expression=current_expr)
        df[[group.by]] = as.factor(current_analysis[[seuratObjectNameToTake]][[group.by, drop=T]])
        if (!is.null(group.by2)) {df[[group.by2]] = as.factor(current_analysis[[seuratObjectNameToTake]][[group.by2, drop=T]])}
        
        # Calculate mean, take 2nd group into account if desired
        if (!is.null(group.by2)) { 
          agr_list=list(df[[group.by]], df[[group.by2]]); names(agr_list) = c(group.by, group.by2) 
        } else {
          agr_list=list(df[[group.by]]); names(agr_list) = c(group.by)
        }
        df_agr = aggregate(x = list(expression=df$expression), by = agr_list, FUN = mean)
        
        # Abbreviate dataset names
        # abbreviations=c('vRooij'='R', 'Hu'='H', 'Teichmann'='T')
        # if (group.by=='annotation_paper_str') {
        #   df[[group.by]]=factor(abbreviations[df[[group.by]]], levels=c('R','H','T'))
        #   df_agr[[group.by]]=factor(abbreviations[df_agr[[group.by]]], levels=c('R','H','T'))
        # }
        
        # Plot violins (simmilar plot as earlier)
        p=ggplot(df, aes_string(x=group.by, y='expression',fill=group.by))+
            geom_violin(scale='width')+
            ggtitle(current_list_name)+theme_bw()+
            theme(legend.position='none', plot.margin = margin(mymargin,mymargin,0,0,'mm'))+
            # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            xlab(element_blank())+give_better_textsize_plot(myfontsize)+
            ylab('Expression')
        p
        p_violin_list[[current_list_name]]=p
        ggsave(filename = paste0(base_dir, 'Rplots/', seuratObjectNameToTake, '_9_customSummaryComposite_', current_list_name,'_for_',group.by,'.pdf'), plot = p,
               height=26.5, width=26.5, units='mm')
        
        # Now plot bars
        p=ggplot(df_agr, mapping=aes_string(x=group.by, y='expression', fill=group.by))+
            geom_bar(data=df_agr, stat='identity')+
            xlab(element_blank())+ylab('Expression')+
            ggtitle(current_list_name)+
            theme_bw()+give_better_textsize_plot(myfontsize)+
            theme(legend.position = 'none', plot.margin = margin(mymargin,mymargin,0,0,'mm'))
            #theme(legend.position = 'none', legend.key.size = unit(3, "mm"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            #ylim(c(min(df$expression),mylims[2]))
            #ylim(c(max(-3,mylims[1]),min(3,mylims[2])))
        p_bar_list[[current_list_name]]=p
        p 
        
        ggsave(filename = paste0(base_dir, 'Rplots/', seuratObjectNameToTake, '_9_customSummaryMean_', current_list_name,'_for_',group.by,'.pdf'), plot = p,
               height=26.5, width=26.5, units='mm')
        
        # If 1st plot, also print legend for reference
        if (current_list_name == names(gene_lists)[1]) {
        p=p+theme(legend.position='right', legend.key.size = unit(3, "mm"))
        #if (group.by==group.by=='annotation_paper_str') {p=p+scale_fill_discrete(name = "Set")}
        if (!is.null(custom_legend)) {p=p+scale_fill_discrete(name = custom_legend)}
        ggsave(filename = paste0(base_dir, 'Rplots/', seuratObjectNameToTake, '_9_customSummaryMean_LEGEND_', current_list_name,'_for_',group.by,'.pdf'), plot = p,
               height=26.5, width=26.5, units='mm')
        
        }
        
        if (!is.null(group.by2)) {
        # Now plot bars
          p=ggplot(df_agr, mapping=aes_string(x=group.by, y='expression'))+
              geom_boxplot(aes_string(fill=group.by), outlier.shape = NA)+
              geom_jitter(size=mypointsize)+
              xlab(element_blank())+ylab('Expression')+
              ggtitle(current_list_name)+
              theme_bw()+give_better_textsize_plot(myfontsize)+
              theme(legend.position = 'none', plot.margin = margin(mymargin,mymargin,0,0,'mm'))
              #theme(legend.position = 'none', legend.key.size = unit(3, "mm"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
              #ylim(c(min(df$expression),mylims[2]))
              #ylim(c(max(-3,mylims[1]),min(3,mylims[2])))
          p_bar_list_g2[[current_list_name]]=p
          ggsave(filename = paste0(base_dir, 'Rplots/', seuratObjectNameToTake, '_9_customSummaryMean_SPLIT2_', current_list_name,'_for_',group.by,'.pdf'), plot = p,
               height=26.5, width=26.5, units='mm')
          
        }
        
    }
    
    return(list(p_violin_list=p_violin_list, p_bar_list=p_bar_list, p_bar_list_g2=p_bar_list_g2))
        
}
    
    




    


