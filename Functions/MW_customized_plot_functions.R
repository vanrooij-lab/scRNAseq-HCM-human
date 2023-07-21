

# For convenience, little fn to change text size ggplot in one go
give_better_textsize_plot <- function(TEXTSIZE, myFamily='Arial',suppressmessage=F){
  
  if (is.null(myFamily)) {
    theme(#legend.position="none",
          text = element_text(size=TEXTSIZE),
          axis.text = element_text(size=TEXTSIZE),
          plot.title = element_text(size=TEXTSIZE))
  } else {
    if(!suppressmessage){print('Set to Arial')}
    theme(#legend.position="none",
          text = element_text(size=TEXTSIZE, family=myFamily),
          axis.text = element_text(size=TEXTSIZE, family=myFamily),
          plot.title = element_text(size=TEXTSIZE, family=myFamily))
  }

}

theme_void_extramw_removeTickText = function() {
  
  theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
                     
}



# Violin plot
# Use lapply and wrap_plots to create a list of plots
shorthand_plotViolinBox_custom = function(myseuratobjectlist, analysis_name, cat_by='seurat_clusters', 
    cat_name = 'category', gene_of_interest=NULL,base_dir, percentile=.02, type='boxplot', 
    myfontsize=8, manual_expr=NULL,manual_expr_featname=NULL,custom_ylab=NULL, custom_title=NULL, custom_ylim=NULL, myfillcolor='#bd001f') {
  # TO DO: rename gene_of_interest to feature_of_interest
    
    # create df for plotting
    #plot_df = data.frame(count=numeric(), cat=character(), gene=character())
    #for (gene in gene_of_interest_fullname) {
        # retrieve current expression
        #current_expr = myseuratobjectlist[[analysis_name]]@assays$RNA@data[gene,]
        #expr_limits=calc_limits(current_expr, percentile = .03)

        #plot_df = rbind(plot_df, data.frame(count=current_expr, cat=myseuratobjectlist[[analysis_name]][[cat_by]], gene=gene))
    #}
    
    # retrieve expression
    if (is.null(manual_expr)) {
      
      # get full gene names
      gene_of_interest_fullname = rownames(myseuratobjectlist[[analysis_name]])[grepl(paste0(paste0(':',gene_of_interest,'$'),collapse='|'),rownames(myseuratobjectlist[[analysis_name]]))]
      # retrieve expression
      print(paste0('Retrieved gene ',gene_of_interest_fullname))
      current_expr = myseuratobjectlist[[analysis_name]]@assays$RNA@data[gene_of_interest_fullname,] 
      
    }else{ print('Using manually supplied expression'); current_expr=manual_expr; gene_of_interest_fullname = 'manual'; gene_of_interest = manual_expr_featname }
    
    
    # create dataframe to plot
    plot_df = data.frame(count=current_expr, cat=myseuratobjectlist[[analysis_name]][[cat_by]], gene=gene_of_interest_fullname)
    colnames(plot_df) = c('count',cat_name)
    
    # decide on limits
    # per group, determine where 2% percentile lies to set that as maximum to the plot
    my_max_limit=max(sapply(unique(plot_df[[cat_name]]), function(x) { calc_limits(current_expr[plot_df[[cat_name]]==x], percentile = .02)[2] } ))

    # start plotting w/ desired options
    p=ggplot(plot_df, aes_string(x=cat_name, y='count'))
    
    if (type=='violin') {p=p+geom_violin(scale='width', fill=myfillcolor)}
    if (type=='boxplot') {p=p+geom_boxplot(fill=myfillcolor)}
    
    p=p+
        #geom_bar(aes_string(y='count', fill=as.factor(cat_name)),stat='summary',fun='mean',fill='grey',color='black')+
        #geom_boxplot(alpha=0)+
        ggtitle(if(is.null(custom_title)){gene_of_interest}else{custom_title})+theme_bw()+
        theme(legend.position='none')+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        xlab(if(is.null(cat_name=='category')){element_blank()}else{cat_name})+
        ylab(if(is.null(custom_ylab)){'UMI count'}else{custom_ylab})+give_better_textsize_plot(myfontsize)
    # p
    
    if (!is.null(custom_ylim)) {p=p+ylim(custom_ylim)}
    
    # save 1
    ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',analysis_name,'_6_',type,'_by-', cat_name, '_',gene_of_interest,'.pdf'), width = 4, height= 6, units='cm', device=cairo_pdf)
    ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',analysis_name,'_6_',type,'_by-', cat_name, '_',gene_of_interest,'_172mm.pdf'), width = 172/3-4, height= 172/3-4, units='mm', device=cairo_pdf)
    
    # add ylim and save 2
    p=p+ylim(c(0,my_max_limit))
    ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',analysis_name,'_6_',type,'_ylim98_by-', cat_name, '_',gene_of_interest,'.pdf'), width = 4, height= 6, units='cm', device=cairo_pdf)
    ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',analysis_name,'_6_',type,'_ylim98_by-', cat_name, '_',gene_of_interest,'_172mm.pdf'), width = 172/3-4, height= 172/3-4, units='mm', device=cairo_pdf)
    
    # return plot
    return(p)
    
}

shorthand_plotFreqPoly_custom = function(myseuratobjectlist, analysis_name, cat_by='seurat_clusters', 
    cat_name = 'category', gene_of_interest=NULL,base_dir, percentile=.02, 
    myfontsize=8, manual_expr=NULL,manual_expr_featname=NULL,custom_ylab=NULL, custom_title=NULL, custom_ylim=NULL, custom_colors = c('#bd0020','#9d9d9c','#575756')) {
  # TO DO: rename gene_of_interest to feature_of_interest
    
    # create df for plotting
    #plot_df = data.frame(count=numeric(), cat=character(), gene=character())
    #for (gene in gene_of_interest_fullname) {
        # retrieve current expression
        #current_expr = myseuratobjectlist[[analysis_name]]@assays$RNA@data[gene,]
        #expr_limits=calc_limits(current_expr, percentile = .03)

        #plot_df = rbind(plot_df, data.frame(count=current_expr, cat=myseuratobjectlist[[analysis_name]][[cat_by]], gene=gene))
    #}
    
    # retrieve expression
    if (is.null(manual_expr)) {
      
      # get full gene names
      gene_of_interest_fullname = rownames(myseuratobjectlist[[analysis_name]])[grepl(paste0(paste0(':',gene_of_interest,'$'),collapse='|'),rownames(myseuratobjectlist[[analysis_name]]))]
      # retrieve expression
      print(paste0('Retrieved gene ',gene_of_interest_fullname))
      current_expr = myseuratobjectlist[[analysis_name]]@assays$RNA@data[gene_of_interest_fullname,] 
      
    }else{ print('Using manually supplied expression'); current_expr=manual_expr; gene_of_interest_fullname = 'manual'; gene_of_interest = manual_expr_featname }
    
    
    # create dataframe to plot
    plot_df = data.frame(count=current_expr, cat=myseuratobjectlist[[analysis_name]][[cat_by]], gene=gene_of_interest_fullname)
    colnames(plot_df) = c('count',cat_name)
    
    # decide on limits
    # per group, determine where 2% percentile lies to set that as maximum to the plot
    my_max_limit=max(sapply(unique(plot_df[[cat_name]]), function(x) { calc_limits(current_expr[plot_df[[cat_name]]==x], percentile = .02)[2] } ))

    # start plotting w/ desired options
    p=ggplot(plot_df, aes(x=count, after_stat(density)))+
        geom_freqpoly(aes_string(color=cat_name), bins=50, )+
        #geom_bar(aes_string(y='count', fill=as.factor(cat_name)),stat='summary',fun='mean',fill='grey',color='black')+
        #geom_boxplot(alpha=0)+
        ggtitle(if(is.null(custom_title)){gene_of_interest}else{custom_title})+theme_bw()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        xlab('UMI count')+ylab('PDF')+
        give_better_textsize_plot(myfontsize)+scale_color_manual(values = custom_colors)
    # p
    
    # # save 1
    # ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',analysis_name,'_6_',type,'_by-', cat_name, '_',gene_of_interest,'.pdf'), width = 4, height= 6, units='cm', device=cairo_pdf)
    # ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',analysis_name,'_6_',type,'_by-', cat_name, '_',gene_of_interest,'_172mm.pdf'), width = 172/3-4, height= 172/3-4, units='mm', device=cairo_pdf)
    # 
    # # add ylim and save 2
    # p=p+ylim(c(0,my_max_limit))
    # ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',analysis_name,'_6_',type,'_ylim98_by-', cat_name, '_',gene_of_interest,'.pdf'), width = 4, height= 6, units='cm', device=cairo_pdf)
    # ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',analysis_name,'_6_',type,'_ylim98_by-', cat_name, '_',gene_of_interest,'_172mm.pdf'), width = 172/3-4, height= 172/3-4, units='mm', device=cairo_pdf)
    
    # return plot
    return(p)
    
}

# myseuratobjectlist=current_analysis
# analysis_name='ALL.SP_btypSel_RID2l_clExtended'
# gene_of_interest = 'SRF'
# cat_by='annotation_paper_beatified'
# cat_name = 'Dataset'
# myfillcolor='#bd001f'
shorthand_plotPercentageCellsHighForGene = function(myseuratobjectlist, analysis_name, gene_of_interest, cat_by, cat_name, myfillcolor='#bd001f', PERCENTILE=.5, myfontsize=8) {
    
    # retrieve expression
    # get full gene names
    gene_of_interest_fullname = rownames(myseuratobjectlist[[analysis_name]])[grepl(paste0(paste0(':',gene_of_interest,'$'),collapse='|'),rownames(myseuratobjectlist[[analysis_name]]))]
    # retrieve expression
    print(paste0('Retrieved gene ',gene_of_interest_fullname))
    current_expr = myseuratobjectlist[[analysis_name]]@assays$RNA@data[gene_of_interest_fullname,] 
  
    # create dataframe to plot
    summary_df = data.frame(count=current_expr, cat=myseuratobjectlist[[analysis_name]][[cat_by]], gene=gene_of_interest_fullname)
    colnames(summary_df) = c('count',cat_name)
    
    cutoff_value = calc_limits(summary_df$count, PERCENTILE)[1] # .5 = median
    
    summary_df$detected = summary_df$count>cutoff_value
    
    
    cat_list = list(item=summary_df[[cat_name]])
    names(cat_list) = cat_name
    plot_df = aggregate(list(perc_above=summary_df$detected), cat_list, FUN=mean)
    
    p1=ggplot(plot_df, aes_string(x=cat_name, y='perc_above')) + 
        geom_bar(stat='identity', fill=myfillcolor)+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ylab('Percentage cells high')+give_better_textsize_plot(myfontsize)+
        ggtitle(paste0(gene_of_interest,'>',round(cutoff_value,2)))
    
    return(p1)

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
# textsize=8; pointsize=1; mypercentile=0.03; custom_title=NULL;mymargin=0.1;zscore=F;myFamily='Arial'; add_box=F
shorthand_seurat_custom_expr = function(seuratObject, gene_of_interest, textsize=8, pointsize=1, mypercentile=0.03, 
                                        custom_title=NULL,mymargin=0.1,zscore=F,myFamily='Arial', add_box=F, parse_annot=T,myColors=NULL) {
    
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
      # gene_of_interest_fullname = shorthand_seurat_fullgenename(seuratObject, gene_of_interest)
      gene_of_interest_fullname = shorthand_seurat_fullgenename_faster(seuratObject, gene_of_interest, return_NA = T)
      if (any(is.na(gene_of_interest_fullname))) {warning('Not all genes found uniquely, filtering those out ..')}
      gene_of_interest_fullname = gene_of_interest_fullname[!is.na(gene_of_interest_fullname)]
    }
    
    if (length(gene_of_interest_fullname)==0) { return(NA) }
    if (is.null(myColors)) {myColors = rainbow_colors}
    
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
    # expr_limits
    
    if (add_box) {boxpad=.1} else {boxpad=0}
    p=ggplot(data.frame(UMAP_1=seuratObject@reductions$umap@cell.embeddings[,1],
                      UMAP_2=seuratObject@reductions$umap@cell.embeddings[,2],
                      expr=current_expr),
            mapping = aes(x=UMAP_1, y=UMAP_2, color=expr))+
        geom_point(size = pointsize, stroke = 0, shape = 16)+
        scale_color_gradientn(colours=myColors, limits=c(0,expr_limits[2]), oob=squish)+
        #+ggtitle(gene_of_interest)+
        annotate("text", -Inf, Inf, label = mytitle, hjust = 0-boxpad, vjust = 1+textsize*boxpad, size=textsize / .pt, family = myFamily, parse=T)+
        theme_void()+
        #give_better_textsize_plot(textsize)+
        theme(legend.position = 'none', plot.margin = margin(mymargin,mymargin,0,0,'mm'))
    
    # p
    
    if (add_box) {
      p=p+theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))
    }
    
    return(p)
        
        
}
# shorthand to get full names
shorthand_seurat_fullgenename = function(seuratObject, gene_names) {
    
    # gene_of_interest_fullname = rownames(seuratObject)[grepl(paste0(paste0(':',gene_names,'$'),collapse = '|'),rownames(seuratObject))]
    gene_of_interest_fullname = 
        sapply(gene_names, function(current_name) {
          hits = grepl(paste0(':',current_name,'$'),rownames(seuratObject))
          cnt = sum(hits)
          if (cnt>1) { stop('Multiple genes could be linked to supplied short name(s)..')}
          if (cnt==0) { stop('Supplied gene name not found..') }
          return(rownames(seuratObject)[hits])
          })
    return(gene_of_interest_fullname)
    
}

shorthand_seurat_fullgenename_faster = function(seuratObject, gene_names, return_NA=F) {
    
  all_short_names = shorthand_cutname(rownames(seuratObject))
  unique_selection = !(all_short_names %in% all_short_names[duplicated(all_short_names)])
  
  unique_names      = rownames(seuratObject)[unique_selection]
  unique_shortnames = all_short_names[unique_selection]
  
  lookupframe=unique_names
  names(lookupframe) = unique_shortnames
  
  if (!return_NA) {
    if (!all(gene_names %in% names(lookupframe))) {
      print(paste0('Supplied list: ',gene_names)); 
      print(paste0('Genes not found: ',toString(gene_names[!(gene_names %in% names(lookupframe))])))
      stop('Couldn\'t find all genes in lookuptable.')}
  }
  
  # Now return the full names
  return(lookupframe[gene_names])        

    
}

################################################################################

shorthand_seurat_custom_umap_groupby = function(seuratObject, textsize=8, pointsize=1, my.group.by = NULL,
                                        custom_title=NULL,mymargin=0.1,myFamily='Arial', add_box=F, parse_annot=T,myColors=NULL,
                                        central_marker_size=2, annotate_clusters=T) {
    
  if(is.null(my.group.by)) {stop('Set my.group.by')}
  
    if (is.null(custom_title)) {
      mytitle = my.group.by
    } else { mytitle = custom_title }
  
    if (is.null(myColors)) {myColors = pals::alphabet2(); names(myColors)=NULL} 
    
    
    df_plot = data.frame(
      umap_1 = seuratObject@reductions$umap@cell.embeddings[,1],
      umap_2 = seuratObject@reductions$umap@cell.embeddings[,2])
    df_plot[[my.group.by]] = as.factor(seuratObject[[my.group.by]][,1])
  
    mylist=list(df_plot[[my.group.by]])
    names(mylist)=my.group.by
    df_annotation = aggregate(x = df_plot[,c('umap_1','umap_2')], by=mylist, FUN=mean)
    
    if (add_box) {boxpad=.1} else {boxpad=0}
    p=ggplot(df_plot,
            mapping = aes_string(x='umap_1', y='umap_2', color=my.group.by))+
        geom_point(size = pointsize, stroke = 0, shape = 16)+
        scale_color_manual(values =myColors)+
        scale_fill_manual(values =myColors)+
        #annotate("text", -Inf, Inf, label = mytitle, hjust = 0-boxpad, vjust = 1+textsize*boxpad, size=textsize / .pt, family = myFamily, parse=T)+
        theme_void()+
        #give_better_textsize_plot(textsize)+
        theme(legend.position = 'none', plot.margin = margin(mymargin,mymargin,0,0,'mm'))
    
    # p
    
    # add markers to indicate center of groups
    if (central_marker_size>0 & annotate_clusters) {
      p=p+geom_point(data=df_annotation, aes_string(fill=my.group.by), color='black', size=3, shape = 22) # color='grey50',
      # p=p+geom_point(data=df_annotation, color='grey50', size=3, shape = 15) 
    }
    
    # code not tested
    if (add_box) {
      p=p+theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))
    }
    
    # Add text annotation last to get order right
    if (annotate_clusters) {
      p=p+geom_text_repel(data=df_annotation, aes_string(label=my.group.by),color='black', 
                        size = textsize/.pt, max.overlaps = Inf, min.segment.length = 0)
    }  
    
    return(p)
        
        
}


################################################################################



# Some extra comparative plots between bulk level data
# 
# load(file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_sorted.Rdata')) # core_regulons_sorted
# current_analysis$ROOIJonly_RID2l$dummyXX = round(runif(length(current_analysis$ROOIJonly_RID2l$annotation_paper_fct))*2)
#
# shorthand_custom_boxplot(seuratObject_list, gene_lists, seuratObjectNameToTake, group.by='annotation_paper_fct', topX=10, mylimits=.01)
shorthand_custom_boxplot = function(seuratObject_list, gene_lists, seuratObjectNameToTake, group.by='annotation_paper_fct', topX=10, mylimits=.01, show=F, SUBDIR='') {
    
    for (current_list_name in names(gene_lists)) {
      
      # current_list_name = names(gene_lists)[1]
        
        current_genes = gene_lists[[current_list_name]]
        current_genes = current_genes[1:(min(length(current_genes), topX))]
        
        # Retrieving full ensembl names if necessary
        print(paste0('Retrieving ',current_genes))
        if ((substr(x = current_genes[1], 1, 3)=='ENS')) {
          gene_of_interest_fullname = current_genes
        } else {
          gene_of_interest_fullname = shorthand_seurat_fullgenename_faster(seuratObject_list[[seuratObjectNameToTake]], current_genes, return_NA = T)
          if (any(is.na(gene_of_interest_fullname))) {warning('Not all genes found uniquely, filtering those out ..')}
          gene_of_interest_fullname = gene_of_interest_fullname[!is.na(gene_of_interest_fullname)]
        }
        
        current_exprs = as.matrix(seuratObject_list[[seuratObjectNameToTake]]@assays$RNA@data[gene_of_interest_fullname,]) # removes dgMatrix class if there
        current_exprs = t(scale(t(current_exprs), center=T))
        rownames(current_exprs) = shorthand_cutname(rownames(current_exprs))
         
        df_ = data.frame(t(current_exprs))
        df_[[group.by]] = as.factor(seuratObject_list[[seuratObjectNameToTake]][[group.by, drop=T]])
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
            labs(fill = element_blank())
        p 
        
        # Save
        ggsave(filename = paste0(base_dir, 'Rplots/', SUBDIR, seuratObjectNameToTake, '_9_customBoxplotGenes_', current_list_name,'.pdf'), plot = p,
               height=42, width=42, units='mm', device=cairo_pdf)
        
        # Save with extra limits
        p=p+ylim(c(max(-3,mylims[1]),min(3,mylims[2])))
        ggsave(filename = paste0(base_dir, 'Rplots/', SUBDIR, seuratObjectNameToTake, '_9_customBoxplotGenes_', current_list_name,'.pdf'), plot = p,
               height=42, width=42, units='mm', device=cairo_pdf)
        
        # If 1st plot, also print legend for reference
        if (current_list_name == names(gene_lists)[1]) {
        p=p+theme(legend.position='right')
        ggsave(filename = paste0(base_dir, 'Rplots/', SUBDIR, seuratObjectNameToTake, '_9_customBoxplotGenes_',current_list_name,'_LEGEND.pdf'), plot = p,
               height=42, width=4.2*min(length(current_genes), topX)+18, units='mm', device=cairo_pdf)      
        }
        
        p=ggplot(df_agr, mapping=aes_string(x='gene', y='expression', fill=group.by))+ # , ymin='ymin', ymax='ymax'
            geom_violin(data=df, position='dodge', alpha=0, lwd=.1, scale='width')+ # alpha=0.2,
            geom_bar(stat='identity',position='dodge', alpha=.9)+
            #geom_errorbar(stat='identity',position='dodge')+
            xlab('Gene')+ylab('Expression')+
            theme_bw()+give_better_textsize_plot(6)+
            theme(legend.position = 'none', legend.key.size = unit(3, "mm"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            #ylim(c(min(df$expression),mylims[2]))
            ylim(c(min(df_agr$expression)-.1,max(df_agr$expression)+.5))
        # p
        # Save
        ggsave(filename = paste0(base_dir, 'Rplots/', SUBDIR, seuratObjectNameToTake, '_9_customBarplotGenes_', current_list_name,'.pdf'), plot = p,
               height=42, width=42, units='mm', device=cairo_pdf)
        
        # 3rd style
        p=ggplot(df, mapping=aes_string(x=group.by, y='expression', fill=group.by))+
            #geom_bar(stat='identity',position='dodge')+
            geom_boxplot(outlier.size = .1, lwd=.25, alpha=.1)+
            geom_jitter(size=.25, aes_string(color=group.by),position=position_dodge(width=1))+
            xlab(element_blank())+ylab('Expression (z-score)')+
            theme_bw()+give_better_textsize_plot(6)+
            theme(legend.position = 'none', legend.key.size = unit(3, "mm"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  panel.spacing = unit(c(0),'mm'))+
            labs(fill = element_blank())+
            facet_grid(cols=vars(gene))
        
        ggsave(filename = paste0(base_dir, 'Rplots/', SUBDIR, seuratObjectNameToTake, '_9_customBoxplotGenes_v3_', current_list_name,'.pdf'), plot = p,
               height=172/6, width=172/3*2/10*min(length(current_genes), topX), units='mm', device=cairo_pdf)
        
        # Now also save the data for stat. analysis
        # save(list = 'df_agr', file=paste0(base_dir, 'Rplots/', seuratObjectNameToTake, '_9_customBoxplotGenes-PerPat-data_', current_list_name,'.Rdata'))
        openxlsx::write.xlsx(x = df, file=paste0(base_dir, 'Rplots/', SUBDIR, seuratObjectNameToTake, '_9_customBoxplotGenes-data_', current_list_name,'.xlsx'), overwrite = T)
        
    }
  
    if (show) {return(p)}
  
}

##########

shorthand_custom_boxplot_perpatient = function(seuratObject_list, gene_lists, seuratObjectNameToTake, group.by='annotation_paper_fct', topX=10, mylimits=.01, aggr.by='annotation_patient_fct',mySize=NULL, myFontSize=6) {
    
    for (current_list_name in names(gene_lists)) {
        
        current_genes = gene_lists[[current_list_name]]
        current_genes = current_genes[1:(min(length(current_genes), topX))]
        
        # Retrieving full ensembl names if necessary
        print(paste0('Retrieving ',current_genes))
        if ((substr(x = current_genes[1], 1, 3)=='ENS')) {
          gene_of_interest_fullname = current_genes
        } else {
          gene_of_interest_fullname = shorthand_seurat_fullgenename_faster(seuratObject_list[[seuratObjectNameToTake]], current_genes, return_NA = T)
          if (any(is.na(gene_of_interest_fullname))) {warning('Not all genes found uniquely, filtering those out ..')}
          gene_of_interest_fullname = gene_of_interest_fullname[!is.na(gene_of_interest_fullname)]
        }
        
        current_exprs = as.matrix(seuratObject_list[[seuratObjectNameToTake]]@assays$RNA@data[gene_of_interest_fullname,]) # removes dgMatrix class if there
        current_exprs = t(scale(t(current_exprs), center=T))
        rownames(current_exprs) = shorthand_cutname(rownames(current_exprs))
         
        df_ = data.frame(t(current_exprs))
        df_[[group.by]]   = as.factor(seuratObject_list[[seuratObjectNameToTake]][[group.by, drop=T]])
        df_[[aggr.by]]  = as.factor(seuratObject_list[[seuratObjectNameToTake]][[aggr.by, drop=T]])
        
        df = reshape2::melt(data = df_, id.vars = list(group.by,aggr.by), value.name = 'expression', variable.name='gene')
        
        agr_list=list(df[[group.by]], df$gene, df[[aggr.by]]); names(agr_list) = c(group.by, 'gene', aggr.by)
        df_agr = aggregate(x = list(expression=df$expression), by = agr_list, FUN = mean)
        #df_agr$sd = aggregate(x = list(sd=df$expression), by = agr_list, FUN = sd)$sd
        #df_agr$ymin=df_agr$expression-df_agr$sd
        #df_agr$ymax=df_agr$expression+df_agr$sd
        
        # mylims=calc_limits(values = df$expression, percentile = mylimits)
        
        p=ggplot(df_agr, mapping=aes_string(x='gene', y='expression', fill=group.by))+
            #geom_bar(data=df_agr, stat='identity',position='dodge')+
            geom_jitter(size=.1, aes_string(color=group.by),position=position_dodge(width=1))+
            geom_boxplot(position=position_dodge(width=1), outlier.size = .1, lwd=.25, alpha=.5)+
            # geom_point(data=df_agr, stat='identity', position=position_dodge(width=1), shape=23, size=.5)+
            xlab('Gene')+ylab('Expression')+
            theme_bw()+give_better_textsize_plot(myFontSize)+
            theme(legend.position = 'none', legend.key.size = unit(3, "mm"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  panel.spacing = unit(c(0),'mm'))+
            #ylim(c(min(df$expression),mylims[2]))
            labs(fill = element_blank())
        #p 
        
        # Save
        ggsave(filename = paste0(base_dir, 'Rplots/', seuratObjectNameToTake, '_9_customBoxplotGenes-PerPat_', current_list_name,'.pdf'), plot = p,
               height=172/6, width=172/3*2/10*length(current_genes), units='mm', device=cairo_pdf)
        
        p=ggplot(df_agr, mapping=aes_string(x=group.by, y='expression', fill=group.by))+
            #geom_bar(stat='identity',position='dodge')+
            geom_boxplot(outlier.size = .1, lwd=.25, alpha=.1)+
            geom_jitter(size=.25, aes_string(color=group.by),position=position_dodge(width=1))+
            xlab(element_blank())+ylab('Expression (z-score)')+
            theme_bw()+give_better_textsize_plot(myFontSize)+
            theme(legend.position = 'none', legend.key.size = unit(3, "mm"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  panel.spacing = unit(c(0),'mm'))+
            labs(fill = element_blank())+
            facet_grid(cols=vars(gene))
        
        ggsave(filename = paste0(base_dir, 'Rplots/', seuratObjectNameToTake, '_9_customBoxplotGenes-PerPat-v2_', current_list_name,'.pdf'), plot = p,
               height=172/6, width=172/3*2/10*length(current_genes), units='mm', device=cairo_pdf)
        if (!is.null(mySize)) {
          ggsave(filename = paste0(base_dir, 'Rplots/', seuratObjectNameToTake, '_9_customBoxplotGenes-PerPat-v2_', current_list_name,'_customSize.pdf'), plot = p,
                 height=mySize[2], width=mySize[1], units='mm', device=cairo_pdf)
        }
          
        # Now also save the data for stat. analysis
        # save(list = 'df_agr', file=paste0(base_dir, 'Rplots/', seuratObjectNameToTake, '_9_customBoxplotGenes-PerPat-data_', current_list_name,'.Rdata'))
        openxlsx::write.xlsx(x = df_agr, file=paste0(base_dir, 'Rplots/', seuratObjectNameToTake, '_9_customBoxplotGenes-PerPat-data_', current_list_name,'.xlsx'), overwrite = T)
    }
}

##########

# zscore=T; myfontsize=6; custom_legend=NULL; mymargin=2; mypointsize=.25
shorthand_custom_compositeplot = function(seuratObject_list, gene_lists, seuratObjectNameToTake, group.by='annotation_paper_fct', group.by2=NULL, zscore=T, myfontsize=6, custom_legend=NULL, mymargin=2, mypointsize=.25, SUBDIR="") {
    
  print('Note: this function assumes input names are long format.')
  
    p_violin_list=list(); p_bar_list = list(); p_bar_list_g2 = list(); p_boxplot_list = list()
    df_agr_list = list()
    for (current_list_name in names(gene_lists)) {
        
        # Overall mean
        current_genes = gene_lists[[current_list_name]]
        
        # Retrieving full ensembl names if necessary
        print(paste0('Retrieving ',current_genes))
        if ((substr(x = current_genes[1], 1, 3)=='ENS')) {
          gene_of_interest_fullname = current_genes
        } else {
          # gene_of_interest_fullname = shorthand_seurat_fullgenename(seuratObject_list[[seuratObjectNameToTake]], current_genes)
          gene_of_interest_fullname = shorthand_seurat_fullgenename_faster(seuratObject_list[[seuratObjectNameToTake]], current_genes, return_NA = T)
          if (any(is.na(gene_of_interest_fullname))) {warning('Not all genes found uniquely, filtering those out ..')}
          gene_of_interest_fullname = gene_of_interest_fullname[!is.na(gene_of_interest_fullname)]
        }
        
        current_exprs = as.matrix(seuratObject_list[[seuratObjectNameToTake]]@assays$RNA@data[gene_of_interest_fullname,]) # removes dgMatrix class if there
        
        # Create composite expression
        zero_rows=apply(current_exprs, 1, function(x) {all(x==0)})
        if (any(zero_rows)) {
          warnings('Rows with all zero values were found and removed.')  
          print('Rows with all zero values were found and removed.')  
          current_exprs=current_exprs[!zero_rows,]
        }
        if (zscore) { 
          scaled_exprs = t(scale(t(current_exprs),center = T, scale = T))
        } else { 
          scaled_exprs = t(scale(t(current_exprs),center = F, scale = T)) 
        }
        
        # current_expr=colSums(scaled_expr)
        current_expr=colMeans(scaled_exprs)
        
        # Plotting data frames
        df = data.frame(expression=current_expr)
        df[[group.by]] = as.factor(seuratObject_list[[seuratObjectNameToTake]][[group.by, drop=T]])
        if (!is.null(group.by2)) {df[[group.by2]] = as.factor(seuratObject_list[[seuratObjectNameToTake]][[group.by2, drop=T]])}
        
        # Calculate mean, take 2nd group into account if desired
        if (!is.null(group.by2)) { 
          agr_list=list(df[[group.by]], df[[group.by2]]); names(agr_list) = c(group.by, group.by2) 
        } else {
          agr_list=list(df[[group.by]]); names(agr_list) = c(group.by)
        }
        df_agr = aggregate(x = list(expression=df$expression), by = agr_list, FUN = mean)
        df_agr$set_name = current_list_name
        df_agr_list[[current_list_name]] = df_agr
        
        # Now export excel file, to do statistics (Note: later will also export joined one)
        openxlsx::write.xlsx(x = df_agr, file = paste0(base_dir, 'Rplots/', SUBDIR, seuratObjectNameToTake, '_9_customSummaryMean_', current_list_name,'_for_',group.by,'_splt2',group.by2,'-data.xlsx'), overwrite = T)

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
        # p
        p_violin_list[[current_list_name]]=p
        ggsave(filename = paste0(base_dir, 'Rplots/', SUBDIR, seuratObjectNameToTake, '_9_customSummaryComposite_', current_list_name,'_for_',group.by,'.pdf'), plot = p,
               height=26.5, width=26.5, units='mm', device=cairo_pdf)
        
        # Plot boxplot (simmilar plot as earlier)
        p=ggplot(df, aes_string(x=group.by, y='expression',fill=group.by))+
            geom_boxplot()+
            ggtitle(current_list_name)+theme_bw()+
            theme(legend.position='none', plot.margin = margin(mymargin,mymargin,0,0,'mm'))+
            # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            xlab(element_blank())+give_better_textsize_plot(myfontsize)+
            ylab('Expression')
        if (nrow(df)<100) {p=p+geom_jitter()}
        # p
        p_boxplot_list[[current_list_name]]=p
        ggsave(filename = paste0(base_dir, 'Rplots/', SUBDIR, seuratObjectNameToTake, '_9_customSummaryCompositeBoxPlot_', current_list_name,'_for_',group.by,'.pdf'), plot = p,
               height=26.5, width=26.5, units='mm', device=cairo_pdf)
        
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
        # p 
        
        ggsave(filename = paste0(base_dir, 'Rplots/', SUBDIR, seuratObjectNameToTake, '_9_customSummaryMean_', current_list_name,'_for_',group.by,'.pdf'), plot = p,
               height=26.5, width=26.5, units='mm', device=cairo_pdf)
        
        # If 1st plot, also print legend for reference
        if (current_list_name == names(gene_lists)[1]) {
        p=p+theme(legend.position='right', legend.key.size = unit(3, "mm"))
        #if (group.by==group.by=='annotation_paper_str') {p=p+scale_fill_discrete(name = "Set")}
        if (!is.null(custom_legend)) {p=p+scale_fill_discrete(name = custom_legend)}
        ggsave(filename = paste0(base_dir, 'Rplots/', SUBDIR, seuratObjectNameToTake, '_9_customSummaryMean_LEGEND_', current_list_name,'_for_',group.by,'.pdf'), plot = p,
               height=26.5, width=26.5, units='mm', device=cairo_pdf)
        
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
          ggsave(filename = paste0(base_dir, 'Rplots/', SUBDIR, seuratObjectNameToTake, '_9_customSummaryMean_SPLIT2_', current_list_name,'_for_',group.by,'.pdf'), plot = p,
               height=26.5, width=26.5, units='mm', device=cairo_pdf)
          
        }
        
    }
    
    # Now export all composite data
    df_agr_combined = Reduce(f = rbind, df_agr_list)
    openxlsx::write.xlsx(x = df_agr_combined, file = paste0(base_dir, 'Rplots/', SUBDIR, seuratObjectNameToTake, '_9_customSummaryMean_ALL-', current_list_name,'-etc_for_',group.by,'_splt2',group.by2,'-data.xlsx'), overwrite = T)
    
    return(list(p_violin_list=p_violin_list, p_bar_list=p_bar_list, p_bar_list_g2=p_bar_list_g2, p_boxplot_list=p_boxplot_list))
        
}
    
    

##########

# Another summary plot tool
# Produces a heatmap that shows average expression per <CLASS1> and <CLASS2>
# For multiple genes of interest
# Also allows for annotation of the genes of interest
# cols = genes of interest
# rows = conditions1 x conditions2
shorthand_heatmap_feature_aggr = function(current_analysis, analysis, feature_list, 
                                          feature_list_annotation, group.by1, group.by2, fontsize=10,
                                          manual_zlims = NULL, only_return_max=F, savepath=NULL) {
    
    # It seems Seurat does something weird when requesting using the
    # current_analysis[[analysis]][["clusters_custom"]]
    # syntax (returns df)
    # So this is workaround
    df_groups_meta = data.frame(group.by1=current_analysis[[analysis]][[group.by1]],group.by2=current_analysis[[analysis]][[group.by2]])
    
    # Convert to df
    df_features_meta =
        data.frame( feature_annot=feature_list_annotation, feature=feature_list )
    
    all_features = shorthand_seurat_fullgenename_faster(seuratObject = current_analysis[[analysis]], gene_names = df_features_meta$feature, return_NA = T)
    
    # Create summary expression matrix
    print(paste0('Creating summary expression matrix -- ', date()))
    df_expr_list =
        lapply(all_features, function(g) {
            print(paste0('Working on gene ',g,' -- ', date()))
            if (is.na(g)) {
                # Create 
                dummy = df_groups_meta
                dummy$expr = NA
                expr1_aggr = aggregate(dummy$expr, df_groups_meta, mean)
                expr1_aggr$Z = NA
                expr1_aggr$f = NA
                print(paste0('Gene ',g,' done -- ', date()))
                return(expr1_aggr)
            } else {
                # collect data
                expr1 = as.vector(current_analysis[[analysis]]@assays$RNA@data[g, ])
                # mean(current_analysis[[analysis]]@assays$RNA@counts[g, ])
                
                # Calculate fraction of cells in which it is expressed
                fraction_aggr = aggregate(1*(expr1>0), df_groups_meta, mean)
                
                # aggregate data
                expr1_aggr = aggregate(expr1, df_groups_meta, mean)
                
                # determine Z-score of aggr
                expr1_aggr$Z = scale(expr1_aggr$x)
                expr1_aggr$f = fraction_aggr$x
                print(paste0('Gene ',g,' done -- ', date()))
                return(expr1_aggr)
            }
        })
    # df_expr = Reduce(rbind, df_expr_list)
    print(paste0('Done -- ', date()))
    
    # Create normalized matrix (Z), and non-normalized matrix (x)
    matrix_expr_Z = sapply(df_expr_list, function(df){df$Z})
    matrix_expr_x = sapply(df_expr_list, function(df){df$x}) 
    matrix_expr_f = sapply(df_expr_list, function(df){df$f})
    
    colnames(matrix_expr_Z) = paste0(df_features_meta$feature, '_', df_features_meta$feature_annot)
    rownames(matrix_expr_Z) = paste0(df_expr_list[[1]][[group.by1]], '_', df_expr_list[[1]][[group.by2]])
    colnames(matrix_expr_x) = paste0(df_features_meta$feature, '_', df_features_meta$feature_annot)
    rownames(matrix_expr_x) = paste0(df_expr_list[[1]][[group.by1]], '_', df_expr_list[[1]][[group.by2]])
    colnames(matrix_expr_f) = paste0(df_features_meta$feature, '_', df_features_meta$feature_annot)
    rownames(matrix_expr_f) = paste0(df_expr_list[[1]][[group.by1]], '_', df_expr_list[[1]][[group.by2]])
    
    annotation_cols = data.frame(annot = df_features_meta$feature_annot)
    rownames(annotation_cols) = paste0(df_features_meta$feature, '_', df_features_meta$feature_annot)
    
    annotation_rows = data.frame(g1=df_expr_list[[1]][[group.by1]], g2=df_expr_list[[1]][[group.by2]])
    rownames(annotation_rows) = paste0(df_expr_list[[1]][[group.by1]], '_', df_expr_list[[1]][[group.by2]])
    names(annotation_rows) = c(group.by1, group.by2)
    
    # collapse group 2
    matrix_expr_Z_g2collapse =
      apply(matrix_expr_Z, 2, function(x) {df_agr=aggregate(x, list(annotation_rows[[group.by2]]), mean); df_agr$x})
    matrix_expr_x_g2collapse =
      apply(matrix_expr_x, 2, function(x) {df_agr=aggregate(x, list(annotation_rows[[group.by2]]), mean); df_agr$x})
    matrix_expr_f_g2collapse =
      apply(matrix_expr_f, 2, function(x) {df_agr=aggregate(x, list(annotation_rows[[group.by2]]), mean); df_agr$x})
    
    str_annot = aggregate(matrix_expr_Z[,1], list(annotation_rows[[group.by2]]), mean)$Group.1
    
    rownames(matrix_expr_Z_g2collapse) = str_annot
    rownames(matrix_expr_x_g2collapse) = str_annot
    rownames(matrix_expr_f_g2collapse) = str_annot
    
    annotation_rows_g2collapse = data.frame(annot=str_annot)
    rownames(annotation_rows_g2collapse) = str_annot
    names(annotation_rows_g2collapse) = group.by2
    print(paste0('Done -- ', date()))
    
    # can be convenient to know maximum of full matrix
    if (only_return_max) {
      print('Returning max values and exiting function.')
      return(list(
        max_val           = max(matrix_expr_x, na.rm=T),
        max_val_collapsed = max(matrix_expr_x_g2collapse, na.rm=T)
      ))
    }
    
    print('Creating heatmap 1')
    # Determine color scale
    if (!is.null(manual_zlims)) {max_val = manual_zlims[1]} else {max_val = max(matrix_expr_x, na.rm=T)}
    mybreaks  = seq(0,max_val,max_val/100)
    # my_colors = viridis_pal()(100)
    
    # Z-score matrix
    p1_Z=pheatmap(matrix_expr_Z, annotation_row = annotation_rows, annotation_col = annotation_cols, 
                  cluster_cols = F, cluster_rows = F, fontsize = fontsize)
    # Expression matrix
    p1_x=pheatmap(matrix_expr_x, annotation_row = annotation_rows, annotation_col = annotation_cols, 
                  cluster_cols = F, cluster_rows = F, fontsize = fontsize, breaks = mybreaks)
    # Fraction of cells detected matrix
    p1_f=pheatmap(matrix_expr_f, annotation_row = annotation_rows, annotation_col = annotation_cols, 
                  cluster_cols = F, cluster_rows = F, fontsize = fontsize, breaks = mybreaks)
    
    # return(p)
    print(paste0('Done -- ', date()))
    
    print('Creating heatmap 2')
    if (!is.null(manual_zlims)) {max_val = manual_zlims[2]} else {max_val = max(matrix_expr_x_g2collapse, na.rm=T)}
    mybreaks  = seq(0,max_val,max_val/100)
    
    p2_Z=pheatmap(matrix_expr_Z_g2collapse, annotation_row = annotation_rows_g2collapse, annotation_col = annotation_cols, 
                  cluster_cols = F, cluster_rows = F, fontsize = fontsize)
    p2_x=pheatmap(matrix_expr_x_g2collapse, annotation_row = annotation_rows_g2collapse, annotation_col = annotation_cols, 
                  cluster_cols = F, cluster_rows = F, fontsize = fontsize, breaks = mybreaks)
    p2_f=pheatmap(matrix_expr_f_g2collapse, annotation_row = annotation_rows_g2collapse, annotation_col = annotation_cols, 
                  cluster_cols = F, cluster_rows = F, fontsize = fontsize, breaks = mybreaks)
    
    print(date())
    print(paste0('Done -- ', date()))
    
    
    # Save the data if a subdir is given
    if (!is.null(savepath)) {
      print(paste0('Saving data -- ', date()))
      
      # To do , redundant, but determine max values again
      max_val           = max(matrix_expr_x, na.rm=T)
      max_val_collapsed = max(matrix_expr_x_g2collapse, na.rm=T)
        
      data_LR_expression =
        list(matrix_expr_x=matrix_expr_x, matrix_expr_Z=matrix_expr_Z, matrix_expr_f=matrix_expr_f, annotation_rows=annotation_rows, annotation_cols=annotation_cols,
              matrix_expr_Z_g2collapse=matrix_expr_Z_g2collapse, matrix_expr_x_g2collapse=matrix_expr_x_g2collapse, matrix_expr_f_g2collapse=matrix_expr_f_g2collapse, annotation_rows_g2collapse=annotation_rows_g2collapse,
             max_val=max_val, max_val_collapsed=max_val_collapsed)
      save(list='data_LR_expression', file = savepath)
      print('Saved to data_LR_expression.')
    }
    
    return(list(p1_x=p1_x, p2_x=p2_x, p1_Z=p1_Z, p2_Z=p2_Z, p1_f = p1_f, p2_f = p2_f))
    
}





# code from earlier function: regulon_overlap_heatmap (in "regulon" script)
geneset_overlap_heatmap = function(gene_sets, myfontsize=8) { #, makeallheatmaps=T, cutree_k=NULL, saveplots=T) {
    
    # Create pairs to compare
    df_compare = tidyr::expand_grid(x=names(gene_sets), y=names(gene_sets))
    
    # Calculate overlaps
    df_compare$overlap = sapply(1:dim(df_compare)[1], function(X) { 
        sum(gene_sets[[df_compare$x[X]]] %in% gene_sets[[df_compare$y[X]]]) /
            min(length(gene_sets[[df_compare$x[X]]]), length(gene_sets[[df_compare$y[X]]]))
        })
    
    # Create matrix
    matrix_compare <- reshape2::acast(df_compare, x~y, value.var="overlap")
    
    # Show heatmap
    p0=pheatmap(matrix_compare, clustering_method = 'ward.D2')
    
    return(list(plot=p0, matrix_compare=matrix_compare))

}
