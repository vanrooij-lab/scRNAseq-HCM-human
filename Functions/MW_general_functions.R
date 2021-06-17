
################################################################################

library(mutoss) # required for something in volcano code
library("limma") # required for Venn

################################################################################
# New volcano stuff




# strip__chrXX = function(gene_query) {
#     return(gsub("__chr(\\d+|[MXY])", '', gene_query))
# }
# 
# get_expression_gene2 = function(expression_matrix, gene_query) {
#     
#     sel_matrix = 
#         expression_matrix[ gsub("__chr(\\d+|[MXY])", '', rownames(expression_matrix)) %in% 
#                             gsub("__chr(\\d+|[MXY])", '', gene_query),]
#     
#     return(sel_matrix)
#     
# }

# Volcano plot of NPPA
#GENE = 'NPPA__chr1'
get_volcano_df3 = function(expression_matrix,my_gene,calc_qvals=T,no_expression_val=0,min_cell_expressed=.1,manual_expression=NULL) {
    
    # get the expression of the gene of interest (goi)
    if (is.null(manual_expression)) {
      #goi_expr = get_expression_gene2(expression_matrix, my_gene)
      goi_expr = expression_matrix[my_gene,]
    } else {
      goi_expr = manual_expression
    }
    if (mean(goi_expr>no_expression_val,na.rm=T)<min_cell_expressed) {
      warning(paste0('Your gene of interest is expressed in only ',round(mean(goi_expr>no_expression_val)*100,1), '% cells (below min treshold))'))
    }
    
    # Determine which genes to take along
    seen_in_frac_cells = apply(expression_matrix>no_expression_val, 1, mean)
    
    # calculate the correlations again, but now for cell_counts vs. genes
    start_time <- Sys.time()
    my_corrs_current = lapply(as.data.frame(t(expression_matrix[seen_in_frac_cells>=min_cell_expressed,])), function(X) {cor.test(X,goi_expr)})
    end_time <- Sys.time(); elapsed_time <- end_time - start_time; elapsed_time
    
    # Put result in dataframe
    my_corrs_df_current = data.frame(corr=unlist(lapply(my_corrs_current,function(X) {X$estimate})),
                  pval=unlist(lapply(my_corrs_current,function(X) {X$p.value})))
    
    # Add some more params
    my_corrs_df_current$pval.adj = p.adjust(my_corrs_df_current$pval, method = 'BH') # BH
    if (calc_qvals) {
      my_corrs_df_current$qval = pval2qval(my_corrs_df_current$pval)[[1]]
    }    
    
    rownames(my_corrs_df_current) = rownames(expression_matrix[seen_in_frac_cells>=min_cell_expressed,])
    my_corrs_df_current$gene_name = rownames(expression_matrix[seen_in_frac_cells>=min_cell_expressed,])
    my_corrs_df_current$gene_name_short = rownames(expression_matrix[seen_in_frac_cells>=min_cell_expressed,])
    
    # add some more parameter to export to excel file
    my_corrs_df_current$nrCellsExpressGene  = apply(expression_matrix[seen_in_frac_cells>=min_cell_expressed,]>no_expression_val,1,sum)
    #my_corrs_df_current$averageExpression   = averageExpression_current
    
    # export correlations to excel file
    #outputDir_sub = paste0(outputDir2, '/Excel/')
    #if (!dir.exists(outputDir_sub)) {dir.create(outputDir_sub)}
    #xlsx::write.xlsx(my_corrs_df_current,   file=paste0(outputDir_sub, 'correlations_',gsub("__chr(\\d+|[MXY])", '',GENE),'.xlsx'), sheetName=paste0(gsub("__chr(\\d+|[MXY])", '',GENE)))
    
    # look at the p-val distribution to get an impression of the FDR
    #hist(my_corrs_df_current$pval,100)
    
    # Order the output
    my_corrs_df_current = my_corrs_df_current[order(my_corrs_df_current$pval.adj),]
    
    return(my_corrs_df_current)
}

# now plot the volcano
plot_volcano3 = function(my_corrs_df_current,mytextsize=15,mycex=3,manual_gene_name=NULL,
        NRLABELED=20,mypvaltreshold=0.01,custom_highlight_group=NULL,NRLABELED_hlgroup=NULL,
                  mypointsize=.1, mylinesize=.25) {

    if (!is.null(manual_gene_name)){
        current_gene=manual_gene_name
    } else {
        current_gene = my_corrs_df_current$gene_name_short[my_corrs_df_current$corr==1&!is.na(my_corrs_df_current$corr)]
    }
    
    my_corrs_df_current = my_corrs_df_current[!my_corrs_df_current$gene_name_short==current_gene,]
    if (sum(my_corrs_df_current$gene_name_short==current_gene)>1) { stop('Error') }
    
    pval_idx = order(my_corrs_df_current$pval.adj, decreasing = F)
    
    p = ggplot(my_corrs_df_current)+
        geom_point(aes(x=corr,y=-log10(pval.adj)),size=mypointsize)+
        geom_point(data= dplyr::filter(my_corrs_df_current,pval.adj<0.01&corr<0),aes(x=corr,y=-log10(pval.adj)),size=mypointsize, color='blue')+
        geom_point(data= dplyr::filter(my_corrs_df_current,pval.adj<0.01&corr>0),aes(x=corr,y=-log10(pval.adj)),size=mypointsize, color='red')+
        xlab('Correlation')+ylab('-log10(p-value)')+
        #geom_text_repel(data=my_corrs_df_current[my_corrs_df_current$pval<.001,],
        #    mapping = aes(x=corr,y=-log10(pval.adj),label=gene_name),cex=4,color='red',segment.color='gray')+
        give_better_textsize_plot(mytextsize)+
        theme_bw()+
        geom_hline(yintercept = -log10(mypvaltreshold), size=mylinesize)+
        #xlim(c(-.6,.6))+
        ggtitle(paste0('Genes correlated with ', current_gene))
    if (is.null(custom_highlight_group)) {
      # add labels to 1st N genes (determined by p-val)
      p=p+ggrepel::geom_text_repel(data=my_corrs_df_current[pval_idx[1:NRLABELED],],
            mapping = aes(x=corr,y=-log10(pval.adj),label=gene_name_short),color='black',cex=mycex,segment.size	=mylinesize)#,segment.color='gray')+#,cex=.3
    } else {
      # select a custom group of genes
      toshow_df = my_corrs_df_current[my_corrs_df_current$gene_name_short %in% custom_highlight_group,]
      # only display 1st N labels if desired
      if (!is.null(NRLABELED_hlgroup)) {toshow_df = toshow_df[1:NRLABELED_hlgroup,]}
      # add labels
      p=p+ggrepel::geom_text_repel(data=toshow_df,
            mapping = aes(x=corr,y=-log10(pval.adj),label=gene_name_short),color='black',cex=mycex)#,segment.color='gray')+#,cex=.3
        
    }
    p
    
}



##########

# Create simple Venn diagram
venn_simple_plot_mw = function(venn_list) {
    
    library("limma")
    venn_genes = union(venn_list[[1]], venn_list[[2]])
    gene_overlap_matrix=as.matrix(data.frame(var1=venn_genes %in% venn_list[[1]], var2=venn_genes %in% venn_list[[2]]))
    colnames(gene_overlap_matrix) = names(venn_list)
    my_ven <- vennCounts(gene_overlap_matrix)
    vennDiagram(my_ven, cex=1) 
}

# plot_Venn_MW_2lists_v3 = function(list1, list2, name1='list1', name2='list2') {
#     
#     # Create overlap matrix
#     all_genes = unique(c(list1, list2))
#     gene_overlap_df=as.matrix(1*data.frame(list1=all_genes %in% list1, 
#                                            list2=all_genes %in% list2))
#     
#     # Now plot this with Venn
#     
#     # library 1
#     # library("limma")
#     my_ven <- vennCounts(gene_overlap_df)
#     vennDiagram(my_ven, cex=1, names = c(name1,name2)) # note: cex is general text scaling parameter
#     
#     # library 2
#     #library(venneuler)
#     #v<-venneuler::venneuler(gene_overlap_df)
#     #v$labels=c(paste0(name1,'\n(',sum(gene_overlap_df[,1]),')'),
#     #           paste0(name2,'\n(',sum(gene_overlap_df[,2]),')'))
#     #par(cex = 1.4) 
#     #p=plot(v, main = paste0("Overlap\n",sum(gene_overlap_df[,1]&gene_overlap_df[,2])), cex.main = 1)
# 
#     #return(p)
#     
# }

 


