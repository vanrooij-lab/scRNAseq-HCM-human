



################################################################################
# Start defining function
MW_determine_regulons_part1 = function(expression_matrix, calculate_p = F, 
                                       outputDir, analysis_name = '1',
                                       min_expression_fraction_genes=.05, 
                                       MYPADJUSTMETHOD='BH') {
    # example settings (but it depends completely on the data):
    # chosen_cutoff_parameter='r'
    # calculate_p = F
    # p_or_r_cutoff = .35
    # connectedness_cutoff = 30
    # hierarchical_cutoff = 1.3

    # note that p_or_r_cutoff is set at the beginning of this document
    # note that r_cutoff_value also needs to be set     
        
    ################################################################################
    # Correlations between genes
    # m.wehrens@hubrecht.eu, 2019-07
    # 
    # This script combines elements from the scripts 
    # - HCM_20190530_HCM_looking_at_R_cutoff_in_detail.R
    # - HCM_20190603_HCM_looking_at_correl_genes_sel.R
    # The idea is that we can determine a cutoff value where we think correlations
    # are significant. Consequently, we can use that information to determine which genes
    # are significantly correlated with multiple other genes. When we look at a subselection
    # of genes that are all correlated with multiple other genes, we can find more structure
    # in the gene-gene correlation matrix.
    #
    
    # Set the desired final p-value selection criterium here
    #p_or_r_cutoff <- 1/1000
    
    ################################################################################
    # Some pre-processing
    
    # Just to make sure
    expression_matrix = as.matrix(expression_matrix) 
    
    # create object where we can store relevant params from this analysis
    regulon_object = list()
    # now store for future use in object
    regulon_object$expression_matrix_original = expression_matrix
    
    # Now do a first gene filter, require min. expression for cells
    # Let's calculate first in how many cells each gene is expressed
    gene_expressed_in_X_cells<-rowSums(1*(expression_matrix>min(expression_matrix)))
    cells_N<-dim(expression_matrix)[2] # number of cells
    sel_genes_05percent_expr<-gene_expressed_in_X_cells>round(cells_N*min_expression_fraction_genes) # selection for genes expressed in >5% cells
    # perform actual filtering
    regulon_object$expression_matrix =
        regulon_object$expression_matrix_original[sel_genes_05percent_expr,]
    
    # 
    print(paste0('Keeping ', sum(sel_genes_05percent_expr), ' genes of ',length(sel_genes_05percent_expr)))
    
    # Generate 100 distinct colors
    n <- 100
    regulon_object$colors_random_100 <- distinctColorPalette(n)
    #pie(rep(1, n), col=colors_random_100)
    
    # Set and generate output dir
    regulon_object$analysis_name = analysis_name
    outputDir_sub = paste0(outputDir, 'analysis_', analysis_name, '/regulon/')
    if (!dir.exists(outputDir_sub)) {dir.create(outputDir_sub)}
    regulon_object$outputDir_sub = outputDir_sub
    
    ################################################################################
    # Determine significance correlations
    
    # First create a scrambled expression matrix that shouldn't have correlations
    # Previously, I created a matrix equal to the size of the expression
    # matrix to do this, but for large expression matrices, this takes quite long..
    # So now I select 1000 genes at random, and calculate a 'scrambled' expression
    # matrix from this
    print('Calculating scrambled expression matrix for 0-hypothesis expectations.')
    start_time <- Sys.time()
    #scrambled_expression_matrix<-expression_matrix[sample(dim(expression_matrix)[1],size = 1000),]
    scrambled_expression_matrix<-expression_matrix
    thedims<-dim(scrambled_expression_matrix)
    for (ii in seq(1,dim(scrambled_expression_matrix)[1])) {
        scrambled_expression_matrix[ii,] <- scrambled_expression_matrix[ii,sample(thedims[2])]
        if ((ii%%(thedims[1]/100))==0) {print(paste0(100*ii/thedims[1],'% done'))}
    }
    end_time <- Sys.time(); print(paste0('Task performed in ', round(end_time - start_time,2), ' ',units((end_time-start_time)),'..'))
        # sanity check:
        # cor(expression_matrix[3,],expression_matrix[2,])
        # cor(scrambled_expression_matrix[3,],scrambled_expression_matrix[2,])
    
    # First determine the correlation matrix for the actual data
    print('Calculating correlation matrix')
    start_time <- Sys.time()
    cor_out <- cor(t(expression_matrix))
    end_time <- Sys.time(); print(paste0('Task performed in ', round(end_time - start_time,2), ' ',units((end_time-start_time)),'..'))
        # this takes a few minutes (4.7 mins on my computer)
    
    # Then determine the correlation matrix for the scrambled data
    print('Calculating correlation matrix for scrambled data')
    start_time <- Sys.time()
    cor_out_scrambled <- cor(t(scrambled_expression_matrix))
    end_time <- Sys.time(); print(paste0('Task performed in ', round(end_time - start_time,2), ' ',units((end_time-start_time)),'..'))
        # this takes a few minutes (also 4.7 mins on my computer)
    
    # Create histograms for both
    # note that the number of bins is rather important here, 
    # peak should be covered by enough bins
    dx=.01
    mybreaks               <- seq(-1,1,dx)
    hist_out_cor           <- hist(as.vector(cor_out[cor_out<1]), breaks=mybreaks, plot=F)
    hist_out_cor_scrambled <- hist(as.vector(cor_out_scrambled[cor_out_scrambled<1]),breaks=mybreaks,plot=F)
    
    ################################################################################
    
    # Now let's plot both distributions and also the cutoff
    p_corr_distr <- ggplot()+
        geom_line(data=data.frame(x=hist_out_cor_scrambled$mids,y=hist_out_cor_scrambled$density),
            aes(x=x,y=y,colour='scrambled')) +
        geom_line(data=data.frame(x=hist_out_cor$mids,y=hist_out_cor$density),     
            aes(x=x,y=y, colour='observed')) +
         #geom_vline(xintercept = r_cutoff_value, colour='blue')+
         #geom_vline(xintercept = -r_cutoff_value, colour='blue')+
         scale_colour_manual("", 
                          breaks = c("scrambled", "observed"),
                          values = c("red", "black"))+
        xlab('Correlation')+ylab('Density')+
        theme_bw()+
        theme(legend.position = c(.8,.8))+
        xlim(-0.5,0.5)+
        give_better_textsize_plot(10)
    print(p_corr_distr)
    
    ggsave(paste0(regulon_object$outputDir_sub,'hist_observed_corrs.pdf'),plot=p_corr_distr,
            units='mm',width=75,height=75,dpi=600)
    
    p_corr_distr_ylog10=p_corr_distr+scale_y_continuous(trans='log10')+theme(legend.position = 'top')
    print(p_corr_distr_ylog10)
    
    ggsave(paste0(regulon_object$outputDir_sub,'hist_observed_corrs_ylog10.pdf'),plot=p_corr_distr_ylog10,
            units='mm',width=75,height=75,dpi=600)

    ################################################################################
    # To determine the cutoff of the correlation, let's look at the estimated
    # false-positive rates:
    
    my_false_positives = hist_out_cor_scrambled$density/hist_out_cor$density
    low_dist = rep(2,length(hist_out_cor_scrambled$density))/(sum(hist_out_cor_scrambled$counts)*dx)/
        hist_out_cor$density
        # two because the matrix has two entries for each pair
    #my_false_positives[my_false_positives>1]=1
    #my_false_positives[my_false_positives<0.05]=0
    delta_0hypothesis  = hist_out_cor$density-hist_out_cor_scrambled$density
    
    # False positive percentages per correlation
    # The dotted line shows the value in case there is only one positive
    # oobservation in the scrambled matrix (helps with filtering out high FP values that 
    # are just result of stochasticity).
    # 
    # This plot should help you decide where to chose your cutoff for significant correlations.
    # It compares the previously plotted distributions of the correlation coefficients
    # observed in your actual data, and the distribution of correlation coefficient in 
    # a scrambled version of the data where no 'real' correlations exist, only 'accidental' ones
    # and calculates which fraction of the observed correlations in your data are probably
    # due to chance (and thus also which ones are real, ie the remaining fraction). 
    # The plot shows three boundaries that might help you deciding the correlation
    # coefficient cutoff. 
    p_cutoff_decision <- ggplot()+
        #geom_line(data=data.frame(x=hist_out_cor_scrambled$mids,y=delta_0hypothesis),
        #    aes(x=x,y=y,color='Difference')) +
        geom_bar(data=data.frame(x=hist_out_cor$mids,y=my_false_positives),     
            aes(x=x,y=y, color='a'), stat='identity', fill='red')+#,linetype='solid') +
        #geom_point(data=data.frame(x=hist_out_cor_scrambled$mids,y=hist_out_cor_scrambled$counts),     
        #    aes(x=x,y=y, color='CNTS')) +
        geom_line(data=data.frame(x=hist_out_cor$mids,y=low_dist),     
            aes(x=x,y=y, color='b'))+#,linetype='longdash') +
        geom_hline(aes(yintercept = 0.01, color='c'))+
        geom_hline(aes(yintercept = 0.03, color='d'))+
        geom_hline(aes(yintercept = 0.05, color='e'))+
        theme_bw()+
         #geom_vline(xintercept = r_cutoff_value, color='blue')+
         #geom_vline(xintercept = -r_cutoff_value, color='blue')+
        scale_colour_manual(  name='Legend',
                              breaks = c('a',    'b', 'c', 'd', 'e'),
                              values = c('red',   'grey',  'black',     'black',     'black'),
                              labels = c('FP',    '1_obs', '1_percent', '3_percent', '5_percent'))+
        scale_x_continuous(limits=c(-1,1), breaks = seq(-1,1,.2), minor_breaks=seq(-1,1,.05))+
        scale_y_continuous(limits=c(0,1),oob = rescale_none)+
        xlab('Correlation')+ylab('Fraction false positives')+
        theme(legend.position = c(.8,.8))+
        #theme(legend.position = 'none')+#c(.8,.8))+
        give_better_textsize_plot(10)
    print(p_cutoff_decision)
    
    ggsave(paste0(regulon_object$outputDir_sub,'Correlations_FPR.pdf'),plot=p_cutoff_decision,
            units='mm',width=100,height=75,dpi=600)
    
    ################################################################################
    # An alternative method to determining significance for multiple testing
    # is adjusting the p-value using established statistical methods
    #
    # So let's do that, and use the p.adjust function
    
    # First, we need a matrix with all p-values
    # See also: https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
    
    # Function to calculate theoretical p-value
    my_freedom <- dim(expression_matrix)[2]-2
    calculate_pvalue <- function(R, df) {
        t1<-sqrt(df) * R/sqrt(1 - R^2)
        p_val <- 2*min(pt(t1,df),pt(t1,df,lower.tail = FALSE))
        return(p_val)
    }
    
    # Let's first also reproduce the above plot, but now with an adjusted axis
    p_vals = -log10(sapply(X = hist_out_cor$mids, FUN = calculate_pvalue, df=my_freedom))
    R_sign = sapply(hist_out_cor$mids>=0, function(x) {if(x) {'positive'} else {'negative'}})
    p_cutoff_decision2 <- ggplot()+
        geom_line(data=data.frame(x=p_vals[!is.infinite(p_vals)],y=my_false_positives[!is.infinite(p_vals)],R_sign=R_sign[!is.infinite(p_vals)]),     
            aes(x=x,y=y), color='red')+#,linetype='solid') + # , stat='identity'
        geom_line(data=data.frame(x=p_vals[!is.infinite(p_vals)],y=low_dist[!is.infinite(p_vals)],R_sign=R_sign[!is.infinite(p_vals)]),     
            aes(x=x,y=y),color='grey')+#,linetype='longdash') +
        geom_hline(yintercept=0.01)+geom_hline(yintercept=0.05)+
        facet_grid(rows = 'R_sign')+
        theme_bw()+xlab('-log10 Theoretical p-value based on R-value')+ylab('Fraction false positives')+
        theme(legend.position = c(.8,.8))+give_better_textsize_plot(10)+
        geom_vline(xintercept = 5)+geom_vline(xintercept = 10)+
        ggtitle('red=FPR, black,FPR=0.01, 0.05\nblack,p=10^-5,10^-10, grey=FPR for 1 obs')
    print(p_cutoff_decision2)
    
    ggsave(paste0(regulon_object$outputDir_sub,'Correlations_FPR_pvals.pdf'),plot=p_cutoff_decision2,
            units='mm',width=100,height=75,dpi=600)
    
    
    # Then calculate the p-values (if desired)
    if (calculate_p) {
        
        print('Calculating p-values')
        start_time <- Sys.time() # just to see how long commands below take..
        
        # Calculate p-values
        p_val_matrix<-matrix(sapply(as.vector(cor_out), calculate_pvalue, df=my_freedom),nrow=dim(cor_out)[1])
        
        # code that is not executed:
        if (F) { 
            # Alternative using multiple cores (this also uses crazy amounts of memory):
            p_val_matrix<-matrix(unlist(mclapply(as.vector(cor_out), calculate_pvalue, df=my_freedom, mc.cores = 10)),nrow=dim(cor_out)[1])
        }
        
        # Then, we produce adjusted the p-values using p.adjust
        p_val_matrix_adjusted<-matrix(p.adjust(as.vector(p_val_matrix), method=MYPADJUSTMETHOD),nrow=dim(p_val_matrix)[2])
        
        end_time <- Sys.time(); print(paste0('Task performed in ', round(end_time - start_time,2), ' ',units((end_time-start_time)),'..')) # timing of operation

        regulon_object$p_val_matrix = p_val_matrix        
        regulon_object$p_val_matrix_adjusted = p_val_matrix_adjusted
    
    }
    ################################################################################
    # This plot just illustrates how the p-value is calculated
    # theoretically from the student-t distribution.
    # But it also shows which r-cutoff you might choose based on 
    # the theoretical calculation of p-values.
    
    if (calculate_p) {
        
        df = dim(expression_matrix)[2]-2
        fn_dtr<-function(r){t=sqrt(df) * r/sqrt(1 - r^2); dt(t, df = df)}
        
        p_tdist = ggplot(data.frame(  pval_adj=as.vector(p_val_matrix_adjusted)[1:1000],
                            pval   = as.vector(p_val_matrix)[1:1000],
                            corr = as.vector(cor_out)[1:1000]))+
            geom_point(aes(x=corr, y=pval_adj, color='adjusted'))+
            geom_point(aes(x=corr, y=pval,     color='non-adjusted'))+
            stat_function(data=data.frame(x=c(-1, 1, .0001)), aes(x=x,color='formula'), fun=fn_dtr, geom="line", size=1) + 
            geom_hline(yintercept = 0.01)+
            scale_x_continuous(limits=c(-.8,.8), breaks = seq(-.8,.8,.2), minor_breaks=seq(-1,1,.05))+
            scale_y_continuous(limits=c(-0.01,.25),oob = rescale_none)+
            theme_bw()+give_better_textsize_plot(10)+
            xlab('Correlation coefficient')+ylab('(Adjusted) p-value')+
            scale_color_manual(name='Legend',breaks = c('adjusted','non-adjusted','formula'), values = c('black','red','slateblue4'))+#, labels = )+
            #theme(legend.position = c(.8,.8))+
            theme(legend.position = 'bottom')+
            ggtitle('Calculating p-val for R (student-t)')
        print(p_tdist)
    
    }
    
    ggsave(paste0(regulon_object$outputDir_sub,'Correlations_Pvalue_Relation.pdf'),plot=p_tdist,
            units='mm',width=75,height=75,dpi=600)

    # Prepare return value and return
    if (calculate_p) {
        regulon_object$p_corr_distr = p_corr_distr
        regulon_object$p_cutoff_decision = p_cutoff_decision
        regulon_object$p_tdist = p_tdist
    }
    
    regulon_object$cor_out = cor_out
    regulon_object$cor_out_scrambled = cor_out_scrambled
    regulon_object$expression_matrix = expression_matrix
    
    return(regulon_object)

}

################################################################################
# now let's actually start looking at the matrix

MW_determine_regulons_part2 = function(regulon_object, chosen_cutoff_parameter='r', 
    p_or_r_cutoff=.35) {
    
    # 'loading' params from object
    cor_out = regulon_object$cor_out
    cor_out_scambled = regulon_object$cor_out_scambled
    p_val_matrix_adjusted = regulon_object$p_val_matrix_adjusted
    
    #if (!is.null(p_or_r_cutoff)) {
    
    # also create a correlation matrix with NA values instead of the diagonal 1 values
    cor_out_NA <- cor_out
    cor_out_NA[row(cor_out_NA)==col(cor_out_NA)]<-NA # remove diagonal ones
    
    # now calculate the amount of significant correlations per gene
    if (chosen_cutoff_parameter == 'r') {
        amount_sign_corrs_per_gene <- 
            apply(1*(cor_out_NA > abs(p_or_r_cutoff) | cor_out_NA < -abs(p_or_r_cutoff)), 2, sum, na.rm=T)
    } else if (chosen_cutoff_parameter == 'p') {
        # remove NA values
        regulon_object$p_val_matrix_adjusted_NA <- p_val_matrix_adjusted
        regulon_object$p_val_matrix_adjusted_NA[row(regulon_object$p_val_matrix_adjusted_NA)==col(regulon_object$p_val_matrix_adjusted_NA)]<-NA # remove diagonal ones
        # calculated significant corrs
        amount_sign_corrs_per_gene <- 
            apply(1*(p_val_matrix_adjusted < p_or_r_cutoff), 2, sum, na.rm=T)
    }
    amount_sign_corrs_per_gene_ordered <- 
            amount_sign_corrs_per_gene[order(amount_sign_corrs_per_gene, decreasing = T)]
    
    # show histogram of how many genes genes are correlated with
    # 
    # This plot should help you to determine the cutoff to take into account genes
    # based on with how many other genes they correlate
    amount_sign_corrs_per_gene
    hist_out_intrx<-hist(amount_sign_corrs_per_gene,100)
    p_connectedness = ggplot(data=data.frame(x=hist_out_intrx$mids,y=hist_out_intrx$density))+
        geom_line(aes(x=x,y=y))+
        xlab('# genes that significantly correlate')+ylab('Density (normalized)')+
        theme_bw()+give_better_textsize_plot(10)+
        ggtitle('Connectedness of genes')
        #xlim(c(0,100))
    print(p_connectedness)
    
    ggsave(paste0(regulon_object$outputDir_sub,'Connectedness.pdf'),plot=p_connectedness,
            units='mm',width=75,height=75,dpi=600)
    
    p_connectedness2 = ggplot(data=data.frame(x=hist_out_intrx$mids,y=sum(hist_out_intrx$counts)-cumsum(hist_out_intrx$counts)))+
        geom_line(aes(x=x,y=y))+
        xlab('Cutoff treshold connectedness')+ylab('Genes that will remain in analysis')+
        theme_bw()+give_better_textsize_plot(10)+
        ggtitle('Connectedness of genes')
    print(p_connectedness2)
    
    ggsave(paste0(regulon_object$outputDir_sub,'Connectedness_vs_treshold.pdf'),plot=p_connectedness2,
            units='mm',width=75,height=75,dpi=600)
    
    regulon_object$p_connectedness  = p_connectedness
    regulon_object$p_connectedness2 = p_connectedness2
    regulon_object$p_or_r_cutoff = p_or_r_cutoff
    regulon_object$cor_out_NA = cor_out_NA
    return(regulon_object)
    
#} else { print('Set p_or_r_cutoff to continue') }
}

MW_determine_regulons_part3 = function(regulon_object, 
                connectedness_cutoff = NULL, 
                max_genes = NULL,
                # min_expression_fraction_genes=.05, # Should be done before calculating p-values!
                show_heatmap=F,
                chosen_cutoff_parameter='r') {
    # min_expression_fraction_genes = in how many cells (as a fraction of total) 
    #       each gene should be expressed to be taken along

#if (!is.null(connectedness_cutoff)) {
    
    # 'loading' some parameters from object
    p_or_r_cutoff = regulon_object$p_or_r_cutoff 
    cor_out = regulon_object$cor_out
    cor_out_NA = regulon_object$cor_out_NA 
    p_val_matrix_adjusted_NA = regulon_object$p_val_matrix_adjusted_NA
    expression_matrix=regulon_object$expression_matrix
    # 'saving' some parameters to object for later use
    regulon_object$connectedness_cutoff = connectedness_cutoff
    regulon_object$max_genes = max_genes
    
    # XXXX
    
    # Now select which genes should be taken into account for "regulon" matrix
    # selection based on 
    # - genes that are expressed in at least X% of cells (see above)
    # - genes that are (a) significantly correlated (cutoff) with 
    # - (b) Y other genes (connectedness_cutoff)
    if (chosen_cutoff_parameter == 'r') {
        count_of_high_corrs_per_gene_stringent <- apply(1*(cor_out_NA > p_or_r_cutoff | cor_out_NA < -p_or_r_cutoff), 2, sum, na.rm=T)
    } else if (chosen_cutoff_parameter == 'p') {
        count_of_high_corrs_per_gene_stringent <- apply(1*(p_val_matrix_adjusted_NA < p_or_r_cutoff), 2, sum, na.rm=T)    
    }
    sel_genes_2<-(count_of_high_corrs_per_gene_stringent>connectedness_cutoff)
    # sel_genes_2<-(count_of_high_corrs_per_gene_stringent>connectedness_cutoff)&sel_genes_05percent_expr # sel_genes_05percent_expr should be done at the start for p-vals
    
    # only take along top X connected genes if desired (in case computational burden too high)
    topX_connected_genes_idxs = order(count_of_high_corrs_per_gene_stringent,decreasing = T)[1:min(max_genes,sum(sel_genes_2))]
    if (!is.null(max_genes)) {
        sel_genes_2[-topX_connected_genes_idxs] = F
    }
    sum(1*sel_genes_2)
    
    # show an updated connectedness plot with selection criteria applied
    hist_out_intrx2 <- hist(count_of_high_corrs_per_gene_stringent,100,plot=F)
    p<-ggplot(data=data.frame(x=hist_out_intrx2$mids,y=hist_out_intrx2$density))+
        geom_line(aes(x=x,y=y))+
        xlab(paste0('# correlations this gene with others, ',chosen_cutoff_parameter,switch(chosen_cutoff_parameter, 'r'='>','p'='<'),round(p_or_r_cutoff,6),''))+ylab('Density (normalized)')+
        ggtitle('Connectedness')+
        # ggtitle(paste0('For genes that are expressed in ',round(min_expression_fraction_genes*100,2),'% of cells'))+
        theme_bw()+give_better_textsize_plot(10)
        #xlim(c(0,100))+
    print(p)
    
    ggsave(paste0(regulon_object$outputDir_sub,'Connectedness_selectedgenes.pdf'),plot=p,
            units='mm',width=75,height=75,dpi=600)
    
    # now select that subsection of the matrix that only contains those genes
    cor_out_selected_2 <- cor_out[sel_genes_2,][,sel_genes_2]
    # and display the heatmap
    if (show_heatmap) {
        pheatmap_out_sel2 <- pheatmap(cor_out_selected_2)
        print(pheatmap_out_sel2)
        #pheatmap(cor_out_selected_2, cluster_rows = F, cluster_cols = F)
        #image(cor_out_selected_2)
    }
    
    # save to object for export
    regulon_object$cor_out_selected_2 = cor_out_selected_2
    regulon_object$chosen_cutoff_parameter = chosen_cutoff_parameter
    # regulon_object$sel_genes_05percent_expr=sel_genes_05percent_expr
    
    return(regulon_object)

#} else {print('Set connectedness_cutoff to continue')}

}
    
################################################################################
# Now let's take a look at clustering the matrix with a bit more
# fine-tuning

MW_determine_regulons_part4 = function(regulon_object) {

# Parameter relevant for this section:
#if (!is.null(hierarchical_cutoff)) {
    
    # 'load' params from to object
    cor_out_selected_2 = regulon_object$cor_out_selected_2
    
    # use the built-in hierarchical clustering algorithm
    # option 1
    # use correlation as distance between points
    #hclust_out      <- hclust(as.dist(1-cor_out_selected_2)) #, method = 'complete') # method='average'
    # option 2
    # calculate distances between rows, which is a more indirect measure
    # (gives a bit clearer blocks though)
    #hclust_out      <- hclust(dist(cor_out_selected_2))#, method='ward.D2') #, method = 'complete') # method='average'
    # option 3
    # take correlations as distance measures between points, but now use method that takes
    # into account multiple points in clusters to calculate cluster-point or cluster-cluster
    # distances; this also gives clearer blocks AND is more intuitive
    hclust_out      <- hclust(as.dist(1-cor_out_selected_2), method = 'ward.D2') # method='average'
    #plot(as.dendrogram(hclust_out))#, hang = -1, cex = 0.6)
    p = ggdendrogram(hclust_out, rotate = TRUE, theme_dendro = FALSE)+
                ggtitle(regulon_object$analysis_name)+
                theme_bw()+give_better_textsize_plot(10)+
                theme(axis.text.y = element_text(size=.5))
        #
        
    # Let's try to systemetically investigate the best cutoff
    clust_dens_out <- plot_clust_join_density(hclust_out)
        # see also comments inside this function; it looks at how
        # the distances between points that are joined develops.
    
    print(clust_dens_out$p2)
    
    ggsave(paste0(regulon_object$outputDir_sub,'regulons_dendrogram_cutoff_analysis1.pdf'),plot=clust_dens_out$p1+give_better_textsize_plot(10),
            units='mm',width=75,height=75,dpi=600)
    ggsave(paste0(regulon_object$outputDir_sub,'regulons_dendrogram_cutoff_analysis2.pdf'),plot=clust_dens_out$p2+give_better_textsize_plot(10),
            units='mm',width=75,height=75,dpi=600)
    
    p=p+geom_hline(yintercept = clust_dens_out$y_cutoff[1], color='blue')+
        geom_hline(yintercept = clust_dens_out$y_cutoff[2], color='red')
        
    ggsave(paste0(regulon_object$outputDir_sub,'regulons_dendrogram.pdf'),plot=p,
            units='mm',width=75,height=150,dpi=600)
    
    print(p)
    
    regulon_object$auto_cutoff1 = clust_dens_out$y_cutoff[1]
    regulon_object$auto_cutoff2 = clust_dens_out$y_cutoff[2]
    regulon_object$hclust_out = hclust_out
    
    return(regulon_object)
}

MW_determine_regulons_part5 = function(regulon_object, hierarchical_cutoff=NULL, KMAX_GAPTEST=20, hierarchical_cutoff_fallback=NULL, BIGB=100) {    
    
    GAPSTATMETHOD= 'Tibs2001SEmax' #'Tibs2001SEmax'
    
    # 'load' params from object
    hclust_out = regulon_object$hclust_out
    cor_out_selected_2 = regulon_object$cor_out_selected_2
    colors_random_100 = regulon_object$colors_random_100
    p_or_r_cutoff = regulon_object$p_or_r_cutoff
    # 'save' params to object
    regulon_object$hierarchical_cutoff = if (!is.null(hierarchical_cutoff)) {hierarchical_cutoff} else {paste0('gap-stat-',GAPSTATMETHOD)}
    
    # Either calculate optimal number of clusters (# regulons) using gap-stat
    if (is.null(hierarchical_cutoff)) {
        
        # It appears an issue can arise in the maxSE function, when not all SE.f >= 0.
        # (for which by the way I don't see a specific reason for occurring ..)
        # So I'll just use a try-catch statement here ..
        return1=tryCatch({
                
            # Using the hierarchical clustering we did previously
            gap_stat <- cluster::clusGap(cor_out_selected_2, FUN = function(matrix, k, mytree) {return(list(cluster=(cutree(mytree, k = k))))}, 
                mytree=hclust_out, K.max = min(KMAX_GAPTEST, nrow(cor_out_selected_2)), B = BIGB)
                # Joep used kmeans before 
                # gap_stat <- cluster::clusGap(cor_out_selected_2, FUN = kmeans, nstart = 10, K.max = 20, B = 10) # note: nstart is #random configs to start, B=bootstrapping, K.max=max. # clusters
            
            nCluster = cluster::maxSE(f=gap_stat$Tab[,'gap'], SE.f = gap_stat$Tab[,'SE.sim'], method = GAPSTATMETHOD)
            p=ggplot(data.frame(gap=gap_stat$Tab[,'gap'],K=1:KMAX_GAPTEST))+
                geom_vline(xintercept = nCluster)+
                geom_line(aes(x=K, y=gap))+geom_point(aes(x=K, y=gap))+theme_bw()+
                ggtitle(paste0('Gap-stat, K=',nCluster,' optimal'))+give_better_textsize_plot(8)
            p
            ggsave(paste0(regulon_object$outputDir_sub,'regulon_K_gap_stat.pdf'),plot=p,units='mm',width=100,height=50,dpi=600)
        
            #gene_clustering_assign_hierarchical <- cutree(hclust_out, k = nCluster)
            #regulon_object$nCluster=nCluster # save for later
            return1 = list(gene_clustering_assign_hierarchical = cutree(hclust_out, k = nCluster),
                            nCluster = nCluster)
            return(return1)
            
        }, error = function(error_condition) {
            # error-handler-code
            
            print('Failed to successfully execute gap-stat strategy.')
            print('Trying once more with less bootstrapping.')
            
            ##########
            
                return2=tryCatch({
                        
                    # ATTEMPT 2, B=10
                    # Using the hierarchical clustering we did previously
                    gap_stat <- cluster::clusGap(cor_out_selected_2, FUN = function(matrix, k, mytree) {return(list(cluster=(cutree(mytree, k = k))))}, 
                        mytree=hclust_out, K.max = min(KMAX_GAPTEST, nrow(cor_out_selected_2)), B = 10)
                        # Joep used kmeans before 
                        # gap_stat <- cluster::clusGap(cor_out_selected_2, FUN = kmeans, nstart = 10, K.max = 20, B = 10) # note: nstart is #random configs to start, B=bootstrapping, K.max=max. # clusters
                    
                    nCluster = cluster::maxSE(f=gap_stat$Tab[,'gap'], SE.f = gap_stat$Tab[,'SE.sim'], method = GAPSTATMETHOD)
                    p=ggplot(data.frame(gap=gap_stat$Tab[,'gap'],K=1:KMAX_GAPTEST))+
                        geom_vline(xintercept = nCluster)+
                        geom_line(aes(x=K, y=gap))+geom_point(aes(x=K, y=gap))+theme_bw()+
                        ggtitle(paste0('Gap-stat, K=',nCluster,' optimal'))+give_better_textsize_plot(8)
                    p
                    ggsave(paste0(regulon_object$outputDir_sub,'regulon_K_gap_stat.pdf'),plot=p,units='mm',width=100,height=50,dpi=600)
                
                    #gene_clustering_assign_hierarchical <- cutree(hclust_out, k = nCluster)
                    #regulon_object$nCluster=nCluster # save for later
                    #regulon_object$errors='2nd attempt'
                    
                    return2 = list(gene_clustering_assign_hierarchical = cutree(hclust_out, k = nCluster),
                                    nCluster = nCluster,
                                    errors = '2nd attempt')
                    return(return2)
                    
                }, error = function(error_condition) {
                    # error-handler-code
                    
                    print('Failed again to execute gap-stat, falling back to hierarchical_cutoff_fallback.')
                    
                    #gene_clustering_assign_hierarchical <- cutree(hclust_out, h=hierarchical_cutoff_fallback)
                    #regulon_object$nCluster=length(unique(gene_clustering_assign_hierarchical))
                    #regulon_object$errors='fallback cutoff used'
                    
                    gene_clustering_assign_hierarchical=cutree(hclust_out, h=hierarchical_cutoff_fallback)
                    return2 = list(gene_clustering_assign_hierarchical = gene_clustering_assign_hierarchical,
                                    nCluster = length(unique(gene_clustering_assign_hierarchical)),
                                    errors = 'fallback cutoff used')
                    return(return2)
                })
            
            return1=return2
            return(return1)
        
            ##########
            
            #gene_clustering_assign_hierarchical <- cutree(hclust_out, h=hierarchical_cutoff_fallback)
            #regulon_object$nCluster=length(unique(gene_clustering_assign_hierarchical))
            
        })
        
        gene_clustering_assign_hierarchical = return1$gene_clustering_assign_hierarchical
        regulon_object$nCluster = return1$nCluster
        regulon_object$errors = return1$errors
        
    # Or use manually decided cutoff 
    } else {
        # Now use informed cutoff to determine clusters:
        # there might be multiple good choices 
        gene_clustering_assign_hierarchical <- cutree(hclust_out, h=hierarchical_cutoff)
    }
    
    gene_clustering_colors <- colors_random_100[gene_clustering_assign_hierarchical]
    
    ##heatmap(cor_out_selected_2)
    #heatmap(cor_out_selected_2)
    
    # shows heatmap but doesn't adjust scales
    #heatmap(cor_out_selected_2, 
    #    ColSideColors = gene_clustering_colors, 
    #    Colv=as.dendrogram(hclust_out), Rowv=as.dendrogram(hclust_out), 
    #        col = viridis_pal(option = "D")(100) )
    
    # 99% interval for color boundaries
    min_val <- cor_out_selected_2[order(as.vector(cor_out_selected_2),decreasing = F)][round(.01*sum(1*(cor_out_selected_2<1)))]
    max_val <- cor_out_selected_2[order(as.vector(cor_out_selected_2),decreasing = F)][round(.99*sum(1*(cor_out_selected_2<1)))]
    #min_val <- min(cor_out_selected_2)
    #max_val <- max(cor_out_selected_2[cor_out_selected_2<1])
    d_val   <- (max_val-min_val)/100
    
    # png(file = paste0(regulon_object$outputDir_sub,'regulons_gene_corrmatrix.png'),
    #     width = 300, height = 200, units = 'mm', res = 600)
    # heatmap.2(cor_out_selected_2, 
    #     ColSideColors = gene_clustering_colors,
    #     RowSideColors = gene_clustering_colors,
    #     Colv=as.dendrogram(hclust_out), Rowv=as.dendrogram(hclust_out), 
    #         col = viridis_pal(option = "D")(100),
    #     breaks = seq(min_val,max_val,d_val),tracecol=NA)
    # #print(p)
    # dev.off()
    
    p=pheatmap(mat = cor_out_selected_2, cluster_rows = hclust_out, cluster_cols = hclust_out, 
        fontsize_row = 2, fontsize_col = 2,color = viridis_pal(option = "D")(100), breaks = seq(min_val,max_val,d_val))
    ggsave(filename = paste0(regulon_object$outputDir_sub,'regulons_gene_corrmatrix.png'), plot=p, units = 'cm',height = dim(cor_out_selected_2)[2]*0.06, width = dim(cor_out_selected_2)[2]*0.06, dpi = 1200)
    
    # binary map
    if (regulon_object$chosen_cutoff_parameter=='r') {
        png(file = paste0(regulon_object$outputDir_sub,'regulons_gene_corrmatrix_blackwhite.png'),
                width = 300, height = 200, units = 'mm', res = 600)
        heatmap.2((cor_out_selected_2>p_or_r_cutoff)*1, 
            ColSideColors = gene_clustering_colors, 
            Colv=as.dendrogram(hclust_out), Rowv=as.dendrogram(hclust_out),tracecol=NA,
            col=c('#FFFFFF','#000000'))
        dev.off()
    }
    #print(p_b)
    
    # 'save' params to object
    #regulon_object$p = p
    regulon_object$gene_clustering_assign_hierarchical = gene_clustering_assign_hierarchical
    
    # Also create lists of which genes per regulon
    max_cl = max(regulon_object$gene_clustering_assign_hierarchical)
    correlated_sel2_gene_names = rownames(regulon_object$cor_out_selected_2)
    regulon_object$the_regulons = list()
    for (cl_idx in 1:max_cl) {
        
        # Get gene names from a cluster
        regulon_object$the_regulons[[cl_idx]] <- correlated_sel2_gene_names[gene_clustering_assign_hierarchical==cl_idx]
    }
    
    return(regulon_object)
    
#} else { print('Set hierarchical_cutoff to continue') }
}

################################################################################
# Now perform GO analysis and export some information to excel

MW_determine_regulons_part6 = function(regulon_object) {
    
    cor_out_selected_2= regulon_object$cor_out_selected_2
    gene_clustering_assign_hierarchical= regulon_object$gene_clustering_assign_hierarchical
    
    correlated_sel2_gene_names <- rownames(cor_out_selected_2)
    gene_names_SCS <- rownames(cor_out_selected_2)
    
    toString(correlated_sel2_gene_names[gene_clustering_assign_hierarchical==1])
    
    # decide on background genes
    genes_background <- rownames(regulon_object$expression_matrix)[regulon_object$sel_genes_05percent_expr]
    
    wb <- openxlsx::createWorkbook("Regulon")
    df_complete_export <- data.frame(Gene_names=character(),GO_terms=character())
    max_cl = max(gene_clustering_assign_hierarchical)
    df_export<-list()
    regulon_list <- list()
    for (cl_idx in 1:max_cl) {
        
        print(paste0("Looking at cluster ", toString(cl_idx),' (of ',max_cl,')'))
        
        # Get gene names from a cluster
        genes_query <- correlated_sel2_gene_names[gene_clustering_assign_hierarchical==cl_idx]
        
        # Perform GO analysis (with old function)
        #goResult <- give_GO_terms(config = cfg, gene_names_SCS = correlated_sel2_gene_names, 
        #                          genes_query = genes_query, 
        #                          genes_background = correlated_sel2_gene_names)
        #    # do we want background genes or not? and which type?
        
        # Perform GO analysis new way
        goResult <- analyzeGeneOntology_MW(config = cfg, 
            all_genes = union(genes_background, genes_query), 
            background_genes = genes_background, 
            genes_query = genes_query, 
            pathwayPCutoff=0.05, 
            GOKegg='GO', includeChildTerms=F)    
            # do we want background genes or not? and which type?
        
        # Create exportable data
        #df_export[[cl_idx]] <- cbind.fill(data.frame(Gene_names=genes_query), data.frame(GO_terms=goResult$Term), fill = "")

        # Create exportable data
        # Note: the function "cbind.fill" was from the "rowr" package, which was removed from CRAN
        # so now I wrote a custom workaround function; this is untested so I hope it works.
        df_export[[cl_idx]] <- bind_cols(pad.df.mw(list(data.frame(Gene_names=genes_query), data.frame(GO_terms=goResult$Term))))
        
        regulon_list[[cl_idx]] <- list()
        regulon_list[[cl_idx]]$regulon_genes          <- strip__chrXX(genes_query)
        regulon_list[[cl_idx]]$regulon_genes__chrX    <- genes_query
        regulon_list[[cl_idx]]$regulon_GO_terms <- goResult$Term
        regulon_list[[cl_idx]]$regulon_full_GO <- goResult
        
        # Export data to Excel sheet
        #if (cl_idx==1) {
            # xlsx::write.xlsx(df_export, file=paste0(outputDir_sub,'GO_analysis_co-correlations_v2.xlsx'), sheetName=paste0('cluster',cl_idx))
        #    openxlsx::write.xlsx(df_export[[cl_idx]], file=paste0(outputDir_sub,'GO_analysis_co-correlations_v2.xlsx'), sheetName=paste0('cluster',cl_idx))
        #} else {
            # xlsx::write.xlsx(df_export, file=paste0(outputDir_sub,'GO_analysis_co-correlations_v2.xlsx'), sheetName=paste0('cluster',cl_idx),append=TRUE)
            openxlsx::addWorksheet(wb = wb, sheetName = paste0(cl_idx))
            openxlsx::writeData(wb, sheet = paste0(cl_idx), df_export[[cl_idx]])
        #}
        
        # Create large data frame with all results (for 2nd export)
        df_complete_export <- rbind(df_complete_export, 
            data.frame(Gene_names=paste0('*** Cluster ',toString(cl_idx),' ***'),GO_terms=''), df_export[[cl_idx]])
    }
    #xlsx::write.xlsx(df_complete_export, file=paste0(outputDir_sub,'GO_analysis_co-correlations_v2.xlsx'), sheetName='all_clusters',append=TRUE)
    
    openxlsx::addWorksheet(wb = wb, sheetName = 'all')
    openxlsx::writeData(wb, sheet = 'all', df_complete_export)
    
    openxlsx::saveWorkbook(wb, file = paste0(regulon_object$outputDir_sub,'GO_analysis_co-correlations_v2.xlsx'), overwrite = TRUE)
    
    regulon_object$regulon_list = regulon_list
    
    return(regulon_object)
    
}

################################################################################
# Now create a function that exports the regulon information

# example: regulon_object = MW_regulon_add_TF_RL_flags(regulon_object, data_container)
MW_regulon_add_TF_RL_flags = function(regulon_object, data_container) {

    regulon_object$regulon_list_extended = list()
     
    for (idx in 1:length(regulon_object$the_regulons)) {
        
        df = data.frame(gene_name = regulon_object$the_regulons[[idx]])
            
        # identify transcription factors
        if (!is.null(data_container$applicable_TF_list)) {
            
            TF_searchresult = give_idxs_on_gene_list(df$gene_name,
                                                        data_container$applicable_TF_list)
            df$is_TF = TF_searchresult$boolean_query
            
        } else { warning('No transcription factor list was set.') }
        
        # identify ligands and receptors (if available)
        if (exists('ligand_receptor_pairs_df')) {
            LR_df = get_LRs(ligand_receptor_pairs_df, df$gene_name)
            df = cbind(df, LR_df)
        }
        
        regulon_object$regulon_list_extended[[idx]] = df
        
    }
    
    return(regulon_object)
       
}

# example: MW_regulon_final_export(regulon_object)
MW_regulon_final_export = function(regulon_object) {
        
    if (is.null(regulon_object$regulon_list_extended)) {
        stop('regulon_list_extended not available in regulon_list_extended, execute MW_regulon_add_TF_RL_flags first.')   
    }
    
    # Initialize data frame
    output_df = regulon_object$regulon_list_extended
    names(output_df) = paste0('regulon',1:length(output_df))
    
    # Also collect settings
    output_df$parameters = data.frame(p_or_r_cutoff=regulon_object$p_or_r_cutoff,
            connectedness_cutoff=regulon_object$connectedness_cutoff, 
            auto_cutoff1=regulon_object$auto_cutoff1, 
            auto_cutoff2=regulon_object$auto_cutoff2, 
            hierarchical_cutoff=regulon_object$hierarchical_cutoff,
            nCluster=regulon_object$nCluster)
    
    # export gnee assignments as one list
    output_df$regulon_assignments = data.frame(gene_name = names(regulon_object$gene_clustering_assign_hierarchical),
                                        assignment = regulon_object$gene_clustering_assign_hierarchical)
        
    # Prepare the dataframe
    openxlsx::write.xlsx(x = output_df, file = paste0(regulon_object$outputDir_sub, 'regulon_analysis_',regulon_object$analysis_name,'.xlsx'))

    # output
    print(paste0('output to: ',regulon_object$outputDir_sub, 'regulon_analysis_',regulon_object$analysis_name,'.xlsx'))
    
}

# 
# example: regulon_object = MW_regulon_reinitilize_outputdir(regulon_object, cfg$outputDir)
MW_regulon_reinitilize_outputdir = function(regulon_object, outputDir) {
    
    # Set and generate output dir
    outputDir_sub = paste0(outputDir, 'analysis_', regulon_object$analysis_name, '/regulon/')
    if (!dir.exists(outputDir_sub)) {dir.create(outputDir_sub)}
    regulon_object$outputDir_sub = outputDir_sub
    
    return(regulon_object)
    
}
    







