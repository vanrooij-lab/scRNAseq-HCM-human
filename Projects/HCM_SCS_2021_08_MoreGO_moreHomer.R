
# GO analysis
# Plus
# Further analysis of HOMER and GO analysis

################################################################################

# Generating background tables for LISA and HOMER, symbol format

# Background treshold: at least 1% of all cells 
# Strict treshold
TRESHOLD_PCT_STRICT = 0.2
background_genes_strict = shorthand_cutname(
        rownames(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data[
            rowMeans(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data>0)>TRESHOLD_PCT_STRICT,]
            ), PART1OR2 = 2)
# Export for additional testing
write.table(x = background_genes_strict, file = paste0(base_dir,'GeneLists/',ANALYSIS_NAME,'_background_table_symbol_strict_',TRESHOLD_PCT_STRICT,'.txt'), 
                row.names = F, col.names = F, quote = F)

TRESHOLD_PCT = 0.05
background_genes = shorthand_cutname(
        rownames(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data[
            rowMeans(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data>0)>TRESHOLD_PCT,]
            ), PART1OR2 = 2)
# Export for additional testing
write.table(x = background_genes, file = paste0(base_dir,'GeneLists/',ANALYSIS_NAME,'_background_table_symbol_',TRESHOLD_PCT,'.txt'), 
                row.names = F, col.names = F, quote = F)


################################################################################

# Run GO analyses
if (F) {

    source('/Users/m.wehrens/Documents/git_repos/SCS_Joep/Functions/functions_mw_copy/MW_GOKEGG_analysis.R')
    # file.edit('/Users/m.wehrens/Documents/git_repos/SCS_Joep/Functions/functions_mw_copy/MW_GOKEGG_analysis.R')
    
    # Load clusters
    ANALYSIS_NAME = 'ROOIJonly_RID2l_clExtended'
    load(paste0(base_dir,'Rdata/enriched_genes_lists_clusters_ROOIJ__ROOIJonly_RID2l_clExtended.Rdata')) # loads enriched_genes_lists_clusters_ROOIJ
    
    # Requires 
    # background_genes, all_genes
    # all_genes is used to generate the conversion table
    all_genes = shorthand_cutname( rownames(current_analysis[[ANALYSIS_NAME]]@assays$RNA@counts), PART1OR2 = 1 )
    # Background treshold: at least 1% of all cells 
    TRESHOLD_PCT = 0.01
    background_genes = shorthand_cutname(
        rownames(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data[
            rowMeans(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data>0)>TRESHOLD_PCT,]
            ), PART1OR2 = 1)
    # Export for additional testing
    write.table(x = background_genes, file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_background_table_',TRESHOLD_PCT,'.txt'), 
                row.names = F, col.names = F, quote = F)
    
    # Perform GO analysis on clusters
    # ==
    
    # Analysis
    config_GO=list(geneIdentifierType='ensembl_id', species='human')
    enriched_genes_lists_clusters_ROOIJ_ = shorthand_cutname_table( enriched_genes_lists_clusters_ROOIJ , PART1OR2 = 1 )
    termGSCnFRAME = generate_termGSCnFrame(config_GO, includeChildTerms=F)
    
    GO_all_clusters = lapply(names(enriched_genes_lists_clusters_ROOIJ_), function(set_name) {
        print(paste0('Analyzing set ',set_name))
        current_genes_query = enriched_genes_lists_clusters_ROOIJ_[[set_name]]
        GO_out = analyzeGeneOntology_MW(config = config_GO, all_genes=all_genes, background_genes=background_genes, 
                            genes_query=current_genes_query, pathwayPCutoff=0.05, GOKegg='GO', 
                            includeChildTerms=F, termGSCnFRAME = termGSCnFRAME)
        return(GO_out)        
    })
    names(GO_all_clusters)=names(enriched_genes_lists_clusters_ROOIJ_)
    
    # Now create a summary parameter that we'll use for plotting
    TOPX=4
    GO_clusters_summary_list = lapply(names(GO_all_clusters), function(n) {x=GO_all_clusters[[n]]; x$regulon=n; return(x[1:TOPX,])})
    GO_clusters_summary = Reduce(f = rbind, x= GO_clusters_summary_list)
    max_str_length = median(nchar(GO_clusters_summary$Term)*2)
    GO_clusters_summary$Term_short =
        sapply(GO_clusters_summary$Term, function(x) {
                if (nchar(x)>max_str_length) {
                    paste0(substring(x, 1, max_str_length),' ...')
                } else {x}
            })
    # put back to get abbreviated names
    GO_clusters_summary_list = split(x = GO_clusters_summary, f = GO_clusters_summary$regulon)
    
    # Plot 1
    ggplot(GO_clusters_summary, aes(x=Term_short, y=-log10(Pvalue), fill=regulon))+
        geom_bar(stat='identity')+
        facet_grid(cols=vars(regulon), scales = 'free_x')+
        theme_bw()+give_better_textsize_plot(8)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none')
    
    # Some colors
    #some_cols = hue_pal()(6)
    #some_cols = brewer.pal(n = 6, name = 'Accent')
    some_cols = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(7)
    
    # Plot 2
    plot_list=
        lapply(1:length(GO_clusters_summary_list), function(i) {
            df=GO_clusters_summary_list[[i]]
            df$Term = factor(df$Term, levels=df[order(df$Pvalue),]$Term)
            df$Term_short = factor(df$Term_short, levels=df[order(df$Pvalue),]$Term_short)
            ggplot(df, aes(x=Term_short, y=-log10(Pvalue)))+
                geom_bar(stat='identity', fill=some_cols[i], color='black', size=.5)+
                theme_bw()+give_better_textsize_plot(4)+
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                      legend.position = 'none',
                      plot.margin = margin(1,1,0,0,'mm'))+
                ylab(element_blank())+xlab(element_blank())#+ylab('-log10(p)')+
        })
    p=wrap_plots(plot_list, nrow=1)
    p
    ggsave(filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_clusters_GO_terms.pdf'), plot = p, 
           width=184.6*2/3-4, height=80, units='mm')
    
    # Plot 3
    plot_list=
        lapply(GO_clusters_summary_list, function(df) {
            df$Term = factor(df$Term, levels=df[order(df$Pvalue),]$Term)
            df$Term_short = factor(df$Term_short, levels=df[order(df$Pvalue),]$Term_short)
            ggplot(df, aes(x=Term_short, y=-log10(Pvalue)))+
                geom_bar(stat='identity')+
                theme_bw()+give_better_textsize_plot(8)+
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none')+
                ylab(element_blank())+xlab(element_blank())+#+ylab('-log10(p)')+
                coord_flip()
        })
    wrap_plots(plot_list, ncol=1)
    
    # Plot 4
    plot_list=
        lapply(1:length(GO_clusters_summary_list), function(i) {
            df=GO_clusters_summary_list[[i]]
            df$Term = factor(df$Term, levels=df[order(df$Pvalue),]$Term)
            #df$Term_short = factor(df$Term_short, levels=rev(df[order(df$Pvalue),]$Term_short))
            ggplot(df, aes(x=factor(Term_short, levels=rev(Term_short)), y=1, size=-log10(Pvalue)))+
                geom_tile(width=1, height=1, fill='white')+#, aes(fill=-log10(Pvalue)))+
                geom_point(shape=21, color='black', fill=some_cols[i], )+
                scale_y_continuous(breaks = 1, expand = c(0,0)) + 
                #scale_x_discrete(expand = c(0,0))+
                theme_bw()+give_better_textsize_plot(6)+
                theme(# axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                      legend.position = 'none',
                      axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                      axis.ticks.y=element_blank(),
                      panel.border = element_blank(),
                      panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                      aspect.ratio = TOPX,
                      plot.margin = margin(1,1,0,0,'mm'))+
                xlab(element_blank())+ylab(element_blank())+
                #xlab(names(GO_clusters_summary_list)[i])+
                scale_size_continuous(range = c(0.25, 2))+
                coord_flip()
                    
                    #+ylab('-log10(p)')
                #coord_flip()+ylim(c(0,2))
        })
    p=wrap_plots(plot_list, ncol=1)
    p
    ggsave(filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_clusters_GO_terms_points.pdf'), plot = p, 
           width=172*2/3-4, height=80, units='mm', device = cairo_pdf)
    
    p=wrap_plots(plot_list, nrow=1)
    p
    ggsave(filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_clusters_GO_terms_points_turned.pdf'), plot = p, 
           width=172*2/3-4, height=50, units='mm', device = cairo_pdf)
    
    
    # ggplot(GO_clusters_summary, aes(x=factor(Term_short, levels=Term_short), y=1, size=-log10(Pvalue), color=regulon))+
    #     geom_point()+
    #     theme_minimal()+
    #     coord_flip()+
    #     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    #                   legend.position = 'none',
    #                   plot.margin = margin(1,1,0,0,'mm'))
    
}
