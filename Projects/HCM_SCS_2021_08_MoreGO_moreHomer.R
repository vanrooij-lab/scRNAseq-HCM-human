
# GO analysis
# Plus
# Further analysis of HOMER and GO analysis

source('/Users/m.wehrens/Documents/git_repos/SCS_Joep/Functions/functions_mw_copy/MW_GOKEGG_analysis.R')
# file.edit('/Users/m.wehrens/Documents/git_repos/SCS_Joep/Functions/functions_mw_copy/MW_GOKEGG_analysis.R')
    

################################################################################

ANALYSIS_NAME='ROOIJonly.sp.bt_RID2l'
if (!exists('current_analysis')) {current_analysis = list()}
current_analysis[[ANALYSIS_NAME]] =
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',ANALYSIS_NAME,'.h5seurat'))
ANALYSIS_NAME = 'ROOIJonly.sp.bt_RID2l_clExtended'
current_analysis[[ANALYSIS_NAME]] =
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',ANALYSIS_NAME,'.h5seurat'))

################################################################################

# Generating background tables for LISA and HOMER, symbol format

ANALYSIS_NAME = 'ROOIJonly.sp.bt_RID2l_clExtended'

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
# For running HOMER and LISA, see:
# - HCM_SCS_2021_08_HOMER.sh
# - HCM_SCS_2021_08_LISA.sh

# Note: MARGE wasn't used in the end.

################################################################################

# Run GO analyses for Rooij clusters
if (F) {
         
    ANALYSIS_NAME = 'ROOIJonly.sp.bt_RID2l_clExtended'
    load(paste0(base_dir,'Rdata/enriched_genes_lists_clusters_ROOIJ__ROOIJonly.sp.bt_RID2l_clExtended.Rdata')) # loads enriched_genes_lists_clusters_ROOIJ
    
    # Sanity check
    # load(file = paste0(base_dir,'Rdata/DE_cluster__','ROOIJonly.sp.bt_RID2l_clExtended','.Rdata')) # DE_cluster
    # enriched1_check_DE = rownames(DE_cluster$ROOIJonly.sp.bt_RID2l_clExtended$'1'[order(DE_cluster$ROOIJonly.sp.bt_RID2l_clExtended$'1'$avg_log2FC, decreasing = T),])[1:10]
    # enriched1_check_en = enriched_genes_lists_clusters_ROOIJ$`1`[1:10]
    # enriched1_check_DE==enriched1_check_en
 
    # Requires 
    # background_genes, all_genes
    # all_genes is used to generate the conversion table
    all_genes = shorthand_cutname( rownames(current_analysis[[ANALYSIS_NAME]]@assays$RNA@counts), PART1OR2 = 1 )
    # Background treshold: at least 1% of all cells 
    TRESHOLD_PCT = 0.05
    background_genes = shorthand_cutname(
        rownames(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data[
            rowMeans(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data>0)>TRESHOLD_PCT,]
            ), PART1OR2 = 1)
    # Export for additional testing
    write.table(x = background_genes, file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_background_table_GO_',TRESHOLD_PCT,'.txt'), 
                row.names = F, col.names = F, quote = F)
    save(list='background_genes', file=paste0(base_dir,'Rplots/GO_analysis__background_genes.Rdata'))
    # load(file=paste0(base_dir,'Rplots/GO_analysis__background_genes.Rdata')) # background_genes
    
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
    
    save(list = 'GO_all_clusters', file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'__GO_all_clusters.Rdata'))
    
    # Now create a summary parameter that we'll use for plotting
    TOPX=25
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
    
    # Export to xlsx
    GO_collected_top_terms_clusters = lapply(GO_clusters_summary_list, function(X) {data.frame(Term=X$Term, pval=paste0('10^-',round(-log10(X$Pvalue),1)))})
    GO_collected_top_terms_clusters$overview = Reduce(cbind, GO_collected_top_terms_clusters)
    # save to xlsx
    openxlsx::write.xlsx(x = GO_collected_top_terms_clusters, file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'__GO_all_clusters_summarized.xlsx'), overwrite = T)
    
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

################################################################################
# GO analysis on the ALL.SP clusters
# Bit redundant code with above, but did it like this because background genes
# should be generated at cluster

# Run GO analyses for ALL.SP clusters

# Part to run @ HPC
if (F) {
         
    # Load necessary data
    # DE_cluster__ALL.SP_btypSel_RID2l_clExtended.Rdata
    # enriched_genes_lists_clusters_ALL.SP__ALL.SP_btypSel_RID2l_clExtended.Rdata
    ANALYSIS_NAME = 'ALL.SP_btypSel_RID2l_clExtended'
    load(paste0(base_dir,'Rdata/enriched_genes_lists_clusters_ALL.SP__',ANALYSIS_NAME,'.Rdata'))
 
    if (!exists('current_analysis')) {current_analysis = list()}
    current_analysis[[ANALYSIS_NAME]] =
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',ANALYSIS_NAME,'.h5seurat'))

    # Requires 
    # background_genes, all_genes
    # all_genes is used to generate the conversion table
    all_genes = shorthand_cutname( rownames(current_analysis[[ANALYSIS_NAME]]@assays$RNA@counts), PART1OR2 = 1 )
    # Background treshold: at least 1% of all cells 
    TRESHOLD_PCT = 0.05
    background_genes = shorthand_cutname(
        rownames(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data[
            rowMeans(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data>0)>TRESHOLD_PCT,]
            ), PART1OR2 = 1)
    # Export for additional testing
    write.table(x = background_genes, file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_background_table_GO_',TRESHOLD_PCT,'.txt'), 
                row.names = F, col.names = F, quote = F)
    save(list='background_genes', file=paste0(base_dir,'Rplots/GO_analysis__',ANALYSIS_NAME,'__background_genes.Rdata'))
}

# Part to run locally
if (F) {
        
    # Load necessary data
    ANALYSIS_NAME = 'ALL.SP_btypSel_RID2l_clExtended'
    load(paste0(base_dir,'Rdata/enriched_genes_lists_clusters_ALL.SP__ALL.SP_RID2l_clExtended.Rdata'))
    load(paste0(base_dir,'Rplots/GO_analysis__',ANALYSIS_NAME,'__background_genes.Rdata')) # background_genes
    
    # Perform GO analysis on clusters
    # ==
    
    # Analysis
    config_GO=list(geneIdentifierType='ensembl_id', species='human')
    enriched_genes_lists_clusters_ALL.SP_clExt_ = shorthand_cutname_table( enriched_genes_lists_clusters_ALL.SP , PART1OR2 = 1 )
    termGSCnFRAME = generate_termGSCnFrame(config_GO, includeChildTerms=F)
    
    GO_all_clusters_ALL.SP = lapply(names(enriched_genes_lists_clusters_ALL.SP_clExt_), function(set_name) {
        print(paste0('Analyzing set ',set_name))
        current_genes_query = enriched_genes_lists_clusters_ALL.SP_clExt_[[set_name]]
        GO_out = analyzeGeneOntology_MW(config = config_GO, all_genes=all_genes, background_genes=background_genes, 
                            genes_query=current_genes_query, pathwayPCutoff=0.05, GOKegg='GO', 
                            includeChildTerms=F, termGSCnFRAME = termGSCnFRAME)
        return(GO_out)        
    })
    names(GO_all_clusters_ALL.SP)=names(enriched_genes_lists_clusters_ALL.SP_clExt_)
    
    save(list = 'GO_all_clusters_ALL.SP', file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'__GO_all_clusters.Rdata'))
    
    # Now create a summary parameter that we'll use for plotting
    TOPX=25
    GO_clusters_summary_list_ALL.SP = lapply(names(GO_all_clusters_ALL.SP), function(n) {x=GO_all_clusters_ALL.SP[[n]]; x$cluster=n; return(x[1:TOPX,])})
    GO_clusters_summary_ALL.SP = Reduce(f = rbind, x= GO_clusters_summary_list_ALL.SP)
    max_str_length = median(nchar(GO_clusters_summary_ALL.SP$Term)*2)
    GO_clusters_summary_ALL.SP$Term_short =
        sapply(GO_clusters_summary_ALL.SP$Term, function(x) {
                if (nchar(x)>max_str_length) {
                    paste0(substring(x, 1, max_str_length),' ...')
                } else {x}
            })
    # put back to get abbreviated names
    GO_clusters_summary_list_ALL.SP = split(x = GO_clusters_summary_ALL.SP, f = GO_clusters_summary_ALL.SP$cluster)
    
    # Export to xlsx
    GO_collected_top_terms_clusters_ALL.SP = lapply(GO_clusters_summary_list_ALL.SP, function(X) {data.frame(Term=X$Term, pval=paste0('10^-',round(-log10(X$Pvalue),1)))})
    GO_collected_top_terms_clusters_ALL.SP$overview = Reduce(cbind, GO_collected_top_terms_clusters_ALL.SP)
    # save to xlsx
    openxlsx::write.xlsx(x = GO_collected_top_terms_clusters_ALL.SP, file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'__GO_all_clusters_summarized.xlsx'), overwrite = T)
    
}

################################################################################
# GO analysis on the SCENIC regulons

# Note that this analysis eventually identifies regulons that are specific
# to the Van Rooij HCM data, but to create input data, technically, 
# the "ALL.SP_btypSel_RID2l_clExtended" dataset is used.

ANALYSIS_NAME_pooled = 'ALL.SP_btypSel_RID2l_clExtended'

if (F) {

    # Requires above to be partially run
    
    load(file = paste0(base_dir,'Rdata/SCENIC_regulons_core_genes_sel.Rdata')) # SCENIC_regulons_core_genes_sel
    load(file=paste0(base_dir,'Rplots/GO_analysis__background_genes.Rdata'))
    
    # Repetition of above
    
    # Analysis
    config_GO=list(geneIdentifierType='ensembl_id', species='human')
    SCENIC_regulons_core_genes_sel__ = lapply(SCENIC_regulons_core_genes_sel, function(X) {
        shorthand_seurat_fullgenename_faster(seuratObject = current_analysis$ROOIJonly.sp.bt_RID2l, gene_names = X)})
    SCENIC_regulons_core_genes_sel_ = shorthand_cutname_table( SCENIC_regulons_core_genes_sel__ , PART1OR2 = 1 )
    termGSCnFRAME = generate_termGSCnFrame(config_GO, includeChildTerms=F)
    
    GO_all_scenic_regs = lapply(names(SCENIC_regulons_core_genes_sel_), function(set_name) {
        print(paste0('Analyzing set ',set_name))
        current_genes_query = SCENIC_regulons_core_genes_sel_[[set_name]]
        GO_out = analyzeGeneOntology_MW(config = config_GO, all_genes=all_genes, background_genes=background_genes, 
                            genes_query=current_genes_query, pathwayPCutoff=0.05, GOKegg='GO', 
                            includeChildTerms=F, termGSCnFRAME = termGSCnFRAME)
        return(GO_out)        
    })
    names(GO_all_scenic_regs)=names(SCENIC_regulons_core_genes_sel)
    
    save(list = 'GO_all_scenic_regs', file = paste0(base_dir,'Rdata/',ANALYSIS_NAME_pooled,'__GO_all_scenic_regs.Rdata'))
        # load(file = paste0(base_dir,'Rdata/',ANALYSIS_NAME_pooled,'__GO_all_scenic_regs.Rdata')) # GO_all_scenic_regs
    
    # Now create a summary parameter that we'll use for plotting
    TOPX=100
    GO_scenic_regs_summary_list = lapply(names(GO_all_scenic_regs), function(n) {x=GO_all_scenic_regs[[n]]; x$regulon=n; return(x[1:min(TOPX,dim(x)[1]),])})
    GO_scenic_regs_summary = Reduce(f = rbind, x= GO_scenic_regs_summary_list)
    max_str_length = median(nchar(GO_scenic_regs_summary$Term)*2)
    GO_scenic_regs_summary$Term_short =
        sapply(GO_scenic_regs_summary$Term, function(x) {
                if (nchar(x)>max_str_length) {
                    paste0(substring(x, 1, max_str_length),' ...')
                } else {x}
            })
    # put back to get abbreviated names
    GO_scenic_regs_summary_list = split(x = GO_scenic_regs_summary, f = GO_scenic_regs_summary$regulon)
    
    # Export to xlsx
    #GO_collected_top_terms_scenic_regs = lapply(GO_scenic_regs_summary_list, function(X) {data.frame(Term=X$Term, pval=paste0('10^-',round(-log10(X$Pvalue),1)))})
    #GO_collected_top_terms_scenic_regs$overview = rbind(rep(c('TR', 'p-val'), length(GO_scenic_regs_summary_list)), GO_collected_top_terms_scenic_regs$overview)
    #colnames(GO_collected_top_terms_scenic_regs$overview)=rep(gsub('\\(\\+\\)','',names(GO_scenic_regs_summary_list)),each=2)
    GO_collected_top_terms_scenic_regs = lapply(GO_scenic_regs_summary_list, function(X) {
        rbind(data.frame(Line=paste0('(p=10^-',round(-log10(X$Pvalue),1),') ',X$Term)),
              data.frame(Line=rep('', TOPX-dim(X)[1]))) # manual padding
        })
    GO_collected_top_terms_scenic_regs$overview = Reduce(cbind, GO_collected_top_terms_scenic_regs)
    colnames(GO_collected_top_terms_scenic_regs$overview)=gsub('\\(\\+\\)','',names(GO_scenic_regs_summary_list))
    # save to xlsx
    openxlsx::write.xlsx(x = GO_collected_top_terms_scenic_regs, file = paste0(base_dir,'Rplots/',ANALYSIS_NAME_pooled,'__GO_all_scenic_regs_summarized.xlsx'), overwrite = T)
    
    library(beepr)
    beepr::beep()
}    


################################################################################
# GO analysis on our own regulons/modules

if (F) {

    ANALYSIS_NAME_ROOIJ = 'ROOIJonly.sp.bt_RID2l'
    
    # Requires above to be partially run
    
    #load(file = paste0(base_dir,'Rdata/SCENIC_regulons_core_genes_sel'))
    load(file=paste0(base_dir,'Rplots/GO_analysis__background_genes.Rdata'))
        # background_genes
    load(file = paste0(base_dir,'Rdata/',ANALYSIS_NAME_ROOIJ,'_core_regulons_sorted.Rdata'))
        # core_regulons_sorted
    
    # Repetition of above
    
    # Analysis
    config_GO=list(geneIdentifierType='ensembl_id', species='human')
    CUSTOM_regulons_core_genes_sel_ = shorthand_cutname_table( core_regulons_sorted , PART1OR2 = 1 )
    termGSCnFRAME = generate_termGSCnFrame(config_GO, includeChildTerms=F)
    
    GO_all_custom_regs = lapply(names(CUSTOM_regulons_core_genes_sel_), function(set_name) {
        print(paste0('Analyzing set ',set_name))
        current_genes_query = CUSTOM_regulons_core_genes_sel_[[set_name]]
        GO_out = analyzeGeneOntology_MW(config = config_GO, all_genes=all_genes, background_genes=background_genes, 
                            genes_query=current_genes_query, pathwayPCutoff=0.05, GOKegg='GO', 
                            includeChildTerms=F, termGSCnFRAME = termGSCnFRAME)
        return(GO_out)        
    })
    names(GO_all_custom_regs)=names(CUSTOM_regulons_core_genes_sel_)
    
    save(list = 'GO_all_custom_regs', file = paste0(base_dir,'Rdata/',ANALYSIS_NAME_ROOIJ,'__GO_all_custom_regs.Rdata')) # GO_all_custom_regs
    
    # Now create a summary parameter that we'll use for plotting
    TOPX=5
    GO_custom_regs_summary_list = lapply(names(GO_all_custom_regs), function(n) {x=GO_all_custom_regs[[n]]; x$regulon=n; return(x[1:TOPX,])})
    GO_custom_regs_summary = Reduce(f = rbind, x= GO_custom_regs_summary_list)
    max_str_length = median(nchar(GO_custom_regs_summary$Term)*2)
    GO_custom_regs_summary$Term_short =
        sapply(GO_custom_regs_summary$Term, function(x) {
                if (nchar(x)>max_str_length) {
                    paste0(substring(x, 1, max_str_length),' ...')
                } else {x}
            })
    # put back to get abbreviated names
    GO_custom_regs_summary_list = split(x = GO_custom_regs_summary, f = GO_custom_regs_summary$regulon)
    
    # Export to xlsx
    #GO_collected_top_terms_custom_regs = lapply(GO_custom_regs_summary_list, function(X) {data.frame(Term=X$Term, pval=paste0('10^-',round(-log10(X$Pvalue),1)))})
    #GO_collected_top_terms_custom_regs$overview = rbind(rep(c('TR', 'p-val'), length(GO_custom_regs_summary_list)), GO_collected_top_terms_custom_regs$overview)
    #colnames(GO_collected_top_terms_custom_regs$overview)=rep(gsub('\\(\\+\\)','',names(GO_custom_regs_summary_list)),each=2)
    GO_collected_top_terms_custom_regs = lapply(GO_custom_regs_summary_list, function(X) {data.frame(Line=paste0('(p=10^-',round(-log10(X$Pvalue),1),') ',X$Term))})
    GO_collected_top_terms_custom_regs$overview = Reduce(cbind, GO_collected_top_terms_custom_regs)
    colnames(GO_collected_top_terms_custom_regs$overview)=gsub('\\(\\+\\)','',names(GO_custom_regs_summary_list))
    # save to xlsx
    openxlsx::write.xlsx(x = GO_collected_top_terms_custom_regs, file = paste0(base_dir,'Rplots/',ANALYSIS_NAME_ROOIJ,'__GO_all_custom_regs_summarized.xlsx'), overwrite = T)
    
    library(beepr)
    beepr::beep()
}    




