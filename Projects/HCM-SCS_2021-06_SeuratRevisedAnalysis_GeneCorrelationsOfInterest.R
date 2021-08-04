

# Generating correlations to genes of interest for each of the datasets
#
#

print('Executing main script to load libraries')

# set script dir (note: overwritten by loading SeuratRevisedAnalysis below)
if (exists('LOCAL')) {
    script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/'
} else {
    script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'
}

desired_command='dummy'
source(paste0(script_dir, 'HCM-SCS_2021-06_SeuratRevisedAnalysis_v2_UmiTools.R'))
    # Also loads already
    # source(paste0(script_dir,'Functions/MW_general_functions.R'))
    # file.edit('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Functions/MW_general_functions.R')
        # Volcano plot functions (other copies of this function exist in other repos)

print('Starting correlation script ..')

library(Matrix) # required to transpose the sparse matrix, dgCMatrix
library(pheatmap)


# For the integrated dataset
# currentSeuratObject_recombined3@assays$integrated@scale.data

################################################################################
# Load the objects

#OBJECTS_TO_ANALYZE = c('ROOIJonly_RID2l', 'ROOIJonly_default_Int1c','HUonly_RID2l', 'HUonly_default_Int1c')
#INTEGRATED_OR_NOT = c('no'              , 'yes'                     ,'no'                 , 'yes')
#names(INTEGRATED_OR_NOT)=OBJECTS_TO_ANALYZE

OBJECTS_TO_ANALYZE = c('ROOIJonly_RID2l','HUonly_RID2l', 'TEICHMANNonly_RID2l')
INTEGRATED_OR_NOT = c('no'              , 'no'                     ,'no')
names(INTEGRATED_OR_NOT)=OBJECTS_TO_ANALYZE

# For testing purposes
#OBJECTS_TO_ANALYZE = c('ROOIJonly_RID2l')
#INTEGRATED_OR_NOT = c('no'              )
#names(INTEGRATED_OR_NOT)=OBJECTS_TO_ANALYZE

current_analysis=list()
for (analysis_name in OBJECTS_TO_ANALYZE) {

    # analysis_name = 'ROOIJonly_RID2l'
    
    print(paste0('Loading ', analysis_name))
    
    current_analysis[[analysis_name]] = 
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',analysis_name,'.h5seurat'))

}

################################################################################
# Some TTN correlations
# To do: automate this such that it just loops over multiple datasets, patients, and the pooled/integrated data

OBJECTS_TO_ANALYZE=OBJECTS_TO_ANALYZE # defined above; this is just a reminder
GENES_OF_INTEREST_ = c('TTN','NPPA')

# let's assume the genes of interest are among the features in each dataset 
GENES_OF_INTEREST = shorthand_seurat_fullgenename(seuratObject = current_analysis[[1]], gene_names = GENES_OF_INTEREST_)

Volcano_df_collection=list()
for (analysis_name in OBJECTS_TO_ANALYZE) {
    
    for (gene_name in GENES_OF_INTEREST) {
        
        # analysis_name='ROOIJonly_RID2l'; gene_name=shorthand_seurat_fullgenename('TTN')
        # analysis_name='ROOIJonly_RID2l'; gene_name=shorthand_seurat_fullgenename('NPPA')
        # analysis_name='ROOIJonly_RID2l'; gene_name=shorthand_seurat_fullgenename('CMYA5')
        # analysis_name='TEICHMANNonly_RID2l'; gene_name='TTN'
        
        # Depending on whether we want to look at integrated data or not, Seurat has stored the
        # expression matrix of interest in a different spot
        assayType = switch(INTEGRATED_OR_NOT[[analysis_name]], 'yes'='integrated', 'no'='RNA')
        
        # Now go over patients
        patients_current_dataset = unique(current_analysis[[analysis_name]]$annotation_patient_str)
        for (current_patient_idx in 0:length(patients_current_dataset)) {
            
            if (current_patient_idx>0) {
                current_patient=patients_current_dataset[current_patient_idx]
                cells_from_patient_selection = (current_analysis[[analysis_name]]$annotation_patient_str == current_patient)
            }else{
                current_patient=paste0(substring(analysis_name, 1,1), 'pooled')
                cells_from_patient_selection=rep(T, length(current_analysis[[analysis_name]]$annotation_patient_str))
            }
            
            # Display progress
            print(paste0('Running analysis for: ',analysis_name,', ',gene_name,', ',current_patient,'..'))
            
            # Get correlation df
            # Note we can take "data" here, since correlation is scaled by definition
            Volcano_df_collection[[analysis_name]][[gene_name]][[current_patient]] = 
                get_volcano_df3(expression_matrix=current_analysis[[analysis_name]]@assays[[assayType]]@data[,cells_from_patient_selection],
                                my_gene=gene_name,calc_qvals=T,no_expression_val=0,min_cell_expressed=.1,manual_expression=NULL)
                
            # Make a plot
            p=plot_volcano3(my_corrs_df_current = Volcano_df_collection[[analysis_name]][[gene_name]][[current_patient]],mycex=3,
                            NRLABELED=20,mypvaltreshold=0.01,manual_gene_name = shorthand_cutname(gene_name),
                              mypointsize=.1, mylinesize=.25, mytextsize=10)+
                ggtitle(paste0('Correlations with ',gene_name,'\n(',analysis_name,'); ',current_patient,''))
            # p
            
            # Save & export it
            ggsave(filename = paste0(base_dir,'Rplots/',analysis_name,'_6_Volcano_',gene_name,'_',current_patient,'.png'), 
                plot = p, height=7.5, width=7.5, units = 'cm')
            
        }    

        # Export xlsx for all patients at once
        openxlsx::write.xlsx(x = Volcano_df_collection[[analysis_name]][[gene_name]], file = paste0(base_dir,'Rplots/',analysis_name,'_6_Volcano_',gene_name,'.xlsx'), overwrite=T)
    }
    
}#; beepr::beep()

save(list = c('Volcano_df_collection'), file = paste0(base_dir,'Rdata/Volcano_df_collection__for_all.Rdata'))

# Prettier plots also
for (analysis_name in OBJECTS_TO_ANALYZE) {
    for (gene_name in GENES_OF_INTEREST) {
        patients_current_dataset = unique(names(Volcano_df_collection[[analysis_name]][[gene_name]]))
        for (current_patient in patients_current_dataset) {
            pvalslog10_=-log10(Volcano_df_collection[[analysis_name]][[gene_name]][[current_patient]]$pval.adj)
            p=
                plot_volcano3(my_corrs_df_current = Volcano_df_collection[[analysis_name]][[gene_name]][[current_patient]],mycex=3,
                            NRLABELED=15,mypvaltreshold=0.01,manual_gene_name = shorthand_cutname(gene_name),
                              mypointsize=.1, mylinesize=.25, mytextsize=8, mylabelsize = 6)+
                              xlab('Correlation coefficient')+ylab('-log10(p)')+
                              ylim(c(-20, 1.1*max(pvalslog10_[is.finite(pvalslog10_)])))+
                                # ggtitle
                                    if (grepl(x = current_patient,pattern = 'pooled')[1]) {    
                                        ggtitle(paste0('Correlations with ',shorthand_cutname(gene_name)))
                                    } else {ggtitle(paste0('Correlations with ',shorthand_cutname(gene_name),'\n(',analysis_name,'); ',current_patient,''))}
            
            p
            ggsave(filename = paste0(base_dir,'Rplots/',analysis_name,'_6_VolcanoPrettier_',gene_name,'_',current_patient,'.pdf'), 
                plot = p, height=50, width=61, units = 'mm')
        }
    }
}

################################################################################
# An additional nice summary plot would be where we show the average and CV
# in a plot

load(paste0(base_dir,'Rdata/Volcano_df_collection__for_all.Rdata'))

for (CURRENT_DATASET in OBJECTS_TO_ANALYZE) {
    for (CURRENT_GENE in c('ENSG00000155657:TTN', 'ENSG00000175206:NPPA')) {

    # CURRENT_DATASET = 'ROOIJonly_RID2l'; CURRENT_GENE = 'ENSG00000155657:TTN'
    # CURRENT_DATASET = 'ROOIJonly_RID2l'; CURRENT_GENE = 'ENSG00000175206:NPPA'
    # CURRENT_DATASET = 'TEICHMANNonly_RID2l'
    # CURRENT_DATASET = 'HUonly_RID2l'
    
    # Volcano_df_collection[[CURRENT_DATASET]][[CURRENT_GENE]]
    
    sel_idx = !grepl('pooled',names(Volcano_df_collection[[CURRENT_DATASET]][[CURRENT_GENE]])  )
    
    # First collect shared names
    all_rownames = lapply(Volcano_df_collection[[CURRENT_DATASET]][[CURRENT_GENE]][sel_idx], rownames)
    shared_genes = Reduce(intersect, all_rownames) 
    shared_genes = shared_genes[!(CURRENT_GENE==shared_genes)]
    
    # Bind data
    current_dfs_donorsNames_selGene =
        lapply(names(Volcano_df_collection[[CURRENT_DATASET]][[CURRENT_GENE]][sel_idx]),
            function(df_name) { df = Volcano_df_collection[[CURRENT_DATASET]][[CURRENT_GENE]][sel_idx][[df_name]][shared_genes,]
                                df$donor = df_name
                                return(df)})
    df_melted = Reduce(rbind, current_dfs_donorsNames_selGene)
    
    df_melted$pval.adj.sign = df_melted$pval.adj<0.05
    
    df_correlations_mean = 
        aggregate(df_melted[,c('corr','pval.adj')], by=list(gene_name=df_melted$gene_name), mean)
    df_correlations_mean$sign.donors = 
        aggregate(df_melted[,c('pval.adj.sign')], by=list(gene_name=df_melted$gene_name), sum)$x
    df_correlations_SE = 
        aggregate(df_melted[,c('corr','pval.adj')], by=list(gene_name=df_melted$gene_name), function(x) {SE = sqrt(var(x)/length(x))})
    
    df_correlations_SE$corr_min = df_correlations_mean$corr-df_correlations_SE$corr
    df_correlations_SE$corr_max = df_correlations_mean$corr+df_correlations_SE$corr
    
    df_correlations_SE$pval_min = df_correlations_mean$pval.adj-df_correlations_SE$pval.adj
    df_correlations_SE$pval_max = df_correlations_mean$pval.adj+df_correlations_SE$pval.adj
    
    df_correlations_SE$corr_mean = df_correlations_mean$corr
    df_correlations_SE$pval_mean = df_correlations_mean$pval.adj
    
    df_correlations_SE_subset = df_correlations_SE[df_correlations_SE$pval_mean<.05, ]
    
    df_correlations_SE_subset$gene=shorthand_cutname(df_correlations_SE_subset$gene_name)
    
    #
    if (nrow(df_correlations_SE_subset)>1) {
    p=ggplot(df_correlations_SE_subset, aes(x=corr_mean, y=-log10(pval_mean))) + 
        geom_errorbarh(aes(xmin = corr_min,xmax = corr_max),color='darkgrey') + 
        geom_errorbar(aes(ymin = -log10(pval_mean),ymax = -log10(pval_max)), color='darkgrey') +
        geom_point()+
        xlim(c(-1,1))+theme_bw()+
        geom_text_repel(data=df_correlations_SE_subset[df_correlations_SE_subset$pval_mean<1e-3, ], aes(label=gene), color='red', max.overlaps = 100, min.segment.length = 0) #+ylim(c())
    } else {p=ggplot(data.frame(x=1,y=1,t='non_found'))+geom_text(aes(x=x,y=y,label=t))}
    # p
    ggsave(filename = paste0(base_dir,'Rplots/',CURRENT_DATASET,'_6_CorrelationsGenes_Means_',CURRENT_GENE,'_.pdf'), 
                plot = p, height=80, width=46*2, units = 'mm') # 185/4
    
    # Slightly adjusted plot (for NPPA, mostly)
    selected_genes=df_correlations_mean[df_correlations_mean$sign.donors>1,]$gene_name
    df_melted_sel = df_melted[df_melted$gene_name %in% selected_genes,]
    df_melted_sel$log10_pval_=-log10(df_melted_sel$pval.adj)
    myymax=max(df_melted_sel$log10_pval_[is.finite(df_melted_sel$log10_pval_)])
    if (any(is.infinite(df_melted_sel$log10_pval_))) {inf_line=T} else {inf_line=F}
    df_melted_sel$log10_pval_[is.infinite(df_melted_sel$log10_pval_)]=myymax*1.1
    # create the plot
    p=ggplot(df_melted_sel, aes(x=corr, y=log10_pval_, color=donor)) 
    if (inf_line) {p=p+geom_hline(yintercept = myymax*1.1, color='black', size=0.5, linetype='dotted')}
    p=p+
        geom_point()+
        #geom_errorbarh(aes(xmin = corr_min,xmax = corr_max)) + 
        #geom_errorbar(aes(ymin = -log10(pval_mean),ymax = -log10(pval_max))) +
        xlim(c(-1,1))+ylim(c(-myylim*.03,myymax*1.2))+theme_bw()+
        geom_text_repel(data=df_melted_sel[order(df_melted_sel$pval),][1:20,],aes(label=gene_name_short), color='red', max.overlaps = 100, min.segment.length = 0)
    # p
    ggsave(filename = paste0(base_dir,'Rplots/',CURRENT_DATASET,'_6_CorrelationsGenesColor_',CURRENT_GENE,'_.pdf'), 
                plot = p, height=80, width=46*2, units = 'mm') # 185/4

    # Heatmap displaying significant correlations
    # determing ordering
    gene_order = df_correlations_mean[df_correlations_mean$gene_name %in% selected_genes,][order(df_correlations_mean[df_correlations_mean$gene_name %in% selected_genes,]$corr, decreasing = T),]$gene
    gene_order_short = shorthand_cutname(gene_order)
    
    ncol_effective=length(unique(df_melted_sel$donor))
    nrow_effective=length(unique(df_melted_sel$gene_name_short))
    # extra selection if necessary (because some plots show MANY genes)
    if (nrow_effective>100) {
        df_melted_sel_sel=df_melted_sel[
            df_melted_sel$gene_name %in% df_correlations_mean[order(df_correlations_mean$sign.donors, decreasing = T),][1:100,]$gene_name,]
        nrow_effective=100
    } else {df_melted_sel_sel=df_melted_sel}
    p=ggplot(df_melted_sel_sel, aes(x=donor, y=factor(gene_name_short, levels=rev(gene_order_short)), fill=corr, color=pval.adj.sign)) +
        geom_tile(size=.5, width=0.7, height=0.7)+
        scale_color_manual(values=c('white','black'))+
        scale_fill_gradientn(colours = rainbow_colors, breaks=seq(-1,1,.5), limits=c(-1,1))+
        geom_text(aes(label=round(corr,2)), color='#666666', size=8/.pt)+
        theme_bw()+ylab(element_blank())+give_better_textsize_plot(8)+
        ggtitle(paste0('Correlations with ',shorthand_cutname(CURRENT_GENE)))+theme(legend.position = 'none')
    # p 
    
    ggsave(filename = paste0(base_dir,'Rplots/',CURRENT_DATASET,'_6_TableSignCorrelations_',CURRENT_GENE,'_.pdf'), 
                    plot = p, height=5.7*nrow_effective, width=ncol_effective*18, units = 'mm') # 185/4 || 46*2
    
    }
}
            


################################################################################

# custom analysis code to compare patients and datasets

# for (gene_name in GENES_OF_INTEREST) {

if (F) {
    
    load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/Rdata/Volcano_df_collection__for_all.Rdata')
    
    current_gene='TTN'

    collected_genes = lapply(OBJECTS_TO_ANALYZE,
        function(current_dataset) {
            patients_current_dataset = names(Volcano_df_collection[[current_dataset]][[current_gene]])
            lapply(patients_current_dataset, function(current_patient) {
                Volcano_df_collection[[current_dataset]][[current_gene]][[current_patient]]$gene_name[order(Volcano_df_collection[[current_dataset]][[current_gene]][[current_patient]]$corr, decreasing = T)[1:20]]
            })
        })
    collected_genes=unique(unlist(collected_genes))
    
    all_correlations_df=data.frame(gene=numeric(),corr=numeric(),patient=character(),dataset=character())
    for (current_dataset in OBJECTS_TO_ANALYZE) {
    
        # current_dataset='ROOIJonly_RID2l'; gene_name='TTN'; current_patient = 'R.P1'
        
        patients_current_dataset = names(Volcano_df_collection[[current_dataset]][[current_gene]])
        for (current_patient in patients_current_dataset) {

            
            current_data=
                data.frame(gene=collected_genes, 
                           corr=Volcano_df_collection[[current_dataset]][[current_gene]][[current_patient]][collected_genes,]$corr,
                           patient=current_patient, dataset=current_dataset)
                
            all_correlations_df=rbind(all_correlations_df,current_data)
            
        }
        
    }
    
    expand.grid()
    
    ggplot(all_correlations_df, aes(x=patient,y=gene,fill=corr))+
        geom_tile()

    mtx <- matrix(NA, nrow=length(collected_genes), ncol=length(unique(all_correlations_df$patient)) )
    dimnames(mtx) <- list( collected_genes, sort(unique(all_correlations_df$patient) ) )
    mtx[cbind(all_correlations_df$gene, all_correlations_df$patient)] <- all_correlations_df$corr
    
    mtx[is.na(mtx)]=0
    p=pheatmap(mtx, fontsize_row = 3)
    p
    ggsave(filename = paste0(base_dir,'Rplots/combined_TTN_corrs.pdf'), plot = p, height=30, width=15, units='cm')
}
################################################################################

# Old code/for manual running

if (F) {
    
    # NPPA Venn
    plot_Venn_MW_2lists_v3(Volcano_df_collection[['ROOIJonly_RID2l']][['NPPA']]$gene_name_short[2:31],
                            Volcano_df_collection[['ROOIJonly_default_Int1c']][['NPPA']]$gene_name_short[2:31],
                            name1='Rooij, RID2',name2='Rooij, Int')
    
    plot_Venn_MW_2lists_v3(Volcano_df_collection[['ROOIJonly_default_Int1c']][['NPPA']]$gene_name_short[2:31],
                            Volcano_df_collection[['HUonly_default_Int1c']][['NPPA']]$gene_name_short[2:31],
                            name1='Rooij, Int',name2='Hu, Int')
    
    Volcano_df_collection[['ROOIJonly_default_Int1c']][['NPPA']]$gene_name_short[2:31][Volcano_df_collection[['ROOIJonly_default_Int1c']][['NPPA']]$gene_name_short[1:30] %in% Volcano_df_collection[['HUonly_default_Int1c']][['NPPA']]$gene_name_short[2:31]]
    
    # TTN Venn
    plot_Venn_MW_2lists_v3(Volcano_df_collection[['ROOIJonly_RID2l']][['TTN']]$gene_name_short[2:31],
                            Volcano_df_collection[['ROOIJonly_default_Int1c']][['TTN']]$gene_name_short[2:31],
                            name1='Rooij, RID2',name2='Rooij, Int')
    
    plot_Venn_MW_2lists_v3(Volcano_df_collection[['ROOIJonly_default_Int1c']][['TTN']]$gene_name_short[2:31],
                            Volcano_df_collection[['HUonly_default_Int1c']][['TTN']]$gene_name_short[2:31],
                            name1='Rooij, Int',name2='Hu, Int')
}

################################################################################






