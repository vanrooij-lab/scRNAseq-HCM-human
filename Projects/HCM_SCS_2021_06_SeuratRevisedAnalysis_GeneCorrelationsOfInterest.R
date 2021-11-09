
################################################################################

# Generating correlations to genes of interest for each of the datasets
#
#

################################################################################
# Loading libraries and setting up parameters

print('Executing main script to load libraries')

# set script dir (note: overwritten by loading SeuratRevisedAnalysis below)
if (exists('LOCAL')) {
    script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Projects/'
} else {
    script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'
}

desired_command='dummy'
source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R'))
    # Also loads already
    # source(paste0(script_dir,'Functions/MW_general_functions.R'))
    # file.edit('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Functions/MW_general_functions.R')
        # Volcano plot functions (other copies of this function exist in other repos)

print('Starting correlation script ..')

library(Matrix) # required to transpose the sparse matrix, dgCMatrix
library(pheatmap)


# Set some info

#OBJECTS_TO_ANALYZE = c('ROOIJonly_RID2l', 'ROOIJonly_default_Int1c','HUonly_RID2l', 'HUonly_default_Int1c')
#INTEGRATED_OR_NOT = c('no'              , 'yes'                     ,'no'                 , 'yes')
#names(INTEGRATED_OR_NOT)=OBJECTS_TO_ANALYZE

OBJECTS_TO_ANALYZE = c('ROOIJonly_RID2l','HUonly_RID2l', 'TEICHMANNonly_RID2l', 'TEICHMANN.SP.only_RID2l')
INTEGRATED_OR_NOT = c('no'              , 'no'                     ,'no'      , 'no' )
names(INTEGRATED_OR_NOT)=OBJECTS_TO_ANALYZE
    # Note:
    # For the integrated dataset
    # currentSeuratObject_recombined3@assays$integrated@scale.data


################################################################################
# Load the objects

# !!DO NOT RUN THIS CODE UNLESS YOU WANT TO REGENERATE ALL CORRELATIONS!!

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
# Actually generate the correlations
# Note: I copy-pasted this into HPC's R command console, and saved the
# result into an Rdata file. (see "save()" below)

# Some TTN correlations
# To do: automate this such that it just loops over multiple datasets, patients, and the pooled/integrated data


# !!DO NOT RUN THIS CODE UNLESS YOU WANT TO REGENERATE ALL CORRELATIONS!!
if (T) {
    
    OBJECTS_TO_ANALYZE=OBJECTS_TO_ANALYZE # defined above; this is just a reminder
    GENES_OF_INTEREST_ = c('TTN','NPPA','CMYA5','XIRP2')
    
    # let's assume the genes of interest are among the features in each dataset 
    GENES_OF_INTEREST = shorthand_seurat_fullgenename(seuratObject = current_analysis$ROOIJonly_RID2l, gene_names = GENES_OF_INTEREST_)
    
    Volcano_df_collection=list()
    
    for (gene_name in GENES_OF_INTEREST) {
        
        # gene_name = 'ENSG00000163092:XIRP2'
        # gene_name = "ENSG00000164309:CMYA5"
        
        for (analysis_name in OBJECTS_TO_ANALYZE) {
            
            # analysis_name='ROOIJonly_RID2l'; gene_name=shorthand_seurat_fullgenename('TTN')
            # analysis_name='ROOIJonly_RID2l'; gene_name=shorthand_seurat_fullgenename(current_analysis$ROOIJonly_RID2l, 'NPPA')
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
                # print(Sys.time()) # corr function lists run time already
                
                # Get correlation df
                # Note we can take "data" here, since correlation is scaled by definition
                Volcano_df_collection[[gene_name]][[analysis_name]][[current_patient]] = 
                    get_volcano_df3(expression_matrix=current_analysis[[analysis_name]]@assays[[assayType]]@data[,cells_from_patient_selection],
                                    my_gene=gene_name,calc_qvals=T,no_expression_val=0,min_cell_expressed=.1,manual_expression=NULL)
                    
                # Make a plot
                p=plot_volcano3(my_corrs_df_current = Volcano_df_collection[[gene_name]][[analysis_name]][[current_patient]],mycex=3,
                                NRLABELED=10,mypvaltreshold=0.01,manual_gene_name = shorthand_cutname(gene_name),
                                  mypointsize=.1, mylinesize=.25, mytextsize=6)+
                    ggtitle(paste0('Correlations with ',gene_name,'\n(',analysis_name,'); ',current_patient,''))
                # p
                
                # Save & export it
                ggsave(filename = paste0(base_dir,'Rplots/',analysis_name,'_6_Volcano_',gene_name,'_',current_patient,'.pdf'), 
                    plot = p, height=172/3-4, width=172/3-4, units = 'mm', device = cairo_pdf)
                p=p+ggtitle(paste0('Correlations with ',shorthand_cutname(gene_name)))
                ggsave(filename = paste0(base_dir,'Rplots/',analysis_name,'_6_Volcano_',gene_name,'_',current_patient,'_v2.pdf'), 
                    plot = p, height=172/3-4, width=172/3-4, units = 'mm', device = cairo_pdf)
                
            }    
    
            # Export xlsx for all patients at once
            openxlsx::write.xlsx(x = Volcano_df_collection[[gene_name]][[analysis_name]], file = paste0(base_dir,'Rplots/',analysis_name,'_6_Volcano_',gene_name,'.xlsx'), overwrite=T)
            
            # export gene corr data per gene, per patient
            Volcano_df_gene_pt = list(); Volcano_df_gene_pt[[gene_name]][[analysis_name]] = Volcano_df_collection[[gene_name]][[analysis_name]]
            save(list = c('Volcano_df_gene_pt'), file = paste0(base_dir,'Rdata/Volcano_df_gene_pt__',gene_name,'_',analysis_name,'.Rdata'))
            rm('Volcano_df_gene_pt')
        }
        
    }#; beepr::beep()
    
    if (F) { # just to prevent accidental overwriting of stuff
        save(list = c('Volcano_df_collection'), file = paste0(base_dir,'Rdata/Volcano_df_collection__for_all.Rdata'))
        # load(file = paste0(base_dir,'Rdata/Volcano_df_collection__for_all.Rdata'))
    }
    
    # In case manual analysis was done separately, construct larger volcano using per-gene data
    if (F) {
        Volcano_df_collection=list()
        for (gene_name in GENES_OF_INTEREST) {
            load(file = paste0(base_dir,'Rdata/Volcano_df_gene_',gene_name,'.Rdata'))
            Volcano_df_collection[[gene_name]] = Volcano_df_gene[[gene_name]]
        }
        rm('Volcano_df_gene')
    }
    
    # Prettier plots also
    for (analysis_name in OBJECTS_TO_ANALYZE) {
        for (gene_name in GENES_OF_INTEREST) {
            patients_current_dataset = unique(names(Volcano_df_collection[[gene_name]][[analysis_name]]))
            for (current_patient in patients_current_dataset) {
                pvalslog10_=-log10(Volcano_df_collection[[gene_name]][[analysis_name]][[current_patient]]$pval.adj)
                p=
                    plot_volcano3(my_corrs_df_current = Volcano_df_collection[[gene_name]][[analysis_name]][[current_patient]],mycex=3,
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
    
}

################################################################################
# An additional nice summary plot would be where we show the average and CV
# in a plot
#
# Run this part to load the data for plotting purposes

SAVEPLOT=F

SP_SWITCH = '.SP.'
OBJECTS_TO_ANALYZE = c("ROOIJonly_RID2l",     "HUonly_RID2l",        paste0("TEICHMANN",SP_SWITCH,"only_RID2l"))

if (F) {
    
    #load(paste0(base_dir,'Rdata/Volcano_df_collection__for_all.Rdata'))
    Volcano_df_collection = list()
    for (CURRENT_DATASET in OBJECTS_TO_ANALYZE) {
        for (CURRENT_GENE in c('ENSG00000155657:TTN', 'ENSG00000175206:NPPA', "ENSG00000163092:XIRP2", "ENSG00000164309:CMYA5")) {
            load(file = paste0(base_dir,'Rdata/Volcano_df_gene_pt__',CURRENT_GENE,'_',CURRENT_DATASET,'.Rdata'))
            Volcano_df_collection[[CURRENT_GENE]][[CURRENT_DATASET]] = Volcano_df_gene_pt[[CURRENT_GENE]][[CURRENT_DATASET]]
        }
    }
    rm('Volcano_df_gene_pt')
    
    # Pretty comparison plots
    df_corr_collection = list()
    for (CURRENT_DATASET in OBJECTS_TO_ANALYZE) {
        for (CURRENT_GENE in c('ENSG00000155657:TTN', 'ENSG00000175206:NPPA', "ENSG00000163092:XIRP2", "ENSG00000164309:CMYA5")) {
    
        # CURRENT_DATASET = 'ROOIJonly_RID2l'; CURRENT_GENE = 'ENSG00000155657:TTN'
        # CURRENT_DATASET = 'ROOIJonly_RID2l'; CURRENT_GENE = 'ENSG00000175206:NPPA'
        # CURRENT_DATASET = 'TEICHMANNonly_RID2l'
        # CURRENT_DATASET = 'HUonly_RID2l'
        
        # Volcano_df_collection[[CURRENT_GENE]][[CURRENT_DATASET]]
        
        sel_idx = !grepl('pooled',names(Volcano_df_collection[[CURRENT_GENE]][[CURRENT_DATASET]])  )
        
        # First collect shared names
        all_rownames = lapply(Volcano_df_collection[[CURRENT_GENE]][[CURRENT_DATASET]][sel_idx], rownames)
        shared_genes = Reduce(intersect, all_rownames) 
        shared_genes = shared_genes[!(CURRENT_GENE==shared_genes)]
        
        # Bind data
        current_dfs_donorsNames_selGene =
            lapply(names(Volcano_df_collection[[CURRENT_GENE]][[CURRENT_DATASET]][sel_idx]),
                function(df_name) { df = Volcano_df_collection[[CURRENT_GENE]][[CURRENT_DATASET]][sel_idx][[df_name]][shared_genes,]
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
        if (SAVEPLOT) {
        ggsave(filename = paste0(base_dir,'Rplots/',CURRENT_DATASET,'_6_CorrelationsGenes_Means_',CURRENT_GENE,'_.pdf'), 
                    plot = p, height=80, width=46*2, units = 'mm', device=cairo_pdf) # 185/4
        }
        
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
            xlim(c(-1,1))+ylim(c(-myymax*.03,myymax*1.2))+theme_bw()+
            geom_text_repel(data=df_melted_sel[order(df_melted_sel$pval),][1:20,],aes(label=gene_name_short), color='red', max.overlaps = 100, min.segment.length = 0)
        # p
        if (SAVEPLOT) {
        ggsave(filename = paste0(base_dir,'Rplots/',CURRENT_DATASET,'_6_CorrelationsGenesColor_',CURRENT_GENE,'_.pdf'), 
                    plot = p, height=80, width=46*2, units = 'mm', device=cairo_pdf) # 185/4
        }
    
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
        if (SAVEPLOT) {
        ggsave(filename = paste0(base_dir,'Rplots/',CURRENT_DATASET,'_6_TableSignCorrelations_',CURRENT_GENE,'_.pdf'), 
                        plot = p, height=5.7*nrow_effective, width=ncol_effective*18, units = 'mm', device=cairo_pdf) # 185/4 || 46*2
        }
        
        # Now save for future use
        df_corr_collection$melt[[CURRENT_GENE]][[CURRENT_DATASET]] = df_melted
        df_corr_collection$melt_sel_sel[[CURRENT_GENE]][[CURRENT_DATASET]] = df_melted_sel_sel
        df_corr_collection$mean[[CURRENT_GENE]][[CURRENT_DATASET]] = df_correlations_mean
        df_corr_collection$SE[[CURRENT_GENE]][[CURRENT_DATASET]] = df_correlations_SE
        
        }
    }
            
}    

################################################################################
# More customized versions of Volcano plots in Rooij data for NPPA and XIRP2
# label:volcanoRooijGOI

# c('ENSG00000155657:TTN', 'ENSG00000175206:NPPA', "ENSG00000163092:XIRP2", "ENSG00000164309:CMYA5")

# Make a plot
TOPX=10
gene_name = 'ENSG00000175206:NPPA'
gene_name = 'ENSG00000163092:XIRP2'
analysis_name='ROOIJonly_RID2l'
current_patient='Rpooled'

for (idx in 1:2) {
    
    gene_name = c('ENSG00000175206:NPPA', 'ENSG00000163092:XIRP2')[idx]
    current_force = c(1000, 100)[idx]
    
    temp_vdf = Volcano_df_collection[[gene_name]][[analysis_name]][[current_patient]]
    temp_vdf = temp_vdf[order(temp_vdf$corr, decreasing = T),]
    custom_top_genes = c(temp_vdf$gene_name_short[2:(TOPX+1)],
        temp_vdf$gene_name_short[(dim(temp_vdf)[1]):(dim(temp_vdf)[1]-TOPX)])
    p=plot_volcano3(my_corrs_df_current = Volcano_df_collection[[gene_name]][[analysis_name]][[current_patient]],mycex=3,
                NRLABELED=10,mypvaltreshold=0.01,manual_gene_name = shorthand_cutname(gene_name),
                  #mypointsize=.1, mylinesize=.1, mytextsize=8, mylabelsize = 5, custom_highlight_group = custom_top_genes, myforce = 33, mydirection = 'y')+ 
                  mypointsize=.1, mylinesize=.1, mytextsize=8, mylabelsize = 5, custom_highlight_group = custom_top_genes, myforce = current_force, mydirection = 'both')+ 
                    ggtitle(paste0('Correlation with ',shorthand_cutname(gene_name)))+xlim(-.8,.8)
                    # ggtitle(paste0('Correlations with ',gene_name,'\n(',analysis_name,'); ',current_patient,''))
    p
    # Save & export it
    ggsave(filename = paste0(base_dir,'Rplots/',analysis_name,'_6_Volcano_',gene_name,'_',current_patient,'-customStyle3.pdf'), 
        plot = p, height=PANEL_WIDTH-4, width=PANEL_WIDTH-4, units = 'mm', device = cairo_pdf)
}

################################################################################
# Another overview plot; expand on what we did above
# This throws all datasets together

SP_SWITCH = '.SP.'
OBJECTS_TO_ANALYZE = c("ROOIJonly_RID2l",     "HUonly_RID2l",        paste0("TEICHMANN",SP_SWITCH,"only_RID2l"))

CUSTOM_PATIENT_ORDER = c('R.P1', 'R.P2', 'R.P3', 'R.P4', 'R.P5', 'H.N1', 'H.N2', 'H.N3', 'H.N4', 'H.N5', 'H.N13', 'H.N14', 'T.D1', 'T.D2', 'T.D3', 'T.D4', 'T.D5', 'T.D6', 'T.D7', 'T.D11', 'T.H2', 'T.H3', 'T.H4', 'T.H5', 'T.H6', 'T.H7')
REFERENCE_DATASET = 'ROOIJonly_RID2l'

# Similar to above, but slightly adjusted

# shorthand_seurat_fullgenename(current_analysis$ROOIJonly_RID2l, c('XIRP2', 'CMYA5'))

if (F) {
        
    gene_lists_customcorrelated = list()
    
    for (CURRENT_GENE in c('ENSG00000155657:TTN', 'ENSG00000175206:NPPA', "ENSG00000163092:XIRP2", "ENSG00000164309:CMYA5")) {
    #for (CURRENT_GENE in c('ENSG00000155657:TTN', "ENSG00000164309:CMYA5")) {
        
        # CURRENT_GENE = 'ENSG00000175206:NPPA'
        # CURRENT_GENE='ENSG00000155657:TTN'
        # CURRENT_GENE="ENSG00000164309:CMYA5"
        # CURRENT_GENE="ENSG00000163092:XIRP2"
        # First collect shared names (between the patients)
        #all_rownames = lapply(Volcano_df_collection[[REFERENCE_DATASET]][[CURRENT_GENE]][sel_idx], rownames)
        #shared_genes = Reduce(intersect, all_rownames) 
        #shared_genes = shared_genes[!(CURRENT_GENE==shared_genes)]
        
        # Bind data
        namesData=c("ROOIJonly_RID2l"="R",     "HUonly_RID2l"="H", "TEICHMANNonly_RID2l"="T", "TEICHMANN.SP.only_RID2l"="T")
        current_dfs_donorsNames_selGene_allData =
            # loops over datasets
            unlist(lapply(OBJECTS_TO_ANALYZE, function(current_dataset) {
                print(paste0('Object: ',current_dataset))
                # This creates a list of dfs that contain the volc dfs per patient
                sel_idx = !grepl('pooled',names(Volcano_df_collection[[CURRENT_GENE]][[current_dataset]])  )
                lapply(names(Volcano_df_collection[[CURRENT_GENE]][[current_dataset]][sel_idx]),
                    function(df_name) { 
                                print(paste0('Looking at df: ', df_name))
                                df = Volcano_df_collection[[CURRENT_GENE]][[current_dataset]][sel_idx][[df_name]]
                                df$donor = df_name
                                df$paper = namesData[current_dataset]
                                return(df)})
            }), recursive = F)
        df_melted = Reduce(rbind, current_dfs_donorsNames_selGene_allData)
        
        df_melted$pval.adj.sign = df_melted$pval.adj<0.05
        
        df_correlations_mean = 
            aggregate(df_melted[,c('corr','pval.adj')], by=list(gene_name=df_melted$gene_name), mean)
        df_correlations_mean$sign.donors = 
            aggregate(df_melted[,c('pval.adj.sign')], by=list(gene_name=df_melted$gene_name), sum)$x
        df_correlations_SE = 
            aggregate(df_melted[,c('corr','pval.adj')], by=list(gene_name=df_melted$gene_name), function(x) {SE = sqrt(var(x)/length(x))})
        
        # Calculate Rooij-specific significance counts and mean corr values
        # pval counts
        pval.adj.sign_temp = df_melted$pval.adj.sign
        pval.adj.sign_temp[!(df_melted$donor %in% paste0('R.P',1:5))]=0
        df_correlations_mean$sign.donors_rooij = 
            aggregate(pval.adj.sign_temp, by=list(gene_name=df_melted$gene_name), sum)$x
        # mean corrs
        corr_temp = df_melted$corr
        corr_temp[!(df_melted$donor %in% paste0('R.P',1:5))]=NA
        df_correlations_mean$corr_rooij = 
            aggregate(corr_temp, by=list(gene_name=df_melted$gene_name), mean, na.rm=T)$x
        
        
        df_correlations_SE$corr_min = df_correlations_mean$corr-df_correlations_SE$corr
        df_correlations_SE$corr_max = df_correlations_mean$corr+df_correlations_SE$corr
        
        df_correlations_SE$pval_min = df_correlations_mean$pval.adj-df_correlations_SE$pval.adj
        df_correlations_SE$pval_max = df_correlations_mean$pval.adj+df_correlations_SE$pval.adj
        
        df_correlations_SE$corr_mean = df_correlations_mean$corr
        df_correlations_SE$pval_mean = df_correlations_mean$pval.adj
        
        df_correlations_SE_subset = df_correlations_SE[df_correlations_SE$pval_mean<.05, ]
        
        df_correlations_SE_subset$gene=shorthand_cutname(df_correlations_SE_subset$gene_name)
        
        # Create gene of interest lists
        TOPX=25
        # Split for pos and neg corr, selection at least 2 sign rooij donors
        df_correlations_mean_Rpos = df_correlations_mean[df_correlations_mean$corr_rooij>0&df_correlations_mean$gene_name!=CURRENT_GENE&
                                                         df_correlations_mean$sign.donors_rooij>2,]
        df_correlations_mean_Rneg = df_correlations_mean[df_correlations_mean$corr_rooij<0&
                                                         df_correlations_mean$sign.donors_rooij>2,]
        # Sort appropriately
        df_correlations_mean_Rpos_sorted = 
            Reduce(rbind, lapply(5:0,
                   function(x) {df_=df_correlations_mean_Rpos[df_correlations_mean_Rpos$sign.donors_rooij==x,]
                                df_[order(df_$corr_rooij, decreasing = T),]}))
        df_correlations_mean_Rneg_sorted = 
            Reduce(rbind, lapply(5:0,
                   function(x) {df_=df_correlations_mean_Rneg[df_correlations_mean_Rneg$sign.donors_rooij==x,]
                                df_[order(abs(df_$corr_rooij), decreasing = T),]}))
        # 
        selected_genes_correlated=list()
        selected_genes_correlated$pos = df_correlations_mean_Rpos_sorted$gene_name[1:min(TOPX,nrow(df_correlations_mean_Rpos_sorted))]
        selected_genes_correlated$neg = df_correlations_mean_Rneg_sorted$gene_name[1:min(TOPX,nrow(df_correlations_mean_Rneg_sorted))]
    
        # Now export to later look at bulk expression levels per patient
        gene_lists_customcorrelated[[CURRENT_GENE]] = selected_genes_correlated
        
        # Mean correlations; note that all datasets are used here, so
        # this is not very representative (e.g. Teichmann is much more data, biases)
        if (nrow(df_correlations_SE_subset)>1) {
        p=ggplot(df_correlations_SE_subset, aes(x=corr_mean, y=-log10(pval_mean))) + 
            geom_errorbarh(aes(xmin = corr_min,xmax = corr_max),color='darkgrey') + 
            geom_errorbar(aes(ymin = -log10(pval_mean),ymax = -log10(pval_max)), color='darkgrey') +
            geom_point()+
            xlim(c(-1,1))+theme_bw()+
            geom_text_repel(data=df_correlations_SE_subset[df_correlations_SE_subset$pval_mean<1e-3, ], aes(label=gene), color='red', max.overlaps = 100, min.segment.length = 0) #+ylim(c())
        } else {p=ggplot(data.frame(x=1,y=1,t='non_found'))+geom_text(aes(x=x,y=y,label=t))}
        # p
        # ggsave(filename = paste0(base_dir,'Rplots/',CURRENT_DATASET,'_6_CorrelationsGenes_Means_',CURRENT_GENE,'_.pdf'), 
        #             plot = p, height=80, width=46*2, units = 'mm') # 185/4
        
        
        # Select genes based on ROOIJ mean plot
        selected_genes=df_corr_collection$mean[[CURRENT_GENE]][[REFERENCE_DATASET]][
            df_corr_collection$mean[[CURRENT_GENE]][[REFERENCE_DATASET]]$sign.donors>1,]$gene_name
        
        # Create log10 pval
        df_melted$log10_pval_=-log10(df_melted$pval.adj)
        myymax=max(df_melted$log10_pval_[is.finite(df_melted$log10_pval_)])
        if (any(is.infinite(df_melted$log10_pval_))) {inf_line=T} else {inf_line=F}
        df_melted$log10_pval_[is.infinite(df_melted$log10_pval_)]=myymax*1.1
        
        
        # Correlations per patients, color coded for paper
        df_melted_sel = df_melted[df_melted$gene_name %in% selected_genes,]
        
        # create the plot
        p=ggplot(df_melted_sel, aes(x=corr, y=log10_pval_, color=paper)) 
        if (inf_line) {p=p+geom_hline(yintercept = myymax*1.1, color='black', size=0.5, linetype='dotted')}
        p=p+
            geom_point()+
            #geom_errorbarh(aes(xmin = corr_min,xmax = corr_max)) + 
            #geom_errorbar(aes(ymin = -log10(pval_mean),ymax = -log10(pval_max))) +
            xlim(c(-1,1))+ylim(c(-myymax*.03,myymax*1.2))+theme_bw()+
            geom_text_repel(data=df_melted_sel[order(df_melted_sel$pval),][1:20,],aes(label=gene_name_short), color='red', max.overlaps = 100, min.segment.length = 0)+
            scale_color_manual(values = col_vector_60)
        # Note: multiple dots per patient: ugly
        # p
        #ggsave(filename = paste0(base_dir,'Rplots/',CURRENT_DATASET,'_6_CorrelationsGenesColor_',CURRENT_GENE,'_.pdf'), 
        #         plot = p, height=172/3-4, width=172/3-4, units = 'mm', device = cairo_pdf) # 185/4
        
        # Heatmap displaying significant correlations, using selection made above, based on 
        # whether they are significant in reference dataset (Rooij) in >1 patient
        # determing ordering
        for (negpos in c('pos','neg')){
        
            # negpos = 'pos'
            # negpos = 'neg'
            
            if(is.na(selected_genes_correlated[[negpos]])) {print(paste0('No genes that are ',negpos,' correlated to ',CURRENT_GENE,' '));next}
            
            # negpos='pos'
            
            # gene_order = df_correlations_mean[df_correlations_mean$gene_name %in% selected_genes,][order(df_correlations_mean[df_correlations_mean$gene_name %in% selected_genes,]$corr, decreasing = T),]$gene
            # gene_order_short = shorthand_cutname(gene_order)
            gene_order_short = shorthand_cutname(selected_genes_correlated[[negpos]])
            
            df_melted_sel2=df_melted[df_melted$gene_name %in% selected_genes_correlated[[negpos]],]
            
            ncol_effective=length(unique(df_melted_sel2$donor))
            nrow_effective=length(unique(df_melted_sel2$gene_name_short))
            
            # extra selection if necessary (because some plots show MANY genes)
            #if (nrow_effective>100) {
            #    df_melted_sel_sel=df_melted_sel[
            #        df_melted_sel$gene_name %in% df_correlations_mean[order(df_correlations_mean$sign.donors_rooij, decreasing = T),][1:100,]$gene_name,]
            #    nrow_effective=length(unique(df_melted_sel_sel$gene_name_short))
            #} else {df_melted_sel_sel=df_melted_sel}
            p=ggplot(df_melted_sel2, aes(x=factor(donor, levels=CUSTOM_PATIENT_ORDER), y=factor(gene_name_short, levels=rev(gene_order_short)), fill=corr, color=pval.adj.sign)) +
                geom_tile(size=.5, width=0.7, height=0.7)+
                scale_color_manual(values=c('white','black'))+
                scale_fill_gradientn(colours = rainbow_colors, breaks=seq(-1,1,.5), limits=c(-1,1))+
                geom_text(aes(label=round(corr,2)), color='#666666', size=6/.pt)+
                theme_bw()+ylab(element_blank())+give_better_textsize_plot(8)+
                ggtitle(paste0('Correlations with ',shorthand_cutname(CURRENT_GENE)))+
                theme(legend.position = 'none', axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
                xlab('donor')
            # p 
            
            ggsave(filename = paste0(base_dir,'Rplots/customALL',SP_SWITCH,'_6_TableSignCorrelations_',CURRENT_GENE,'_Corr',negpos,'_.pdf'), 
                             plot = p, height=2.5*nrow_effective, width=172, units = 'mm', device=cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5
            
            ##########
            # label:XIRP2corrHeatm
            # Slightly adjusted heatmap in case of more information
            df_melted_sel2$donor_beautified = df_melted_sel2$donor 
            df_melted_sel2$donor_beautified = gsub('^R\\.','HCM.',df_melted_sel2$donor_beautified); df_melted_sel2$donor_beautified=gsub('^H\\.','Ctrl1.',df_melted_sel2$donor_beautified); df_melted_sel2$donor_beautified=gsub('^T\\.','Ctrl2.',df_melted_sel2$donor_beautified)
            CUSTOM_PATIENT_ORDER_beautified = CUSTOM_PATIENT_ORDER
            CUSTOM_PATIENT_ORDER_beautified = gsub('^R\\.','HCM.',CUSTOM_PATIENT_ORDER_beautified); CUSTOM_PATIENT_ORDER_beautified=gsub('^H\\.','Ctrl1.',CUSTOM_PATIENT_ORDER_beautified); CUSTOM_PATIENT_ORDER_beautified=gsub('^T\\.','Ctrl2.',CUSTOM_PATIENT_ORDER_beautified)
            
            p=ggplot(df_melted_sel2, aes(x=factor(donor_beautified, levels=CUSTOM_PATIENT_ORDER_beautified), y=factor(gene_name_short, levels=rev(gene_order_short)), fill=corr)) +
                geom_tile(size=.5, width=1, height=1)+
                geom_point(data=df_melted_sel2[df_melted_sel2$pval.adj.sign, ])+
                scale_color_manual(values=c('white','black'))+
                scale_fill_gradientn(colours = rainbow_colors, breaks=seq(-1,1,.5), limits=c(-1,1))+
                #geom_text(aes(label=round(corr,2)), color='#666666', size=6/.pt)+
                theme_bw()+ylab(element_blank())+give_better_textsize_plot(8)+
                ggtitle(paste0('Correlations with ',shorthand_cutname(CURRENT_GENE)))+
                theme(legend.position = 'none', legend.key.height = unit(2,"mm"),
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
                xlab('donor')
            # p
            ggsave(filename = paste0(base_dir,'Rplots/customALL',SP_SWITCH,'_6_TableSignCorrelations-style2_',CURRENT_GENE,'_Corr',negpos,'_.pdf'), 
                             plot = p, height=min(172,3*nrow_effective+20), width=2/3*172-4, units = 'mm', device=cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5
            # Now with legend
            p=p+theme(legend.position='bottom') # p
            ggsave(filename = paste0(base_dir,'Rplots/customALL',SP_SWITCH,'_6_TableSignCorrelations-style2_LEGEND_',CURRENT_GENE,'_.pdf'), 
                             plot = p, height=min(172,3.5*nrow_effective+5), width=2/3*172-4, units = 'mm', device=cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5
            
            ##########
            # label:NPPAcorrHeatm
            # Now same with only Rooij data displayed (could have been made earlier, but OK)
            # Slightly adjusted heatmap in case of more information
            df_melted_sel2$donor_short = gsub('^.\\.','',df_melted_sel2$donor)
            CUSTOM_PATIENT_ORDER_short = gsub('^.\\.','',CUSTOM_PATIENT_ORDER)
            p=ggplot(df_melted_sel2[df_melted_sel2$paper=='R',], aes(x=factor(donor_short, levels=CUSTOM_PATIENT_ORDER_short), y=factor(gene_name_short, levels=rev(gene_order_short)), fill=corr)) +
                geom_tile(size=.5, width=1, height=1)+
                geom_point(data=df_melted_sel2[df_melted_sel2$pval.adj.sign&df_melted_sel2$paper=='R', ])+
                scale_color_manual(values=c('white','black'))+
                scale_fill_gradientn(colours = rainbow_colors, breaks=seq(-1,1,.5), limits=c(-1,1))+
                #geom_text(aes(label=round(corr,2)), color='#666666', size=6/.pt)+
                theme_bw()+ylab(element_blank())+give_better_textsize_plot(8)+
                ggtitle(paste0('Correlations with ',shorthand_cutname(CURRENT_GENE)))+
                theme(legend.position = 'none', legend.key.height = unit(2,"mm"),
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
                xlab('donor')
            # p
            ggsave(filename = paste0(base_dir,'Rplots/','ROOIJonly_RID2l','_6_TableSignCorrelations-style2_',CURRENT_GENE,'_Corr',negpos,'_.pdf'), 
                             plot = p, height=min(172,3*nrow_effective+20), width=1/3*172-4, units = 'mm', device=cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5
            # Now with legend
            p=p+theme(legend.position='bottom') # p
            ggsave(filename = paste0(base_dir,'Rplots/','ROOIJonly_RID2l','_6_TableSignCorrelations-style2_LEGEND_',CURRENT_GENE,'_.pdf'), 
                             plot = p, height=min(172,3.5*nrow_effective+5), width=2/3*172-4, units = 'mm', device=cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5
            
            
            # # Less sophisticated heatmap
            # # Note: this isn't 100 rows, because not all top-100 mean genes are also sign. in >1 patient
            # current_genes = gene_order_short[gene_order_short %in% df_melted_sel_sel$gene_name_short]
            # mtx <- matrix(NA, nrow=length(current_genes), ncol=length(unique(df_melted_sel_sel$donor)) )
            # dimnames(mtx) <- list( current_genes,  CUSTOM_PATIENT_ORDER[CUSTOM_PATIENT_ORDER %in% df_melted_sel_sel$donor] )
            # mtx[cbind(df_melted_sel_sel$gene_name_short, df_melted_sel_sel$donor)] <- df_melted_sel_sel$corr
            # 
            # mtx[is.na(mtx)]=0
            # 
            # p=pheatmap(mtx, fontsize_row = 3, cluster_rows = F, cluster_cols = F)
            # p
        }

    }
    
    save(list = 'gene_lists_customcorrelated', file = paste0(base_dir,'Rdata/gene_lists_customcorrelated__Rooijbased.Rdata'))
}

################################################################################

# custom analysis code to compare patients and datasets

# for (gene_name in GENES_OF_INTEREST) {

if (F) {
    
    load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/Rdata/Volcano_df_collection__for_all.Rdata')
    
    # current_gene = 'ENSG00000155657:TTN'
    # current_gene = 'ENSG00000175206:NPPA'
    current_gene="ENSG00000164309:CMYA5"
    # current_gene="ENSG00000163092:XIRP2"
        

    # Retrieves top 20s
    collected_genes = lapply(OBJECTS_TO_ANALYZE,
        function(current_dataset) {
            patients_current_dataset = names(Volcano_df_collection[[CURRENT_GENE]][[CURRENT_DATASET]])
            lapply(patients_current_dataset, function(current_patient) {
                Volcano_df_collection[[CURRENT_GENE]][[CURRENT_DATASET]][[current_patient]]$gene_name[order(Volcano_df_collection[[CURRENT_GENE]][[CURRENT_DATASET]][[current_patient]]$corr, decreasing = T)[1:20]]
            })
        })
    collected_genes=unique(unlist(collected_genes))
    
    # Retrieves correlations
    all_correlations_df=data.frame(gene=numeric(),corr=numeric(),patient=character(),dataset=character())
    for (current_dataset in OBJECTS_TO_ANALYZE) {
    
        # current_dataset='ROOIJonly_RID2l'; gene_name='TTN'; current_patient = 'R.P1'
        
        patients_current_dataset = names(Volcano_df_collection[[CURRENT_GENE]][[CURRENT_DATASET]])
        for (current_patient in patients_current_dataset) {

            
            current_data=
                data.frame(gene=collected_genes, 
                           corr=Volcano_df_collection[[CURRENT_GENE]][[CURRENT_DATASET]][[current_patient]][collected_genes,]$corr,
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
    ggsave(filename = paste0(base_dir,'Rplots/combined_TTN_corrs.pdf'), plot = p, height=30, width=15, units='cm', device=cairo_pdf)
}
################################################################################

# Old code/for manual running

if (F) {
    
    # NPPA Venn
    plot_Venn_MW_2lists_v3(Volcano_df_collection[['NPPA']][['ROOIJonly_RID2l']]$gene_name_short[2:31],
                            Volcano_df_collection[['NPPA']][['ROOIJonly_default_Int1c']]$gene_name_short[2:31],
                            name1='Rooij, RID2',name2='Rooij, Int')
    
    plot_Venn_MW_2lists_v3(Volcano_df_collection[['NPPA']][['ROOIJonly_default_Int1c']]$gene_name_short[2:31],
                            Volcano_df_collection[['NPPA']][['HUonly_default_Int1c']]$gene_name_short[2:31],
                            name1='Rooij, Int',name2='Hu, Int')
    
    Volcano_df_collection[['NPPA']][['ROOIJonly_default_Int1c']]$gene_name_short[2:31][Volcano_df_collection[['ROOIJonly_default_Int1c']][['NPPA']]$gene_name_short[1:30] %in% Volcano_df_collection[['HUonly_default_Int1c']][['NPPA']]$gene_name_short[2:31]]
    
    # TTN Venn
    plot_Venn_MW_2lists_v3(Volcano_df_collection[['TTN']][['ROOIJonly_RID2l']]$gene_name_short[2:31],
                            Volcano_df_collection[['TTN']][['ROOIJonly_default_Int1c']]$gene_name_short[2:31],
                            name1='Rooij, RID2',name2='Rooij, Int')
    
    plot_Venn_MW_2lists_v3(Volcano_df_collection[['TTN']][['ROOIJonly_default_Int1c']]$gene_name_short[2:31],
                            Volcano_df_collection[['TTN']][['HUonly_default_Int1c']]$gene_name_short[2:31],
                            name1='Rooij, Int',name2='Hu, Int')
}

################################################################################






