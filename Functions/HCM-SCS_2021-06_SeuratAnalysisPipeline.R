################################################################################
# FUNCTION DEFINITIONS TO RUN ANALYSES
# Now perform the analysis

# Let's always define these functions

# To test:
# mySeuratObject = RHL_SeuratObject_merged_noMito_sel

mySeuratAnalysis = function(mySeuratObject, run_name,
    normalization.method='LogNormalize', scale.factor=10000,
    do.scale=T,do.center=T,scale.max=10,features_to_use_choice='variable',
    remove_genes=T, cluster_resolution=1) {
    # type_features_to_use_PCA: variable or all
    
    # First remove genes that are only in <5 cells
    if (remove_genes) {
        gene_in_cell_count = rowSums(mySeuratObject@assays$RNA@counts>0)
        sel_genes = names(gene_in_cell_count)[gene_in_cell_count>4]
        mySeuratObject = subset(mySeuratObject, features= sel_genes)
        print(paste0('Genes removed: ', sum(gene_in_cell_count<5))) 
        print(paste0('Genes kept: ', sum(gene_in_cell_count>4))) 
    }

    # Collect some info
    all_genes = rownames(mySeuratObject)
    
    # Normalize data
    mySeuratObject <- NormalizeData(mySeuratObject, normalization.method = normalization.method, scale.factor = scale.factor) # values are defaults
    
    # Find variable features
    # We won't necessarily use them for all of the analyses
    mySeuratObject <- FindVariableFeatures(mySeuratObject, selection.method = "vst", nfeatures = 2000)

    # Identify the 10 most highly variable genes
    top30 <- head(VariableFeatures(mySeuratObject), 30)
    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(mySeuratObject)
    ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_1_VariableFeatures.png'), plot = plot1)
    plot2 <- LabelPoints(plot = plot1, points = top30, repel = TRUE, xnudge=0, ynudge=0)
    ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_1_VariableFeatures_labeled.png'), plot = plot2, width=15, height=15)

    # Determine which features to use
    features_to_use = if (features_to_use_choice=='all') { all_genes
                            } else {VariableFeatures(mySeuratObject)}

    # scale the data
    mySeuratObject <- ScaleData(mySeuratObject, features = features_to_use, do.scale=do.scale, do.center=do.center, scale.max=scale.max) # features = all_variable_features is default
    mySeuratObject <- RunPCA(object=mySeuratObject, npcs = 30, verbose = FALSE, features=features_to_use)
    mySeuratObject <- RunUMAP(mySeuratObject, reduction = "pca", dims = 1:30)
    mySeuratObject <- FindNeighbors(mySeuratObject, reduction = "pca", dims = 1:30)
    mySeuratObject <- FindClusters(mySeuratObject, resolution = cluster_resolution)
    
    # save(list = c('mySeuratObject'), file = paste0(base_dir,'Rdata/WRL_object_bigbig_sel_AnaDone.Rdata'))
    
    return(mySeuratObject)
    
}

mySeuratAnalysis_verybasic_part2only = function(mySeuratObject,
    do.scale=T,do.center=T,scale.max=10,features_to_use_choice='all', cluster_resolution=1, skip_scaling=F) {
    
    assayType = if ('integrated' %in% names(mySeuratObject@assays)) {'integrated'} else {'RNA'}
    
    features_to_use = if (features_to_use_choice=='all') { rownames(mySeuratObject@assays$RNA@data)
        } else {VariableFeatures(mySeuratObject)}
    
    if (!skip_scaling) {
        mySeuratObject <- ScaleData(mySeuratObject, features = features_to_use, do.scale=do.scale, do.center=do.center, scale.max=scale.max) # features = all_variable_features is default
    }
    mySeuratObject <- RunPCA(object=mySeuratObject, npcs = 30, verbose = FALSE, features=features_to_use)
    mySeuratObject <- RunUMAP(mySeuratObject, reduction = "pca", dims = 1:30)
    mySeuratObject <- FindNeighbors(mySeuratObject, reduction = "pca", dims = 1:30)
    mySeuratObject <- FindClusters(mySeuratObject, resolution = cluster_resolution)
    return(mySeuratObject)
}


# e.g. mySeurat_genenames(current_analysis$HUonly, c('TTN','NPPA'))
mySeurat_genenames = function(mySeuratObject, gene_names, complete=T) {
    
    if (complete) {
        search_pattern = paste0(':',gene_names, '$')
    } else { search_pattern = gene_names }
    
    extended_names = sapply(search_pattern, function(current_search) { 
        hits = grepl(pattern = current_search, x = rownames(mySeuratObject))
        if (sum(hits)==1) { return(rownames(mySeuratObject)[hits]) } else { return(NA) }
    } ) 
    
    return(extended_names)
    
}

mySeuratCommonPlots = function(mySeuratObject, 
    run_name,
    mymarkers = c('MALAT1', 'TTN', 'MYH7', 'MYH6', 'NPPA', 'NPPB', 'ACTA1','MYL2','SORBS2','CSRP3','NDUFA4','CRYAB','HSPB1', 'KCNQ1OT1'),
    mypointsize=NULL,
    add_custom_fields=NULL) {    
    
    mymarkers_ext = mySeurat_genenames(mySeuratObject, mymarkers)
    
    print(paste0('Markers not found: ', paste0(mymarkers[is.na(mymarkers_ext)], collapse = ', ')))
    
    clusterplus1mapping = 1:length(levels(Idents(mySeuratObject)))
    names(clusterplus1mapping) = as.character(levels(Idents(mySeuratObject)))
    mySeuratObject[['Seurat_Clusters_plus1']] = plyr::revalue(Idents(mySeuratObject), clusterplus1mapping)
    
    # Show umap with annotations
    for (current_annotation in c('annotation_sample_str','annotation_patient_str','annotation_paper_str','annotation_region_str','ident','Seurat_Clusters_plus1', add_custom_fields)) {
        
        # create labeled and unlabeled version
        p=DimPlot(mySeuratObject, group.by = current_annotation, cols = rep(col_vector_60,2), label = T, repel = T, label.size = 7/.pt, pt.size = mypointsize)+
            theme_void()+ggtitle(element_blank())+theme(legend.position = 'none')
        # p
        p_nl=DimPlot(mySeuratObject, group.by = current_annotation, cols = rep(col_vector_60,2))
        # remove legend if it's too large (prevents errors)
        if(!(current_annotation=='ident')) {if (dim(unique(mySeuratObject[[current_annotation]]))[1]>30) {
         p=p+theme(legend.position = 'none') # now redundant
         p_nl=p_nl+theme(legend.position = 'none')
        }}
        
        # test
        ggsave(plot = p, filename = paste0(base_dir,'Rplots/',run_name,'_2_umapLabeled_by_',current_annotation,'_style3.pdf'), height=172/3-4, width=172/3-4, units='mm')
        
        # Save files
        ggsave(plot = p, filename = paste0(base_dir,'Rplots/',run_name,'_2_umapLabeled_by_',current_annotation,'.pdf'), height=5.5, width=5.5, units='cm')
        ggsave(plot = p_nl, filename = paste0(base_dir,'Rplots/',run_name,'_2_umap_by_',current_annotation,'.pdf'), height=7, width=7, units='cm')
        
        # More customized/stylized version, specific size
        if(!(current_annotation=='ident')) {current_annotation='Seurat_Clusters_plus1'}
        p=DimPlot(mySeuratObject, group.by = current_annotation, cols = rep(col_vector_60,2), label = T, repel = T, label.size = 6/.pt, pt.size = 1, label.box=T)+
                theme_void()+ggtitle(element_blank())+theme(legend.position = 'none')
        # p
        ggsave(plot = p, filename = paste0(base_dir,'Rplots/',run_name,'_2_umapLabeled_by_',current_annotation,'_morecustomized.pdf'), height=172/3-4, width=172/3-4, units='mm', device = cairo_pdf)
        
    }
    
    for (marker in mymarkers_ext[!is.na(mymarkers_ext)]) {
        
        # marker = 'KCNQ1OT1'
        # marker = 'ACTC1'
        # marker = 'ENSG00000120937:NPPB'
        
        # marker projection on umaps
        p_cm = FeaturePlot(mySeuratObject, features = marker, cols = rainbow_colors)+theme_void()+ggtitle(element_blank())+theme(legend.position = 'none')
        ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_2_umap_markers_',marker,'.png'), plot = p_cm, height=5, width=5)
        
        # Violins
        pViol_m = VlnPlot(object = mySeuratObject, features = marker, group.by = 'annotation_paper_str') #, group.by = 'from_paper')
        ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_4_Violin_markers_',marker,'.png'), plot = pViol_m, height=7.5, width=7.5)

        # Also create histograms
        pHist=ggplot(data.frame(expression=mySeuratObject@assays[[mySeuratObject@active.assay]]@data[marker,], 
                                source=mySeuratObject$annotation_paper_fct))+
            geom_histogram(aes(x=expression, fill=source, after_stat(density)))+
            facet_grid(rows='source')+theme_bw()+theme(legend.position = 'none', )+ggtitle(marker)+give_better_textsize_plot(10)
        ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_4_histogram_markers_',marker,'.png'), plot = pHist, height=10, width=7.5, units='cm')
        
        # More complicated histogram, split, x-lims @98% of data, such that shape of curve is visible in case of outliers
        currentlims=list()
        for (source in unique(mySeuratObject$annotation_paper_fct)) {
            currentlims[[source]]=calc_limits(mySeuratObject@assays[[mySeuratObject@active.assay]]@data[marker,mySeuratObject$annotation_paper_fct==source], percentile = .02)
        }
        # then calculate breaks
        max_val = max(sapply(currentlims, function(x) {x[2]}))
        currentbreaks = seq(from=-(max_val/29)/2,by=(max_val/29),to=max_val+(max_val/29)/2) 
        # And make histogram, use density as y-value
        pHist=ggplot(data.frame(expression=mySeuratObject@assays[[mySeuratObject@active.assay]]@data[marker,],
                           source=mySeuratObject$annotation_paper_fct),
                     )+
            geom_histogram(aes(x=expression, y=..density.., fill=source), breaks=currentbreaks)+
            give_better_textsize_plot(10)+
            facet_grid(rows='source')+theme_bw()+theme(legend.position = 'none')+ggtitle(marker)
        ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_4_histogram98Lim_markers_',marker,'.png'), plot = pHist, height=10, width=7.5, units='cm')

        
    }
    # Add one plot that has legend
    p_cm = p_cm+theme(legend.position = 'bottom')
    ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_2_umap_markers_',marker,'_LEGEND.png'), plot = p_cm, height=5, width=5)
    
    # Custom plot that shows distribution of patients over clusters
    # mymaxy=1.5*max(table(Idents(mySeuratObject)))
    p=ggplot(data.frame( cluster = Idents(mySeuratObject)+1,
                    Donor = mySeuratObject$annotation_patient_fct))+
        geom_bar(aes(x=cluster, fill=Donor))+theme_bw()+
        xlab('Cluster')+ylab('Number of cells')+
        give_better_textsize_plot(8)+
        theme(legend.position = 'right', legend.key.size = unit(3, "mm"))
        #theme(legend.position=c(.99,.99), legend.justification = c(1,1), legend.key.size = unit(3, "mm"))
        # ylim(c(0,mymaxy)); p
    # p
    ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_5_Barplot_PatientCluster_distr.pdf'), 
        plot = p, height=172/3-4, width=172/3-4, units='mm', device = cairo_pdf)
    
    # Custom plot that shows distribution of patients over clusters
    # mymaxy=1.5*max(table(Idents(mySeuratObject)))
    Donor2=gsub(pattern = '^R\\.|^H\\.|^T\\.', replacement = '', mySeuratObject$annotation_patient_str)
    mySeuratObject$annotation_patient_fct_NS = factor(Donor2, sort(unique(Donor2)))
    p=ggplot(data.frame( cluster = mySeuratObject$Seurat_Clusters_plus1,
                    Donor = mySeuratObject$annotation_patient_fct_NS))+
        geom_bar(aes(x=cluster, fill=Donor))+theme_bw()+
        xlab('Cluster')+ylab('Number of cells')+
        give_better_textsize_plot(8)+
        theme(legend.position = 'right', legend.key.size = unit(3, "mm"))
        #theme(legend.position=c(.99,.99), legend.justification = c(1,1), legend.key.size = unit(3, "mm"))
        # ylim(c(0,mymaxy)); p
    # p
    ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_5_Barplot_PatientCluster_distr_DonorNoSet.pdf'), 
        plot = p, height=172/3-4, width=172/3-4, units='mm', device = cairo_pdf)
    
    
    
}

diff_express_clusters = function(mySeuratObject, mc.cores=8, custom_ident=NULL) {
    
    # set custom clustering gruop if desired
    if (!is.null(custom_ident)) {
        Idents(mySeuratObject) <- custom_ident
    }
    
    all_markers = 
        mclapply(X =  levels(Idents(mySeuratObject)), 
            FUN = function(cluster_index) {
                
                print(paste0('Initiating analysis of cl ',cluster_index))
                
                # find all markers of cluster X
                current_markers <- FindMarkers(mySeuratObject, ident.1 = cluster_index, min.pct = .05)
                
                # done
                print(paste0("Cluster ",cluster_index," DE done"))
                return(current_markers)
            },
            mc.cores = mc.cores)
    names(all_markers) = levels(Idents(mySeuratObject))
    
    return(all_markers)
    
}

# note things are saved like
# load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/Rdata/DE_cluster_ROOIJonly_RID2l.Rdata')
# all_markers=DE_cluster[[ANALYSIS_NAME]]
diff_express_clusters_save_results = function(all_markers, run_name, base_dir, topX=10, easy_names=T, save=T, pval_cutoff=0.05, extendedOutput=F, FC_cutoff=1.2) {
     
    # ****
    
    # Now collect the top-10 (log2fold) genes for each
    topHitsPerCluster=
        data.frame(lapply(all_markers, function(x) {
            x_sel = x[x$p_val_adj<pval_cutoff&x$avg_log2FC>0,]
            c(  rownames(x_sel[order(x_sel$avg_log2FC, decreasing = T),][1:min(length(x_sel$avg_log2FC), topX),]), 
                rep(NA, topX-min(length(x_sel$avg_log2FC), topX)))
            } ))
    colnames(topHitsPerCluster)=names(all_markers)
    if (easy_names) {topHitsPerCluster=shorthand_cutname_table(topHitsPerCluster)}
    # and save
    if (save) {openxlsx::write.xlsx(x = as.data.frame(topHitsPerCluster), file = paste0(base_dir,'Rplots/ClusterTopHits_',run_name,'.xlsx'), overwrite=T)}
    
    # ****
    
    if (extendedOutput) {
        
        # SAME table as above BUT INCLUDING FC
        # Also return an overview that includes FC
        # Now collect the top-10 (log2fold) genes for each
        topHitsPerCluster_FC=
            Reduce(cbind, lapply(all_markers, function(x) {
                x_sel = x[x$p_val_adj<pval_cutoff&x$avg_log2FC>0,]
                data.frame(
                    GENE=c(  shorthand_cutname (  rownames(x_sel[order(x_sel$avg_log2FC, decreasing = T),][1:min(length(x_sel$avg_log2FC), topX),])   ), 
                        rep(NA, topX-min(length(x_sel$avg_log2FC), topX))),
                    FC=c(  sprintf("%.2f",    round(2^(x_sel[order(x_sel$avg_log2FC, decreasing = T),][1:min(length(x_sel$avg_log2FC), topX),]$avg_log2FC),2)     )     , 
                        rep(NA, topX-min(length(x_sel$avg_log2FC), topX)))
                )
                } ))
        colnames(topHitsPerCluster_FC)=rep(names(all_markers), each=2)
        # and save
        if (save) {openxlsx::write.xlsx(x = as.data.frame(topHitsPerCluster_FC), file = paste0(base_dir,'Rplots/ClusterTopHits_',run_name,'_plusFC.xlsx'), overwrite=T)}
        
        # SAME table as above BUT INCLUDING FC, AND MOST NEGATIVE GENES
        # Also return an overview that includes FC
        # Now collect the top-10 (log2fold) genes for each
        topHitsPerCluster_FC_neg=
            Reduce(cbind, lapply(all_markers, function(x) {
                x_sel = x[x$p_val_adj<pval_cutoff&x$avg_log2FC<0,] # !!
                data.frame(
                    GENE=c(  shorthand_cutname (  rownames(x_sel[order(x_sel$avg_log2FC, decreasing = F),][1:min(length(x_sel$avg_log2FC), topX),])   ), 
                        rep(NA, topX-min(length(x_sel$avg_log2FC), topX))),
                    FC=c(  sprintf("%.2f",    round(2^(x_sel[order(x_sel$avg_log2FC, decreasing = F),][1:min(length(x_sel$avg_log2FC), topX),]$avg_log2FC),2)     )     , 
                        rep(NA, topX-min(length(x_sel$avg_log2FC), topX)))
                )
                } ))
        colnames(topHitsPerCluster_FC_neg)=rep(names(all_markers), each=2)
        # and save
        if (save) {openxlsx::write.xlsx(x = as.data.frame(topHitsPerCluster_FC_neg), file = paste0(base_dir,'Rplots/ClusterTopHits_',run_name,'_plusFC_DOWNgenes.xlsx'), overwrite=T)}
        
        
        # ****
        
        # Enriched genes list based on FC_cutoff and pval_cutoff
        enriched_genes_lists=list(); enriched_genes_lists_down=list()
        for (subset_name in names(all_markers)) { 
            
            marker_df = all_markers[[subset_name]]
            
            # Determine selection for genes up and down resp.
            x_sel      = marker_df[marker_df$p_val_adj<pval_cutoff & 2^(marker_df$avg_log2FC)>FC_cutoff,]
            x_sel_down = marker_df[marker_df$p_val_adj<pval_cutoff & 2^(marker_df$avg_log2FC)<(1/FC_cutoff),]
            
            current_genes      = rownames(x_sel[order(x_sel$avg_log2FC, decreasing = T),])
            current_genes_down = rownames(x_sel_down[order(x_sel_down$avg_log2FC, decreasing = F),])
            
            print(paste0('Genes for cl.',subset_name,':  ',length(current_genes)))
            print(paste0('Genes down for cl.',subset_name,':  ',length(current_genes_down)))
                    
            write.table(x=data.frame(gene=shorthand_cutname(current_genes)), file=paste0(base_dir,'GeneLists/ClusterHits_',run_name,'_table_cl',subset_name,'.txt'), 
                        quote = F, row.names = F, col.names = F)
            write.table(x=data.frame(gene=shorthand_cutname(current_genes_down)), file=paste0(base_dir,'GeneLists/ClusterHitsDown_',run_name,'_table_cl',subset_name,'.txt'), 
                        quote = F, row.names = F, col.names = F)
            
            enriched_genes_lists[[subset_name]] = current_genes
            enriched_genes_lists_down[[subset_name]] = current_genes_down
        }
        
        # Now also produce lists of genes of interest
    }
    
    
    # Also save the whole thing to Excel
    all_markers = lapply(all_markers, function(df) {df$gene = rownames(df); df$gene_short = shorthand_cutname(df$gene); df = df[order(df$avg_log2FC,decreasing = T),]; return(df)})
    if (save) {openxlsx::write.xlsx(x = all_markers, file = paste0(base_dir,'Rplots/Cluster_DE_full_',run_name,'.xlsx'), overwrite=T)}
    
    if (extendedOutput) {
        return(list(topHitsPerCluster=topHitsPerCluster, enriched_genes_lists=enriched_genes_lists, enriched_genes_lists_down=enriched_genes_lists_down))
    } else {
        return(topHitsPerCluster)
    }
}


# This function produces a Batch-corrected Seurat object according to our 
# settings.
# 
# (Note to self: See also dev folder, this was tested preliminary under the name "C1b";
# C1b, now including all genes, but correcting only for the patients.)
#
# This follows the tutorial, section "Performing integration on datasets normalized with SCTransform"
# https://satijalab.org/seurat/articles/integration_introduction.html
# Introduction to scRNA-seq integration
# Compiled: June 14, 2021
# Source: vignettes/integration_introduction.Rmd
#
Seurat_batch_correction_1c = function(set_name, base_dir, MinCellsRequired=50) {

    # Load uncorrected pre-processed analysis object 
    current_analysis=list()
    current_analysis[[set_name]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',set_name,'.h5seurat'))
    
    # Split object and apply SCT transformation
    currentSeuratObject_list3 <- SplitObject(current_analysis[[set_name]], split.by = "annotation_patient_str")
    
    # Count number of cells in each, and select only ones with enough cells
    nrCellsPerPatient = sapply(currentSeuratObject_list3, ncol)
    currentSeuratObject_list3 = currentSeuratObject_list3[nrCellsPerPatient>MinCellsRequired]
    if (sum(!(nrCellsPerPatient>MinCellsRequired)) > 0) {
        warning(paste0('Throwing out the following patients, because too few cells: ', 
            paste0(names(currentSeuratObject_list3)[!(nrCellsPerPatient>MinCellsRequired)], collapse = '; ')))
    }
    
    # SCTransform
    currentSeuratObject_list3 <- lapply(X = currentSeuratObject_list3, FUN = SCTransform); beepr::beep()
    
    all_genenames = Reduce(intersect, lapply(currentSeuratObject_list3, rownames))
    
    #currentFeatures3 <- SelectIntegrationFeatures(object.list = currentSeuratObject_list3, nfeatures = 3000)
    currentFeatures3 <- all_genenames
    
    currentSeuratObject_list3 <- PrepSCTIntegration(object.list = currentSeuratObject_list3, anchor.features = currentFeatures3); beepr::beep()
    currentAnchors3 <- FindIntegrationAnchors(object.list = currentSeuratObject_list3, 
        normalization.method = 'SCT', anchor.features = currentFeatures3, dims=1:30); beepr::beep() 
        # dims=1:10, k.anchor = ..

    # Instead of    
    #currentSeuratObject_recombined3 <- IntegrateData(anchorset = currentAnchors3, normalization.method = 'SCT'); beepr::beep()
    #currentSeuratObject_recombined3 <- IntegrateData(anchorset = currentAnchors3, normalization.method = 'SCT', k.weight=10); beepr::beep()
    # Automate these two options by trying first one first:
    currentSeuratObject_recombined3 = tryCatch({
                    IntegrateData(anchorset = currentAnchors3, normalization.method = 'SCT')
                }, error = function(e) {
                    print('Setting k.weight=10 because default k.weight failed.')
                    IntegrateData(anchorset = currentAnchors3, normalization.method = 'SCT', k.weight=10)
                }) ; beepr::beep()
    
    ############################################################
    
    # Run analysis again
    currentSeuratObject_recombined3 = mySeuratAnalysis_verybasic_part2only(mySeuratObject = currentSeuratObject_recombined3)
    
    ############################################################
    
    print(paste0('Saving object to: ', base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',paste0(set_name, '_Int1c'),'.h5seurat'))
    SaveH5Seurat(object = currentSeuratObject_recombined3,
                filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',paste0(set_name, '_Int1c'),'.h5seurat'))
    
    return(currentSeuratObject_recombined3)
    
}






################################################################################