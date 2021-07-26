################################################################################
# FUNCTION DEFINITIONS TO RUN ANALYSES
# Now perform the analysis

# Let's always define these functions

# To test:
# mySeuratObject = RHL_SeuratObject_merged_noMito_sel

mySeuratAnalysis = function(mySeuratObject, run_name,
    normalization.method='LogNormalize', scale.factor=10000,
    do.scale=T,do.center=T,scale.max=10,features_to_use_choice='variable',
    remove_genes=T) {
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
    mySeuratObject <- FindClusters(mySeuratObject, resolution = 0.5)
    
    # save(list = c('mySeuratObject'), file = paste0(base_dir,'Rdata/WRL_object_bigbig_sel_AnaDone.Rdata'))
    
    return(mySeuratObject)
    
}

mySeuratAnalysis_verybasic_part2only = function(mySeuratObject,
    do.scale=T,do.center=T,scale.max=10,features_to_use_choice='all') {
    
    assayType = if ('integrated' %in% names(mySeuratObject@assays)) {'integrated'} else {'RNA'}
    
    features_to_use = if (features_to_use_choice=='all') { rownames(mySeuratObject@assays$RNA@data)
        } else {VariableFeatures(mySeuratObject)}
    
    mySeuratObject <- ScaleData(mySeuratObject, features = features_to_use, do.scale=do.scale, do.center=do.center, scale.max=scale.max) # features = all_variable_features is default
    mySeuratObject <- RunPCA(object=mySeuratObject, npcs = 30, verbose = FALSE, features=features_to_use)
    mySeuratObject <- RunUMAP(mySeuratObject, reduction = "pca", dims = 1:30)
    mySeuratObject <- FindNeighbors(mySeuratObject, reduction = "pca", dims = 1:30)
    mySeuratObject <- FindClusters(mySeuratObject, resolution = 0.5)
    return(mySeuratObject)
}

mySeuratCommonPlots = function(mySeuratObject, 
    run_name,
    mymarkers = c('MALAT1', 'TTN', 'MYH7', 'MYH6', 'NPPA', 'NPPB', 'ACTA1','MYL2','SORBS2','CSRP3','NDUFA4','CRYAB','HSPB1', 'KCNQ1OT1')) {    
    
    mymarkers_found = mymarkers[mymarkers %in% rownames(mySeuratObject)]
    mymarkers_notfound = mymarkers[!(mymarkers %in% rownames(mySeuratObject))]
    print(paste0('Markers not found: ', paste0(mymarkers_notfound, collapse = ', ')))
    
    # Show umap with annotations
    for (current_annotation in c('annotation_sample_str','annotation_patient_str','annotation_paper_str','annotation_region_str','ident')) {
        # create labeled and unlabeled version
        p=DimPlot(mySeuratObject, group.by = current_annotation, cols = rep(col_vector_60,2), label = T, repel = T, label.size = 7)
        p_nl=DimPlot(mySeuratObject, group.by = current_annotation, cols = rep(col_vector_60,2))
        # remove legend if it's too large (prevents errors)
        if(!(current_annotation=='ident')) {if (dim(unique(mySeuratObject[[current_annotation]]))[1]>30) {
         p=p+theme(legend.position = 'none')   
         p_nl=p_nl+theme(legend.position = 'none')
        }}
        # Save files
        ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_2_umapLabeled_by_',current_annotation,'.png'), plot = p, height=7, width=7)
        ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_2_umap_by_',current_annotation,'.png'), plot = p_nl, height=7, width=7)
    }
    
    for (marker in mymarkers_found) {
        
        # marker = 'KCNQ1OT1'
        # marker = 'ACTC1'
        
        # marker projection on umaps
        p_cm = FeaturePlot(mySeuratObject, features = marker, cols = rainbow_colors)
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
    
    # Custom plot that shows distribution of patients over clusters
    p=ggplot(data.frame( cluster = Idents(mySeuratObject),
                    Donor = mySeuratObject$annotation_patient_fct))+
        geom_bar(aes(x=cluster, fill=Donor))+theme_bw()+
        xlab('Cluster')+ylab('Number of cells')+give_better_textsize_plot(10)
    ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_5_Barplot_PatientCluster_distr.png'), 
        plot = p, height=7.5, width=7.5, units='cm')
    
}

diff_express_clusters = function(mySeuratObject, mc.cores=8) {
    
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

diff_express_clusters_save_results = function(all_markers, run_name, base_dir, topX=10) {
     
    # Now collect the top-10 (log2fold) genes for each
    topHitsPerCluster=
        data.frame(lapply(all_markers, function(x) {rownames(x[order(x$avg_log2FC, decreasing = T),][1:topX,])} ))
    names(topHitsPerCluster)=names(all_markers)
    # and save
    openxlsx::write.xlsx(x = topHitsPerCluster, file = paste0(base_dir,'Rplots/ClusterTopHits_',run_name,'.xlsx'))
    
    # Also save the whole thing to Excel
    all_markers = lapply(all_markers, function(df) {df$gene = rownames(df); df = df[order(df$avg_log2FC,decreasing = T),]; return(df)})
    openxlsx::write.xlsx(x = all_markers, file = paste0(base_dir,'Rplots/Cluster_DE_full_',run_name,'.xlsx'))
       
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