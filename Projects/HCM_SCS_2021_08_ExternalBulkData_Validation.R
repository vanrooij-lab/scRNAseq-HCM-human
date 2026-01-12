
################################################################################

# Data from:
# 
# - Hemerich, D., Pei, J., Harakalova, M., Van Setten, J., Boymans, S., Boukens, B. J., … Asselbergs, F. W. (2019). Integrative Functional Annotation of 52 Genetic Loci Influencing Myocardial Mass Identifies Candidate Regulatory Variants and Target Genes. Circulation: Genomic and Precision Medicine, 12(2), 76–83. https://doi.org/10.1161/CIRCGEN.118.002328
# - Pei, J., Schuldt, M., Nagyova, E., Gu, Z., el Bouhaddani, S., Yiangou, L., … Harakalova, M. (2021). Multi-omics integration identifies key upstream regulators of pathomechanisms in hypertrophic cardiomyopathy due to truncating MYBPC3 mutations. Clinical Epigenetics, 13(1), 1–20. https://doi.org/10.1186/s13148-021-01043-3

################################################################################
# Pre-processing the bulk data and conversion to Seurat dataframe

if (F) {
    
    # Load external data set
    
    louk_data_ = read.table(paste0(base_dir, 'External_data/Raw_read_counts.txt'), header = 1)
    
    load(file = paste0(base_dir,'Rdata/ens_to_sym_conv_table.Rdata'))
    
    louk_data = louk_data_[,-1]
    rownames(louk_data) = 
        paste0(louk_data_$GeneID, ':', ens_to_sym_conv_table[louk_data_$GeneID])
    
    
    # Load it into Seurat
    
    Seurat_Object_BulkData = CreateSeuratObject(counts = louk_data, project = 'bulkData')
    
    rownames(Seurat_Object_BulkData)[grepl(':MT-', rownames(Seurat_Object_BulkData))]
    
    # Determine mitochondrial percentage
    Seurat_Object_BulkData[["percent.mt"]] <- PercentageFeatureSet(Seurat_Object_BulkData, pattern = ":MT-")
    
    # Remove mitochondrial reads
    mito_genes = rownames(Seurat_Object_BulkData)[grepl(':MT-',rownames(Seurat_Object_BulkData))]
    Seurat_Object_BulkData@misc$desired_genes_excl_mito = 
                rownames(Seurat_Object_BulkData)[!(rownames(Seurat_Object_BulkData) %in% mito_genes)]
    Seurat_Object_BulkData_noMito =
            subset(Seurat_Object_BulkData, features= Seurat_Object_BulkData@misc$desired_genes_excl_mito)
        
    # Determine total counts and features 
    CalcN_out = Seurat:::CalcN(Seurat_Object_BulkData_noMito)
    Seurat_Object_BulkData_noMito[['nCount.nMT']]   = CalcN_out$nCount
    Seurat_Object_BulkData_noMito[['nFeature.nMT']] = CalcN_out$nFeature
    
    # Run Seurat analysis (I only need the scaling, actually)
    
    #Seurat_Object_BulkData_noMito = mySeuratAnalysis(Seurat_Object_BulkData_noMito,
    #            run_name='BulkData_noMito',
    #            normalization.method='RC', scale.factor=median(Seurat_Object_BulkData_noMito$nCount.nMT),
    #            do.scale=F,do.center=F,scale.max=Inf, features_to_use_choice = 'variable') # variable because otherwise too large calculation ..
    
    # Perform scaling
    Seurat_Object_BulkData_noMito <- NormalizeData(Seurat_Object_BulkData_noMito, normalization.method = 'RC', scale.factor = median(Seurat_Object_BulkData_noMito$nCount.nMT))
    
    # Annotate samples
    Seurat_Object_BulkData_noMito$annotation_str = NA
    Seurat_Object_BulkData_noMito$annotation_str[grepl('^HCM', colnames(Seurat_Object_BulkData_noMito))] = 'HCM'
    Seurat_Object_BulkData_noMito$annotation_str[grepl('^CONTROL', colnames(Seurat_Object_BulkData_noMito))] = 'Ctrl'
    Seurat_Object_BulkData_noMito$annotation_fct=factor(Seurat_Object_BulkData_noMito$annotation_str, levels=c('HCM', 'Ctrl'))

    SaveH5Seurat(object = Seurat_Object_BulkData_noMito, overwrite = T,
            filename = paste0(base_dir,'Rdata/H5_Seurat_Object_BulkData_noMito.h5seurat'))
    

    
}

################################################################################
# Loading data

CURRENT_RUNNAME='BulkData_noMito'
current_analysis[[CURRENT_RUNNAME]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_Seurat_Object_BulkData_noMito.h5seurat'))

################################################################################
# Running analysis for regulon genes

# Box plots
shorthand_custom_boxplot(seuratObject_list=current_analysis, 
                         gene_lists=core_regulons_sorted, 
                         seuratObjectNameToTake=CURRENT_RUNNAME, 
                         group.by='annotation_fct', 
                         topX=10, mylimits=.01) 

# Summary composite plots
p_lists=shorthand_custom_compositeplot(seuratObject_list=current_analysis, 
                         gene_lists=core_regulons_sorted, 
                         seuratObjectNameToTake=CURRENT_RUNNAME, 
                         group.by='annotation_fct', 
                         #group.by2='annotation_patient_fct',
                         zscore=F) 

p=wrap_plots(p_lists$p_boxplot_list, nrow=1)
p
ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customCOMPOSITE_REG_customCOEXP_boxplot.pdf'), plot = p,
       height=(172-4)/5, width=172-4, units='mm', device=cairo_pdf)

# p=wrap_plots(p_lists$p_violin_list, nrow=1)
# p
# ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customCOMPOSITE_REG.pdf'), plot = p,
#        height=26.5, width=min(184.6-4, 26.5*6), units='mm', device=cairo_pdf)
# p=wrap_plots(p_lists$p_bar_list_g2, nrow=1)
# ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customCOMPOSITE_REG_g2.pdf'), plot = p,
#        height=26.5, width=min(184.6-4, 26.5*6), units='mm', device=cairo_pdf)

# Also do SCENIC regulons
# ===
load(file=paste0(base_dir, 'Rdata/SCENIC_regulons_core_genes_sel.Rdata'))

# Box plots
shorthand_custom_boxplot(seuratObject_list=current_analysis, 
                         gene_lists=SCENIC_regulons_core_genes_sel, 
                         seuratObjectNameToTake=CURRENT_RUNNAME, 
                         group.by='annotation_paper_oneletter_fct', 
                         topX=10, mylimits=.01) 

# Summary plots
p_lists=shorthand_custom_compositeplot(seuratObject_list=current_analysis, 
                         gene_lists=SCENIC_regulons_core_genes_sel, 
                         seuratObjectNameToTake=CURRENT_RUNNAME, 
                         group.by='annotation_fct', 
                         #group.by2='annotation_patient_fct',
                         zscore=T) 

p=wrap_plots(p_lists$p_boxplot_list, nrow=4)
p
ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customCOMPOSITE_REG_SCENIC_boxplot.pdf'), plot = p,
       height=(172-4), width=172-4, units='mm', device=cairo_pdf)
# p=wrap_plots(p_lists$p_violin_list, nrow=4)
# ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customCOMPOSITE_REG_SCENIC.pdf'), plot = p,
#        height=26.5*4, width=min(184.6-4, 26.5*4), units='mm', device=cairo_pdf)
# p=wrap_plots(p_lists$p_bar_list_g2, nrow=4)
# ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customCOMPOSITE_REG_SCENIC_g2.pdf'), plot = p,
#        height=26.5*4, width=min(184.6-4, 26.5*4), units='mm', device=cairo_pdf)

# UMAP
p_list = lapply(names(SCENIC_regulons_core_genes_sel), function(reg_name) {
            shorthand_seurat_custom_expr(seuratObject = Seurat_Object_BulkData_noMito, 
                             gene_of_interest = SCENIC_regulons_core_genes_sel[[reg_name]],
                             textsize = 6, pointsize = .5, custom_title = reg_name, mymargin = .1, zscore = T) 
    })
p=wrap_plots(p_list, nrow=4)
ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_REG_SCENIC.pdf'), plot = p,
       height=4*30, width=min(184.6-4, 30*4), units='mm', device=cairo_pdf)      

################################################################################
# Analysis for the custom-correlation genes (correlated to NPPA, XIRP2, etc)


# Load and re-organize the gene lists from the correlation-analysis
load(paste0(base_dir,'Rdata/gene_lists_customcorrelated__Rooijbased.Rdata'))
# Re-organize lists
gene_lists_customcorrelated_reorganized=list()
for (gene in names(gene_lists_customcorrelated)) {
    if (!is.na(gene_lists_customcorrelated[[gene]]$pos)) {
        gene_lists_customcorrelated_reorganized[[paste0('posCorrWith_',gene)]] = gene_lists_customcorrelated[[gene]]$pos
    }
    if (!is.na(gene_lists_customcorrelated[[gene]]$neg)) {
        gene_lists_customcorrelated_reorganized[[paste0('negCorrWith_',gene)]] = gene_lists_customcorrelated[[gene]]$neg
    }
}

# Box plots
shorthand_custom_boxplot(seuratObject_list=current_analysis, 
                         gene_lists=gene_lists_customcorrelated_reorganized, 
                         seuratObjectNameToTake=CURRENT_RUNNAME, 
                         group.by='annotation_fct', 
                         topX=10, mylimits=.01)

################################################################################
# Let's also quickly look at the genes that were lower expressed in cluster 5

DE_genes_Cl5_ordered_rev = DE_cluster$ROOIJonly_RID2l_clExtended$`5`[order(DE_cluster$ROOIJonly_RID2l_clExtended$`5`$avg_log2FC),]
    # View(DE_genes_Cl5_ordered_rev)

Cl5_list = list(Cl5_down = rownames(DE_genes_Cl5_ordered_rev)[1:10])

shorthand_custom_boxplot(seuratObject_list=current_analysis, 
                         gene_lists=Cl5_list, 
                         seuratObjectNameToTake=CURRENT_RUNNAME, 
                         group.by='annotation_fct', 
                         topX=10, mylimits=.01)

################################################################################
# SCENIC TFs

shorthand_custom_boxplot(seuratObject_list=current_analysis, 
                         gene_lists=list(SCENIC_TFs=gsub('\\(\\+\\)','',names(SCENIC_regulons_core_genes))), 
                         seuratObjectNameToTake=CURRENT_RUNNAME, 
                         group.by='annotation_fct', 
                         topX=10, mylimits=.01)

# See main script for same across all datasets (@ "# Now do the SCENIC TFs")



