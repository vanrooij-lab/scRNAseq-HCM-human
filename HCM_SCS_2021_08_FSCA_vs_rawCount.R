

# Comparing FSCA data with RNA read counts
# Note that I use the raw data here, since I want to include (i) non-filtered cells and (ii) mitochondrial counts



# Loading raw data
current_analysis$RHL_SeuratObject_Rooij_raw = 
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_Rooij_raw.h5seurat'))

# Loading FSCA data
load(file = paste0(base_dir,'Rdata/FSCA__indexData_MW.Rdata')) # indexData_MW, colnames match cell names
    # e.g. use indexData_MW[colnames(current_analysis$ROOIJonly_RID2l_clExtended),]$FSC_A

# Insert into the raw data
current_analysis$RHL_SeuratObject_Rooij_raw$FSCA = indexData_MW[colnames(current_analysis$RHL_SeuratObject_Rooij_raw),]$FSC_A

gene_names = rownames(current_analysis$RHL_SeuratObject_Rooij_raw)
current_mt_genes = gene_names[grepl(':MT-',gene_names)]

current_analysis$RHL_SeuratObject_Rooij_raw$mt_totalCounts = 
    colSums(current_analysis$RHL_SeuratObject_Rooij_raw@assays$RNA[current_mt_genes,])

# Filtered data
# 
# current_analysis$ROOIJonly_RID2l_clExtended_FSCA =
#     LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_','ROOIJonly_RID2l_clExtended_FSCA','.h5seurat'))
#     # ROOIJonly_RID2l_clExtended_FSCA

df_toplot = 
    data.frame(
    RNA_count = current_analysis$RHL_SeuratObject_Rooij_raw$nCount_RNA,
    #RNA_count_nMT = current_analysis$ROOIJonly_RID2l_clExtended_FSCA$nCount.nMT,
    FSCA = current_analysis$RHL_SeuratObject_Rooij_raw$FSCA,
    patient_ = current_analysis$RHL_SeuratObject_Rooij_raw$annotation_patient_str,
    pct_mito = current_analysis$RHL_SeuratObject_Rooij_raw$mt_totalCounts/current_analysis$RHL_SeuratObject_Rooij_raw$nCount_RNA)
df_toplot$patient = as.factor(gsub('^R\\.','',df_toplot$patient_))

# RNA violin
df_toplot_ptSel = df_toplot[df_toplot$patient %in% c('P4','P5'),]
p=ggplot(df_toplot_ptSel, aes(x=patient, y=log10(.1+RNA_count)))+
    geom_jitter(size=.05, shape=22)+ # 20
    geom_violin(fill='gray',alpha=.5, lwd=.25)+
    geom_hline(yintercept = log10(.1+500), size=.25)+
    geom_hline(yintercept = log10(.1+30000), size=.25)+
    theme_bw()+give_better_textsize_plot(8)
p
ggsave(filename = paste0(base_dir,'Rplots/FSCA_raw_Violin-Count_inclMT.pdf'), plot = p, 
           width=PANEL_WIDTH-4, height=(PANEL_WIDTH-4), units='mm', device = cairo_pdf)

# FSCA violin
ggplot(df_toplot, aes(x=patient, y=FSCA))+
    geom_violin()+
    theme_bw()+give_better_textsize_plot(8)

# RNA 
ggplot(df_toplot, aes(x=log10(.1+RNA_count)))+
    geom_freqpoly()+
    theme_bw()+give_better_textsize_plot(8)

#ggplot(df_toplot, aes(x=log10(.1+RNA_count), y=log10(.1+RNA_count_nMT)))+
#    geom_point()+
#    theme_bw()+give_better_textsize_plot(8)

# RNA vs. FSCA
df_toplot_sel1 = df_toplot[df_toplot$RNA_count>500&df_toplot$RNA_count<30000&df_toplot$patient%in%c('P4','P5'),]
ggplot(df_toplot_sel1, aes(x=FSCA, y=RNA_count))+ # , color=pct_mito
    geom_point()+
    theme_bw()+give_better_textsize_plot(8)+
    facet_grid(cols = vars(patient))

# now per patient
MIN_COUNT = 500
# MIN_COUNT = 100
p_list = lapply(c('P4','P5'), function(current_pt) {
        # current_pt='P4'
        df_toplot_sel1 = df_toplot[df_toplot$RNA_count>MIN_COUNT&df_toplot$RNA_count<30000&df_toplot$patient==current_pt,]
        
        lm_out = lm(formula = RNA_count ~ FSCA, data = df_toplot_sel1) 
        cor.test_out = cor.test(df_toplot_sel1$RNA_count, df_toplot_sel1$FSCA)
        R=cor.test_out$estimate; p=cor.test_out$p.value

        p_string = if (p >= .001) {paste0('p=',round(p,3))} else {paste0('p=10^',round(log10(p),1))}
        
        ggplot(df_toplot_sel1, aes(x=FSCA, y=RNA_count))+ # , color=pct_mito
            geom_point(size=.25)+
            theme_bw()+give_better_textsize_plot(8)+
            geom_abline(slope = lm_out$coefficients[2], intercept = lm_out$coefficients[1], size=.25)+
            ggtitle(paste0('Patient ',substring(current_pt, 2,2), ', R=',round(R,2), ', ', p_string))
    })

p=wrap_plots(p_list)
p

ggsave(filename = paste0(base_dir,'Rplots/FSCA_raw_Scatter-Count-FSCA_inclMT.pdf'), plot = p, 
           width=PANEL_WIDTH*2-4, height=(PANEL_WIDTH*2-4)/2, units='mm', device = cairo_pdf)
    




