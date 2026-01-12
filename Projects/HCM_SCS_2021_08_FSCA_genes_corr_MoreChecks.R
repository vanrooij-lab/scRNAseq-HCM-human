


current_analysis$ROOIJonly_RID2l_clExtended$sel14 = Idents(current_analysis$ROOIJonly_RID2l_clExtended) %in% 1:4
DimPlot(current_analysis$ROOIJonly_RID2l_clExtended, group.by='sel14')
celnames14 = colnames(current_analysis$ROOIJonly_RID2l_clExtended)[current_analysis$ROOIJonly_RID2l_clExtended$sel14]

current_analysis$ROOIJonly_RID2l_clExtended$sel5 = Idents(current_analysis$ROOIJonly_RID2l_clExtended) %in% 5
celnames5 = colnames(current_analysis$ROOIJonly_RID2l_clExtended)[current_analysis$ROOIJonly_RID2l_clExtended$sel5]
DimPlot(current_analysis$ROOIJonly_RID2l_clExtended, group.by='sel5')


#######

MANUAL_SELECTIONS = list(all=T,
                     c14=colnames(current_analysis$ROOIJonly_RID2l_FSCA) %in% celnames14, 
                     c5=colnames(current_analysis$ROOIJonly_RID2l_FSCA) %in% celnames5)

correlations_FSCA_per_patient_list_selected=list()
for (set_name in names(MANUAL_SELECTIONS)) {
    
    MANUAL_SELECTION=MANUAL_SELECTIONS[[set_name]]
    the_patients = unique(current_analysis$ROOIJonly_RID2l_FSCA$annotation_patient_str)
    correlations_FSCA_per_patient_list_selected[[set_name]] = 
        lapply(the_patients, function(current_patient) {
            
            print(paste0('Analyzing set ',set_name,', patient ',current_patient))
            
            cell_sel2 = (current_analysis$ROOIJonly_RID2l_FSCA$annotation_patient_str == current_patient) & MANUAL_SELECTION
            genes = rownames(current_analysis$ROOIJonly_RID2l_FSCA@assays$RNA@data)[rowMeans(current_analysis$ROOIJonly_RID2l_FSCA@assays$RNA@data[,cell_sel2]>0)>GENE_MIN_PERCENTAGE]
            
            # Perform test
            correlations_FSCA = as.data.frame(t(data.frame(lapply(genes, function(g) {
                    cor_out = cor.test(current_analysis$ROOIJonly_RID2l_FSCA@assays$RNA@data[g,cell_sel2], indexData_MW$FSC_A[cell_sel2])        
                    return(c(cor_out$estimate, p=cor_out$p.value))
                }))))
            rownames(correlations_FSCA) = genes
            correlations_FSCA$fdr=p.adjust(correlations_FSCA$p, method='fdr')
            correlations_FSCA$p.adjust=p.adjust(correlations_FSCA$p, method='BH')
            
            return(correlations_FSCA)
            
        })
    names(correlations_FSCA_per_patient_list_selected[[set_name]]) = the_patients
}

################################################################################

correlations_FSCA_per_patient_combined_forSel=list()
p_list=list()
for (set_name in names(MANUAL_SELECTIONS)) {
    
    # View(correlations_FSCA_per_patient_list_selected[[set_name]]$R.P5)
    
    # Combine the two frames
    #correlations_FSCA_per_patient_list_selected[[set_name]]$R.P4$patient = 'R.P4'
    #correlations_FSCA_per_patient_list_selected[[set_name]]$R.P5$patient = 'R.P5'
    genes_overlap = intersect(rownames(correlations_FSCA_per_patient_list_selected[[set_name]]$R.P4), rownames(correlations_FSCA_per_patient_list_selected[[set_name]]$R.P5))
    if (!any(grepl('P4_',colnames(correlations_FSCA_per_patient_list_selected[[set_name]]$R.P4)))) {colnames(correlations_FSCA_per_patient_list_selected[[set_name]]$R.P4) = paste0('P4_',colnames(correlations_FSCA_per_patient_list_selected[[set_name]]$R.P4))}
    if (!any(grepl('P5_',colnames(correlations_FSCA_per_patient_list_selected[[set_name]]$R.P5)))) {colnames(correlations_FSCA_per_patient_list_selected[[set_name]]$R.P5) = paste0('P5_',colnames(correlations_FSCA_per_patient_list_selected[[set_name]]$R.P5))}
    correlations_FSCA_per_patient_combined_forSel[[set_name]] = cbind(correlations_FSCA_per_patient_list_selected[[set_name]]$R.P4[genes_overlap,], correlations_FSCA_per_patient_list_selected[[set_name]]$R.P5[genes_overlap,])
    correlations_FSCA_per_patient_combined_forSel[[set_name]]$gene_symbol = shorthand_cutname( rownames(correlations_FSCA_per_patient_combined_forSel[[set_name]]) )
    correlations_FSCA_per_patient_combined_forSel[[set_name]]$gene = rownames(correlations_FSCA_per_patient_combined_forSel[[set_name]]) 
    
    tresholds_p4=calc_limits(correlations_FSCA_per_patient_combined_forSel[[set_name]]$P4_cor, 0.15)[2]
    tresholds_p4
    tresholds_p5=calc_limits(correlations_FSCA_per_patient_combined_forSel[[set_name]]$P5_cor, 0.15)[2]
    tresholds_p5
    
    selected_genes=correlations_FSCA_per_patient_combined_forSel[[set_name]]$P4_cor>.2&correlations_FSCA_per_patient_combined_forSel[[set_name]]$P5_cor>.1
    selected_genes=correlations_FSCA_per_patient_combined_forSel[[set_name]]$P4_cor>tresholds_p4&correlations_FSCA_per_patient_combined_forSel[[set_name]]$P5_cor>tresholds_p5
    shorthand_cutname(rownames(correlations_FSCA_per_patient_combined_forSel[[set_name]][selected_genes,]))
    
    selected_genes = correlations_FSCA_per_patient_combined_forSel[[set_name]]$gene_symbol=='MYL2'
    
    p_list[[set_name]]=ggplot(correlations_FSCA_per_patient_combined_forSel[[set_name]], aes(x=P4_cor, y=P5_cor))+
        geom_point(size=.5, shape=1, color='grey')+
        geom_smooth(method = 'lm')+
        geom_point(size=.5, shape=1, data=correlations_FSCA_per_patient_combined_forSel[[set_name]][selected_genes,], color='red')+
        geom_text_repel(data=correlations_FSCA_per_patient_combined_forSel[[set_name]][selected_genes,], color='red',
                        mapping=aes(label=gene_symbol), 
                        max.overlaps = Inf, min.segment.length = 0, force = 4, size=6/.pt, segment.size=.25)+
        theme_bw()+give_better_textsize_plot(6)+
        xlim(c(min(correlations_FSCA_per_patient_combined_forSel[[set_name]]$P4_cor)*1.01,max(correlations_FSCA_per_patient_combined_forSel[[set_name]]$P4_cor)*1.01))+
        ylim(c(min(correlations_FSCA_per_patient_combined_forSel[[set_name]]$P5_cor)*1.01,max(correlations_FSCA_per_patient_combined_forSel[[set_name]]$P5_cor)*1.01))+
        xlab('Gene-cell size correlation patient 4')+ylab('Gene-cell size correlation patient 5')+
        ggtitle(set_name)
    p_list[[set_name]]

}

p=wrap_plots(p_list)
p
ggsave(plot = p,filename = paste0(base_dir, 'Rplots/ROOIJonly_RID2l_clExtended_9_FSCA_3way.pdf'), 
       width = 172-4, height= 172/3, units='mm', device = cairo_pdf)
#ggsave(plot = p,filename = paste0(base_dir, 'Rplots/ROOIJonly_RID2l_clExtended_9_FSCA_3way_MYL2.pdf'), 
#       width = 172-4, height= 172/3, units='mm', device = cairo_pdf)










