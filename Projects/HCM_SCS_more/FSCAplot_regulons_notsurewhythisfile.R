




current_highlight_genes = current_regulon_set[[reg_name]]

correlations_FSCA_per_patient_combined$regulon_membership = 'no'
correlations_FSCA_per_patient_combined$regulon_membership[correlations_FSCA_per_patient_combined$gene_symbol %in% current_highlight_genes] = 'yes'
correlations_FSCA_per_patient_combined$regulon_membership = as.factor(correlations_FSCA_per_patient_combined$regulon_membership)

correlations_FSCA_per_patient_combined_sel = correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$regulon_membership=='yes',]
toplist_p4=correlations_FSCA_per_patient_combined_sel$gene_symbol[order(correlations_FSCA_per_patient_combined_sel$P4_cor, decreasing = T)][1:20]
toplist_p5=correlations_FSCA_per_patient_combined_sel$gene_symbol[order(correlations_FSCA_per_patient_combined_sel$P5_cor, decreasing = T)][1:20]
top_both = toplist_p4[toplist_p4 %in% toplist_p5]
toString(top_both)

lims_p4=c(min(c(-.11, correlations_FSCA_per_patient_combined_sel$P4_cor)), max(c(.21, correlations_FSCA_per_patient_combined_sel$P4_cor)))
lims_p5=c(min(c(-.11, correlations_FSCA_per_patient_combined_sel$P5_cor)), max(c(.21, correlations_FSCA_per_patient_combined_sel$P5_cor)))
        
p=ggplot(correlations_FSCA_per_patient_combined %>% subset(regulon_membership=='yes'), aes(x=P4_cor, y=P5_cor))+
            geom_hline(yintercept = 0)+
            geom_vline(xintercept = 0)+
            geom_point()+
            geom_point(data=correlations_FSCA_per_patient_combined %>% subset(gene_symbol%in%top_both), color='red')+
            geom_text_repel(aes(label=gene_symbol), size=8/.pt, max.overlaps = Inf)+
            theme_bw()+give_better_textsize_plot(8)+ggtitle(gsub('\\(\\+\\)','',reg_name))+
            theme(legend.position = 'none')+
            scale_x_continuous(breaks = seq(-1, 1, by = .2), limits = lims_p4)+
            scale_y_continuous(breaks = seq(-1, 1, by = .2), limits = lims_p5)
p
ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',DATASET_NAME,'_9_custom_size-correlations_module2.pdf'), 
                   width = (172)-4, height= (172)-4, units='mm', device = cairo_pdf)


p=ggplot(correlations_FSCA_per_patient_combined %>% subset(regulon_membership=='yes'), aes(x=P4_cor, y=P5_cor))+
            geom_point()+
            geom_text_repel(aes(label=gene_symbol))+
            theme_bw()+give_better_textsize_plot(6)+ggtitle(gsub('\\(\\+\\)','',reg_name))+
            theme(legend.position = 'none')+
            scale_x_continuous(breaks = seq(-1, 1, by = .2), limits = lims_p4)+
            scale_y_continuous(breaks = seq(-1, 1, by = .2), limits = lims_p5)
p
