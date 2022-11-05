
# Note that tfs and ligands should be available; check new repo for this

human_TFs_combined_temp = unique(human_TFs, TF_human_trrust_table$V1)

load(file = paste0(base_dir,'Rdata/FSCA__correlations_FSCA_per_patient_combined.Rdata')) # correlations_FSCA_per_patient_combined

correlations_FSCA_per_patient_combined$ens = 
    shorthand_cutname(correlations_FSCA_per_patient_combined$gene, PART1OR2 = 1)

# First let's look at TFs
FONTSIZE=10
ggplot(correlations_FSCA_per_patient_combined, aes(x=P4_cor, y=P5_cor))+ #  %>% arrange(regulon_membership)
            geom_point(size=.1, color='blue')+      #shape=1,alpha=.1)+
            geom_point(data=correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$gene_symbol %in% human_TFs_combined_temp,],
                       size=.5, color='red')+
            geom_text_repel(data=correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$gene_symbol %in% human_TFs_combined_temp,],
                            aes(label=gene_symbol), size=FONTSIZE/.pt, max.overlaps = Inf, force=20)+
            geom_rect(aes(xmin = calc_limits(correlations_FSCA_per_patient_combined$P4_cor,percentile = .25)[2], xmax = max(correlations_FSCA_per_patient_combined$P4_cor), 
                          ymin = calc_limits(correlations_FSCA_per_patient_combined$P5_cor,percentile = .25)[2], ymax = max(correlations_FSCA_per_patient_combined$P5_cor)), 
                            fill=NA, color='black')+
            theme_bw()+give_better_textsize_plot(FONTSIZE)

            
# Then let's look at ligands
FONTSIZE=10
ggplot(correlations_FSCA_per_patient_combined, aes(x=P4_cor, y=P5_cor))+ #  %>% arrange(regulon_membership)
            geom_point(size=.1, color='blue')+      #shape=1,alpha=.1)+
            geom_point(data=correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$ens %in% ligand_list_unique_ENS,],
                       size=.5, color='red')+
            geom_text_repel(data=correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$ens %in% ligand_list_unique_ENS,],
                            aes(label=gene_symbol), size=FONTSIZE/.pt, max.overlaps = Inf, force=20)+
            geom_rect(aes(xmin = calc_limits(correlations_FSCA_per_patient_combined$P4_cor,percentile = .25)[2], xmax = max(correlations_FSCA_per_patient_combined$P4_cor), 
                          ymin = calc_limits(correlations_FSCA_per_patient_combined$P5_cor,percentile = .25)[2], ymax = max(correlations_FSCA_per_patient_combined$P5_cor)), 
                            fill=NA, color='black')+
            geom_hline(yintercept = calc_limits(correlations_FSCA_per_patient_combined$P5_cor,percentile = .25)[2], linetype='dashed')+
            theme_bw()+give_better_textsize_plot(FONTSIZE)

correlations_FSCA_per_patient_combined_LigSel =
    correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$P5_cor>calc_limits(correlations_FSCA_per_patient_combined$P5_cor,percentile = .25)[2]&
                                           correlations_FSCA_per_patient_combined$ens %in% ligand_list_unique_ENS,]

View(correlations_FSCA_per_patient_combined_LigSel)

calc_limits(correlations_FSCA_per_patient_combined$P4_cor,percentile = .1)[2]     
calc_limits(correlations_FSCA_per_patient_combined$P5_cor,percentile = .1)[2]

            
ggplot(correlations_FSCA_per_patient_combined)+
    geom_point(aes())