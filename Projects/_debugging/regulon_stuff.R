library(randomcoloR)
library(ggplot2)
library(scales)
library(gridExtra)

current_matrix=current_matrix[1:3000,1:50]

run_name='test'

EXPRESSION_ZERO_CUTOFF = 0
dir.create(paste0(base_dir,'regulon/analysis_',run_name,'/regulon/'), recursive = T)

genes_expressed_in_cells = rowSums(current_matrix>EXPRESSION_ZERO_CUTOFF)
genes_expr_in_cells_fraction = genes_expressed_in_cells/dim(current_matrix)[2]
gene_selection = genes_expr_in_cells_fraction>0.1
print(paste0('Genes in run: ',sum(gene_selection)))

regulon_object = MW_determine_regulons_part1(
        expression_matrix = current_matrix[gene_selection,], 
        calculate_p = T,
        outputDir = paste0(base_dir,'regulon/'),
        analysis_name = run_name)
regulon_object = MW_determine_regulons_part2(regulon_object=regulon_object, 
    chosen_cutoff_parameter='p',
    p_or_r_cutoff=0.00001) 
regulon_object = MW_determine_regulons_part3(regulon_object, 
                    connectedness_cutoff = 1, 
                    min_expression_fraction_genes=.1,
                    show_heatmap=F,
                    chosen_cutoff_parameter = 'p')
regulon_object = MW_determine_regulons_part4(regulon_object = regulon_object)
regulon_object = MW_determine_regulons_part5(regulon_object = regulon_object,
        hierarchical_cutoff = NULL) # current_cutoff)
