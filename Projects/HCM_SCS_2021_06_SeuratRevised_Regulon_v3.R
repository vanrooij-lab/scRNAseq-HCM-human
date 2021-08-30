
# Applying the regulon analysis to all of the different datasets separately,
# separately to each donor/patient.

# OK, convenient to again use the already made Seurat objects,
# this time, I will use the ones that aren't corrected, since
# the regulon analysis was appied per patient earlier,
# and i can split it out again

# For reference, this was what we got from the previous analysis:
# load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/Previous_analysis_for_reference/JoepAnalysis_Regulons.Rdata')
# See regulons_README_objects_saved
# See exporatory analyses for comparison earlier and current data analyses

# Note, for HPC custom analyses, this script can be sourced as
# desired_command_regulon='dummy'; source('/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/HCM_SCS_2021_06_SeuratRevised_Regulon_v3.R')

######################################################################
# libraries

print('Executing main script to load libraries')

myargs = commandArgs(trailingOnly = T)
print(paste0('myargs=',myargs))

# set script dir (note: overwritten by loading SeuratRevisedAnalysis below)
if (exists('LOCAL')) {
    script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/'
    desired_command='dummy'
    source(paste0(script_dir, 'Projects/HCM-SCS_2021-06_SeuratRevisedAnalysis_v2_UmiTools.R'))
    rm('desired_command')
} else {
    script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'
    desired_command='dummy'
    source(paste0(script_dir, 'HCM-SCS_2021-06_SeuratRevisedAnalysis_v2_UmiTools.R'))
    rm('desired_command')
}

# install.packages(c('ggdendro','pheatmap', 'randomcoloR', 'UpSetR'))

library(ggdendro)
library(gplots) # heatmap.2
library(pheatmap)

library(randomcoloR)
library(gridExtra)

library(parallel)

library(openxlsx)

library(UpSetR)

source(paste0(script_dir, 'Functions/MW_regulon_analysis_COPY.R'))
source(paste0(script_dir, 'Functions/MW_regulon_analysis_supportFns_COPIES.R'))
    # file.edit(paste0(script_dir, 'Functions/MW_general_functions.R'))

########################################################################

if (!exists('desired_command_regulon')) {
    cmd=command_interpreter_mw(myargs)
    desired_command_regulon=cmd$desired_command
    config=cmd$config
}

########################################################################
# To be able to run this script from the command line on the cluster,
# let's make sure we can pass arguments to it
#
# Note that some sections are executed together.

# stuff below is now handled by SeuratRevisedAnalysis
# args = commandArgs(trailingOnly = T)
# 
# library(stringr)
# if (!exists('desired_command')) {
#     if (length(args)==1) {
#         desired_command = unlist(str_split(string = args, pattern = '-'))
#         print(paste0('Sections to execute: ',paste0(desired_command,collapse = ', ')))
#         if(any(grepl(x = desired_command,pattern = 'CORES'))) {MYMCCORES=str_split(string = desired_command[grepl(x = desired_command,pattern = 'CORES')], pattern = '=')[[1]][2];print(paste0('MYMCCORES set to ',MYMCCORES))}
#     } else {
#         warning('Please pass 1 string, in the form \'arg1-arg2-arg3\' to give what section(s) to execute.')
#         print('Will proceed (this loads functions).')
#         desired_command='dummy'
#     }
# }

######################################################################
######################################################################
# Functions

# For testing purposes:
# get_GO_terms=F; required_minCellFraction=.05; connectedness_cutoff=10; chosen_cutoff_parameter='p'; p_or_r_cutoff=10^-3
#
# Using Tomo code!
giveMeRegulons_SeuratVersion = function(run_name, base_dir, current_matrix, 
    get_GO_terms=F, strip_object=F, MAX_GENES=NULL,
    required_minCellFraction=.05,
    connectedness_cutoff=10, 
    chosen_cutoff_parameter='p', p_or_r_cutoff=10^-3, BIGB=100) {
    
    EXPRESSION_ZERO_CUTOFF = 0 # 0.1
    run_name # run_name = 'groupedSCS_Patient1Mod'
    # run_name = 'ROOIJonly_RID2l_test_J_intersect' # DEBUG
    
    dir.create(paste0(base_dir,'regulon/analysis_',run_name,'/regulon/'), recursive = T)
    
    # Note that we can either take @data or @scaled.data, the latter only contains the 
    
    # current_analysis_temp
    
    genes_expressed_in_cells = rowSums(current_matrix>EXPRESSION_ZERO_CUTOFF)
    genes_expr_in_cells_fraction = genes_expressed_in_cells/dim(current_matrix)[2]
    gene_selection = genes_expr_in_cells_fraction>required_minCellFraction
    print(paste0('Genes in run: ',sum(gene_selection)))
    
    # Select genes first based on how many 
    # current_analysis_temp@misc$genes_expr_in_cells =
    #     rowSums(current_analysis_temp@assays$RNA@data>0)
    # current_analysis_temp@misc$genes_expr_in_cells_fraction =
    #     current_analysis_temp@misc$genes_expr_in_cells/
    #         dim(current_analysis_temp@assays$RNA@data)[2]
        # length(current_analysis_temp@misc$genes_expr_in_cells_fraction)
    #gene_selection = current_analysis_temp@misc$genes_expr_in_cells_fraction>0.1
    # print(paste0('Genes in run: ',sum(gene_selection)))
    
    regulon_object = MW_determine_regulons_part1(
        expression_matrix = current_matrix[gene_selection,], 
        calculate_p = T,
        outputDir = paste0(base_dir,'regulon/'),
        analysis_name = run_name, 
        MYPADJUSTMETHOD = 'fdr', # Joep's choice, default is BH in my code
        min_expression_fraction_genes = required_minCellFraction)
    
    # Next, we can assess which 'connectedness' cutoff we want.
    # With connectedness, I mean with how many other genes a gene is significantly
    # correlated. 
    # In principle, any gene that has a few connections might be of interest to
    # take along in the analysis. However, to narrow down the analysis to
    # genes that are most interesting, it can be usefull to only take along
    # genes that have many connections to other genes. 
    # This function provides plots with statistics on with how many other
    # genes genes are significantely correlated.
    regulon_object = MW_determine_regulons_part2(regulon_object=regulon_object, 
        chosen_cutoff_parameter=chosen_cutoff_parameter,
        # this parameter should be 'r' or 'p', respectively selection of 
        # genes based on r-value or p-value cutoff (r being the correlation 
        # coefficient value itself)
        #p_or_r_cutoff=.25
        p_or_r_cutoff=p_or_r_cutoff) # 1/100000, value Joep, has pretty high FPR though
        # this parameter provides the cutoff you have chosen earlier,
        # e.g. consider r-values (correlation coefficients) >.35 or <-.35 
        # significant.
    
    # Now, we can already create a correlation heatmap, which shows the structure
    # of correlation between genes. 
    # Per default, I don't show the heatmap here, because we can later determine
    # some parameters and set a cutoff to perform clustering on this correlation
    # matrix.
    regulon_object = MW_determine_regulons_part3(regulon_object, 
                    connectedness_cutoff = connectedness_cutoff, # this is rather arbitrary so let's set to 25 instead of 40
                    max_genes = MAX_GENES,
                    #min_expression_fraction_genes=.1,
                    #min_expression_fraction_genes=.1,
                    show_heatmap=F,
                    chosen_cutoff_parameter = 'p')
                    # Note that when using MAX_GENES, top connected genes are chosen;
                    # this might be biasing towards larger regulons
                    # So one might be able to improve this algorithm by e.g.
                    # picking genes whose average correlation to the X genes 
                    # it's most correlated with is highest or something.
                    #   
                    # genes in analysis: 
                    # dim(regulon_object$cor_out_selected_2)
    
    # To further quantify the pattern of correlations, we can sort the 
    # correlation data using hierarchical clustering, and also classify clusters 
    # based on this method.
    # To create the clustering, a cutoff should be determined. This can be done
    # based on the observed distances between points that are joined during the
    # hierarchical clustering procedure. The plot that is shown by this function
    # looks for a sudden change in length scales in the points which are joined,
    # and tries to suggest a good cutoff
    # (Given to the user in regulon_object$auto_cutoff1 and 
    # regulon_object$auto_cutoff2.) Also the dendrogram is shown, such that you
    # can inspect whether the suggested cutoff makes sense.
    regulon_object = MW_determine_regulons_part4(regulon_object = regulon_object)
    
    # Finally, we can create a plot which shows the correlation coefficient matrix
    # and also the hierarchical clustering dendrogram and cluster classification
    # performed on it.
    # With hierarchical_cutoff = NULL, gap_stat will be used OR Choose a cutoff value, e.g.:
    # - current_cutoff = regulon_object$auto_cutoff2
    regulon_object = MW_determine_regulons_part5(regulon_object = regulon_object,
        hierarchical_cutoff = NULL, KMAX_GAPTEST = 75, hierarchical_cutoff_fallback = regulon_object$auto_cutoff2, BIGB=BIGB) # current_cutoff)
        # hierarchical_cutoff = NULL; KMAX_GAPTEST = 75; hierarchical_cutoff_fallback = regulon_object$auto_cutoff2
    # regulon_object2 = MW_determine_regulons_part5(regulon_object = regulon_object, hierarchical_cutoff = NULL, KMAX_GAPTEST = 75, hierarchical_cutoff_fallback = regulon_object$auto_cutoff2, BIGB=BIGB) # current_cutoff)
        
    
    # Now we can add some information about the genes (are they TF, ligand, receptor?)
    data_container=list() # Change this to get TF ligand stuff (should be easily doable)
    print('See above: add TF/ligand stuff')
    regulon_object = MW_regulon_add_TF_RL_flags(regulon_object, data_container)
    # And we can export the regulon analysis to an excel file:
    MW_regulon_final_export(regulon_object)
    
    # GO analysis
    # Additionally, we can export the classification of genes to clusters to an
    # excel file.
    # Additionally, we can perform a GO analysis on these genes, with all genes
    # taken along in the matrix as background (TODO: perhaps reconsider which
    # genes to use as background panel??-- perhaps all genes in the expression
    # matrix would make more sense). 
    # This function also exports the GO term for each of the clusters to an 
    # excel file.
    if (get_GO_terms) {
        regulon_object = MW_determine_regulons_part6(regulon_object = regulon_object)
    }
    
    if (strip_object) {
        regulon_object = remove_data_regulon_object(regulon_object)
    }
    
    return(regulon_object)
}

remove_data_regulon_object = function(reg_object) {
        
    reg_object$rownames_matrices = rownames(reg_object$expression_matrix)
    reg_object$expression_matrix = NULL
    reg_object$p_val_matrix = NULL
    reg_object$p_val_matrix_adjusted = NULL
    reg_object$p_val_matrix_adjusted_NA = NULL
    reg_object$cor_out = NULL
    reg_object$cor_out_scrambled = NULL
    reg_object$cor_out_NA = NULL
    reg_object$rownames_cor_out_selected_2 = rownames(reg_object$cor_out_selected_2)
    reg_object$cor_out_selected_2 = NULL
    
    # remove plots that are data-heavy
    reg_object$p_corr_distr = NULL
    reg_object$p_cutoff_decision = NULL
    reg_object$p_tdist = NULL
    reg_object$p_connectedness = NULL
    reg_object$p_connectedness2 = NULL

        
    return(reg_object)
}


regulon_overlap_heatmap = function(pooled_regulons, base_dir, run_name, MYTREECUTTINGHEIGHT=2, myfontsize=8, makeallheatmaps=T) {
    
    # Create pairs to compare
    df_compare = tidyr::expand_grid(x=names(pooled_regulons), y=names(pooled_regulons))
    
    # Calculate overlaps
    df_compare$overlap = sapply(1:dim(df_compare)[1], function(X) { 
        sum(pooled_regulons[[df_compare$x[X]]] %in% pooled_regulons[[df_compare$y[X]]]) /
            min(length(pooled_regulons[[df_compare$x[X]]]), length(pooled_regulons[[df_compare$y[X]]]))
        })
    
    # Create matrix
    matrix_compare <- reshape2::acast(df_compare, x~y, value.var="overlap")
    
    # Show heatmap
    if (makeallheatmaps) {
        p0=pheatmap(matrix_compare, clustering_method = 'ward.D2')
        p0
        #ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_7_Regulons_noCats.png'), plot = p0, width=nrow(matrix_compare)*.4, height=nrow(matrix_compare)*.4, units='cm')
        ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_7_Regulons_noCats.pdf'), plot = p0, width=nrow(matrix_compare)*.44, height=nrow(matrix_compare)*.44, units='cm')
        #
        # Without legend
        p0=pheatmap(matrix_compare, clustering_method = 'ward.D2', legend = F)
        ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_7_Regulons_noCats_noLegend.pdf'), plot = p0, width=nrow(matrix_compare)*.44, height=nrow(matrix_compare)*.44, units='cm')
    }
    
    # Check out how they group using hclust
    hclust_out = hclust(dist(matrix_compare), method='ward.D2')
    
    # Plot dendrogram
    p=ggdendrogram(hclust_out)+give_better_textsize_plot(6)+theme_bw()+scale_y_continuous(minor_breaks = seq(0, max(hclust_out$height), .1)) # p
    ggsave(plot = p, filename = paste0(base_dir,'Rplots/',run_name,'_7_Regulons_Dendrogram_',run_name,'.png'), width = 80, height=80, units='mm')
    
    # 
    pcjd_out = plot_clust_join_density(hclust_out)
    p=pcjd_out$p2+theme_bw()+give_better_textsize_plot(8)
    ggsave(plot = p, filename = paste0(base_dir,'Rplots/',run_name,'_7_Regulons_DendrogramCutoffAnalysis_',run_name,'.png'), width = 80, height=80, units='mm')
    
    # Now again use gap-stat to test proper # of "core" regulons
    # That didn't work so well; line is monotically increasing, but by eye there
    # are clearly a well-defined amount of clusters
    # gap_stat <- cluster::clusGap(matrix_compare, FUN = function(matrix, k, mytree) {return(list(cluster=(cutree(mytree, k = k))))}, 
    #     mytree=hclust_out, K.max = min(KMAX_GAPTEST, nrow(matrix_compare)-1))
    # #
    # nCluster = cluster::maxSE(f=gap_stat$Tab[,'gap'], SE.f = gap_stat$Tab[,'SE.sim'], method = 'firstSEmax')#'Tibs2001SEmax')
    # p=ggplot(data.frame(gap=gap_stat$Tab[,'gap'],K=1:(min(KMAX_GAPTEST, nrow(matrix_compare)-1))))+
    # geom_vline(xintercept = nCluster)+
    # geom_line(aes(x=K, y=gap))+geom_point(aes(x=K, y=gap))+theme_bw()+
    # ggtitle(paste0('Gap-stat, K=',nCluster,' optimal'))+give_better_textsize_plot(8)
    # p

    cutree_out = cutree(hclust_out, h = MYTREECUTTINGHEIGHT) # MYTREECUTTINGHEIGHT=2.3
    # cutree_out = cutree(hclust_out, k = nCluster)
    cutree_df  = as.data.frame(as.factor(cutree_out)); colnames(cutree_df) = c('group')
    #annotation_colors = col_Dark2[1:max(cutree_out)]
    # Create a little heatmap
    annotation_colors = col_vector_60[1:max(cutree_out)]
    names(annotation_colors) = unique(cutree_out)
    annotation_colors=list(group=annotation_colors)
    
    # Heatmap Version 1
    p=pheatmap(matrix_compare, cluster_rows = hclust_out,cluster_cols = hclust_out, 
        annotation_col = cutree_df, annotation_row = cutree_df, annotation_colors = annotation_colors, 
        fontsize = myfontsize, fontsize_col = myfontsize, fontsize_row = myfontsize, treeheight_row = 0)
        #annotation_colors = list(colors=col_Dark2[1:max(cutree_out)])))
    print(p)
    ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_7_Regulons.pdf'), 
        plot = p, width=184.6*2/3-4, height=184.6*2/3-4, units='mm')
    
    ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_7_Regulons-L.pdf'), 
        plot = p, width=length(pooled_regulons)*2.5, height=length(pooled_regulons)*2.5, units='mm', limitsize = F)
    
    # Version 1b: no annotation, clustered but no tree shown
    if (makeallheatmaps) {
        
        p=pheatmap(matrix_compare, cluster_rows = hclust_out,cluster_cols = hclust_out, 
            fontsize = myfontsize, fontsize_col = myfontsize, fontsize_row = myfontsize, treeheight_row = 0, treeheight_col = 0, legend = F, limits=c(0,1), border_color = NA)
            #annotation_colors = list(colors=col_Dark2[1:max(cutree_out)])))
        print(p)
        ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_7_Regulons-v1b-L.pdf'), 
            plot = p, width=length(pooled_regulons)*myfontsize/.pt*1.1+15, height=length(pooled_regulons)*myfontsize/.pt*1.1+15, units='mm', limitsize = F)
    }
    
    # Version 2 (smaller)
    if (makeallheatmaps) {
        p=pheatmap(matrix_compare, cluster_rows = hclust_out,cluster_cols = hclust_out, 
            annotation_col = cutree_df, annotation_row = cutree_df, annotation_colors = annotation_colors, 
            fontsize = 5, fontsize_col = 5, fontsize_row = 5, treeheight_col = 0, legend=F, annotation_legend=F)
        ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_7_Regulons-v2.pdf'), 
            plot = p, width=184.6*2/3-4-12, height=184.6/2-4, units='mm')
    
        # Now with legend
        p=pheatmap(matrix_compare, cluster_rows = hclust_out,cluster_cols = hclust_out, 
            annotation_col = cutree_df, annotation_row = cutree_df, annotation_colors = annotation_colors, 
            fontsize = 5, fontsize_col = 5, fontsize_row = 5, treeheight_col = 0, legend=T)
        ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_7_Regulons_LEGEND.pdf'), 
            plot = p, width=184.6/2-4, height=184.6/2-4, units='mm')
    }
    
    return(list(cutree_df=cutree_df))
    
    #View(df_compare)
    #ggplot(df_compare)+
    #    geom_tile(aes(x=x, y=y, fill=overlap))
    
    #pheatmap(matrix_compare, cluster_rows = F, cluster_cols = F)
    
}

######################################################################
######################################################################

################################################################################
# Run the actual analysis

# Execute for all patients in a specific analysis

# unique(current_analysis[['ROOIJonly_default']]$annotation_patient_str)
# ANALYSIS_NAME = 'ROOIJonly_RID2l' #  c('ROOIJonly_RID2l', 'HUonly_RID2l')
# ANALYSIS_NAME = 'HUonly_RID2l' #  c('ROOIJonly_RID2l', 'HUonly_RID2l')
# ANALYSIS_NAME = 'TEICHMANN.SP.only_RID2l'

MAX_GENES=1000 # previously 500

if ('run_regulon_step1' %in% desired_command_regulon) {
    
    ANALYSIS_NAME=config$dataset
    
    current_analysis=list()
    current_analysis[[ANALYSIS_NAME]] = 
            LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',ANALYSIS_NAME,'.h5seurat'))
    
    all_patients  = unique(current_analysis[[ANALYSIS_NAME]]$annotation_patient_str)
        
    # Version using loop ..
    # collected_regulon_objects = list()
    # for (current_patient in all_patients) {
    # 
    #     print(paste0('Starting patient ',current_patient))
    #     
    #     current_analysis[[paste0(ANALYSIS_NAME, current_patient)]]=
    #         subset(current_analysis[[ANALYSIS_NAME]], annotation_patient_str == current_patient)
    #     
    #     current_matrix = 
    #         current_analysis[[paste0(ANALYSIS_NAME, current_patient)]]@assays$RNA@data
    #         # expression matrix now in
    #         # current_anaysis_temp@assays$RNA@data
    # 
    #     collected_regulon_objects[[ANALYSIS_NAME]][[current_patient]] = 
    #         giveMeRegulons_SeuratVersion(run_name=paste0(ANALYSIS_NAME,'_',current_patient),base_dir=base_dir,current_matrix=current_matrix)
    # }
    
    # Only take along genes that are present in 20% of all pooled cells
    current_matrix_gene_filter = rowMeans(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data>0)>.2
    genes_to_take_along = rownames(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data)[current_matrix_gene_filter]
    whole_SeuratData_filtered = subset(current_analysis[[ANALYSIS_NAME]], features = genes_to_take_along)
    print(paste0('Taking along ',sum(current_matrix_gene_filter), ' genes for analysis (20% cutoff).'))
    
    # Version that is parallelized
    collected_regulon_objects = list()
        collected_regulon_objects[[ANALYSIS_NAME]] =
        mclapply(X = all_patients, FUN = function(current_patient) {
        # for (current_patient in all_patients) {
                
                print(paste0('Starting patient ',current_patient))
                
                current_analysis_for_patient=
                    subset(whole_SeuratData_filtered, annotation_patient_str == current_patient)
                
                # if (!is.null(regulon_cell_selection)) {
                #     if (is.null(current_analysis_for_patient$regulon_cell_selection)) {
                #         
                #     }
                # }
                
                current_matrix = 
                    current_analysis_for_patient@assays$RNA@data
                    # expression matrix now in
                    # current_anaysis_temp@assays$RNA@data
                
                # collected_regulon_objects[[ANALYSIS_NAME]][[current_patient]] = giveMeRegulons_SeuratVersion(run_name=paste0(ANALYSIS_NAME,'_',current_patient), base_dir=base_dir,current_matrix=current_matrix, strip_object = T, MAX_GENES = MAX_GENES, BIGB=30) }
                
                return( 
                    giveMeRegulons_SeuratVersion(run_name=paste0(ANALYSIS_NAME,'_',current_patient),
                        base_dir=base_dir,current_matrix=current_matrix, 
                        strip_object = T, MAX_GENES = MAX_GENES, BIGB=10)
                        # strip_object = T; BIGB=30; run_name=paste0(ANALYSIS_NAME,'_',current_patient)
                )
        
            }, mc.cores = MYMCCORES)
    names(collected_regulon_objects[[ANALYSIS_NAME]]) = all_patients

    # We can't save the whole thing, because it also contains the expression matrices
    # So we need to create a "light" version without those, perhaps keeping the rownames ..
    # --> not applicable since I built in light version into giveMeRegulons_SeuratVersion
    # # Remove data
    # collected_regulon_objects_small=list()
    # collected_regulon_objects_small[[ANALYSIS_NAME]] =
    #     lapply(collected_regulon_objects[[ANALYSIS_NAME]], remove_data_regulon_object)
    # rm('collected_regulon_objects')
    
    # Save analysis outcome
    save(list='collected_regulon_objects', file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_regulons_per_patient.Rdata'))
    # load(file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_regulons_per_patient.Rdata'))
    
}    
    

################################################################################

# ANALYSIS_NAME='ROOIJonly_RID2l'
# ANALYSIS_NAME='HUonly_RID2l'
# ANALYSIS_NAME='TEICHMANNonly_RID2l'
# ANALYSIS_NAME='TEICHMANN.SP.only_RID2l'
#
# ANALYSIS_NAME='ROOIJonly_RID2l_test'
# ANALYSIS_NAME='original_HCM_SCS_data'
# ANALYSIS_NAME = 'original_HCM_SCS_data_sel'

CUTS_PER_DATASET = c(ROOIJonly_RID2l=2.3,TEICHMANNonly_RID2l=3.423,HUonly_RID2l=3, TEICHMANN.SP.only_RID2l=3)
    # TODO (remove this) TEICHMANNonly_RID2l and HUonly_RID2l not calibrated yet

MIN_GENES=100

if ('XXXXXXX' %in% desired_command_regulon) {
    
    load(file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_regulons_per_patient.Rdata'))
    
    # Create a structure that holds the gene names that were assigned to regulons
    # Make it "flat" for patients, i.e. it will have entries R.P1.R.2 = Rooij patient 1 regulon 2.
    regulon_gene_names=
        do.call(c, 
            lapply(collected_regulon_objects[[ANALYSIS_NAME]], function(x) {
                if (class(x)!='list') {return(NA)}
                if (length(x$rownames_cor_out_selected_2)<MIN_GENES) {return(NA)}
                names(x$the_regulons) = paste0('R.',1:length(x$the_regulons)); x$the_regulons}
                )
        )
        # sizes of regulons:
        # sapply(regulon_gene_names, length)
    regulon_gene_names=regulon_gene_names[!is.na(regulon_gene_names)]
    
    save(list='regulon_gene_names', file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_regulons_per_patient_geneNamesOnly.Rdata'))
    
    # Now make the little comparison again
    
    overlap_output=
        regulon_overlap_heatmap(pooled_regulons = regulon_gene_names, base_dir=base_dir, run_name =ANALYSIS_NAME, 
            MYTREECUTTINGHEIGHT = CUTS_PER_DATASET[ANALYSIS_NAME], myfontsize = 7)
        #regulon_overlap_heatmap(pooled_regulons = regulon_gene_names, base_dir=base_dir, run_name =ANALYSIS_NAME, 
        #    MYTREECUTTINGHEIGHT = CUTS_PER_DATASET[ANALYSIS_NAME], myfontsize = 3)

    save(list='overlap_output', file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_regulons_per_patient_overlapData.Rdata'))
    
}    

################################################################################
# Just plot all dendrograms for the separate patients/donors

if ('XXXXXXX' %in% desired_command_regulon) {
    
    for (ANALYSIS_NAME in c('HUonly_RID2l', 'TEICHMANNonly_RID2l','ROOIJonly_RID2l')) {
    
        load(file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_regulons_per_patient.Rdata'))
        
        donors=names(collected_regulon_objects[[ANALYSIS_NAME]])
        
        for (DONOR in donors) {
            
            current_hclust = collected_regulon_objects[[ANALYSIS_NAME]][[DONOR]]$hclust_out
            chosen_cutoff = collected_regulon_objects[[ANALYSIS_NAME]][[DONOR]]$nCluster
            cutree_out = cutree(collected_regulon_objects[[ANALYSIS_NAME]][[DONOR]]$hclust_out, k = chosen_cutoff)
            df_assignements = as.data.frame(cutree_out)
            names(df_assignements)='assignment'
        
            p=pheatmap(df_assignements, cluster_rows = current_hclust, cluster_cols = F, 
                fontsize_row = 1, color = col_vector_60, cellwidth = 10, legend = F, labels_col=element_blank())
        
            nr_rows=length(collected_regulon_objects[[ANALYSIS_NAME]][[DONOR]]$rownames_cor_out_selected_2)
            ggsave(filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_Regulons_DendrogramAss_',DONOR,'.pdf'), plot = p, width=50, height=nr_rows/2, units='mm')
        }
        
    }
}
    
################################################################################
# Analyze overlap in regulons between the different patients

# ANALYSIS_NAME = "ROOIJonly_RID2l"
# ANALYSIS_NAME = "HUonly_RID2l"
# ANALYSIS_NAME = "TEICHMANNonly_RID2l"
# ANALYSIS_NAME = "TEICHMANN.SP.only_RID2l"

if ('XXXXXXX' %in% desired_command_regulon) {
    
    # load overlap_output
    load(file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_regulons_per_patient_overlapData.Rdata'))
    load(file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_regulons_per_patient_geneNamesOnly.Rdata'))
    
    core_regulons=list(); core_regulons_extended=list()
    core_regulons_gene_info_table=list()
    
    cutree_df = overlap_output$cutree_df
    pooled_regulons=regulon_gene_names
    for (group_idx in unique(cutree_df$group)) {
        
        # group_idx=4
        
        # rownames(cutree_df[cutree_df$group==group_idx,,drop=F])
        
        # Now calculate how often each gene occurs in the joined regulons
        
        regulon_group = rownames(cutree_df[cutree_df$group==group_idx,,drop=F])
        
        genes_in_current_group = unique(unlist(pooled_regulons[regulon_group]))
        gene_table_regulon_group = data.frame(sapply(regulon_group, function(x) {1*(genes_in_current_group %in% pooled_regulons[[x]])}))
        rownames(gene_table_regulon_group) = genes_in_current_group
        
        gene_table_regulon_group$total_occurence = apply(gene_table_regulon_group, 1, sum)
        gene_table_regulon_group = gene_table_regulon_group[order(gene_table_regulon_group$total_occurence,decreasing = T),]
        
        # Check in how many of the contributing regulons the gene is found (this is equal to # patients, since
        # regulons from the same patient cannot have the same gene)
        core_regulons[[paste0('s.R.',group_idx)]] = rownames(gene_table_regulon_group)[gene_table_regulon_group$total_occurence>=3]
        core_regulons_gene_info_table[[paste0('s.R.',group_idx)]] = gene_table_regulon_group
        # add extended regulons (all genes belonging to that group)
        core_regulons_extended[[paste0('s.R.',group_idx)]] = rownames(gene_table_regulon_group)
        
        print(paste0('s.R.',group_idx,': ',toString(core_regulons[[paste0('s.R.',group_idx)]])))
        
    }
    
    # show lengths
    core_regulons_length = sapply(core_regulons, length)
    core_regulons_length

    save(list='core_regulons_extended', file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_extended.Rdata'))
    
    # Create output
    #matrix_core_regulons_padded = t(plyr::ldply(core_regulons, rbind))
    #df_core_regulons_padded = data.frame(matrix_core_regulons_padded[-1,])
    #colnames(df_core_regulons_padded)=matrix_core_regulons_padded[1,]
    #openxlsx::write.xlsx(x= df_core_regulons_padded, file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons.xlsx'))
    
    #df_core_regulons_padded_shortname=shorthand_cutname_table(df_core_regulons_padded)
    #openxlsx::write.xlsx(x= df_core_regulons_padded_shortname, file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_shortname.xlsx'))
    
    #####
    ## Part 2 ## 
    # We'd also like to more clearly define which genes are typical per regulon
    
    # This involves LOADING THE ANALYSIS because some correlations need to be re-calculated
    # (I don't save correlations because that's too much data.)
    if (!exists('current_analysis')) {current_analysis = list()}
    current_analysis[[ANALYSIS_NAME]] =
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',ANALYSIS_NAME,'.h5seurat'))
    
    # Note this currently requires the previous section to be run also ..
    
    # So perhaps basically repeat regulon analyis at the shared regulon level
    # such that we can determine which genes are most correlated to the others
    # I SHOULD MAKE THIS INTO A FUNCTION
    df_core_regulon_overview_list=list()
    for (core_reg_name in names(core_regulons)) {
        
        if (length(core_regulons[[core_reg_name]])<1) { break }
        
        # core_reg_name='s.R.4'
        
        all_patients = unique(current_analysis[[ANALYSIS_NAME]]$annotation_patient_str)
        #all_patients = names(collected_regulon_objects$TEICHMANNonly_RID2l)
        df_reg_genescores = data.frame(gene=character(), corr=numeric())
        for (current_patient in all_patients) {
        
            print(paste0('Calculating ',core_reg_name, ' - ',current_patient))
            
            patient_selection = current_analysis[[ANALYSIS_NAME]]$annotation_patient_str==current_patient
            
            # current_analysis[[ANALYSIS_NAME]]@assays$RNA@counts[core_regulons[[core_reg_name]], ]
            # Let's do this per patient again
            cor_out = cor(t(as.matrix(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data[core_regulons[[core_reg_name]], patient_selection])))
            cor_reg_avg = colMeans(cor_out)
            #View(cor_reg_avg)
            
            df_reg_genescores=rbind(df_reg_genescores, data.frame(gene=names(cor_reg_avg), corr=cor_reg_avg))
            
        }
        df_reg_genescores_sum       = aggregate(x = list(median_corr=df_reg_genescores$corr), by = list(gene=df_reg_genescores$gene), median)
        df_reg_genescores_sum$nr_patients = core_regulons_gene_info_table[[core_reg_name]][df_reg_genescores_sum$gene,]$total_occurence
        rownames(df_reg_genescores_sum) = df_reg_genescores_sum$gene
        
        # Sort first by patients, then by correlation (median)
        df_reg_genescores_sum$nr_patients <- ordered(df_reg_genescores_sum$nr_patients, levels = sort(unique(df_reg_genescores_sum$nr_patients), decreasing = T))
        df_core_regulon_overview_list[[core_reg_name]]=
            do.call(rbind, lapply(split(df_reg_genescores_sum, df_reg_genescores_sum$nr_patients), function(x) x[order(x$median_corr,decreasing = T),]))
        
    }
    
    # Create output
    core_regulons_sorted = lapply(df_core_regulon_overview_list, function(x) {x$gene})
    matrix_core_regulons_padded = t(plyr::ldply(core_regulons_sorted, rbind))
    df_core_regulons_padded = data.frame(matrix_core_regulons_padded[-1,])
    colnames(df_core_regulons_padded)=matrix_core_regulons_padded[1,]
    openxlsx::write.xlsx(x= df_core_regulons_padded, file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons.xlsx'), overwrite = T)
    
    df_core_regulons_padded_shortname=as.data.frame(shorthand_cutname_table(df_core_regulons_padded))
    openxlsx::write.xlsx(x= df_core_regulons_padded_shortname, file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_shortname.xlsx'), overwrite = T)

    # More extensive
    openxlsx::write.xlsx(x= df_core_regulon_overview_list, file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_withInfo.xlsx'), overwrite = T)
    
    # Rdata output
    save(list = 'core_regulons_sorted', file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_sorted.Rdata'))
        # load(file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_sorted.Rdata'))
    core_regulons_sorted_shortname = shorthand_cutname_table(core_regulons_sorted)
    save(list = 'core_regulons_sorted_shortname', file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_sorted_shortname.Rdata'))
    
    # And export for HOMER
    for (reg_name in names(core_regulons_sorted_shortname)) { 
        write.table(x = core_regulons_sorted_shortname[[reg_name]], file = paste0(base_dir,'Homer/',ANALYSIS_NAME,'_HOMER_core_regulons_',reg_name,'.txt'), quote = F, row.names = F, col.names = F)
    }
    # Add background genes
    background_genes = 
        # Treshold: at least 5% of all cells
        shorthand_cutname(rownames(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data[rowMeans(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data>0)>.05,]))
        # Treshold: At least in 2 cells
        # shorthand_cutname(rownames(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data[rowSums(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data>0)>2,]))
    write.table(x = background_genes, file = paste0(base_dir,'Homer/',ANALYSIS_NAME,'_HOMER_backgroundGenes.txt'), quote = F, row.names = F, col.names = F)    
    save(list = 'background_genes', file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_background_genes.Rdata'))
}    
    
################################################################################

# Calculate regulon expression (composite)

if (F) {
    
    # X=shorthand_cutname(core_regulons$s.R.4)
    
    core_regulons_ = lapply(core_regulons[sapply(core_regulons,length)>0],shorthand_cutname)
    plotlist=lapply(1:length(core_regulons_),
        function(X) {shorthand_seurat_custom_expr(seuratObject = current_analysis[[ANALYSIS_NAME]], 
                                                    gene_of_interest = core_regulons_[[X]], textsize=6, pointsize=.5, 
                                                    custom_title = names(core_regulons_)[X], mymargin = .5)})
    p=wrap_plots(plotlist, nrow = 2)
    
    ggsave(filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_Regulons_UMAP_compositeExpr.pdf'), 
        plot = p, width=30*3, height=30*2, units='mm') # 184.6/3*2-4
    
    # 2nd version of plots
    p=wrap_plots(plotlist, nrow = 1)
    ggsave(filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_Regulons_UMAP_compositeExpr.pdf'), 
        plot = p, width=184/3*2-4, height=20*1, units='mm') # 184.6/3*2-4
    
    
}

################################################################################
# Further analyse the Homer results


if (F) {
    homerData_dflist=list()
    for (idx in 1:6) {
    
        homerData_dflist[[idx]] =
            read.csv(paste0(base_dir, 'Homer/output_s.R.',idx,'/knownResults.txt'), sep = '\t', header = 1)
        colnames(homerData_dflist[[idx]]) = gsub(pattern = '[0-9]*',replacement = '',x = colnames(homerData_dflist[[idx]]))
        homerData_dflist[[idx]]$regulon = paste0('s.R.',idx)
    
    }    
    homerData_df = Reduce(x = homerData_dflist, f = rbind)
    
    # Add brief gene names
    homerData_df$Motif.Name_short = 
        gsub(pattern = '\\(.*$', replacement = '', x = homerData_df$Motif.Name)
    # Better way
    homerData_df$Motif.Name_short2 = 
        gsub(pattern = '^[^/]*/[^-]*-', replacement = '', x = homerData_df$Motif.Name)
    homerData_df$Motif.Name_short2 = 
        gsub(pattern = '[-/\\.(].*$', replacement = '', x = homerData_df$Motif.Name_short2)
    homerData_df$Motif.Name_short2 = toupper( homerData_df$Motif.Name_short2 )
        
    homerData_df_sel = homerData_df[homerData_df$P.value<.05&homerData_df$q.value..Benjamini.<.9,]
    #homerData_df_sel = homerData_df[homerData_df$q.value..Benjamini.<.05,]
    #View(homerData_df_sel)
    
    # Now tabulate them
    lapply(split(x = homerData_df_sel, f = homerData_df_sel$regulon), function(x) {x$Motif.Name_short2})
    
    lapply(split(x = homerData_df_sel, f = homerData_df_sel$regulon), function(x) {x[,c('Motif.Name','P.value','q.value..Benjamini.')]})

    homerData_df_sel_list = split(x=homerData_df_sel, f=homerData_df_sel$regulon)
    
    # Plot of TFs
    p_list=lapply(homerData_df_sel_list, function(df) {
        ggplot(df, aes(y=1, x=Motif.Name_short2))+
            #geom_tile(aes(fill=q.value..Benjamini.))+
            geom_tile(width=1, height=1, aes(fill=q.value..Benjamini.))+
            geom_point(aes(size=-log10(P.value)))+
            theme_bw()+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  aspect.ratio = 1/dim(df)[1], legend.position='none',
                  axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                  panel.border = element_blank(),
                  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                  plot.margin = margin(1,1,0,0,'mm'))+
            give_better_textsize_plot(6)+
            scale_fill_gradientn(colors=c('blue','white','red'),limits=c(0,1))+#rainbow_colors)
            xlab(element_blank())+
            scale_size_continuous(range = c(0.25, 2))
        })
    p=wrap_plots(p_list, nrow=1)
    p
    ggsave(filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_Regulons_HOMER_TFs.pdf'), plot = p, 
           width=184.6*2/3-4, height=25, units='mm')
    
    some_cols = colorRampPalette(brewer.pal(n = 7, name ="Accent"))(7)
    p_list=lapply(1:length(homerData_df_sel_list), function(i) {
        df=homerData_df_sel_list[[i]]
        ggplot(df, aes(y=1, x=Motif.Name_short2))+
            geom_tile(width=1, height=1, fill=some_cols[i])+
            #geom_point(aes(size=-log10(P.value)))+
            theme_bw()+
            theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  aspect.ratio = .33/dim(df)[1], legend.position='none',
                  axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                  axis.ticks.x=element_blank(),axis.text.x=element_blank(),
                  panel.border = element_blank(),
                  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                  plot.margin = margin(1,1,0,0,'mm'))+
            #scale_fill_gradientn(colors=c('white','red'))+#rainbow_colors)
            xlab(element_blank())
        })
    p=wrap_plots(p_list, nrow=1)
    p
    ggsave(filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_Regulons_HOMER_TFs_annotation.pdf'), plot = p, 
           width=184.6*2/3-4, height=10, units='mm')
    
    
    +#, aes(fill=-log10(Pvalue)))+
        geom_point(aes(size=-log10(P.value), fill=q.value..Benjamini.), shape=21)+
        scale_y_continuous(breaks = 1, expand = c(0,0)) + 
        #scale_x_discrete(expand = c(0,0))+
        theme_bw()+give_better_textsize_plot(6)+
        theme(# axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
              legend.position = 'none',
              axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
              axis.ticks.y=element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              aspect.ratio = TOPX,
              plot.margin = margin(1,1,0,0,'mm'))+
        xlab(element_blank())+ylab(element_blank())+
        #xlab(names(GO_regulons_summary_list)[i])+
        scale_size_continuous(range = c(0.25, 2))#+
        coord_flip()
        
    
}

################################################################################

# Run GO analyses

if (F) {

    source('/Users/m.wehrens/Documents/git_repos/SCS_Joep/Functions/functions_mw_copy/MW_GOKEGG_analysis.R')
    # file.edit('/Users/m.wehrens/Documents/git_repos/SCS_Joep/Functions/functions_mw_copy/MW_GOKEGG_analysis.R')
    
    ANALYSIS_NAME = 'ROOIJonly_RID2l'
    
    # Requires 
    # background_genes, all_genes
    # all_genes is used to generate the conversion table
    all_genes = shorthand_cutname( rownames(current_analysis[[ANALYSIS_NAME]]@assays$RNA@counts), PART1OR2 = 1 )
    # Background treshold: at least 5% of all cells 
    background_genes = shorthand_cutname(
        rownames(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data[
            rowMeans(current_analysis[[ANALYSIS_NAME]]@assays$RNA@data>0)>.01,]
            ), PART1OR2 = 1)
    # Export for additional testing
    write.table(x = background_genes, file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_background_table.txt'), 
                row.names = F, col.names = F, quote = F)
    
    # Perform GO analysis on regulons
    # ==
    
    # Load 'm first
    load(paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_sorted.Rdata'))    
    
    # Analysis
    config_GO=list(geneIdentifierType='ensembl_id', species='human')
    core_regulons__ = shorthand_cutname_table( core_regulons[sapply(core_regulons, length)>0] , PART1OR2 = 1 )
    termGSC_regulons = generate_termGSC(config_GO, includeChildTerms=F)
    GO_all_regulons = lapply(names(core_regulons__), function(regulon_name) {
        print(paste0('Analyzing regulon ',regulon_name))
        current_genes_query = core_regulons__[[regulon_name]]
        GO_out = analyzeGeneOntology_MW(config = config_GO, all_genes=all_genes, background_genes=background_genes, 
                            genes_query=current_genes_query, pathwayPCutoff=0.05, GOKegg='GO', 
                            includeChildTerms=F, termGSC = termGSC_regulons)
        return(GO_out)        
    })
    names(GO_all_regulons)=names(core_regulons__)
    
    # Now create a summary parameter that we'll use for plotting
    TOPX=4
    GO_regulons_summary_list = lapply(names(GO_all_regulons), function(n) {x=GO_all_regulons[[n]]; x$regulon=n; return(x[1:TOPX,])})
    GO_regulons_summary = Reduce(f = rbind, x= GO_regulons_summary_list)
    max_str_length = median(nchar(GO_regulons_summary$Term)*2)
    GO_regulons_summary$Term_short =
        sapply(GO_regulons_summary$Term, function(x) {
                if (nchar(x)>max_str_length) {
                    paste0(substring(x, 1, max_str_length),' ...')
                } else {x}
            })
    # put back to get abbreviated names
    GO_regulons_summary_list = split(x = GO_regulons_summary, f = GO_regulons_summary$regulon)
    
    # Plot 1
    ggplot(GO_regulons_summary, aes(x=Term_short, y=-log10(Pvalue), fill=regulon))+
        geom_bar(stat='identity')+
        facet_grid(cols=vars(regulon), scales = 'free_x')+
        theme_bw()+give_better_textsize_plot(8)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none')
    
    # Some colors
    #some_cols = hue_pal()(6)
    #some_cols = brewer.pal(n = 6, name = 'Accent')
    some_cols = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(7)
    
    # Plot 2
    plot_list=
        lapply(1:length(GO_regulons_summary_list), function(i) {
            df=GO_regulons_summary_list[[i]]
            df$Term = factor(df$Term, levels=df[order(df$Pvalue),]$Term)
            df$Term_short = factor(df$Term_short, levels=df[order(df$Pvalue),]$Term_short)
            ggplot(df, aes(x=Term_short, y=-log10(Pvalue)))+
                geom_bar(stat='identity', fill=some_cols[i], color='black', size=.5)+
                theme_bw()+give_better_textsize_plot(4)+
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                      legend.position = 'none',
                      plot.margin = margin(1,1,0,0,'mm'))+
                ylab(element_blank())+xlab(element_blank())#+ylab('-log10(p)')+
        })
    p=wrap_plots(plot_list, nrow=1)
    p
    ggsave(filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_Regulons_GO_terms.pdf'), plot = p, 
           width=184.6*2/3-4, height=80, units='mm')
    
    # Plot 3
    plot_list=
        lapply(GO_regulons_summary_list, function(df) {
            df$Term = factor(df$Term, levels=df[order(df$Pvalue),]$Term)
            df$Term_short = factor(df$Term_short, levels=df[order(df$Pvalue),]$Term_short)
            ggplot(df, aes(x=Term_short, y=-log10(Pvalue)))+
                geom_bar(stat='identity')+
                theme_bw()+give_better_textsize_plot(8)+
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none')+
                ylab(element_blank())+xlab(element_blank())+#+ylab('-log10(p)')+
                coord_flip()
        })
    wrap_plots(plot_list, ncol=1)
    
    # Plot 4
    plot_list=
        lapply(1:length(GO_regulons_summary_list), function(i) {
            df=GO_regulons_summary_list[[i]]
            df$Term = factor(df$Term, levels=df[order(df$Pvalue),]$Term)
            #df$Term_short = factor(df$Term_short, levels=rev(df[order(df$Pvalue),]$Term_short))
            ggplot(df, aes(x=factor(Term_short, levels=rev(Term_short)), y=1, size=-log10(Pvalue)))+
                geom_tile(width=1, height=1, fill='white')+#, aes(fill=-log10(Pvalue)))+
                geom_point(shape=21, color='black', fill=some_cols[i], )+
                scale_y_continuous(breaks = 1, expand = c(0,0)) + 
                #scale_x_discrete(expand = c(0,0))+
                theme_bw()+give_better_textsize_plot(6)+
                theme(# axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                      legend.position = 'none',
                      axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                      axis.ticks.y=element_blank(),
                      panel.border = element_blank(),
                      panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                      aspect.ratio = TOPX,
                      plot.margin = margin(1,1,0,0,'mm'))+
                xlab(element_blank())+ylab(element_blank())+
                #xlab(names(GO_regulons_summary_list)[i])+
                scale_size_continuous(range = c(0.25, 2))+
                coord_flip()
                    
                    #+ylab('-log10(p)')
                #coord_flip()+ylim(c(0,2))
        })
    p=wrap_plots(plot_list, ncol=1)
    p
    ggsave(filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_Regulons_GO_terms_points.pdf'), plot = p, 
           width=184.6*2/3-4, height=80, units='mm')
    
    
    # ggplot(GO_regulons_summary, aes(x=factor(Term_short, levels=Term_short), y=1, size=-log10(Pvalue), color=regulon))+
    #     geom_point()+
    #     theme_minimal()+
    #     coord_flip()+
    #     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    #                   legend.position = 'none',
    #                   plot.margin = margin(1,1,0,0,'mm'))
    
}

################################################################################
# GO analysis by Homer (of regulons)
# Not used.

if (F) {
 
    # Get all data in one frame
    homerGOData_dflist=list()
    for (idx in 1:6) {
    
        homerGOData_dflist[[idx]] =
            read.csv(paste0(base_dir, 'Homer/output_s.R.',idx,'/biological_process.txt'), sep = '\t', header = 1)
        homerGOData_dflist[[idx]]$regulon = paste0('s.R.',idx)
    
    }    
    homerGOData_df = Reduce(x = homerGOData_dflist, f = rbind)
    
    homerGOData_df_selTop10 = Reduce(x = lapply(homerGOData_dflist,function(x){x[1:10,]}), f = rbind)
       
    # Select some data of interest
    ggplot(homerGOData_df_selTop10, aes(x=Term, y=Enrichment))+
        geom_bar(stat='identity')+
        facet_grid(cols = vars(regulon))
    # lapply(homerGOData_dflist, function(x) {})
}


################################################################################
# dataset to dataset comparison of the regulons
# approach 1: compare the core regulons

name_conversion = c(ROOIJonly_RID2l='R', TEICHMANNonly_RID2l='T', TEICHMANN.SP.only_RID2l='TS', HUonly_RID2l='H')
SP_SWITCH='.SP.'
# SP_SWITCH=''

if (F) {
    
    sorted_regulon_collection=list()
    for (ANALYSIS_NAME in c('ROOIJonly_RID2l', paste0('TEICHMANN',SP_SWITCH,'only_RID2l'))) { #'TEICHMANNonly_RID2l')) {
        
        # ANALYSIS_NAME='TEICHMANN.SP.only_RID2l'
        
        load(file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_sorted.Rdata'))
        sorted_regulon_collection[[ANALYSIS_NAME]] = core_regulons_sorted
        names(sorted_regulon_collection[[ANALYSIS_NAME]]) = paste0(name_conversion[[ANALYSIS_NAME]], '.', names(sorted_regulon_collection[[ANALYSIS_NAME]]))
        
    }
    
    # now pool 'm
    newnames = unlist(sapply(sorted_regulon_collection, names)) # unlist screws up names, so store 'm
    sorted_regulon_pool = unlist(sorted_regulon_collection, recursive = F)
    names(sorted_regulon_pool) = newnames
    
    # Upset plot
    upset(fromList(sorted_regulon_pool), order.by = "freq", empty.intersections = "on")
    
    regulon_overlap_heatmap(pooled_regulons = sorted_regulon_pool, base_dir = base_dir, run_name = paste0('pool.R.T',SP_SWITCH), MYTREECUTTINGHEIGHT = 1.66, myfontsize = 6)
    
    if (SP_SWITCH=='') { REG_COMBI = c('R.s.R.2','T.s.R.4')
        } else {REG_COMBI = c('R.s.R.2','TS.s.R.4') }
    
    # sorted_regulon_pool$T.s.R.4 and sorted_regulon_pool$R.s.R.2 seem similar
    vlist=list(sorted_regulon_pool[[REG_COMBI[1]]], sorted_regulon_pool[[REG_COMBI[2]]]);names(vlist)=REG_COMBI
    venn_simple_plot_mw(vlist)
    
    # What are top genes?
    shorthand_cutname(sorted_regulon_pool[[REG_COMBI[1]]][1:20])
    shorthand_cutname(sorted_regulon_pool[[REG_COMBI[2]]][1:20])
    sum(shorthand_cutname(sorted_regulon_pool$T.s.R.4[1:20]) %in% shorthand_cutname(sorted_regulon_pool$R.s.R.2[1:20]))
        # 4/10 overlap in top 10, 9/20 in top 20
    shorthand_cutname(sorted_regulon_pool$R.s.R.2[(sorted_regulon_pool$R.s.R.2 %in% sorted_regulon_pool$T.s.R.4)][1:10])
    shorthand_cutname(sorted_regulon_pool$T.s.R.4[!(sorted_regulon_pool$T.s.R.4 %in% sorted_regulon_pool$R.s.R.2)][1:10])
    shorthand_cutname(sorted_regulon_pool$R.s.R.2[!(sorted_regulon_pool$R.s.R.2 %in% sorted_regulon_pool$T.s.R.4)][1:10])
    
    # Overview, so print top 20 of both, show whether either both or specific
    symbol_map = c('','X')
    df_R = data.frame(genes_R=sorted_regulon_pool[[REG_COMBI[1]]][1:20], genes_R_short=shorthand_cutname(sorted_regulon_pool[[REG_COMBI[1]]][1:20]), overlapping_R=symbol_map[1+sorted_regulon_pool[[REG_COMBI[1]]][1:20] %in% sorted_regulon_pool[[REG_COMBI[2]]][1:20]*1],
                      genes_T=sorted_regulon_pool[[REG_COMBI[2]]][1:20], genes_T_short=shorthand_cutname(sorted_regulon_pool[[REG_COMBI[2]]][1:20]), overlapping_T=symbol_map[1+sorted_regulon_pool[[REG_COMBI[2]]][1:20] %in% sorted_regulon_pool[[REG_COMBI[1]]][1:20]*1])
    openxlsx::write.xlsx(x= df_R, file = paste0(base_dir,'Rplots/customALL',SP_SWITCH,'_overlap_',paste0(REG_COMBI,collapse='-'),'.xlsx'), overwrite = T)
    
    # For reference, a somewhat lower overlap
    vlist=list(sorted_regulon_pool[[REG_COMBI[1]]], sorted_regulon_pool[[REG_COMBI[2]]]); names(vlist)=REG_COMBI
    venn_simple_plot_mw(list(R.s.R.5=sorted_regulon_pool$R.s.R.5, T.s.R.7=sorted_regulon_pool$T.s.R.7))
    
    
}

################################################################################
# dataset to dataset comparison of the regulons
# approach 2: a comparison could also be made at the single patient level

if (F) {
 
    
    collected_regulon_collectionDatasets_flat=list()
    for (ANALYSIS_NAME in c('ROOIJonly_RID2l', 'TEICHMANNonly_RID2l')) {
        
        # ANALYSIS_NAME='TEICHMANNonly_RID2l'
        
        load(file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_regulons_per_patient.Rdata'))

        for (donor_name in names(collected_regulon_objects[[ANALYSIS_NAME]])) {
            for (reg_idx in 1:length(collected_regulon_objects[[ANALYSIS_NAME]][[donor_name]]$the_regulons)) {
                collected_regulon_collectionDatasets_flat[[paste0(donor_name, '.R.',reg_idx)]] = 
                    collected_regulon_objects[[ANALYSIS_NAME]][[donor_name]]$the_regulons[[reg_idx]]
            }
        }
    }
    
    # Now make the little comparison again
    overlap_output=
        regulon_overlap_heatmap(pooled_regulons = collected_regulon_collectionDatasets_flat, base_dir=base_dir, run_name ='detailpool_R.T', 
            MYTREECUTTINGHEIGHT = 2, myfontsize = 3)
       
    # It appears there is only one group of regulons that are clearly mixed from the Rooij and Teichmann datasets;
    # This is consistent with the core regulons only showing two core regulons overlapping from respective datasets
    # --> now just check whether indeed that mixed group found here matches those two core regulons
    #
    # I did this manually, so T.s.R.4 = 
    Tsep_Rnames = c('T.D2.R.4', 'T.D5.R.4', 'T.D1.R.4', 'T.D4.R.4', 'T.D11.R.6', 'T.D6.R.6', 'T.D7.R.2', 'T.H5.R.10', 'T.H2.R.15', 'T.H3.R.9', 'T.H4.R.9')
    # R.s.R.2 =
    Rsep_Rnames = c('R.P3.R.7', 'R.P4.R.9', 'R.P5.R.6', 'R.P1.R.2', 'R.P2.R.2', 'R.P3.R.4', 'R.P3.R.3', 'R.P5.R.2')
    #
    # The shared group identified here is:
    TR_Rnames = c('T.D7.R.2', 'T.D2.R.4', 'T.D5.R.4', 'T.D6.R.6', 'T.D4.R.4', 'T.D1.R.4', 'T.D11.R.6', 'T.H5.R.10', 'T.H2.R.15', 'T.H3.R.9', 'T.H4.R.9', 'R.P4.R.9', 'R.P5.R.6', 'R.P3.R.7', 'R.P1.R.2', 'R.P2.R.2')
    
    venn_simple_plot_mw(list(grouped=TR_Rnames, T_separate=Tsep_Rnames))
    venn_simple_plot_mw(list(grouped=TR_Rnames, R_separate=Rsep_Rnames))
}

################################################################################

if (F) {
    
    # Comparison with earlier version of this analysis
    # For reference, this was what we got from the previous analysis:
    # load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/Previous_analysis_for_reference/JoepAnalysis_Regulons.Rdata')
    # See regulons_README_objects_saved
    # See exporatory analyses for comparison earlier and current data analyses
    
    
    # Sizes of old regulons (ie comparison with old data)
    old_regulon_sizes = lapply(joined_regulons, length)
    ggplot(data.frame(old_regulon_size=unlist(old_regulon_sizes), name=names(old_regulon_sizes)))+
        geom_bar(aes(x=name,y=old_regulon_size), stat='identity')+theme_bw()+
        coord_flip()
    
    # Sizes of current regulons
    new_regulon_sizes = lapply(regulon_gene_names, length)
    ggplot(data.frame(old_regulon_size=unlist(new_regulon_sizes), name=names(new_regulon_sizes)))+
        geom_bar(aes(x=name,y=old_regulon_size), stat='identity')+theme_bw()+
        coord_flip()
    
    ggplot(data.frame(old_regulon_size=c(unlist(new_regulon_sizes),unlist(old_regulon_sizes)), 
                      name=c(names(new_regulon_sizes),names(old_regulon_sizes)),
                      version=rep(c('new','old'),times=c(length(new_regulon_sizes), length(old_regulon_sizes)))
                ))+
        geom_bar(aes(x=name,y=old_regulon_size, fill=version), stat='identity')+theme_bw()+
        coord_flip()
    
    
    #######
    
    # comparison core regulons
    
    core_regulons_length 
    old_core_regulons_length = sapply(shared_regulon_core_list, length)
    
    ggplot(data.frame(old_regulon_size=c(unlist(core_regulons_length),unlist(old_core_regulons_length)), 
                      name=c(names(core_regulons_length),names(old_core_regulons_length)),
                      version=rep(c('new','old'),times=c(length(core_regulons_length), length(old_core_regulons_length)))
                ))+
        geom_bar(aes(x=name,y=old_regulon_size, fill=version), stat='identity')+theme_bw()+
        coord_flip()
    
    
    # All genes old analysis
    unique(unlist(joined_regulons)) # ±1000 genes
    unique(unlist(regulon_gene_names)) # ±771 genes

}

################################################################################

print('Done; end of regulon script.')

################################################################################



