
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

######################################################################
# libraries

# local base dir
if (exists('LOCAL')) {
    base_dir = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/'
    base_dir_secondary = '/Volumes/workdrive_m.wehrens_hubrecht/R-sessions/2021_SCS_HCM_Seurat/'
    script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/'
} else {
    base_dir = '/hpc/hub_oudenaarden/mwehrens/data/Wang/'
    base_dir_secondary = base_dir # change if desired to put some less important and/or big files
    script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'
}
dir.create(paste0(base_dir_secondary,'Rdata/'))
    
MYMCCORES=8 # Can be overwritten by using CORES=X argument

library(Seurat)
library(SeuratDisk)

library(ggdendro)
library(gplots) # heatmap.2
library(pheatmap)

library(randomcoloR)
library(ggplot2)
library(scales)
library(gridExtra)

library(parallel)

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector_60 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

source(paste0(script_dir, 'Functions/MW_regulon_analysis_COPY.R'))
source(paste0(script_dir, 'Functions/MW_general_functions.R'))
source(paste0(script_dir, 'Functions/Load-Pool-Scale_Simple_MW.R'))
source(paste0(script_dir, 'Functions/MW_regulon_analysis_supportFns_COPIES.R'))
    # file.edit(paste0(script_dir, 'Functions/MW_general_functions.R'))

# Not necessary
# source('/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Functions/MW_analysis_functions.R')

# Redundant, libs required are now loaded above
if (F) {
    cfg=list()
    cfg$SCRIPT_LOCATION_MW = '/Users/m.wehrens/Documents/git_repos/Tomoseq_mw/'
    cfg$SCRIPT_LOCATION_MW_2 = '/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/'
    cfg$species='human'
    source('/Users/m.wehrens/Documents/git_repos/Tomoseq_mw/sub_functions/MW_load_libraries.R') # doesn't work completely ?
        # file.edit('/Users/m.wehrens/Documents/git_repos/Tomoseq_mw/sub_functions/MW_load_libraries.R')
    source('/Users/m.wehrens/Documents/git_repos/Tomoseq_mw/sub_functions/MW_regulon_analysis.R')
        # file.edit('/Users/m.wehrens/Documents/git_repos/Tomoseq_mw/sub_functions/MW_regulon_analysis.R')
}

########################################################################
# To be able to run this script from the command line on the cluster,
# let's make sure we can pass arguments to it
#
# Note that some sections are executed together.

args = commandArgs(trailingOnly = T)

library(stringr)
if (!exists('desired_command')) {
    if (length(args)==1) {
        desired_command = unlist(str_split(string = args, pattern = '-'))
        print(paste0('Sections to execute: ',paste0(desired_command,collapse = ', ')))
        if(any(grepl(x = desired_command,pattern = 'CORES'))) {MYMCCORES=str_split(string = desired_command[grepl(x = desired_command,pattern = 'CORES')], pattern = '=')[[1]][2];print(paste0('MYMCCORES set to ',MYMCCORES))}
    } else {stop('Please pass 1 string, in the form \'arg1-arg2-arg3\' to give what section(s) to execute.')}
}

######################################################################
######################################################################

# Using Tomo code!
giveMeRegulons_SeuratVersion = function(run_name, base_dir, current_matrix, 
    get_GO_terms=F, strip_object=F, MAX_GENES=NULL,
    required_minCellFraction=.05) {
    
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
        analysis_name = run_name)
    
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
        chosen_cutoff_parameter='p',
        # this parameter should be 'r' or 'p', respectively selection of 
        # genes based on r-value or p-value cutoff (r being the correlation 
        # coefficient value itself)
        #p_or_r_cutoff=.25
        p_or_r_cutoff=0.00001) # 1/100000
        # this parameter provides the cutoff you have chosen earlier,
        # e.g. consider r-values (correlation coefficients) >.35 or <-.35 
        # significant.
    
    # Now, we can already create a correlation heatmap, which shows the structure
    # of correlation between genes. 
    # Per default, I don't show the heatmap here, because we can later determine
    # some parameters and set a cutoff to perform clustering on this correlation
    # matrix.
    regulon_object = MW_determine_regulons_part3(regulon_object, 
                    connectedness_cutoff = 40, 
                    max_genes = MAX_GENES,
                    min_expression_fraction_genes=.1,
                    show_heatmap=F,
                    chosen_cutoff_parameter = 'p')
    
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
    # Choose a cutoff value, e.g.:
    # current_cutoff = regulon_object$auto_cutoff2
    regulon_object = MW_determine_regulons_part5(regulon_object = regulon_object,
        hierarchical_cutoff = NULL) # current_cutoff)
        # With hierarchical_cutoff = NULL, gap_stat will be used
    
    # Now we can add some information about the genes (are they TF, ligand, receptor?)
    data_container=list() # Change this to get TF ligand stuff (should be easily doable)
    warning('See above: add TF/ligand stuff')
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


regulon_overlap_heatmap = function(pooled_regulons, base_dir, run_name, MYTREECUTTINGHEIGHT=2) {
    
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
    p0=pheatmap(matrix_compare, clustering_method = 'ward.D2')
    p0
    ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_7_Regulons_noCats.png'), plot = p0, width=nrow(matrix_compare)*.4, height=nrow(matrix_compare)*.4, units='cm')
    
    hclust_out = hclust(dist(matrix_compare), method='ward.D2')
    plot(as.dendrogram(hclust_out))
    cutree_out = cutree(hclust_out, h = MYTREECUTTINGHEIGHT)
    cutree_df  = as.data.frame(as.factor(cutree_out)); colnames(cutree_df) = c('group')
    #annotation_colors = col_Dark2[1:max(cutree_out)]
    # Create a little heatmap
    annotation_colors = col_vector_60[1:max(cutree_out)]
    names(annotation_colors) = unique(cutree_out)
    annotation_colors=list(group=annotation_colors)
    p=pheatmap(matrix_compare, cluster_rows = hclust_out,cluster_cols = hclust_out, 
        annotation_col = cutree_df, annotation_row = cutree_df, annotation_colors = annotation_colors)
        #annotation_colors = list(colors=col_Dark2[1:max(cutree_out)])))
    print(p)
    ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_7_Regulons.png'), plot = p, width=nrow(matrix_compare)*.4, height=nrow(matrix_compare)*.4, units='cm')
    
    return(list(cutree_df=cutree_df))
    
    #View(df_compare)
    #ggplot(df_compare)+
    #    geom_tile(aes(x=x, y=y, fill=overlap))
    
    #pheatmap(matrix_compare, cluster_rows = F, cluster_cols = F)
    
}

######################################################################
######################################################################

################################################################################

# Execute for all patients in a specific analysis

# unique(current_analysis[['ROOIJonly_default']]$annotation_patient_str)
# ANALYSIS_NAME = 'ROOIJonly_RID2l' #  c('ROOIJonly_RID2l', 'HUonly_RID2l')
# ANALYSIS_NAME = 'HUonly_RID2l' #  c('ROOIJonly_RID2l', 'HUonly_RID2l')

MAX_GENES=1000 # previously 500

if ('run_regulon_step1' %in% desired_command) {
    
    ANALYSIS_NAME=desired_command[2]
    
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
    
    # Version that is parallelized
    collected_regulon_objects = list()
    collected_regulon_objects[[ANALYSIS_NAME]] =
        mclapply(X = all_patients, FUN = function(current_patient) {
                
                print(paste0('Starting patient ',current_patient))
                
                current_analysis_for_patient=
                    subset(current_analysis[[ANALYSIS_NAME]], annotation_patient_str == current_patient)
                
                current_matrix = 
                    current_analysis_for_patient@assays$RNA@data
                    # expression matrix now in
                    # current_anaysis_temp@assays$RNA@data
            
                return( 
                    giveMeRegulons_SeuratVersion(run_name=paste0(ANALYSIS_NAME,'_',current_patient),
                        base_dir=base_dir,current_matrix=current_matrix, strip_object = T, MAX_GENES = MAX_GENES)
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

if ('XXXXXXX' %in% desired_command) {
    
    load(file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_regulons_per_patient.Rdata'))
    
    # Create a structure that holds the gene names that were assigned to regulons
    # Make it "flat" for patients, i.e. it will have entries R.P1.R.2 = Rooij patient 1 regulon 2.
    regulon_gene_names=
        do.call(c, 
            lapply(collected_regulon_objects[[ANALYSIS_NAME]], function(x) {names(x$the_regulons) = paste0('R.',1:length(x$the_regulons)); x$the_regulons})
        )
    
    save(list='regulon_gene_names', file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_regulons_per_patient_geneNamesOnly.Rdata'))
    
    # Now make the little comparison again
    
    overlap_output=
        regulon_overlap_heatmap(regulon_gene_names, base_dir=base_dir, run_name =ANALYSIS_NAME)

    save(list='overlap_output', file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_regulons_per_patient_overlapData.Rdata'))
    
}    
    
################################################################################

if ('XXXXXXX' %in% desired_command) {
    
    # load overlap_output
    load(file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_regulons_per_patient_overlapData.Rdata'))
    load(file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_regulons_per_patient_geneNamesOnly.Rdata'))
    
    core_regulons=list()
    
    cutree_df = overlap_output$cutree_df
    pooled_regulons=regulon_gene_names
    for (group_idx in unique(cutree_df$group)) {
        # rownames(cutree_df[cutree_df$group==group_idx,,drop=F])
        
        # Now calculate how often each gene occurs in the joined regulons
        
        regulon_group = rownames(cutree_df[cutree_df$group==group_idx,,drop=F])
        
        genes_in_current_group = unique(unlist(pooled_regulons[regulon_group]))
        gene_table_regulon_group = data.frame(sapply(regulon_group, function(x) {1*(genes_in_current_group %in% pooled_regulons[[x]])}))
        rownames(gene_table_regulon_group) = genes_in_current_group
        
        gene_table_regulon_group$total_occurence = apply(gene_table_regulon_group, 1, sum)
        
        core_regulons[[paste0('s.R.',group_idx)]] = rownames(gene_table_regulon_group)[gene_table_regulon_group$total_occurence>=3]
    
        print(paste0('s.R.',group_idx,': ',toString(core_regulons[[paste0('s.R.',group_idx)]])))
        
    }
    
    # show lengths
    core_regulons_length = sapply(core_regulons, length)
    
    # Create output
    matrix_core_regulons_padded = t(plyr::ldply(core_regulons, rbind))
    df_core_regulons_padded = data.frame(matrix_core_regulons_padded[-1,])
    colnames(df_core_regulons_padded)=matrix_core_regulons_padded[1,]
    openxlsx::write.xlsx(x= df_core_regulons_padded, file = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons.xlsx'))
    
}

################################################################################

if (F) {
    
    # Comparison with earlier version of this analysis
    # For reference, this was what we got from the previous analysis:
    # load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/Previous_analysis_for_reference/JoepAnalysis_Regulons.Rdata')
    # See regulons_README_objects_saved
    # See exporatory analyses for comparison earlier and current data analyses
    
    
    # Sizes of old regulons 
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
    unique(unlist(regulon_gene_names)) # ±375 genes

}




