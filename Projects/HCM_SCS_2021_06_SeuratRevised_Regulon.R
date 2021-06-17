


# Applying the regulon analysis to all of the different datasets separately,
# separately to each donor/patient.

# Maybe I should use Joep's functions?? --> at least the ones used in the 
# submitted manuscript..



# OK, convenient to again use the already made Seurat objects,
# this time, I will use the ones that aren't corrected, since
# the regulon analysis was appied per patient earlier,
# and i can split it out again


################################################################################
# Load the objects

OBJECTS_TO_ANALYZE = c('ROOIJonly_default', 'HUonly_default')

current_analysis=list()
for (analysis_name in OBJECTS_TO_ANALYZE) {

    # analysis_name = 'HUonly_RID2l'
    
    current_analysis[[analysis_name]] = 
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',analysis_name,'.h5seurat'))

}

################################################################################

# Try first for one patient

# unique(current_analysis[['ROOIJonly_default']]$annotation_patient_str)
CURRENT_PATIENT = 'R.P1'

current_analysis[[paste0('ROOIJonly_default', CURRENT_PATIENT)]]=
    subset(current_analysis[['ROOIJonly_default']], annotation_patient_str == CURRENT_PATIENT)

current_anaysis_temp = 
    current_analysis[[paste0('ROOIJonly_default', CURRENT_PATIENT)]]
    # expression matrix now in
    # current_anaysis_temp@assays$RNA@data

CorrelationMatrix =  
    generateCorrelationMatrix(config=NULL, groupedSCS, groupNames, excludeOutlierCells=T, minCellFraction=0, minCellExpression=0.1, desiredPValue=0.00001, adjustP=T, saveMatrices=F, overrideDir=F, outputMode='pdf', which_genes_to_select=T, filename_precursor='')

# Test run of regulon analysis, using my Tomo seq code

######################################################################

# local base dir
base_dir = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/'

library(ggdendro)
library(gplots) # heatmap.2

source('/Users/m.wehrens/Documents/git_repos/Tomoseq_mw/sub_functions/MW_regulon_analysis.R')
# source('/Users/m.wehrens/Documents/git_repos/Tomoseq_mw/sub_functions/MW_load_libraries.R') # doesn't work completely ?
source('/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/Functions/MW_analysis_functions.R')

cfg=list()
cfg$outputDir='/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/regulon/'

# Calibrate some settings (Rooij)
CURRENT_RUNNAME = paste0('ROOIJ','only_','RID2l')
cell_selection = current_SeuratObject[[CURRENT_RUNNAME]]$annotation_patient_str=='R.P1'
GROUP_NAME = 'test_Rooij_pt1' # choose from names(data_container$groups)
dir.create(paste0(cfg$outputDir,'analysis_',GROUP_NAME,'/regulon/'), recursive = T)

# Calibrate some settings (Hu)
CURRENT_RUNNAME = paste0('HU','only_','RID2l')
cell_selection = current_SeuratObject[[CURRENT_RUNNAME]]$annotation_patient_str=='H.N1'
GROUP_NAME = 'test_Hu_pt1' # choose from names(data_container$groups)
dir.create(paste0(cfg$outputDir,'analysis_',GROUP_NAME,'/regulon/'), recursive = T)

######################################################################

current_SeuratObject=list()
current_SeuratObject[[CURRENT_RUNNAME]] =
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',CURRENT_RUNNAME,'.h5seurat'))

    # Note that we can either take @data or @scaled.data, the latter only contains the 
#

# Select genes first based on how many 
current_SeuratObject[[CURRENT_RUNNAME]]@misc$genes_expr_in_cells =
    rowSums(current_SeuratObject[[CURRENT_RUNNAME]]@assays$RNA@data>0)
current_SeuratObject[[CURRENT_RUNNAME]]@misc$genes_expr_in_cells_fraction =
    current_SeuratObject[[CURRENT_RUNNAME]]@misc$genes_expr_in_cells/
        dim(current_SeuratObject[[CURRENT_RUNNAME]]@assays$RNA@data)[2]
    # length(current_SeuratObject[[CURRENT_RUNNAME]]@misc$genes_expr_in_cells_fraction)
print(paste0('Genes in run: ',sum(current_SeuratObject[[CURRENT_RUNNAME]]@misc$genes_expr_in_cells_fraction>0.1)))
gene_selection = current_SeuratObject[[CURRENT_RUNNAME]]@misc$genes_expr_in_cells_fraction>0.05
#gene_selection = current_SeuratObject[[CURRENT_RUNNAME]]@misc$genes_expr_in_cells_fraction>0.1

regulon_object = MW_determine_regulons_part1(
    expression_matrix = current_SeuratObject[[CURRENT_RUNNAME]]@assays$RNA@data[gene_selection,cell_selection], 
    calculate_p = T,
    outputDir = cfg$outputDir,
    analysis_name = GROUP_NAME)

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
    p_or_r_cutoff=0.00001)
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
# cfg$regulon_dendogram_cutoff[[GROUP_NAME]] = regulon_object$auto_cutoff2
# cfg$regulon_dendogram_cutoff[[GROUP_NAME]] = 3
regulon_object = MW_determine_regulons_part5(regulon_object = regulon_object,
    hierarchical_cutoff = cfg$regulon_dendogram_cutoff[[GROUP_NAME]])

# Now we can add some information about the genes (are they TF, ligand, receptor?)
data_container=list() # Change this to get TF ligand stuff (should be easily doable)
warning('See above: add TF/ligand stuff')
regulon_object = MW_regulon_add_TF_RL_flags(regulon_object, data_container)
# And we can export the regulon analysis to an excel file:
MW_regulon_final_export(regulon_object)

# Additionally, we can export the classification of genes to clusters to an
# excel file.
# Additionally, we can perform a GO analysis on these genes, with all genes
# taken along in the matrix as background (TODO: perhaps reconsider which
# genes to use as background panel??-- perhaps all genes in the expression
# matrix would make more sense). 
# This function also exports the GO term for each of the clusters to an 
# excel file.
regulon_object = MW_determine_regulons_part6(regulon_object = regulon_object)
