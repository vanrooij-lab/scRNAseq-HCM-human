
# Checking out the litvinuka2020 data
# Note that the original h5ad couldn't be loaded using my installation of R
# 
# I used python to load in those h5ad files, and converted them to the "loom" format, 
# and also saved them again as h5ad (there appears to be some issue with versions of both
# R and the hdf5r library that result in a non-descript error).
# 
# The h5ad files that I saved again resulted in memory issues after loading.
# The loom file seems to work.
# See https://satijalab.org/loomR/loomR_tutorial.html for a tutorial on loomR.
#
# Note that I needed both the latest version of R and Python to make this work.
# (Personal note: I used Anaconda environments to make use of the latest version of python
# -- see word document for more info on that.)
#
# Data downloaded from 
# https://www.heartcellatlas.org/
#
# Data described in paper: 
# Litviňuková, M., Talavera-López, C., Maatz, H., Reichart, D., Worth, C. L., Lindberg, E. L., … Teichmann, S. A. (2020). Cells of the adult human heart. Nature, (September). https://doi.org/10.1038/s41586-020-2797-4

library(Seurat)
library(hdf5r)
library(loomR)
library(ggplot2)

# Load some of my standard stuff
cfg=list(); cfg$SCRIPT_LOCATION_MW = '/Users/m.wehrens/Documents/git_repos/Tomoseq_mw/'
cfg$SCRIPT_LOCATION_MW_2 = '/Users/m.wehrens/Documents/git_repos/SCS_RaceID3/'
cfg$species="human"
source('/Users/m.wehrens/Documents/git_repos/Tomoseq_mw/sub_functions/MW_load_libraries.R')

# Attempting to load original data from https://www.heartcellatlas.org/
# This didn't work
# Data_Litvinukova2020 = ReadH5AD('/Users/m.wehrens/Data/2020_12_Litvinukova2020/global_raw.h5ad')
# Data_Litvinukova2020_V_CM = ReadH5AD('/Users/m.wehrens/Data/2020_12_Litvinukova2020/hca_heart_ventricular_CM_raw.h5ad')

# Loading dataset that was saved again after loading original file in python did seem to work
# But R ran out of memory
# Data_Litvinukova2020 = ReadH5AD('/Volumes/workdrive_m.wehrens_hubrecht/data/2020_12_Litvinukova2020/h5ad_vMW/Litvinukova2020_re-exported.h5ad')
# Loading dataset in loom format (from aformentioned python session)
#Data_Litvinukova2020_loom = loomR::connect(filename = '/Volumes/workdrive_m.wehrens_hubrecht/data/2020_12_Litvinukova2020/loom/Litvinukova2020.loom', mode = 'r')
Data_Litvinukova2020_loom = loomR::connect(filename = '/Volumes/workdrive_m.wehrens_hubrecht/data/2020_12_Litvinukova2020/loom/Litvinukova2020_work.loom', mode = 'r+')
    # ! Note that it is relevant to create a working copy of the loom file because it can become corrupted when working on it.
    # - also, loomR seems to be in beta stage of development, some things don't seem to work

# Just looking around in data
Data_Litvinukova2020_loom[['col_attrs']]
Data_Litvinukova2020_loom[['row_attrs']]
Data_Litvinukova2020_loom[['col_attrs']][['cell_type']][1:100]
Data_Litvinukova2020_loom[['row_attrs']][['gene_ids-Sanger-Cells']][1:100]
Data_Litvinukova2020_loom[['row_attrs']][['var_names']][1:1000]

# Loading some relevant data (note: the full dataset contains 
# 486134 cells x 33538 genes)
umap1=Data_Litvinukova2020_loom[['col_attrs']][['X_umap']][1,1:100000]
umap2=Data_Litvinukova2020_loom[['col_attrs']][['X_umap']][2,1:100000]
celltype=Data_Litvinukova2020_loom[['col_attrs']][['cell_type']][1:100000]

# Make a plot of the data
ggplot(data.frame(umap1=umap1, umap2=umap2, celltype=celltype))+
    geom_point(aes(x=umap1,y=umap2, color=celltype),size=.1)+theme_bw()+ 
    guides(colour = guide_legend(override.aes = list(size=4)))+
    scale_color_manual(values = col_vector_60)

# Testing the apply function
Data_Litvinukova2020_loom$apply(FUN = median, MARGIN = 1, display.progress = T, 
                                name = 'row_attrs/test', dataset.use = "matrix", chunk.size = 500) # index.use

# Select  celltypes of interest
all_celltypes= celltype=Data_Litvinukova2020_loom[['col_attrs']][['cell_type']][]
subset(x=Data_Litvinukova2020_loom, filename='/Volumes/workdrive_m.wehrens_hubrecht/data/2020_12_Litvinukova2020/loom/Data_Litvinukova2020_loom_sel_A_CM.loom', m=which(all_celltypes=='Atrial_Cardiomyocyte'), )
    # this unfortunately results in this error:
    # https://github.com/mojaveazure/loomR/issues/54
    #
    # So I think I should download (and convert) the datasets separately to get "typical cell type gene expression patterns"



