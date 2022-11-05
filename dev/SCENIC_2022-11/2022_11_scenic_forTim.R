########################################################################
# SCENIC for Tim's dataset


# Note the following websites are convenient:
# https://scenic.aertslab.org/
# Running SCENIC in R: http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html
# Python, command line interface (CLI): https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_SCENIC-protocol-CLI.ipynb
    # Note that they first do all kinds of non-SCENIC related stuff

########################################################################
# Libs

library(loomR)

########################################################################
# Forget about this

# For convenience, on HPC, this script can be sourced by:
# script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')
#
# Local:
# LOCAL=1; script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Projects/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

########################################################################
# Set directory and load data
 
base_dir = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2022_Data_Tim/'

# Data from tim (122 samples)
dataTim_cellInfo = readRDS('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2022_Data_Tim/cellInfo.Rds') # meta data
dataTim_colVars = readRDS('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2022_Data_Tim/colVars.Rds') # colors
dataTim_exprMat = readRDS('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2022_Data_Tim/exprMat.Rds') # expr matrix
    # dim(dataTim_exprMat)

########################################################################
# Export the loom file

# Updating cellInfo
dataTim_cellInfo_matched = dataTim_cellInfo[colnames(dataTim_exprMat), ]
    
# Creating a loom file
lfile <- loomR::create(filename = paste0(base_dir, 'Ldata/', 'dataTim_exprMat', '.loom'), 
                       data = dataTim_exprMat, overwrite = T, 
                       cell.attrs = dataTim_cellInfo_matched) 
                       #feature.attrs = list(gene_info=), 
                       #cell.attrs = list(cell=colnames(current_analysis$R.P1RID2l@assays$RNA@data)))

########################################################################
# Now continue with the shell script..


# (..)






