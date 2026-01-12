
# 2024-10-24
# Checking consistency between raw data and analyzed (Seurat object) data
#
# See also ~/Documents/git_repos/SCS_More_analyses/Projects/howtousedata_HCMSCS/example-load-plot-data.R
# for convenient Rds files regarding this data.

# Load the analyzed data
# LOCAL=1; script_dir = '/Users/m.wehrens/Documents/git_repos/_Hubrecht-n-before/SCS_More_analyses/Projects/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

DATASET_NAME='ROOIJonly.sp.bt_RID2l'
ANALYSIS_NAME=DATASET_NAME
ANALYSIS_NAME_clExtended = paste0(DATASET_NAME, '_clExtended')

final_filename = paste0('H5_RHL_SeuratObject_nM_sel_',ANALYSIS_NAME_clExtended,'.h5seurat')
if (!exists('current_analysis')) {current_analysis = list()}
current_analysis[[ANALYSIS_NAME_clExtended]] =
    LoadH5Seurat(file = paste0(base_dir,'Rdata/',final_filename))

# Convert to raw data matrix
matrix_HCM_raw_big = as.matrix(current_analysis$ROOIJonly.sp.bt_RID2l_clExtended@assays$RNA@counts)
View(matrix_HCM_raw_big)

# Now load one of the raw data files
# The data files names are stored in the following parameter:
# (defined in the script ~/Documents/git_repos/SCS_More_analyses/Projects/HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R)
HCM_SCS_ids
dataset_list_paths_HCM
dataset_list_paths_HCM_filenameonly = paste0(HCM_SCS_ids, '_pT.nonRibo_E99_Aligned.out.counts.table.tsv')
    # data.frame(HCM_SCS_ids)

# Set main dir with data
main_data_dir_sanitycheck = '/Users/m.wehrens/Data_Hubrecht/_2019_02_HCM_SCS/2021_UmiTools_CountTables_Wang_HCMSCS/exclMulti_Filtered/ROOIJ/'
# main_data_dir_sanitycheck = '/Users/m.wehrens/Data_Hubrecht/_2019_02_HCM_SCS/2024_datasharing/counttables_ROOIJ_per-entry/'

# Loop over tables and perform consistency checks
for (IDX in 1:length(HCM_SCS_ids)) {
    # IDX=3
    
    print(names(HCM_SCS_ids)[IDX])
    file_name = paste0(main_data_dir_sanitycheck, dataset_list_paths_HCM_filenameonly[IDX])
    data = read.table(file_name, header=TRUE, row.names=1, sep=' ')
    
    # make colnames consistent with seurat object data
    colnames(data) = paste0(names(HCM_SCS_ids)[IDX], '_', colnames(data))
    # make rownames consistent 
    # quick sanity check, count number of underscores per name
    rownames_data = rownames(data)
    all(sapply(rownames_data, str_count, '_')==1)
    # now replace the underscore by a :
    rownames(data) = gsub('_', ':', rownames_data)
    
    # now for comparison, select the same cols and rows from both matrices
    # rownames
    rownames_rawfile = rownames(data)
    rownames_seurat = rownames(matrix_HCM_raw_big)
    # colnames
    colnames_rawfile = colnames(data)
    colnames_seurat = colnames(matrix_HCM_raw_big)
    # intersect col and rownames
    rownames_shared = intersect(rownames_rawfile, rownames_seurat)
    colnames_shared = intersect(colnames_rawfile, colnames_seurat)
    
    # now create the intersect subset
    data_shared = data[rownames_shared, colnames_shared]
    matrix_HCM_raw_big_shared = matrix_HCM_raw_big[rownames_shared, colnames_shared]
    
    # now check whether they are the same
    rowcount_shared = dim(data_shared)[1]
    rowcount_seurat = dim(matrix_HCM_raw_big)[1]
    all(data_shared == matrix_HCM_raw_big_shared)
    if (all(data_shared == matrix_HCM_raw_big_shared)) {
        print(paste0('check passed, nr of genes:', rowcount_shared,'/',rowcount_seurat)) }
    else {
        print('check failed') }

}

colSums(data_shared[,1:10])
colSums(matrix_HCM_raw_big_shared[,1:10])

# Now also check another location with raw data, whether that one is the same
main_data_dir_sanitycheck_loc2 = '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_UmiTools_CountTables_Wang_HCMSCS/exclMulti_Filtered/ROOIJ/'
file_name_loc2 = paste0(main_data_dir_sanitycheck_loc2, dataset_list_paths_HCM_filenameonly[IDX])
data_loc2 = read.table(file_name_loc2, header=TRUE, row.names=1, sep=' ')
# make colnames consistent with seurat object data
colnames(data_loc2) = paste0(names(HCM_SCS_ids)[IDX], '_', colnames(data_loc2))
# make rownames consistent
rownames_data_loc2 = rownames(data_loc2)
all(sapply(rownames_data_loc2, str_count, '_')==1)
rownames(data_loc2) = gsub('_', ':', rownames_data_loc2)
# now select the shared data
data_loc2_shared = data_loc2[rownames_shared, colnames_shared]

all(data_loc2_shared == matrix_HCM_raw_big_shared)
dim(data_loc2_shared)
dim(matrix_HCM_raw_big_shared)

# This leads to NA values, as not all data is present at loc2
any(!(rownames(data_loc2_shared) %in% rownames_shared))
any(is.na(!(rownames(data_loc2_shared) %in% rownames_shared)))
any(rownames(data_loc2_shared)=='NA')
rownames(data_loc2_shared)[!(rownames(data_loc2_shared) %in% rownames_shared)]
rownames_shared[!(rownames_shared %in% rownames(data_loc2_shared))]


#####################
# Checking consistency seurat objects that were saved at different times
suerat_Rd1_may = readRDS( '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3b/Rdata/H5_RHL_SeuratObject_nM_sel_ROOIJonly.sp.bt_RID2l_clExtended.Rds' )
seurat_Rd2_jan = readRDS( '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3b/Rdata/2024/H5_RHL_SeuratObject_nM_sel_ROOIJonly.sp.bt_RID2l_clExtended.Rds' )

all(suerat_Rd1_may@assays$RNA@scale.data == seurat_Rd2_jan@assays$RNA@scale.data)
all(current_analysis$ROOIJonly.sp.bt_RID2l_clExtended@assays$RNA@scale.data == seurat_Rd2_jan@assays$RNA@scale.data) 

# all correct



#######################
# Addition 2025-02

# Conversion of joep ID to plate IDs
Rooij_conversion_pts.plt = c('JE01'='R.P1.plt1', 'JE02'='R.P1.plt2', 'JE03'='R.P1.plt3', 'JE04'='R.P1.plt4',
                        'JE05'='R.P2.plt1', 'JE06'='R.P2.plt2', 'JE07'='R.P2.plt3', 'JE08'='R.P2.plt4',
                        'MW05'='R.P3.plt1', 'MW06'='R.P3.plt2', 'MW07'='R.P3.plt3', 'MW08'='R.P3.plt4',
                        'JE10'='R.P4.plt1', 'JE11'='R.P4.plt2',
                        'AL01'='R.P5.plt1', 'AL02'='R.P5.plt2')
current_analysis$ROOIJonly.sp.bt_RID2l_clExtended$patient.plate_fct =
    factor(Rooij_conversion_pts.plt[current_analysis$ROOIJonly.sp.bt_RID2l_clExtended$orig.ident], 
           levels=Rooij_conversion_pts.plt)

# Another sanity check to confirm above mapping
# md5sum was used to confirm.

# the original Joep files with plate X annotation (downloaded from GEO) are in:
# /Users/m.wehrens/Data_Hubrecht/_2019_02_HCM_SCS/_countdata_downloadedGEO-2025

# the original Joep files with IDs are in:
# /Users/m.wehrens/Data_Hubrecht/_2019_02_HCM_SCS/_countdata






