

# base_dir = '/Users/m.wehrens/Data/_2020_03_Wang/'
base_dir = '/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/'
# base_dir = '/hpc/hub_oudenaarden/mwehrens/data/Wang/'

library(plyr)

# The Wang data is split into two batches



# First batch of cells, patients N1-N12
# Conveniently, there's also RSA numbers in Wang_MetaData_GSE109816
# (though not in other file)
Wang_MetaData_GSE109816 = read.csv(paste0(base_dir,'Metadata/GSE109816_metadata_barcodes_9994cells.txt'), sep='\t')
colnames(Wang_MetaData_GSE109816)[1] = 'ID' # rename rows such that they're consisent with the other tables
colnames(Wang_MetaData_GSE109816)[colnames(Wang_MetaData_GSE109816)=='characteristics..individual'] = 'Individual'
dim(Wang_MetaData_GSE109816) # 9994   40
length(unique(Wang_MetaData_GSE109816$ID)) # checking whether IDs are unique (l=9994, yes)
sort(unique(Wang_MetaData_GSE109816$characteristics..individual)) # "N1"  "N10" "N11" "N12" "N2"  "N3"  "N4"  "N5"  "N6"  "N7"  "N8"  "N9" 
colnames(Wang_MetaData_GSE109816)


# Second batch of cells, HF and N13 and N14
Wang_Info_GSE121893  = read.csv(paste0(base_dir,'Metadata/MW_GSE121893_human_heart_sc_info.txt'), sep='\t', comment.char = '#') 
dim(Wang_Info_GSE121893) # 4933   35
length(unique(Wang_Info_GSE121893$ID)) # checking whether IDs are unique (l=4933, yes)
sort(unique(Wang_Info_GSE121893$Individual)) # "C1"  "C2"  "D1"  "D2"  "D4"  "D5"  "N13" "N14"


# Let's merge the two tables into one metadata table, which I'll call combined_metadata_Wang
###
# Are the colnames the same for Wang_Info_GSE121893 and Wang_MetaData_GSE109816?
names(Wang_Info_GSE121893)
names(Wang_MetaData_GSE109816)
names(Wang_Info_GSE121893) %in% names(Wang_MetaData_GSE109816) # mostly yes
# so we can merge these two meta-data dataframes.
combined_metadata_Wang = rbind.fill(Wang_Info_GSE121893, Wang_MetaData_GSE109816)
# merge(Wang_Info_GSE121893,Wang_MetaData_GSE109816,by.x=0,by.y=0,all=TRUE)
dim(combined_metadata_Wang) # 14927 44
length(unique(combined_metadata_Wang$ID)) # IDs are still unique, also n=14927

# Information on the clustering that they performed, 
# relating to all cells
# Here's also information about the celltypes that 
# they determined, 
# we'll be looking for LV-CMs + LA-CMs, 
# which after filtering should be 3894 cells
Wang_ClusterInfo_GSE121893 = read.csv(paste0(base_dir,'Metadata/GSE121893_all_heart_cell_cluster_info.txt'), sep='\t')
dim(Wang_ClusterInfo_GSE121893) # 11377     8
length(unique(Wang_ClusterInfo_GSE121893$ID)) # checking whether IDs are unique
sort(unique(Wang_ClusterInfo_GSE121893$sample)) # "C1"  "C2"  "D1"  "D2"  "D4"  "D5"  "N1"  "N10" "N11" "N12" "N13" "N14" "N2"  "N3"  "N4"  "N5"  "N6"  "N7"  "N8"  "N9" 
colnames(Wang_ClusterInfo_GSE121893)
View(Wang_ClusterInfo_GSE121893)

# Let's look at those cell types
sort(unique(Wang_ClusterInfo_GSE121893$ident))
    # see e.g. fig 2e, LV1-LV7 are our samples of interest, those are LV-CMs
    # LV-CMs + LA-CMs should be 3894 cells
sum(grepl('LV|LA',Wang_ClusterInfo_GSE121893$ident))
sort(unique(Wang_ClusterInfo_GSE121893$ident))
sum(grepl('^LV[0-9]$|^LA[0-9]$',Wang_ClusterInfo_GSE121893$ident))
    # these are more, which is a bit weird, given 
    # that presumably ident was determined 
    # after applying selection criteria
# Unfortunately, BC information can't be found in this table, so i'll have to look it up
# in another table
metadata_Wang_full_table = 
    dplyr::full_join(Wang_ClusterInfo_GSE121893, combined_metadata_Wang, by='ID')
    # note that the number of rows don't match, which is probably because 
    # metadata tables are not QC filtered ..
# Checking whether it worked, using SC_92563_0_17 as example
Wang_ClusterInfo_GSE121893[Wang_ClusterInfo_GSE121893$ID=='SC_92563_0_17',]
metadata_Wang_full_table[metadata_Wang_full_table$ID=='SC_92563_0_17',]

# Anyways, let's see if we can get the names of the cells that we want
desired_celltypes = grepl('^LV[0-9]*$',metadata_Wang_full_table$ident)
desired_individuals = grepl('^N',metadata_Wang_full_table$Individual)
desired_selection = desired_celltypes&desired_individuals
sum((grepl('^N',unique(metadata_Wang_full_table$Individual)))) # just checking
sum(desired_celltypes)
sum(desired_celltypes&desired_individuals) # 1893, a bit more than expected 1864, but that might be because of LA samples classified as LV (e.g. N6 is found below), and vice versa
# View(metadata_Wang_full_table[desired_celltypes&desired_individuals,])
desired_cellnames = metadata_Wang_full_table[desired_celltypes&desired_individuals,]$ID
# now compose a name that we can use for selection
sum(is.na(metadata_Wang_full_table$Barcode)) # just checking
desired_cells_mwName = paste0(metadata_Wang_full_table[desired_selection,]$Individual, '-', metadata_Wang_full_table[desired_selection,]$Barcode)


save(list = c('desired_cells_mwName'),file = paste0(base_dir,'Rdata/desired_cells_mwName.Rdata'))

# Done with stuff that's also needed @HPC
########################################################################
########################################################################
########################################################################


########################################################################
# double checking barcodes

# assmeble barcodes that were used
unique(Wang_Info_GSE121893$Barcode)
unique(Wang_MetaData_GSE109816$Barcode)
unique(c(Wang_Info_GSE121893$Barcode,Wang_MetaData_GSE109816$Barcode))
barcodes_used_wang = unique(c(Wang_Info_GSE121893$Barcode,Wang_MetaData_GSE109816$Barcode))

# let's cross-ref with their list of barcodes
takara_barcodes = read.table('/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/takarabio_protocol/config/barcode_map_tcr_v1.csv', sep=',', header=1)
thetakdir='/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/takarabio_protocol/config/'
thetakfile='barcode_map_preprinted_v2.csv'
takara_barcodes = read.table(paste0(thetakdir,thetakfile), sep=',', header=1)
View(takara_barcodes)

# check whether all barcodes used are in fact in my ref file
all(barcodes_used_wang %in% takara_barcodes$Barcode)
    # yes they are ..
    # also, length(barcodes_used_wang) = 4877
    # length(takara_barcodes$Barcode) = 5184
# Let's create a barcode file from this
df_barcodes_takara = data.frame(BC=takara_barcodes$Barcode, BC_name = paste0('R',takara_barcodes$Row,'.C',takara_barcodes$Column,'.',takara_barcodes$Barcode))

write.table(file='/Users/m.wehrens/Documents/git_repos/mapping_aa_private_vasa/customized_scripts/bc_takara_map_tcr_v1.tsv',x = df_barcodes_takara,sep='\t',quote = F, col.names = F, row.names = F)

rm('Wang_Data1')
Wang_Data1 = read.table(paste0(base_dir,'Original_counttables/GSE109816_normal_heart_umi_matrix.csv'),row.names = 1, header=1, sep=',')
View(Wang_Data1[1:100,1:100])

# rm('Wang_Table1')

########################################################################

