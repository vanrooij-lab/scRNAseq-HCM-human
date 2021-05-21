

# Conveniently, there's also RSA numbers here
# So, GSE109816 is not so complicated
Wang_MetaData_GSE109816 = read.csv('/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/GSE109816_metadata_barcodes_9994cells.txt', sep='\t')
View(Wang_MetaData_GSE109816)

unique(Wang_MetaData_GSE109816$Barcode)
length(unique(Wang_MetaData_GSE109816$SRR.run.accession))
dim(Wang_MetaData_GSE109816)

Wang_ClusterInfo_GSE121893 = read.csv('/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/GSE121893_all_heart_cell_cluster_info.txt', sep='\t')
View(Wang_ClusterInfo_GSE121893)

Wang_Info_GSE121893  = read.csv('/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/MW_GSE121893_human_heart_sc_info.txt', sep='\t', comment.char = '#') 
View(Wang_Info_GSE121893)

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
Wang_Data1 = read.table('/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/GSE109816_normal_heart_umi_matrix.csv',row.names = 1, header=1, sep=',')
View(Wang_Data1[1:100,1:100])

# rm('Wang_Table1')



