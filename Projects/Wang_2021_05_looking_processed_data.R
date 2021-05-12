

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

rm('Wang_Data1')
Wang_Data1 = read.csv('/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/GSE109816_normal_heart_umi_matrix.csv')
View(Wang_Data1[1:100,1:100])

# rm('Wang_Table1')