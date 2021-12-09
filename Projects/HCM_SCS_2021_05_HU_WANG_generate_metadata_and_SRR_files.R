

######################################################################

# - Wang et al. dataset
#       Wang, L., Yu, P., Zhou, B., Song, J., Li, Z., Zhang, M., … Hu, S. (2020). Single-cell reconstruction of the adult human heart during heart failure and recovery reveals the cellular landscape underlying cardiac function. Nature Cell Biology, 22(1), 108–119. https://doi.org/10.1038/s41556-019-0446-7
#       Metadata was obtained via 
#       https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121893 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109816 
#       Raw data was downloaded using fasterq-dump and mapped with custom pipeline to make it more comparable to our own data

######################################################################

# base_dir = '/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/'
# base_dir = '/hpc/hub_oudenaarden/mwehrens/data/Wang/'

library(plyr)

######################################################################
# Collection of important Wang metadata into one file

# The Wang data is split into two batches

# First batch of cells, patients N1-N12
# Conveniently, there's also RSA numbers in Wang_MetaData_GSE109816
# (though not in other file)
Wang_MetaData_GSE109816 = read.csv(paste0(base_dir,'Wang_metadata/GSE109816_metadata_barcodes_9994cells.txt'), sep='\t')
colnames(Wang_MetaData_GSE109816)[1] = 'ID' # rename rows such that they're consisent with the other tables
colnames(Wang_MetaData_GSE109816)[colnames(Wang_MetaData_GSE109816)=='characteristics..individual'] = 'Individual'
dim(Wang_MetaData_GSE109816) # 9994   40
length(unique(Wang_MetaData_GSE109816$ID)) # checking whether IDs are unique (l=9994, yes)
sort(unique(Wang_MetaData_GSE109816$characteristics..individual)) # "N1"  "N10" "N11" "N12" "N2"  "N3"  "N4"  "N5"  "N6"  "N7"  "N8"  "N9" 
colnames(Wang_MetaData_GSE109816)

# Second batch of cells, HF and N13 and N14
Wang_Info_GSE121893  = read.csv(paste0(base_dir,'Wang_metadata/MW_GSE121893_human_heart_sc_info.txt'), sep='\t', comment.char = '#') 
dim(Wang_Info_GSE121893) # 4933   35
length(unique(Wang_Info_GSE121893$ID)) # checking whether IDs are unique (l=4933, yes)
sort(unique(Wang_Info_GSE121893$Individual)) # "C1"  "C2"  "D1"  "D2"  "D4"  "D5"  "N13" "N14"

# so for the Wang_MetaData_GSE109816 part there's already SRR info, but now we'd also
# like to get the SRR info for the Wang_Info_GSE121893 part.
# This can be found in the file (couples SRR to {plate,x,y}
GSE121893_files = read.table(paste0(base_dir, 'Wang_metadata/GSE121893_files.txt'))
# Let's modify this data
GSE121893_SRR_ = matrix(data=GSE121893_files$V1[!GSE121893_files$V1=='END'], ncol =3, byrow = T)
GSE121893_SRR = data.frame(SRR=gsub(pattern = ':',replacement = '',x = GSE121893_SRR_[,1]), 
                            ID=gsub(pattern = '_r1.fq.gz',replacement = '',x = GSE121893_SRR_[,2]),
                            row.names = gsub(pattern = '_r1.fq.gz',replacement = '',x = GSE121893_SRR_[,2]))
# So now we can add this to the previous table
Wang_Info_GSE121893$SRR.run.accession = GSE121893_SRR[Wang_Info_GSE121893$ID,]$SRR


# Let's merge the two tables into one metadata table, which I'll call combined_metadata_Wang
###
# Are the colnames the same for Wang_Info_GSE121893 and Wang_MetaData_GSE109816?
names(Wang_Info_GSE121893)
names(Wang_MetaData_GSE109816)
names(Wang_Info_GSE121893) %in% names(Wang_MetaData_GSE109816) # mostly yes
# so we can merge these two meta-data dataframes.
# first annotate origin
Wang_Info_GSE121893$GSE='GSE121893'
Wang_MetaData_GSE109816$GSE='GSE109816'
# then merge
combined_metadata_Wang = rbind.fill(Wang_Info_GSE121893, Wang_MetaData_GSE109816)
# merge(Wang_Info_GSE121893,Wang_MetaData_GSE109816,by.x=0,by.y=0,all=TRUE)
dim(combined_metadata_Wang) # 14927 44
length(unique(combined_metadata_Wang$ID)) # IDs are still unique, also n=14927





# Information on the clustering that they performed, 
# !!**relating to all cells**!!
# Here's also information about the celltypes that 
# they determined, 
# we'll be looking for LV-CMs + LA-CMs, 
# which after filtering should be 3894 cells
Wang_ClusterInfo_GSE121893 = read.csv(paste0(base_dir,'Wang_metadata/GSE121893_all_heart_cell_cluster_info.txt'), sep='\t')
dim(Wang_ClusterInfo_GSE121893) # 11377     8
length(unique(Wang_ClusterInfo_GSE121893$ID)) # checking whether IDs are unique
sort(unique(Wang_ClusterInfo_GSE121893$sample)) # "C1"  "C2"  "D1"  "D2"  "D4"  "D5"  "N1"  "N10" "N11" "N12" "N13" "N14" "N2"  "N3"  "N4"  "N5"  "N6"  "N7"  "N8"  "N9" 
colnames(Wang_ClusterInfo_GSE121893)
# View(Wang_ClusterInfo_GSE121893)

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
# still add a field that identifies the [plate/]run that they used
metadata_Wang_full_table$plate_nr = sapply(str_split(string = metadata_Wang_full_table$ID, pattern = '_'), function(x) {x[[2]]})
# now also add an unique ID (also as a sanity check)
metadata_Wang_full_table$ID_MW = paste0(metadata_Wang_full_table$sample, '-', metadata_Wang_full_table$plate_nr, '-', metadata_Wang_full_table$Barcode)
# change rownames
rownames(metadata_Wang_full_table) = metadata_Wang_full_table$ID

# Anyways, let's see if we can get the names of the cells that we want
desired_celltypes = grepl('^LV[0-9]*$',metadata_Wang_full_table$ident)
desired_individuals = grepl('^N',metadata_Wang_full_table$Individual)
desired_selection = desired_celltypes&desired_individuals
sum((grepl('^N',unique(metadata_Wang_full_table$Individual)))) # just checking
sum(desired_celltypes)
sum(desired_celltypes&desired_individuals) # 1893, a bit more than expected 1864, but that might be because of LA samples classified as LV (e.g. N6 is found below), and vice versa
# View(metadata_Wang_full_table[desired_celltypes&desired_individuals,])

# So now we know which cells we're interested in, let's create files
# that allow us to gather their raw data
metadata_Wang_full_table_selection = metadata_Wang_full_table[desired_celltypes&desired_individuals,]


######################################################################

# let's create the SRR files
# these will be per patient, per plate now
# the bar codes will identify the cells
all_samples_list = c()
SRR_export_summary_info=list()
current_SRR_of_interest_list=list()
for (current_patient in c('N1','N2','N3','N4','N5','N13','N14')) {
    
    # current_patient = 'N2'
    
    # obtain plates that belong to this patient
    plates_for_this_patient = unique(metadata_Wang_full_table_selection$plate_nr[metadata_Wang_full_table_selection$sample==current_patient])
    
    # now gather SRR numbers that belong to this patient, to this plate
    for (current_plate in plates_for_this_patient) {
        
        print(paste0('=======','PATIENT ',current_patient, ', PLATE ',current_plate,'======='))
        
        # current_plate = '97452'
        current_SRR_of_interest = metadata_Wang_full_table_selection[metadata_Wang_full_table_selection$sample==current_patient&
                                                                     metadata_Wang_full_table_selection$plate_nr==current_plate,]$SRR.run.accession
        # Just another sanity check
        current_SRR_of_interest_list[[length(current_SRR_of_interest_list)+1]] = current_SRR_of_interest
        
        # Also collect current barcodes, this is just to do another sanity check 
        # whether these are unique
        current_BCs_of_interest = metadata_Wang_full_table_selection[metadata_Wang_full_table_selection$sample==current_patient&
                                                                     metadata_Wang_full_table_selection$plate_nr==current_plate,]$Barcode
        if (length(current_BCs_of_interest) == length(unique(current_BCs_of_interest))) {
            print('Barcodes are unnique.')    
        } else { warning ('Barcodes are NOT unique'); print('Barcodes are NOT unique') }

        
        # now put those in two separate files
        # (This is just to speed things up)
        halfway=ceiling(length(current_SRR_of_interest)/2)
        current_SRR_Part1 = current_SRR_of_interest[1:halfway]
        current_SRR_Part2 = current_SRR_of_interest[(halfway+1):length(current_SRR_of_interest)]
        
        # create a name for this "sample"
        sample_name_part1 = paste0('p.',current_patient,'.plate.',current_plate, '.part.1')
        sample_name_part2 = paste0('p.',current_patient,'.plate.',current_plate, '.part.2')
        
        # Create directories 
        current_directory_part1 = paste0(base_dir, 'SRR/', sample_name_part1, '/')
        current_directory_part2 = paste0(base_dir, 'SRR/', sample_name_part2, '/')
        dir.create(current_directory_part1, recursive = T)
        dir.create(current_directory_part2, recursive = T)
        # Write tables
        if (!exists('NOWRITESRR')){ # set NOWRITESRR to not write these files but still complete loop
            write.table(x = current_SRR_Part1, file = paste0(current_directory_part1,sample_name_part1,'_SRR_Acc_List.txt'), quote = F, row.names = F, col.names = F)
            write.table(x = current_SRR_Part2, file = paste0(current_directory_part2,sample_name_part2,'_SRR_Acc_List.txt'), quote = F, row.names = F, col.names = F)
        }
        
        # User info
        print(paste0('Exporting ',sample_name_part1, ' with ', length(current_SRR_Part1),' cells SRR ..'))
        print(paste0('Exporting ',sample_name_part2, ' with ', length(current_SRR_Part2),' cells SRR ..'))
        SRR_export_summary_info[[sample_name_part1]] = length(current_SRR_Part1)
        SRR_export_summary_info[[sample_name_part2]] = length(current_SRR_Part2)
        
        # keep track of samples
        all_samples_list = append(all_samples_list, c(sample_name_part1, sample_name_part2))
        
    }
        
}
# export sample list
write.table(data.frame(as.list(all_samples_list)), file = paste0(base_dir, 'SRR/','sample_list.txt'), quote = F, row.names = F, col.names = F)

# Another sanity check
print('Are the collected SRRs now unique?')
length(unique(unlist(current_SRR_of_interest_list)))==length(unlist(current_SRR_of_interest_list))

# length(SRR_N3_97438)
# length(unique(SRR_N3_97438))

# notes:
# [1] "Exporting p.N1.plate.97474.part.1 with 130 cells SRR .."
# [1] "Exporting p.N1.plate.97474.part.2 with 129 cells SRR .."
# [1] "Exporting p.N2.plate.97452.part.1 with 128 cells SRR .."
# [1] "Exporting p.N2.plate.97452.part.2 with 128 cells SRR .."
# [1] "Exporting p.N2.plate.97493.part.1 with 21 cells SRR .."
# [1] "Exporting p.N2.plate.97493.part.2 with 20 cells SRR .."
# [1] "Exporting p.N3.plate.97438.part.1 with 180 cells SRR .."
# [1] "Exporting p.N3.plate.97438.part.2 with 180 cells SRR .."
# [1] "Exporting p.N4.plate.97461.part.1 with 34 cells SRR .."
# [1] "Exporting p.N4.plate.97461.part.2 with 34 cells SRR .."
# [1] "Exporting p.N5.plate.97458.part.1 with 169 cells SRR .."
# [1] "Exporting p.N5.plate.97458.part.2 with 169 cells SRR .."
# [1] "Exporting p.N13.plate.100355.part.1 with 37 cells SRR .."
# [1] "Exporting p.N13.plate.100355.part.2 with 36 cells SRR .."
# [1] "Exporting p.N14.plate.104720.part.1 with 14 cells SRR .."
# [1] "Exporting p.N14.plate.104720.part.2 with 14 cells SRR .."

################################################################################
# I've noted that the numbers don't quite add up
# This is due to some striking features of the Wang data

# Problem 1: the table with selected cells contains 1893 cells; but only
# 1423 SRR values are exported, why are there cells missing?
# --> 
sum(unlist(SRR_export_summary_info)) # 1423 SRR values
# Turns out that those "missing cells" belong to patients from who no LV was collected.
dim(metadata_Wang_full_table_selection[metadata_Wang_full_table_selection$sample %in% c('N6', 'N7', 'N8', 'N9', 'N10', 'N11', 'N12'),])
dim(metadata_Wang_full_table_selection[grepl('_LA',metadata_Wang_full_table_selection$source.name) & metadata_Wang_full_table_selection$sample %in% c('N6', 'N7', 'N8', 'N9', 'N10', 'N11', 'N12'),])
    # So, there's 470 cells that are classified as LV, whilst they are from LA source

# Also double check this with "original" tables
# Can be seen from "cluster info" table directly 
dim(Wang_ClusterInfo_GSE121893[grepl('_LA',Wang_ClusterInfo_GSE121893$condition) & Wang_ClusterInfo_GSE121893$sample %in% c('N6', 'N7', 'N8', 'N9', 'N10', 'N11', 'N12') &
                                grepl('^LV',Wang_ClusterInfo_GSE121893$ident),])

# Let's summarize this table
Wang_ClusterInfo_GSE121893_selectedFields = Wang_ClusterInfo_GSE121893[,c('condition','group','sample')]
Wang_ClusterInfo_GSE121893_selectedFields$identClass = gsub('[0-9]*$','',Wang_ClusterInfo_GSE121893$ident)
Wang_ClusterInfo_GSE121893_selectedFields_summaryTable_ = as.data.frame(table(Wang_ClusterInfo_GSE121893_selectedFields))
Wang_ClusterInfo_GSE121893_selectedFields_summaryTable = Wang_ClusterInfo_GSE121893_selectedFields_summaryTable_[Wang_ClusterInfo_GSE121893_selectedFields_summaryTable_$Freq>0,]
Wang_ClusterInfo_GSE121893_selectedFields_summaryTable_customSel = Wang_ClusterInfo_GSE121893_selectedFields_summaryTable[grepl('^N_',Wang_ClusterInfo_GSE121893_selectedFields_summaryTable$condition),]
Wang_ClusterInfo_GSE121893_selectedFields_summaryTable_customSel2 = Wang_ClusterInfo_GSE121893_selectedFields_summaryTable[grepl('^N_',Wang_ClusterInfo_GSE121893_selectedFields_summaryTable$condition) & 
                                                                                                                         Wang_ClusterInfo_GSE121893_selectedFields_summaryTable$identClass=='LV',]

# custom order
mapping_temp = paste0('N',c(paste0('0',1:9),10:20))
names(mapping_temp) = paste0('N',1:20)
Wang_ClusterInfo_GSE121893_selectedFields_summaryTable_customSel2_Sorted = Wang_ClusterInfo_GSE121893_selectedFields_summaryTable_customSel2[order(mapping_temp[as.character(Wang_ClusterInfo_GSE121893_selectedFields_summaryTable_customSel2$sample)]),]
openxlsx::write.xlsx(x = as.data.frame(Wang_ClusterInfo_GSE121893_selectedFields_summaryTable_customSel2_Sorted), 
                     file = paste0(base_dir,'Rplots/mapped_Wang_cells_overview_stats2_healthy-IdentLV.xlsx'), overwrite=T)


# Problem 2: there are 1864 cells listed as CMs in the excel table, coming from LV
# source. But now I only have 1423 cells. Are these then classified as LA?
sum(grepl('^LA',metadata_Wang_full_table[metadata_Wang_full_table$sample %in% c('N1','N2','N3','N4','N5','N13','N14'),]$ident))
    # indeed, 500 of those "LV" cells end up being classified as LA cells
sum(!grepl('^LV',metadata_Wang_full_table[metadata_Wang_full_table$sample %in% c('N1','N2','N3','N4','N5','N13','N14'),]$ident))
unique(metadata_Wang_full_table[metadata_Wang_full_table$sample %in% c('N1','N2','N3','N4','N5','N13','N14'),]$ident)
dim(metadata_Wang_full_table[metadata_Wang_full_table$sample %in% c('N1','N2','N3','N4','N5','N13','N14') &
                             metadata_Wang_full_table$ident %in% c("LV2","LV4","LV7","LV5","LV6","LA3","LA4","LA2","LV1","LA1","LA6","LA5","LV3"),])


################################################################################

# Additionally, I'd like to make a mapping that maps the original names to 
# my new names, so 
# GSE109816_SC_92563_2_44 <--> N3_97438_a_R15.C69.ATTGGAGATTG
# This should be possible with the table

# New way
extended_names =  paste0(  metadata_Wang_full_table_selection$sample, '.',
                            metadata_Wang_full_table_selection$plate_nr, '_',
                            'R',metadata_Wang_full_table_selection$Chip.Row.ID,'.C',metadata_Wang_full_table_selection$Chip.Column.ID,'.',
                            metadata_Wang_full_table_selection$Barcode)
rownames(metadata_Wang_full_table_selection) = extended_names
# metadata_Wang_full_table_selection[cellname]$ID now gives the old name (or other info)
save(list = c('metadata_Wang_full_table_selection'),file = paste0(base_dir,'Rdata/metadata_Wang_full_table_selection.Rdata'))

# Also save the full table
save(list = c('metadata_Wang_full_table'),file = paste0(base_dir,'Rdata/metadata_Wang_full_table.Rdata'))
# load(paste0(base_dir,'Rdata/metadata_Wang_full_table.Rdata')) # metadata_Wang_full_table


################################################################################

# Old alternative for name conversion

simplified_names = paste0(  metadata_Wang_full_table_selection$sample, '_',
                            metadata_Wang_full_table_selection$plate_nr, '_',
                            metadata_Wang_full_table_selection$Barcode)
conversion_table_wang_names = metadata_Wang_full_table_selection$ID
names(conversion_table_wang_names) = simplified_names
wang_convert_name_new_to_old = function(newname, add_GSE=F) {
    newname_short =
        gsub(pattern = '_[ab]_R[0-9]+.C[0-9]+.', replacement='_', x = newname)
    if (add_GSE) { 
        if (sum(metadata_Wang_full_table_selection$ID==conversion_table_wang_names[newname_short])!=1) {stop(paste0('irregularity w ',newname))}
        GSE=metadata_Wang_full_table_selection$GSE[metadata_Wang_full_table_selection$ID==conversion_table_wang_names[newname_short]]
        return(paste0(GSE, '_', conversion_table_wang_names[newname_short])) 
    } else { 
        return(conversion_table_wang_names[newname_short]) 
    }
}

################################################################################
# I previously made a mistake, ie not realizing they did run multiple plates,
# which led to inadvertedly joining files of multiple cells.
# The code below still has that mistake.

# desired_cellnames = metadata_Wang_full_table[desired_celltypes&desired_individuals,]$ID
# # now compose a name that we can use for selection
# sum(is.na(metadata_Wang_full_table$Barcode)) # just checking
# desired_cells_mwName = paste0(metadata_Wang_full_table[desired_selection,]$Individual, '-', metadata_Wang_full_table[desired_selection,]$Barcode)
# 
#     # This is assuming that for each patient, only one plate was sequenced, such that patient-barcode combinations are unique
#     # IS THIS THE CASE ??
#     # --> how can we check this --> use $ID, and per ID check number of unique patient-BC combinations there are
#     # metadata_Wang_full_table$MW_ID = paste0(metadata_Wang_full_table$Individual,'-',metadata_Wang_full_table$Barcode)
#     # table_MW_ID = table(metadata_Wang_full_table$MW_ID)
#     # table_MW_ID[order(table_MW_ID, decreasing=T)][1:100]
#     # N2-CCGTCATTGCG
#     # length(table_MW_ID[table_MW_ID>1])
#     # [1] 417
#     # So there's ±417 cells that are double, whilst they're merged in my data ..
#     #
#     # So the question has become: does the first number in the ID part make the cells unique? 
#     # (and does it correspond to a plate?)
#     
# 
# metadata_Wang_full_table$IDpart1=sapply(str_split(metadata_Wang_full_table$ID, pattern = '_'), function(x) {x[2]})
# metadata_Wang_full_table$MW_ID2 = paste0(metadata_Wang_full_table$Individual,'-',metadata_Wang_full_table$Barcode,'-',metadata_Wang_full_table$IDpart1)
# 
# library(dplyr)
# Wang_overview_cellcounts = dplyr::count(metadata_Wang_full_table, Individual, IDpart1)
# rename(count(df, y, m), Freq = n)
# 
# # Answer:
# # yes the number makes the cell IDs unique, 
# # and SRR correspondence for the second data sets can be found in GSE121893_files.txt
# 
# save(list = c('desired_cells_mwName'),file = paste0(base_dir,'Rdata/desired_cells_mwName.Rdata'))
# save(list = c('metadata_Wang_full_table'),file = paste0(base_dir,'Rdata/metadata_Wang_full_table.Rdata'))



######################################################################
# We're done; some random other stuff:

########################################################################
# Double checking whether the barcodes I used for mapping are indeed
# the Takara Bio barcodes I think they've used.

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

