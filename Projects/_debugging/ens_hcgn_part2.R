
# Regarding

# Result what I did originally
# hcgn were original row annotation, sum(hcgn_names3=="") == 0
# hcgn_names3  = sapply(strsplit(rownames(RHL_SeuratObject_list$TEICH), split = ':'), function(x){if(length(x)>1){x[[2]]}else{""}})
# ens_names3  = sapply(strsplit(rownames(RHL_SeuratObject_list$TEICH), split = ':'), function(x){x[[1]]})

# Load list of raw objects
load(paste0(base_dir,'Rdata/RHL_SeuratObject_list.Rdata'))

# Load some other info (code redundant with main script)
RHL_SeuratObject_list$TEICH <- LoadH5Seurat('/hpc/hub_oudenaarden/mwehrens/data/Teichmann/hca_heart_ventricular_CM_raw.h5seurat')
        # Alt. way
        # Teichmann <- LoadH5Seurat('/hpc/hub_oudenaarden/mwehrens/data/Teichmann/hca_heart_ventricular_CM_raw.h5seurat')
        # RHL_SeuratObject_list$TEICH = Teichmann
        # rm('Teichmann')

# Retreive ensembl names from the object
ensembl_IDs = RHL_SeuratObject_list$TEICH@assays$RNA@meta.features$'gene_ids-Harvard-Nuclei'
names(ensembl_IDs) = rownames(RHL_SeuratObject_list$TEICH@assays$RNA@meta.features)
    # Note that all(RHL_SeuratObject_list$TEICH@assays$RNA@meta.features$'gene_ids-Harvard-Nuclei' == RHL_SeuratObject_list$TEICH@assays$RNA@meta.features$'gene_ids-Sanger-Nuclei')

ensembl_IDs_check = ensembl_IDs[rownames(RHL_SeuratObject_list$TEICH)]

# How I converted the names: 
TEICH_original <- LoadH5Seurat('/hpc/hub_oudenaarden/mwehrens/data/Teichmann/hca_heart_ventricular_CM_raw.h5seurat')
genes_Teich_original = rownames(TEICH_original)
ens_redoOriginalConv = ensembl_IDs[genes_Teich_original]
#all(ens_names3==ens_redoOriginalConv)
    # length(ensembl_IDs)==length(unique(ensembl_IDs))
    # length(names(ensembl_IDs))==length(unique(names(ensembl_IDs)))
    # Note: this conversion is unique

# Regarding our data
genes_Rsample = rownames(RHL_SeuratObject_list$JE10)
genes_Rsample_hcgn = sapply(strsplit(genes_Rsample, split = ':'), function(x){if(length(x)>1){x[[2]]}else{""}})
genes_Rsample_ens = sapply(strsplit(genes_Rsample, split = ':'), function(x){x[[1]]})
length(genes_Rsample_hcgn)
sum(genes_Rsample_hcgn=='')
sum(genes_Rsample_hcgn=='')/length(genes_Rsample_hcgn)
genes_Rsample[genes_Rsample_hcgn=='']

# Is there an issue?
# Are there genes from Rooij and Teichmann with equal ensembl annotation, that don't match in hcgn?

# Sanity check: is our hcgn table unique?
# ===
load(file = paste0(base_dir,'Rdata/ens_to_sym_conv_table__','93','_v2.Rdata'))  # ens_to_sym_conv_table__XX_v2
ens_to_sym_conv_table__XX_v2[1:10]
length(ens_to_sym_conv_table__XX_v2)==length(unique(ens_to_sym_conv_table__XX_v2))
length(ens_to_sym_conv_table__XX_v2)
length(unique(ens_to_sym_conv_table__XX_v2))
View(table(ens_to_sym_conv_table__XX_v2))
# An hcgn names can be linked to multiple ensemble IDs
# ENSG00000129535 ENSG00000285493 
#          "NRL"           "NRL"
# 
# length(names(ens_to_sym_conv_table__XX_v2))==length(names(unique(ens_to_sym_conv_table__XX_v2)))

# Conversely, a ens name can give multiple gene names
View(table(names(ens_to_sym_conv_table__XX_v2)))
# ENSG00000274917 ENSG00000274917 ENSG00000274917 ENSG00000274917 ENSG00000274917 
#    "RNA5-8SN5"     "RNA5-8SN4"     "RNA5-8SN3"     "RNA5-8SN2"     "RNA5-8SN1"
length(names(ens_to_sym_conv_table__XX_v2))==length(unique(names(ens_to_sym_conv_table__XX_v2)))
ens_to_sym_conv_table__XX_v2[names(ens_to_sym_conv_table__XX_v2)=='ENSG00000274917']
#
# Conclusion: it is not
#
# But let's check which names are double:
table_ens_duplicity = table(names(ens_to_sym_conv_table__XX_v2))
ens_with_multiple = names(table_ens_duplicity)[table_ens_duplicity>1]
ens_to_sym_conv_table__XX_v2[names(ens_to_sym_conv_table__XX_v2) %in% ens_with_multiple] 
View(data.frame(sym = ens_to_sym_conv_table__XX_v2[names(ens_to_sym_conv_table__XX_v2) %in% ens_with_multiple],
           name = names(ens_to_sym_conv_table__XX_v2)[names(ens_to_sym_conv_table__XX_v2) %in% ens_with_multiple]))

# First check whether the names between rooij-hcgn names and teichmann hcgn-names are equal for those 
# where we found a name
ensembl_IDs_inv = names(ensembl_IDs)
names(ensembl_IDs_inv) = ensembl_IDs
genes_Rsample_hcgn_asTeich = ensembl_IDs_inv[genes_Rsample_ens]
all(genes_Rsample_hcgn==genes_Rsample_hcgn_asTeich)
genes_Rsample_hcgn[genes_Rsample_hcgn!=genes_Rsample_hcgn_asTeich]
sum(genes_Rsample_hcgn=="")
sum(genes_Rsample_hcgn_asTeich=="")
# what names are different?
df_comp_RTnaming_all = data.frame(
    name_R=genes_Rsample_hcgn[(genes_Rsample_hcgn!=genes_Rsample_hcgn_asTeich)],
    name_RbyT=genes_Rsample_hcgn_asTeich[(genes_Rsample_hcgn!=genes_Rsample_hcgn_asTeich)])

df_comp_RTnaming = data.frame(
    name_R=genes_Rsample_hcgn[(genes_Rsample_hcgn!=genes_Rsample_hcgn_asTeich)&(genes_Rsample_hcgn!="")],
    name_RbyT=genes_Rsample_hcgn_asTeich[(genes_Rsample_hcgn!=genes_Rsample_hcgn_asTeich)&(genes_Rsample_hcgn!="")])
df_comp_RTnaming[1:10,]
# what names did our nameless genes get?
genes_Rsample_hcgn_asTeich[(genes_Rsample_hcgn=="")]
 # "AC093323.1"    "AC020659.1"    "AC007601.1"    "AL358781.1"    "AL450998.2" # etc
length(genes_Rsample_hcgn[(genes_Rsample_hcgn!=genes_Rsample_hcgn_asTeich)&(genes_Rsample_hcgn!="")])



# What if we were to use org.Hs.eg.db
library("org.Hs.eg.db")
# mapIds(org.Hs.eg.db, keys = 'ENSG00000047578', keytype = "ENSEMBL", column="SYMBOL") #  
genes_Rsample_hcgn_by_org.Hs <- mapIds(org.Hs.eg.db, keys = genes_Rsample_ens, keytype = "ENSEMBL", column="SYMBOL") #  
all(genes_Rsample_hcgn_by_org.Hs==genes_Rsample_hcgn_asTeich)
sum(is.na(genes_Rsample_hcgn_by_org.Hs))
    # 508 genes not found
df_comp_Ths_naming = data.frame(
    Hs = genes_Rsample_hcgn_by_org.Hs[genes_Rsample_hcgn_by_org.Hs!=genes_Rsample_hcgn_asTeich],
    Tei = genes_Rsample_hcgn_asTeich[genes_Rsample_hcgn_by_org.Hs!=genes_Rsample_hcgn_asTeich], 
    Roo = genes_Rsample_hcgn[genes_Rsample_hcgn_by_org.Hs!=genes_Rsample_hcgn_asTeich])

df_comp_Ths_naming[df_comp_Ths_naming$Hs!=df_comp_Ths_naming$Tei,][1:20,]

# Check Rooij compliancy


