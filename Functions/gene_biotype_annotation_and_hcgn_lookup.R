
########################################################################

# Lookup table for biotypes of genes
# Lookup table for hcgn symbol for ens IDs
# 
# m.wehrens@hubrecht.eu

# Re-written 22-12-2021

########################################################################

library(biomaRt)
        
########################################################################
# BIOPTYPE table

# mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version="93")

ENSEMBL_VERSION_BIOTYPES = 93

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   version = paste0(ENSEMBL_VERSION_BIOTYPES))

biotype_table <- getBM(
  attributes=c("ensembl_gene_id", 'gene_biotype'), # "hgnc_symbol",
  mart = mart 
  )
head(biotype_table)

gene_biotypes_XX_v2 = biotype_table$gene_biotype
names(gene_biotypes_XX_v2) = biotype_table$ensembl_gene_id
  
# Now save as v2, since another script perviously generated this file
save(list='gene_biotypes_XX_v2', file = paste0(base_dir,'Rdata/gene_biotypes_',ENSEMBL_VERSION_BIOTYPES,'_v2.Rdata')) 
# save(list='gene_biotypes_XX', file = paste0(base_dir,'Rdata/gene_biotypes_',ENSEMBL_VERSION_BIOTYPES,'.Rdata')) 

# Sanity check to check whether version are equal
# load(file = paste0(base_dir,'Rdata/gene_biotypes_',ENSEMBL_VERSION_BIOTYPES,'.Rdata')) # v1 
# all(names(gene_biotypes_XX_v2)==names(gene_biotypes_XX)) # TRUE
# all(gene_biotypes_XX_v2==gene_biotypes_XX) # TRUE
# Conclusion: files v1 and v2 are exactly equal, and can be used interchangibly.

########################################################################


ENSEMBL_VERSION_BIOTYPES = 93

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   version = paste0(ENSEMBL_VERSION_BIOTYPES))

ens_to_sym_conv_table_ <- getBM(
  attributes=c("ensembl_gene_id", 'hgnc_symbol'), # "hgnc_symbol",
  mart = mart 
  )
head(ens_to_sym_conv_table_)

ens_to_sym_conv_table__XX_v2 = ens_to_sym_conv_table_$hgnc_symbol
names(ens_to_sym_conv_table__XX_v2) = ens_to_sym_conv_table_$ensembl_gene_id
  
# Now save as v2, since another script perviously generated this file
save(list='ens_to_sym_conv_table__XX_v2', file = paste0(base_dir,'Rdata/ens_to_sym_conv_table__',ENSEMBL_VERSION_BIOTYPES,'_v2.Rdata')) 
# load(file = paste0(base_dir,'Rdata/ens_to_sym_conv_table__','93','_v2.Rdata'))  # ens_to_sym_conv_table__XX_v2

########################################################################

# Sanity check to check whether version are equal
# load(paste0(base_dir,'Rdata/ens_to_sym_conv_table.Rdata')) # v1
# all(names(ens_to_sym_conv_table)==names(ens_to_sym_conv_table__XX_v2)) # TRUE
# all(ens_to_sym_conv_table==ens_to_sym_conv_table__XX_v2) # TRUE
# Conclusion: files v1 and v2 are NOT exactly equal
# (This is probably because either ENS not linked to hcgn in v2, or becuase other synonym was now chosen)
# 

View(ens_to_sym_conv_table)
View(ens_to_sym_conv_table__XX_v2)

# Some genes are not in one but are in the other and vice versa
all(names(ens_to_sym_conv_table) %in% names(ens_to_sym_conv_table__XX_v2))
all(names(ens_to_sym_conv_table__XX_v2) %in% names(ens_to_sym_conv_table))

shared_genes = names(ens_to_sym_conv_table__XX_v2)[names(ens_to_sym_conv_table__XX_v2) %in% names(ens_to_sym_conv_table)]

# But is the lookup the same?
all(ens_to_sym_conv_table[shared_genes] == ens_to_sym_conv_table__XX_v2[shared_genes])
unequal_info = shared_genes[ens_to_sym_conv_table[shared_genes] != ens_to_sym_conv_table__XX_v2[shared_genes]]

length(unequal_info)

ens_to_sym_conv_table[unequal_info]
ens_to_sym_conv_table__XX_v2[unequal_info]
sum(ens_to_sym_conv_table__XX_v2[unequal_info]=="")
sum(ens_to_sym_conv_table__XX_v2[unequal_info]!="")

ens_to_sym_conv_table["ENSG00000166398"]
ens_to_sym_conv_table__XX_v2["ENSG00000166398"]


########################################################################

if (F) {
  
  # Let's do for latest version ()
  # Which is 105 currently speaking
  
  listEnsemblArchives()
  
  mart_latest <- useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl") 
                     # version = paste0(ENSEMBL_VERSION_BIOTYPES))
  
  ens_to_sym_conv_table_latest_ <- getBM(
    attributes=c("ensembl_gene_id", 'hgnc_symbol'), # "hgnc_symbol",
    mart = mart_latest 
    )
  
  ens_to_sym_conv_table__latest_v2 = ens_to_sym_conv_table_latest_$hgnc_symbol
  names(ens_to_sym_conv_table__latest_v2) = ens_to_sym_conv_table_latest_$ensembl_gene_id
    
  
  ens_to_sym_conv_table__XX_v2['ENSG00000047578']
  # ENSG00000047578 
  #      "KIAA0556" 
  ens_to_sym_conv_table__latest_v2['ENSG00000047578']
  # ENSG00000047578 
  #        "KATNIP"
  
  # And NIBAN/FAM
  ens_to_sym_conv_table__latest_v2['ENSG00000135842']
  # ENSG00000135842 
  #        "NIBAN1" 
  ens_to_sym_conv_table__XX_v2['ENSG00000135842']
  # ENSG00000135842 
  #       "FAM129A"
  
  # So this is probably where the name discrepancy (a problem that is now fixed) came from ..

}


