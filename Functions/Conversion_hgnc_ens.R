########################################################################

# Conversion ensembl names to hgnc symbols and vice versa,
# using biomaRt
# 
# 

# Note: one can use
# View(listAttributes(mart=mart))
# to see which value you can get from mart

# m.wehrens@hubrecht.eu

########################################################################

library(biomaRt)
        
# change this to your favorite directory
SAVEDIR = "/Users/m.wehrens/Data/__resources/"

########################################################################

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
conversion_table <- getBM(
  attributes=c("hgnc_symbol","ensembl_gene_id"),
  mart = mart)
head(conversion_table)

# Look at table
# View(conversion_table)

# The above code takes a while, better to save the table to Rdata file
# and load it when you need it, see below

########################################################################

# Create a conversion list
conversion_ens2hgnc = 
    conversion_table$hgnc_symbol
names(conversion_ens2hgnc) = conversion_table$ensembl_gene_id

# example usage
conversion_ens2hgnc[c('ENSG00000210049', 'ENSG00000209082')]

########################################################################

# Other way around (using is.na beaause some values are missing)
conversion_hgnc2ens = 
    conversion_table$ensembl_gene_id[conversion_table$hgnc_symbol != ""]
names(conversion_hgnc2ens) = conversion_table$hgnc_symbol[conversion_table$hgnc_symbol != ""]

# example usage
conversion_hgnc2ens[c('MT-CO2','MYH7')]

########################################################################

# save for later usage
save(list = c('conversion_ens2hgnc', 'conversion_hgnc2ens'),file=paste0(SAVEDIR,'conversion_hgnc2ens_ens2hgnc.Rdata'))
save(list=c('conversion_table'),file=paste0(SAVEDIR,'conversion_hgnc_ens_table.Rdata'))
    # note that the conversion_table is not really necessary any more
    
# then later use 
load(file=paste0(SAVEDIR,'conversion_hgnc2ens_ens2hgnc.Rdata'))
    # which will load the conversion vectors conversion_ens2hgnc and conversion_hgnc2ens

########################################################################



