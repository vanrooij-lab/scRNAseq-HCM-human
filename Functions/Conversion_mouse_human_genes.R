# Basic function to convert mouse to human gene names
require("biomaRt")

# rnorvegicus_gene_ensembl
convertMouseGeneList <- function(x, the_other_animal = 'mouse', make_unique=T, table_dump=F) { # , VERSION='93'){

    if (the_other_animal == 'mouse') {
        the_ensembl = 'mmusculus_gene_ensembl'
        the_symbol = 'mgi_symbol'
    }
    if (the_other_animal == 'rat') {
        the_ensembl = 'rnorvegicus_gene_ensembl'
        the_symbol = 'rgd_symbol'
    }
    
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    other_animal = useMart("ensembl", dataset = the_ensembl)
    
    
    #human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version=VERSION)
    #mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
    
    
    genesV2 = getLDS(attributes = c(the_symbol), filters = the_symbol, values = x , 
                     mart = other_animal, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    
    if (table_dump) {
        return(genesV2)
    }
    
    humanx = genesV2[,2]
    
    if (make_unique) {
        humanx <- unique(humanx)
    }
    
    # Print the first 6 genes found to the screen
    print(head(humanx))
    return(humanx)
    
}


