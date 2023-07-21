

################################################################################
# I think this is the easiest way
#
# This is based on data downloaded from http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt


mouse_human_genes_df_reOrganized = 
    read.table(paste0('/Users/m.wehrens/Documents/git_repos/Resources_bioinf_RNAseq/data_tsv/','mouse_human_genes_df_reOrganized__custom_rerganized_jax.tsv'))
    
conversion_human_to_mouse_symbols = mouse_human_genes_df_reOrganized$Symbol_mouse
names(conversion_human_to_mouse_symbols) = mouse_human_genes_df_reOrganized$Symbol_human
conversion_mouse_to_human_symbols = mouse_human_genes_df_reOrganized$Symbol_human
names(conversion_mouse_to_human_symbols) = mouse_human_genes_df_reOrganized$Symbol_mouse


if (F) {
    
    mouse_human_genes_df = read.csv("/Users/m.wehrens/Data/__resources/HOM_MouseHumanSequence.rpt.txt",sep="\t")
    
    mouse_human_genes_df_mouseOnly = mouse_human_genes_df[grepl('mouse', mouse_human_genes_df$Common.Organism.Name),]
    mouse_human_genes_df_humanOnly = mouse_human_genes_df[grepl('human', mouse_human_genes_df$Common.Organism.Name),]
    
    View(mouse_human_genes_df_mouse)
    
    mouse_human_genes_df_mouseOnly_compact = mouse_human_genes_df_mouseOnly[,c('DB.Class.Key','Symbol')]
    mouse_human_genes_df_humanOnly_compact = mouse_human_genes_df_humanOnly[,c('DB.Class.Key','Symbol')]
    
    mouse_human_genes_df_reOrganized = 
        merge(mouse_human_genes_df_mouseOnly_compact, mouse_human_genes_df_humanOnly_compact, by='DB.Class.Key', suffixes=c('_mouse','_human'))
    
    write.table(x = mouse_human_genes_df_reOrganized, file = 
                    paste0('/Users/m.wehrens/Documents/git_repos/Resources_bioinf_RNAseq/data_tsv/','mouse_human_genes_df_reOrganized__custom_rerganized_jax.tsv'))

}

################################################################################

# Files manually downloaded at http://www.ensembl.org/biomart/martview
# Dataset: Ensembl genes 107, Human genes (GRCh38.p13)
#
# This links ensembl IDs of human with those of mice

# From tutorial at https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
# 1. Go to: http://www.ensembl.org/biomart/martview
# 2. Choose “Ensembl 52”
# 3. Choose “Homo sapiens genes (NCBI36)”
# 4. Click on “Filters” in the left menu
# 5. Unfold the “MULTI SPECIES COMPARISONS” box, tick the “Homolog filters” option and choose 
# “Orthologous Mouse Genes” from the drop-down menu.
# 6. Click on “Attributes” in the left menu
# 7. Click on “Homologs”
# 8. Unfold the “MOUSE ORTHOLOGS” box and select the data you want to get (most probably the gene ID and maybe the orthology type as well).
# 9. Click on the “Results” button (top left)
# 10. Choose your favorite output



mouse_human_conversion_ensmanual_df =
    read.table('/Users/m.wehrens/Documents/git_repos/Resources_bioinf_RNAseq/data_tsv/ensembl_biomart_manual_conversion__v.107__GRCh38.p13.tsv', sep = '\t', header = 1)

convert_enshuman_to_ensmouse_lookup = mouse_human_conversion_ensmanual_df$Mouse.gene.stable.ID
names(convert_enshuman_to_ensmouse_lookup) = mouse_human_conversion_ensmanual_df$Gene.stable.ID

# ps. note that there is also annotation for human symbols from the same source, but it's probably
# easier to use earlier generated tables for this..


if (F) {
    conversion_ens_sym_107manual_df =
        read.table('/Users/m.wehrens/Desktop/gene_conversion/ensembl_biomart_manual_conversion_allgenes__v.107__GRCh38.p13.tsv', sep = '\t', header = 1)
    
    conversion_ens_sym_107manual_lookup = conversion_ens_sym_107manual_df$Gene.name
    names(conversion_ens_sym_107manual_lookup) = conversion_ens_sym_107manual_df$Gene.stable.ID
}
################################################################################


# Basic function to convert mouse to human gene names
require("biomaRt")

# rnorvegicus_gene_ensembl
# Mostly copied code from https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
convertMouseGeneList <- function(x, the_other_animal = 'mouse', make_unique=T, table_dump=F) { # , VERSION='93'){

    if (the_other_animal == 'mouse') {
        the_ensembl = 'mmusculus_gene_ensembl'
        the_symbol = 'mgi_symbol'
    }
    if (the_other_animal == 'rat') {
        the_ensembl = 'rnorvegicus_gene_ensembl'
        the_symbol = 'rgd_symbol'
    }
    
    human        = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    other_animal = useMart("ensembl", dataset = the_ensembl)
    
    # human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version=VERSION)
    # mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
    
    
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

#### 


convertMouseGeneList_usingJAX = function(x) {
    
    # Download it:
    # mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
    # Or use local:
    mouse_human_genes = read.csv("/Users/m.wehrens/Data/__resources/HOM_MouseHumanSequence.rpt.txt",sep="\t")
    
    #sapply(x,
    #function(g) {
    #    mouse_human_genes[mouse_human_genes$Common.Organism.Name == 'human' & mouse_human_genes$symbol == g,]$DB.Class.Key}
    #)
    
    result1 =
    mouse_human_genes[mouse_human_genes$DB.Class.Key %in%
        mouse_human_genes[mouse_human_genes$Symbol %in% x,1:7]$DB.Class.Key,1:7]
    
    result = sort(result1[result1$Common.Organism.Name=='mouse, laboratory',]$Symbol)
    
    return(result)
    
}









