# Let's obtain the nomura data
Nomura_data = openxlsx::read.xlsx('/Users/m.wehrens/Data/_2019_02_HCM_SCS/Nomura_data/41467_2018_6639_MOESM7_ESM.xlsx')
genes_in_Nomura_hypertr_modules = Nomura_data[Nomura_data$assigned.module %in% c('M1', 'M2', 'M5', 'M11', 'M16'),]$gene.name

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){

require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
humanx <- unique(genesV2[, 2])

# Print the first 6 genes found to the screen
print(head(humanx))
return(humanx)
}

human_eq_genes_in_Nomura_hypertr_modules = convertMouseGeneList(genes_in_Nomura_hypertr_modules)



# Now load the SCENIC regulons and custom modules



# now look at some stats
sum(shared_regulon_genes_list[[1]] %in% human_eq_genes_in_Nomura_hypertr_modules)
length(shared_regulon_genes_list[[1]])
sum(shared_regulon_genes_list[[1]] %in% human_eq_genes_in_Nomura_hypertr_modules)/length(shared_regulon_genes_list[[1]])
sum(shared_regulon_genes_list[[2]] %in% human_eq_genes_in_Nomura_hypertr_modules)/length(shared_regulon_genes_list[[2]])
sum(shared_regulon_genes_list[[3]] %in% human_eq_genes_in_Nomura_hypertr_modules)/length(shared_regulon_genes_list[[3]])
sum(shared_regulon_genes_list[[4]] %in% human_eq_genes_in_Nomura_hypertr_modules)/length(shared_regulon_genes_list[[4]])

sum(shared_regulon_core_list$SharedRegulon1 %in% human_eq_genes_in_Nomura_hypertr_modules)/length(shared_regulon_core_list$SharedRegulon1)
sum(shared_regulon_core_list$SharedRegulon2 %in% human_eq_genes_in_Nomura_hypertr_modules)/length(shared_regulon_core_list$SharedRegulon2)
sum(shared_regulon_core_list$SharedRegulon3 %in% human_eq_genes_in_Nomura_hypertr_modules)/length(shared_regulon_core_list$SharedRegulon3)
sum(shared_regulon_core_list$SharedRegulon4 %in% human_eq_genes_in_Nomura_hypertr_modules)/length(shared_regulon_core_list$SharedRegulon4)





