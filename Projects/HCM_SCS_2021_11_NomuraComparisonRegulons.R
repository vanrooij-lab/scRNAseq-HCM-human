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
# SCENIC regulons
# load(paste0(base_dir, 'Rdata/SCENIC_regulons_core_genes_sel.Rdata')) # SCENIC_regulons_core_genes_sel
load(paste0(base_dir, 'Rdata/SCENIC_reg_top_genes_sorted_full.Rdata')) # SCENIC_reg_top_genes_sorted_full
# Custom modules
ANALYSIS_NAME = "ROOIJonly_RID2l"
load(paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_sorted_shortname.Rdata')) # core_regulons_sorted_shortname


# now look at some stats
sum(core_regulons_sorted_shortname[[1]] %in% human_eq_genes_in_Nomura_hypertr_modules)
length(shared_regulon_genes_list[[1]])
module_nomura_overlap = sapply(1:5,function(X){
    sum(core_regulons_sorted_shortname[[X]] %in% human_eq_genes_in_Nomura_hypertr_modules)/length(core_regulons_sorted_shortname[[X]])})
names(module_nomura_overlap) = names(core_regulons_sorted_shortname)
module_nomura_overlap
round(module_nomura_overlap*100,0)

SCENIC_nomura_overlap = sapply(1:25,function(X){
    sum(SCENIC_reg_top_genes_sorted_full[[X]] %in% human_eq_genes_in_Nomura_hypertr_modules)/length(SCENIC_reg_top_genes_sorted_full[[X]])})
names(SCENIC_nomura_overlap)=names(SCENIC_reg_top_genes_sorted_full)
SCENIC_nomura_overlap
sapply(SCENIC_reg_top_genes_sorted_full, length)
round(sort(SCENIC_nomura_overlap, decreasing = T)*100,0)

ggplot(data.frame(regulon=names(module_nomura_overlap), overlap=module_nomura_overlap*100), 
       aes(x=regulon, y=overlap))+
    geom_bar(stat='identity')+theme_bw()+coord_flip()

ggplot(data.frame(regulon=names(SCENIC_nomura_overlap), overlap=SCENIC_nomura_overlap*100), 
       aes(x=regulon, y=overlap))+
    geom_bar(stat='identity')+theme_bw()+coord_flip()

################################################################################

# Create Nomura overlap sets of interest

# Set 1: module 2 - nomura overlapping genes
genelist_module2_nomura_overlap = 
    core_regulons_sorted_shortname[[2]][core_regulons_sorted_shortname[[2]] %in% human_eq_genes_in_Nomura_hypertr_modules]
genelist_module2_nomura_overlap

save(list='genelist_module2_nomura_overlap', file = paste0(base_dir,'Rdata/zcustom__genelist_module2_nomura_overlap.Rdata'))
# load(paste0(base_dir,'Rdata/zcustom__genelist_module2_nomura_overlap.Rdata')) # genelist_module2_nomura_overlap

# Set 2: SCENIC top overlapping - nomura overlapping genes
genelist_SCENIC.FOXN3_nomura_overlap = 
    SCENIC_reg_top_genes_sorted_full[['FOXN3']][SCENIC_reg_top_genes_sorted_full[['FOXN3']] %in% human_eq_genes_in_Nomura_hypertr_modules]
genelist_SCENIC.FOXN3_nomura_overlap

save(list='genelist_SCENIC.FOXN3_nomura_overlap', file = paste0(base_dir,'Rdata/zcustom__genelist_SCENIC.FOXN3_nomura_overlap.Rdata'))
# load(paste0(base_dir,'Rdata/zcustom__genelist_SCENIC.FOXN3_nomura_overlap.Rdata')) # genelist_SCENIC.FOXN3_nomura_overlap







