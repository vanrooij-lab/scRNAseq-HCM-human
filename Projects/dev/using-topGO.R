# BiocManager::install("topGO")
# BiocManager::install('ALL')

library(topGO)

# To use topGO, we need to gene -> GO mappings
Terms = toTable(org.Hs.egGO2EG)
EnsIDs = toTable(org.Hs.egENSEMBL)
Terms$ensembl_id = EnsIDs[Terms$gene_id,]$ensembl_id

# 
# id2go_custom = lapply(unique(Terms$gene_id), function(x) {Terms[Terms$gene_id==x,]$go_id})
if (file.exists(paste0(base_dir,'Rdata/id2go_custom.Rdata'))) {
    load(file=paste0(base_dir,'Rdata/id2go_custom.Rdata'))
} else {
    id2go_custom = lapply(unique(Terms$ensembl_id), function(x) {unique(Terms[Terms$ensembl_id==x,]$go_id)})
    names(id2go_custom) = unique(Terms$ensembl_id)
    save(list='id2go_custom', file=paste0(base_dir,'Rdata/id2go_custom.Rdata'))
}

# For testing purposes
background_genes  = Terms$ensembl_id[1:100]
interesting_genes = sample(background_genes, length(background_genes) / 10)
    # interesting_genes = current_genes_query

# Create gene list of all genes, with factor indicating gois
geneNames = background_genes
geneList <- factor(as.integer(geneNames %in% interesting_genes))
names(geneList) <- geneNames

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
    annot = annFUN.gene2GO, gene2GO = id2go_custom)

gene_scores <- geneScore(GOdata, whichGenes = interesting_genes)


# GSEA test = Kolmogorov-Smirnov = ..
test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
resultKS <- getSigGroups(GOdata, test.stat)

KS_res_table = GenTable(GOdata, #classic = resultFis, 
                KS = resultKS, 
                #weight = resultWeight,
                # orderBy = "weight",
                orderBy = "score", 
                #ranksOf = "classic",
                ranksOf = "KS", 
                topNodes = 20)

View(KS_res_table)
 
# Some other tests
resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)

allRes <- GenTable(GOdata, classic = resultFis, 
                KS = resultKS, 
                weight = resultWeight,
                orderBy = "weight",
                #orderBy = "score", 
                ranksOf = "classic",
                #ranksOf = "KS", 
                topNodes = 20)


