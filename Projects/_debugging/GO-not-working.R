




# source('/Users/m.wehrens/Documents/git_repos/SCS_Joep/Functions/functions_mw_copy/MW_GOKEGG_analysis.R')

if (F) {
    install.packages('BiocManager') # V
    BiocManager::install('GO.db', force=T) # V
    BiocManager::install('org.Hs.eg.db', force=T) # V
    BiocManager::install('GSEABase', force=T) # V
    BiocManager::install('Category', force=T) # P
}

library(GO.db)
library(org.Hs.eg.db)
library(GSEABase)
library(Category)

Terms = toTable(org.Hs.egGO2EG)
includeChildTerms = 'WithoutChildren'

termFrameData = data.frame(Terms$go_id, Terms$Evidence, Terms$gene_id)
termFrame = GOFrame(termFrameData,organism='human')
termFrame = GOAllFrame(termFrame)
termGSC = GeneSetCollection(termFrame, setType = GOCollection())

# For testing
genes_query=c("ENSG00000084754", "ENSG00000084764", "ENSG00000085063", "ENSG00000085117", "ENSG00000085224", "ENSG00000085231")
backgroundGenes=genes_query

params <- GSEAGOHyperGParams(
            name="GO GSEA",
            geneSetCollection=termGSC,
            geneIds = genes_query,
            universeGeneIds = backgroundGenes,
            ontology = 'BP',
            pvalueCutoff = 0.05,
            conditional = T,
            testDirection = 'over'
          )


# S4Vectors
