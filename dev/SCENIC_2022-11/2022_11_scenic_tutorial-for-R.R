### Load data
loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom")
library(SCopeLoomR)
loom <- open_loom(loomPath)
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

### Initialize settings
library(SCENIC)
scenicOptions <- initializeScenic(org="mgi", dbDir="cisTarget_databases", nCores=10)
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

# Optional: Binarize activity
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
# savedSelections <- shiny::runApp(aucellApp)
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

# Export:
# saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo")) # Temporary, to add to loom
export2loom(scenicOptions, exprMat)

# To save the current status, or any changes in settings, save the object again:
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Exploring output 
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox8"]
viewMotifs(tableSubset) 

# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 

# Cell-type specific regulators (RSS): 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)