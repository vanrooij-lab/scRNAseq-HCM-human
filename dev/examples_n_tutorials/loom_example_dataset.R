
# Loom example dataset

loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom")
Open the loom file and load the expression matrix (and cell annotation if available)

library(SCopeLoomR)
loom <- open_loom(loomPath)
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

dim(exprMat)