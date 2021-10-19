
# Apply command gives an issue:
#
# gene_in_cell_count = apply(RHL_SeuratObject_merged_noMito_sel@assays$RNA@counts>0, 1, sum)
#  *** caught segfault ***
# address 0x2b19aee3205c, cause 'memory not mapped'

library(Seurat)
library(SeuratDisk)

base_dir = '/hpc/hub_oudenaarden/mwehrens/data/Wang/'


########################################################################
# Code that gives the error:

# Load stuff
library(Seurat)
base_dir = '/hpc/hub_oudenaarden/mwehrens/data/Wang/'
load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_merged.Rdata'))

# Actuall error line
print('here we go')
gene_in_cell_count = apply(RHL_SeuratObject_merged@assays$RNA@counts>0, 1, sum)

# Trying to pinpoint
testMat = (RHL_SeuratObject_merged@assays$RNA@counts, "matrix")



########################################################################
# test case (works fine)

fakematrix = matrix(c(rep(0,3),1:9),nrow = 3)
rownames(fakematrix) = c('gene1','gene2','gene3')
colnames(fakematrix) = c('cell1','cell2','cell3','cell4')
testObject = CreateSeuratObject(fakematrix)

# no problem
gene_in_cell_count = apply(testObject@assays$RNA@counts>0, 1, sum)

# Save as H5
SaveH5Seurat(object = testObject, overwrite = T,
    filename = paste0(base_dir,'Rdata/testObject.h5seurat'))

# Load the H5 and try apply command
testObject_H5 = LoadH5Seurat(file = paste0(base_dir,'Rdata/testObject.h5seurat'))
gene_in_cell_count = apply(testObject_H5@assays$RNA@counts>0, 1, sum)

########################################################################
# Test case w/ earlier version of the object

# Load separate objects
load(file = paste0(base_dir,'Rdata/RHL_SeuratObject_list.Rdata'))
gene_in_cell_count = apply(RHL_SeuratObject_list$JE10@assays$RNA@counts>0, 1, sum)  # Works
gene_in_cell_count = apply(RHL_SeuratObject_list$TEICH@assays$RNA@counts>0, 1, sum)  # Fails

# I think the problem lies in sparse matrices that are large
#
# testMat = as.matrix(RHL_SeuratObject_merged@assays$RNA@counts)
# Error in asMethod(object) : 
#   Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 102
# 
# as.matrix also gives an error, but specifically throws an error also
# (see below)

# So the solution is to use "rowSums" and "colSums" since
# those are written specifically for sparse matrices.

gene_in_cell_count = rowSums(RHL_SeuratObject_list$TEICH@assays$RNA@counts>0)  # works

# PS.
# https://stackoverflow.com/questions/49190251/caught-segfault-memory-not-mapped-error-in-r
# suggested it was a dependency issue (already compiled deps not properly doing stuff)
# but I don't think that's the case here.

########################################################################
# Installed package "Matrix" @ HPC, 

# testMat = as.matrix(RHL_SeuratObject_merged@assays$RNA@counts)
# Error in asMethod(object) : 
#   Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 102

instPkgPlusDeps <- function(pkg, install = FALSE,
                            which = c("Depends", "Imports", "LinkingTo"),
                            inc.pkg = TRUE) {
  stopifnot(require("tools")) ## load tools
  ap <- available.packages() ## takes a minute on first use
  ## get dependencies for pkg recursively through all dependencies
  deps <- package_dependencies(pkg, db = ap, which = which, recursive = TRUE)
  ## the next line can generate warnings; I think these are harmless
  ## returns the Priority field. `NA` indicates not Base or Recommended
  pri <- sapply(deps[[1]], packageDescription, fields = "Priority")
  ## filter out Base & Recommended pkgs - we want the `NA` entries
  deps <- deps[[1]][is.na(pri)]
  ## install pkg too?
  if (inc.pkg) {
    deps = c(pkg, deps)
  }
  ## are we installing?
  if (install) {
    install.packages(deps)
  }
  deps ## return dependencies
}

########################################################################
# I did use the following to reinstall seurat plus all deps

# instPkgPlusDeps('Seurat', install=T)
# 
#   [1] "Seurat"          "cowplot"         "fitdistrplus"    "future"         
#   [5] "future.apply"    "ggplot2"         "ggrepel"         "ggridges"       
#   [9] "httr"            "ica"             "igraph"          "irlba"          
#  [13] "jsonlite"        "leiden"          "lmtest"          "matrixStats"    
#  [17] "miniUI"          "patchwork"       "pbapply"         "plotly"         
#  [21] "png"             "RANN"            "RColorBrewer"    "Rcpp"           
#  [25] "RcppAnnoy"       "reticulate"      "rlang"           "ROCR"           
#  [29] "Rtsne"           "scales"          "scattermore"     "sctransform"    
#  [33] "SeuratObject"    "shiny"           "spatstat.core"   "spatstat.geom"  
#  [37] "tibble"          "uwot"            "RcppEigen"       "RcppProgress"   
#  [41] "gtable"          "digest"          "globals"         "listenv"        
#  [45] "parallelly"      "glue"            "isoband"         "withr"          
#  [49] "plyr"            "curl"            "mime"            "openssl"        
#  [53] "R6"              "magrittr"        "pkgconfig"       "zoo"            
#  [57] "htmltools"       "viridisLite"     "base64enc"       "htmlwidgets"    
#  [61] "tidyr"           "dplyr"           "vctrs"           "lazyeval"       
#  [65] "crosstalk"       "purrr"           "data.table"      "promises"       
#  [69] "rappdirs"        "gplots"          "farver"          "labeling"       
#  [73] "lifecycle"       "munsell"         "reshape2"        "gridExtra"      
#  [77] "RcppArmadillo"   "httpuv"          "xtable"          "sourcetools"    
#  [81] "later"           "crayon"          "fastmap"         "commonmark"     
#  [85] "bslib"           "cachem"          "ellipsis"        "spatstat.data"  
#  [89] "spatstat.utils"  "spatstat.sparse" "abind"           "tensor"         
#  [93] "goftest"         "deldir"          "polyclip"        "fansi"          
#  [97] "pillar"          "FNN"             "RSpectra"        "dqrng"          
# [101] "sass"            "jquerylib"       "generics"        "tidyselect"     
# [105] "BH"              "sitmo"           "gtools"          "caTools"        
# [109] "yaml"            "colorspace"      "askpass"         "cli"            
# [113] "utf8"            "stringr"         "cpp11"           "sys"            
# [117] "bitops"          "fs"              "stringi"    
# Warning messages:
# 1: In install.packages(deps) :
#   installation of package ‘fansi’ had non-zero exit status
# 2: In install.packages(deps) :
#   installation of package ‘Rtsne’ had non-zero exit status






