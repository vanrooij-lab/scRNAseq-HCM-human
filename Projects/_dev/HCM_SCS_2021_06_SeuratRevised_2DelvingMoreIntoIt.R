
################################################################################
# I'll have to look into multiple things:

# CHECKS
# v What about these MT pseudogenes that pop up since the re-mapping?  ---> are also in Teichmann
# _ how does clsuter/patient overlap look for the "default" seurat settings?
#       _ project default clustering on the rid2 one
# _ check what the input requirements are for DE analysis (FC values seem off)

# "REAL" ANALYSES
#       (basically, question is, can we reproduce what we saw before)
# v DE genes 
# _ TTN correlated genes
# _ Perform analyses on Hu/Teichmann separately
# _ GO analysis 

# More optional currently
# _ Jaccard score?
# _ QC plots?

################################################################################
# PSEUDOGENES

# Let's check whether the MT pseudogenes were also mapped in the previous data,
# and whether they're present in the data mapped by the Teichmann group

# Check our previous mapping
load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_sessions/HCM2020may_re-analysis-v4_groupedSCS-object-only.Rdata')
myrownames_SCS=rownames(groupedSCS$patientAllMod@ndata)
myrownames_SCS[grepl('^MTND1P23_',myrownames_SCS)] # pseudo gene
myrownames_SCS[grepl('^MTND4P12_',myrownames_SCS)] # pseudo gene
myrownames_SCS[grepl('^MTATP6P1',myrownames_SCS)] # pseudo gene
myrownames_SCS[grepl('^MTRNR2L1_',myrownames_SCS)] # protein coding

# These genes are all detected
myrownames_TEICH = rownames(current_analysis[[seuratObject_name]])
myrownames_TEICH[grepl('^MTND1P23$',myrownames_TEICH)] # pseudo gene
myrownames_TEICH[grepl('^MTND4P12$',myrownames_TEICH)] # pseudo gene
myrownames_TEICH[grepl('^MTATP6P1$',myrownames_TEICH)] # pseudo gene




