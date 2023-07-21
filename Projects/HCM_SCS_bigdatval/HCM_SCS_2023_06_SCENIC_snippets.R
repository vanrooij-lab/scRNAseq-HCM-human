

# Little overview of which TFs are members of other regulons etc.




###

# Fetch module genes
load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3b/Rdata/ROOIJonly.sp.bt_RID2l_core_regulons_sorted.Rdata')

###

ourTFs=c('ANKRD1', 'NFE2L1','MAFK', 'NRAP', 'ZEB1' ,'FOXN3','MEF2A','MEF2B','MEF2C','MEF2D')
# ourTFs=names(SCENIC_regulons_core_genes_sel)

membersofthisregulon =
    lapply(ourTFs, function(X) {ourTFs[ourTFs %in% SCENIC_regulons_core_genes_sel[[X]]]})
names(membersofthisregulon)=ourTFs

regulon_sizes=sapply(ourTFs, function(X){length(SCENIC_regulons_core_genes_sel[[X]])})
names(regulon_sizes) = names(SCENIC_regulons_core_genes_sel)

membersofthisregulon
regulon_sizes



####

ourothergenesinterest = c('ANKRD1','CD44','XIRP2')
membersofthisregulon2 =
    lapply(ourTFs, function(X) {ourothergenesinterest[ourothergenesinterest %in% SCENIC_regulons_core_genes_sel[[X]]]})
names(membersofthisregulon2)=ourTFs
membersofthisregulon2




ourothergenesinterest = c('ACTN2','CTNNA3','TP53INP2','CCND1','AHNAK')
membersofthisregulon3 =
    lapply(ourTFs, function(X) {ourothergenesinterest[ourothergenesinterest %in% SCENIC_regulons_core_genes_sel[[X]]]})
names(membersofthisregulon3)=ourTFs
membersofthisregulon3

