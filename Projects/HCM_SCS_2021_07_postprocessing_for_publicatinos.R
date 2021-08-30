

# Some custom post-processing

cluster_toptable=openxlsx::read.xlsx('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/Rplots/ClusterTopHits_ROOIJonly_RID2l.xlsx')

gene_symbol_table = 
    sapply(cluster_toptable, function(X) {
        sapply(str_split(X, pattern = ':'), function(Y) {Y[[2]]})
        })

openxlsx::write.xlsx(x = gene_symbol_table,
                     file= '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/Rplots/ClusterTopHits_ROOIJonly_RID2l_symbols.xlsx')
