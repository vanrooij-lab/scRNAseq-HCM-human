

homerData_dflist=list()

for (idx in 1:6) {

    # /Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/Homer/CLUSTERS/output_cl.2/knownResults.txt
    # /Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/Homer/CLUSTERS/output_cl.2/homerResults.html
    
    #homerData_dflist[[idx]] =
    #    read.csv(paste0(base_dir, 'Homer/CLUSTERS/output_cl.2/knownResults.txt'), sep = '\t', header = 1)
    
    # Read the known results
    homer_colnames =
        read.table(paste0(base_dir, 'Homer/CLUSTERS/output_cl.',idx,'/knownResults.txt'), sep = '\t', nrow = 1, comment.char = '')
    homerData_dflist[[idx]] =
        read.csv(paste0(base_dir, 'Homer/CLUSTERS/output_cl.',idx,'/knownResults.txt'), sep = '\t', skip = 1)
    colnames(homerData_dflist[[idx]])  = homer_colnames
    # 
    colnames(homerData_dflist[[idx]]) = gsub(pattern = '\\(of [0-9]*\\)',replacement = '',x = colnames(homerData_dflist[[idx]]))
    colnames(homerData_dflist[[idx]]) = gsub(pattern = ' ',replacement = '\\.',x = colnames(homerData_dflist[[idx]]))
    # 
    homerData_dflist[[idx]]$cluster = paste0('cl.',idx)

}    
homerData_df = Reduce(x = homerData_dflist, f = rbind)

# Add brief gene names
homerData_df$Motif.Name_short2 = 
    gsub(pattern = '^[^/]*/[^-]*-', replacement = '', x = homerData_df$Motif.Name)
homerData_df$Motif.Name_short2 = 
    gsub(pattern = '[-/\\.(].*$', replacement = '', x = homerData_df$Motif.Name_short2)
homerData_df$Motif.Name_short2 = toupper( homerData_df$Motif.Name_short2 )
    


homerData_df_sel = homerData_df[homerData_df$`P-value`<.05&homerData_df$`q-value.(Benjamini)`<.9,]
homerData_df_sel_list = split(homerData_df_sel, f=homerData_df_sel$cluster)
homer_simple_export_list = sapply(homerData_df_sel_list, function(X) {len=length(X$Motif.Name_short2); X$Motif.Name_short2[1:min(5, len)]})

################################################################################

# Export known results
homer_export_df_simple =
    data.frame(cluster_TFs=sapply(homer_simple_export_list, toString), cluster=names(homer_simple_export_list))


openxlsx::write.xlsx(x=homer_export_df_simple, file=paste0(base_dir,'Rplots/Homer_clusters_short_summary.xlsx'), overwrite = T)

################################################################################




