

for (DATA_TYPE in c('CLUSTERS', 'REGULONS')) {
    
    # DATA_TYPE = 'REGULONS'
    # DATA_TYPE = 'CLUSTERS'
    
    homerData_dflist=list()
    
    idxlist = list('CLUSTERS'=1:6,'REGULONS'=1:5)[[DATA_TYPE]]
    for (idx in idxlist) {
    
        # idx=1
        
        # /Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/Homer/CLUSTERS/output_cl.2/knownResults.txt
        # /Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/Homer/CLUSTERS/output_cl.2/homerResults.html
        
        #homerData_dflist[[idx]] =
        #    read.csv(paste0(base_dir, 'Homer/CLUSTERS/output_cl.2/knownResults.txt'), sep = '\t', header = 1)
        
        # Read the known results
        string1 = list('CLUSTERS'='cl','REGULONS'='s.R')[[DATA_TYPE]]
        homer_colnames =
            read.table(paste0(base_dir, 'Homer/', DATA_TYPE, '/output_',string1,'.',idx,'/knownResults.txt'), sep = '\t', nrow = 1, comment.char = '')
        homerData_dflist[[idx]] =
            read.csv(paste0(base_dir, 'Homer/', DATA_TYPE, '/output_',string1,'.',idx,'/knownResults.txt'), sep = '\t', skip = 1, header = F)
        colnames(homerData_dflist[[idx]])  = homer_colnames
        # 
        colnames(homerData_dflist[[idx]]) = gsub(pattern = '\\(of [0-9]*\\)',replacement = '',x = colnames(homerData_dflist[[idx]]))
        colnames(homerData_dflist[[idx]]) = gsub(pattern = ' ',replacement = '\\.',x = colnames(homerData_dflist[[idx]]))
        # 
        homerData_dflist[[idx]]$cluster = paste0(string1, idx)
    
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
    homer_simple_export_list = lapply(homerData_df_sel_list, function(X) {len=length(X$Motif.Name_short2); X$Motif.Name_short2[1:min(5, len)]})
    
    ###
    
    # Export known results
    homer_export_df_simple =
        data.frame(sapply(homer_simple_export_list, toString), names(homer_simple_export_list))
    colnames(homer_export_df_simple) = paste0(list('CLUSTERS'='cluster','REGULONS'='regulon')[[DATA_TYPE]],c('_TF',''))
    
    
    openxlsx::write.xlsx(x=homer_export_df_simple, file=paste0(base_dir,'Rplots/Homer_', DATA_TYPE, '_short_summary.xlsx'), overwrite = T)

}




