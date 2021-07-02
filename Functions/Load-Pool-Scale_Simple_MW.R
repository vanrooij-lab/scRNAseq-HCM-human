
library(stringr)
library(parallel)

# For convenience, little fn to change text size ggplot in one go
give_better_textsize_plot <- function(TEXTSIZE){
  theme(#legend.position="none",
        text = element_text(size=TEXTSIZE),
        axis.text = element_text(size=TEXTSIZE),
        plot.title = element_text(size=TEXTSIZE))
}

loadData_MW = function(dataset_list_paths, toPool = NULL, prefix=T) {
    # input: 
    # - dataset_list_paths: c(sample_name=path1, sample2_name=path2, ..)
    # - toPool: list of vectors that list datasets that should be pooled, e.g.
    # list(group1=c('sample1', 'sample2'), ..)
    # output:
    # groupedData[[groupName]]

    sample_names = names(dataset_list_paths)
    output = list()
    
    for (idx in 1:length(dataset_list_paths)) {
    
        output[[sample_names[idx]]] =
            read.table(dataset_list_paths[[idx]], row.names = 1, header=1)
        if (prefix) { # prefix colnames with sample name
          colnames(output[[sample_names[idx]]]) = 
              paste0(sample_names[idx],'.',colnames(output[[sample_names[idx]]]))
        }
        
        print(paste0('Loaded: ', sample_names[idx]))
        
    }
    
    return(output)
    
}

loadData_MW_parallel = function(dataset_list_paths, mc.cores=1, prefix=T) {
    # input: 
    # - dataset_list_paths: c(sample_name=path1, sample2_name=path2, ..)
    # - toPool: list of vectors that list datasets that should be pooled, e.g.
    # list(group1=c('sample1', 'sample2'), ..)
    # output:
    # groupedData[[groupName]]

    sample_names = names(dataset_list_paths)
    
    output = mclapply(X = 1:length(dataset_list_paths), FUN=function(idx) {
    
        sample_output = 
            read.table(dataset_list_paths[[idx]], row.names = 1, header=1)
        if (prefix) { # prefix colnames with sample name
          colnames(sample_output) = 
              paste0(sample_names[idx],'.',colnames(sample_output))
        }
          
        print(paste0('Loaded: ', sample_names[idx]))
        return(sample_output)
        
    })
    names(output) = sample_names
    
    return(output)
    
}

# multiple_dfs = SCS_df_list_data
pool_df_mw = function(multiple_dfs) {
  # multiple_dfs should be a list of dataframes
    
    # rnames_all = unique(unlist(lapply(multiple_dfs, rownames)))
    pooled_df = Reduce(
        function(x,y) {
          x <- merge(x,y,by.x=0,by.y=0,all=TRUE)
          rownames(x) <- x[,1]
          x[,"Row.names"] <- NULL
          x[is.na(x)] = 0
          print('.')
          return(x)
        },
        multiple_dfs
      )
    
    return(pooled_df)

}

# scaling
manual_scale_table = function(countTable) {
    
    readsPerWell      = apply(countTable,2,sum)
    median_readsPerWell = median(readsPerWell)
    GenesHowManyCells = apply(countTable>0,1,sum)
    #gene_old_selection = sum(GenesHowManyCells>3)
    
    countTable_scaled = 
        sapply(1:dim(countTable)[2], function(X) {
            countTable[,X]/readsPerWell[X]*median_readsPerWell})
    rownames(countTable_scaled) = rownames(countTable)
    colnames(countTable_scaled) = colnames(countTable)
    
    # apply(countTable_scaled,2,sum)
    return(list(countTable_scaled=countTable_scaled, readsPerWell=readsPerWell,
            median_readsPerWell=median_readsPerWell, GenesHowManyCells=GenesHowManyCells ))
}


# So we'll only take along uniquely mapped genes for now ..
# And we'll convert those names to only the symbol
preprocess_convertAAnames_toSymbol = function(df, revert_to_hgnc=F, script_dir=NULL) {
  # set revert_to_hgnc and script_dir to convert names back to hgnc names
  
    # df = SCS_df_list_data$AL1
    # Take only uniquely identified genes
    the_duplicity=(sapply(rownames(df), str_count, pattern='-')+1)
    df=df[the_duplicity==1,]
    # Now take the symbol name
    symbol_names=sapply(strsplit(rownames(df),'_'), function(x) {x[[2]]})
    
    # For now we'll just sum the duplicates that remain
    # (Note: those will have e.g. a different ENSG ID, 
    # probably due to nomenclature issues.)
    # (Note2: presumably, their mapping also was different, so 
    # they will be different molecules, and their UMIs can be
    # added to each other.)
    # (Note3: the issue might also arise due to different mappings
    # in different samples, which is another argument for adding them)
    table_symbol_names_freq = table(symbol_names)
    symbol_names_plus_duplicity = data.frame(name=symbol_names, freq=as.vector(table_symbol_names_freq[symbol_names]))

    # first already update the names that have duplicity one
    new_df_names = rownames(df)
    new_df_names[symbol_names_plus_duplicity$freq==1] = symbol_names[symbol_names_plus_duplicity$freq==1]
    rownames(df) = new_df_names
    
    # then correct the ones that are present dupl>1
    for (name_to_correct in unique(symbol_names_plus_duplicity$name[symbol_names_plus_duplicity$freq>1])) {
        
        current_rows = grepl(paste0('_',name_to_correct,'_'), rownames(df))
        if (sum(current_rows)<1) { stop('Error, ',name_to_correct,' not found or <1..') }
        
        summed_data = apply(df[current_rows,,drop=F],2,sum)
        df=df[!current_rows,]
        df[name_to_correct,]=summed_data
        print(paste0('Correct: ',name_to_correct,' (count=',sum(summed_data),')'))
        
    }
    
    if (revert_to_hgnc) {
      
      lookup_table = name_conversion_anna(script_dir=script_dir)
      newnames = lookup_table[rownames(df)]
      if (any(is.na(newnames))) {
        print(paste0('Keeping unknown names: ',toString(rownames(df)[is.na(newnames)]))) 
        newnames[is.na(newnames)] = rownames(df)[is.na(newnames)]
      }
      rownames(df) = newnames
      
    }
    
    return(df)
    
}



name_conversion_anna = function(script_dir) {
    
    if (file.exists(paste0(script_dir, 'Resources/Anna_name_conversion.Rdata'))) {
        
        print('Using saved name conversion lib ..')
        load(file=paste0(script_dir, 'Resources/Anna_name_conversion.Rdata'))
        
    } else {
        
        print('Generating name conversion lib using biomart ..')
        
        library(biomaRt)
        
        mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        all_gene_symbols <- getBM(
          attributes=c("hgnc_symbol"),
          mart = mart)
        head(all_gene_symbols)
        
        length(all_gene_symbols$hgnc_symbol)
        length(unique(all_gene_symbols$hgnc_symbol))
        
        write.table(file='/Users/m.wehrens/Data/__resources/biomaRt_all_hgnc_symbol.csv', x=all_gene_symbols, quote = F, col.names = F, row.names=F)
        
        #Anna_name_conversion = data.frame(hgnc_symbol=all_gene_symbols$hgnc_symbol, anna_name=gsub(pattern='-','\\.',x = all_gene_symbols$hgnc_symbol))
        #rownames(Anna_name_conversion) = Anna_name_conversion$anna_name
        
        Anna_name_conversion = all_gene_symbols$hgnc_symbol
        names(Anna_name_conversion) = anna_name=gsub(pattern='-','\\.',x = all_gene_symbols$hgnc_symbol)
        
        
        save(file=paste0(script_dir, 'Resources/Anna_name_conversion.Rdata'), list =c('Anna_name_conversion'))
    }
    
    return(Anna_name_conversion)
}



# pool_df_mw = function(multiple_dfs) {
# 
#     print('update this with above!!')
#     
#     rnames_all = unique(unlist(lapply(multiple_dfs, rownames)))
#     pooled_df =
#         # do.call(cbind, lapply(multiple_dfs, function(df) {df[rnames_all,]}))
#         bind_cols(lapply(multiple_dfs, function(df) {df[rnames_all,]}))
#     rownames(pooled_df) = rnames_all
#     
#     return(pooled_df)
# 
# }
# 
# # some minimal processing
# manual_scale_table = function(countTable) {
#     
#     readsPerWell      = apply(countTable,2,sum)
#     median_readsPerWell = median(readsPerWell)
#     GenesHowManyCells = apply(countTable>0,1,sum)
#     #gene_old_selection = sum(GenesHowManyCells>3)
#     
#     countTable_scaled = 
#         sapply(1:dim(countTable)[2], function(X) {
#             countTable[,X]/readsPerWell[X]*median_readsPerWell})
#     rownames(countTable_scaled) = rownames(countTable)
#     colnames(countTable_scaled) = colnames(countTable)
#     
#     # apply(countTable_scaled,2,sum)
#     return(list(countTable_scaled=countTable_scaled, readsPerWell=readsPerWell,
#             median_readsPerWell=median_readsPerWell, GenesHowManyCells=GenesHowManyCells ))
# }

