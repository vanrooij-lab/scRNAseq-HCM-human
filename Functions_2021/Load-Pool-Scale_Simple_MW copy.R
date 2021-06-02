
pool_df_mw = function(multiple_dfs) {
  # multiple_dfs should be a list of dataframes
    
    rnames_all = unique(unlist(lapply(multiple_dfs, rownames)))
    pooled_df = Reduce(
        function(x,y) {
          x <- merge(x,y,by.x=0,by.y=0,all=TRUE)
          rownames(x) <- x[,1]
          x[,"Row.names"] <- NULL
          x[is.na(x)] = 0
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

