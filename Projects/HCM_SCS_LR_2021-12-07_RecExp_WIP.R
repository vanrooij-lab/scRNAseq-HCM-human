


load(paste0(base_dir,'LR_analysis/LR_of_interest.Rdata'))
     # ligands_of_interest
     # receptors_of_interest

test_rec_list = receptors_of_interest$cl.1

feature_list = unlist(test_rec_list)
feature_list_annotation = rep(names(test_rec_list), times=sapply(test_rec_list, length))

# Gather averages per cell type, per patient
group.by1='annotation_patient_str'
group.by2='clusters_custom'
analysis='ROOIJonly_RID2l_clExtended'

shorthand_heatmap_feature_aggr(current_analysis=current_analysis, analysis='ROOIJonly_RID2l_clExtended', 
                               feature_list=feature_list, feature_list_annotation=feature_list_annotation, 
                               group.by1=group.by1, group.by2=group.by2)

shorthand_heatmap_feature_aggr = function(current_analysis, analysis, feature_list, feature_list_annotation, group.by1, group.by2) {
    
    # It seems Seurat does something weird when requesting using the
    # current_analysis[[analysis]][["clusters_custom"]]
    # syntax (returns df)
    # So this is workaround
    df_groups_meta = data.frame(group.by1=current_analysis[[analysis]][[group.by1]],group.by2=current_analysis[[analysis]][[group.by2]])
    
    # Convert to df
    df_features_meta =
        data.frame( feature_annot=feature_list_annotation, feature=feature_list )
    
    all_features = shorthand_seurat_fullgenename_faster(seuratObject = current_analysis[[analysis]], gene_names = df_features_meta$feature, return_NA = T)
    
    # Create summary expression matrix
    df_expr_list =
        lapply(all_features, function(g) {
            if (is.na(g)) {
                # Create 
                dummy = df_groups_meta
                dummy$expr = NA
                expr1_aggr = aggregate(expr1, df_groups_meta, mean)
                expr1_aggr$Z = NA
                return(expr1_aggr)
            } else {
                # collect data
                expr1 = as.vector(current_analysis[[analysis]]@assays$RNA@data[g, ])
                # mean(current_analysis[[analysis]]@assays$RNA@counts[g, ])
                
                # aggregate data
                expr1_aggr = aggregate(expr1, df_groups_meta, mean)
                
                # determine Z-score of aggr
                expr1_aggr$Z = scale(expr1_aggr$x)
                return(expr1_aggr)
            }
        })
    # df_expr = Reduce(rbind, df_expr_list)
    
    matrix_expr = 
        sapply(df_expr_list, function(df){df$Z})
    colnames(matrix_expr) = paste0(df_features_meta$feature, '_', df_features_meta$feature_annot)
    rownames(matrix_expr) = paste0(df_expr_list[[1]][[group.by1]], '_', df_expr_list[[1]][[group.by2]])
    
    annotation_cols = data.frame(annot = df_features_meta$feature_annot)
    rownames(annotation_cols) = paste0(df_features_meta$feature, '_', df_features_meta$feature_annot)
    
    annotation_rows = data.frame(g1=df_expr_list[[1]][[group.by1]], g2=df_expr_list[[1]][[group.by2]])
    rownames(annotation_rows) = paste0(df_expr_list[[1]][[group.by1]], '_', df_expr_list[[1]][[group.by2]])
    names(annotation_rows) = c(group.by1, group.by2)
    
    library(pheatmap)
    p=pheatmap(matrix_expr, annotation_row = annotation_rows, annotation_col = annotation_cols, cluster_cols = F, cluster_rows = F)
    return(p)

}

# heatmap(matrix_expr, Rowv = NA, Colv = NA)
# 
# ggplot(melt(matrix_expr), aes(x=Var1, y=Var2, fill=value))+
#     geom_tile()
# 
# library(ggplot2)
# ggplot(data.frame(x=c(1,2), y=c(4,5)), aes(x=x, y=y))+
#     geom_point()


