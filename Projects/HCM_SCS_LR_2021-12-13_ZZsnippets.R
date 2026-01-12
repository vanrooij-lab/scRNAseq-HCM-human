

'B2M__' %in% 
    colnames(data_LR_expression$matrix_expr_Z_g2collapse[,max_R_Z_forEach>1.5&max_R_fraction_forEach>.4&colsNoNa3])
'MFGE8__' %in% 
    colnames(data_LR_expression$matrix_expr_Z_g2collapse[,max_R_Z_forEach>1.5&max_R_fraction_forEach>.4&colsNoNa3])
'CALR__' %in% 
    colnames(data_LR_expression$matrix_expr_Z_g2collapse[,max_R_Z_forEach>1.5&max_R_fraction_forEach>.4&colsNoNa3])


'HFE__' %in% 
    colnames(data_LR_expression$matrix_expr_Z_g2collapse[,max_R_Z_forEach>1.5&max_R_fraction_forEach>.4&colsNoNa3])
'TFRC__' %in% 
    colnames(data_LR_expression$matrix_expr_Z_g2collapse[,max_R_Z_forEach>1.5&max_R_fraction_forEach>.4&colsNoNa3])


pheatmap(data_LR_expression$matrix_expr_Z_g2collapse[,max_R_Z_forEach>1.5&max_R_fraction_forEach>.4&colsNoNa3], 
                               annotation_col = annot_col3[max_R_Z_forEach>1.5&max_R_fraction_forEach>.4&colsNoNa3,], 
                               fontsize = 7)

pheatmap(data_LR_expression$matrix_expr_Z_g2collapse[,c('TFRC__','HFE__')], 
                               annotation_col = annot_col3[c('TFRC__','HFE__'),], 
                               fontsize = 7)


pheatmap(data_LR_expression$matrix_expr_Z_g2collapse[,c('B2M__'),drop=F], 
                               annotation_col = annot_col3[c('B2M__'),,drop=F], 
                               fontsize = 7, cluster_cols = F, cluster_rows = F)
