

custom_list_GNAS_R = paste0(receptors_of_interest$HCM_up$GNAS,'__')
         
nonNA = !(apply(is.na(data_LR_expression$matrix_expr_Z_g2collapse[,custom_list_GNAS_R]),2,any))

 pheatmap(data_LR_expression$matrix_expr_Z_g2collapse[,custom_list_GNAS_R][,nonNA], 
                       annotation_col = annot_col3[custom_list_GNAS_R,][nonNA,], 
                       fontsize = 7)
         
  pheatmap(data_LR_expression$matrix_expr_Z_g2collapse[,custom_list_GNAS_R][,nonNA], 
                       annotation_col = annot_col3[custom_list_GNAS_R,][nonNA,], 
                       fontsize = 7)

# VTN
pheatmap(data_LR_expression$matrix_expr_Z_g2collapse[,'GNAS__',drop=F], 
                       annotation_col = annot_col3['GNAS__',,drop=F], 
                       fontsize = 7, cluster_rows = F, cluster_cols = F) 
    
# L VTN w/ R ITGA5  
  
custom_list_VTN_R = paste0(receptors_of_interest$HCM_up$VTN,'__')
toString(receptors_of_interest$HCM_up$VTN)      
  
nonNA = !(apply(is.na(data_LR_expression$matrix_expr_Z_g2collapse[,custom_list_VTN_R]),2,any))

 pheatmap(data_LR_expression$matrix_expr_Z_g2collapse[,custom_list_VTN_R][,nonNA], 
                       annotation_col = annot_col3[custom_list_VTN_R,][nonNA,], 
                       fontsize = 7)
           
# VTN
pheatmap(data_LR_expression$matrix_expr_Z_g2collapse[,'VTN__',drop=F], 
                       annotation_col = annot_col3['VTN__',,drop=F], 
                       fontsize = 7, cluster_rows = F, cluster_cols = F) 






# DDR2

L_for_DDR2 = paste0(quicky_give_L_for_R_OI('DDR2'),'__')

# L for DDR2
pheatmap(data_LR_expression$matrix_expr_Z_g2collapse[,L_for_DDR2,drop=F], 
                       annotation_col = annot_col3[L_for_DDR2,,drop=F], 
                       fontsize = 7, cluster_rows = F, cluster_cols = F) 

# L for DDR2
pheatmap(data_LR_expression$matrix_expr_Z_g2collapse[,'DDR2__',drop=F], 
                       annotation_col = annot_col3['DDR2__',,drop=F], 
                       fontsize = 7, cluster_rows = F, cluster_cols = F)



custom_list_DDR2_R = paste0(receptors_of_interest$HCM_up$COL1A1,'__')
toString(receptors_of_interest$HCM_up$)      
  
nonNA = !(apply(is.na(data_LR_expression$matrix_expr_Z_g2collapse[,custom_list_VTN_R]),2,any))

 pheatmap(data_LR_expression$matrix_expr_Z_g2collapse[,custom_list_VTN_R][,nonNA], 
                       annotation_col = annot_col3[custom_list_VTN_R,][nonNA,], 
                       fontsize = 7)
