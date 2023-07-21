


df_ICD_genes_estigoy = 
    read.csv('/Users/m.wehrens/Documents/Naslag/articles_supp_mics/Estigoy2009/12551_2008_7_MOESM1_ESM_rowswcategorielabelsremoved.csv')

ICD_gene_list =
    df_ICD_genes_estigoy$PROTEIN..ExPASy..Gene.


sapply(ICD_gene_list, function(X) {str_match(X, "\\(\\s*(.*?)\\s*\\)")[,2]})



pattern <- "\\(\\s*(.*?)\\s*\\)"
myresult_=
    sapply(ICD_gene_list, function(X) {
        str_brackets = regmatches(X, regexec(pattern, X))[[1]][2]
        gene_name = str_split(string = str_brackets, pattern = ', ')
        gene_name[[1]][2]}
        )
names(myresult_) = NULL
myresult=myresult_[!is.na(myresult_)]

ICD_gene_list = myresult


genes_module4[genes_module4 %in% ICD_gene_list]


'MAFK' %in% genes_module4
which('MAFK' == genes_module4)
