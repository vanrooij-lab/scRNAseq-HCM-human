

# For reference, this was what we got from the previous analysis:
load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/Previous_analysis_for_reference/JoepAnalysis_Regulons.Rdata')
# regulons_README_objects_saved



# Now compare

regulon_object$the_regulons
View(regulon_object$regulon_list_extended)


View(regulons$patient1Mod$regulons)
    #regulons_README_objects_saved
    


# Comparing two

old_regulons = regulons$patient1Mod$regulons
old_regulons=lapply(old_regulons, strip__chrXX)
names(old_regulons) = paste0('o.R.',1:length(old_regulons))
new_regulons = regulon_object$the_regulons
new_regulons=lapply(new_regulons, strip__chrXX)
names(new_regulons) = paste0('n.R.',1:length(new_regulons))
    
pooled_regulons = c(old_regulons, new_regulons)

#df_compare = data.frame(t(combn(names(pooled_regulons), 2)))
df_compare = tidyr::expand_grid(x=names(pooled_regulons), y=names(pooled_regulons))

df_compare$overlap = sapply(1:dim(df_compare)[1], function(X) { 
    sum(pooled_regulons[[df_compare$x[X]]] %in% pooled_regulons[[df_compare$y[X]]]) /
        min(length(pooled_regulons[[df_compare$x[X]]]), length(pooled_regulons[[df_compare$y[X]]]))
    })
View(df_compare)

ggplot(df_compare)+
    geom_tile(aes(x=x, y=y, fill=overlap))

matrix_compare <- reshape2::acast(df_compare, x~y, value.var="overlap")

pheatmap(matrix_compare, cluster_rows = F, cluster_cols = F)
pheatmap(matrix_compare, clustering_method = 'ward.D2')


#########
# Now compare from same source

# regulon_object_MW_groupedSCS_patient1mod=regulon_object




old_regulons = regulons$patient1Mod$regulons
old_regulons=lapply(old_regulons, strip__chrXX)
names(old_regulons) = paste0('o.R.',1:length(old_regulons))
new_regulons = regulon_object$the_regulons
names(new_regulons) = paste0('n.R.',1:length(new_regulons))
    
pooled_regulons = c(old_regulons, new_regulons)

#df_compare = data.frame(t(combn(names(pooled_regulons), 2)))
df_compare = tidyr::expand_grid(x=names(pooled_regulons), y=names(pooled_regulons))

df_compare$overlap = sapply(1:dim(df_compare)[1], function(X) { 
    sum(pooled_regulons[[df_compare$x[X]]] %in% pooled_regulons[[df_compare$y[X]]]) /
        min(length(pooled_regulons[[df_compare$x[X]]]), length(pooled_regulons[[df_compare$y[X]]]))
    })
View(df_compare)

ggplot(df_compare)+
    geom_tile(aes(x=x, y=y, fill=overlap))

matrix_compare <- reshape2::acast(df_compare, x~y, value.var="overlap")

pheatmap(matrix_compare, cluster_rows = F, cluster_cols = F)
pheatmap(matrix_compare)

######################################################################

# The two data sources:
# groupedSCS_patient1mod@ndata
# current_analysis[['ROOIJonly_RID2l']]@assays$RNA@data

# selection of genes doesn't really matter though, since most important ones are there
selJ = rowSums(groupedSCS_patient1mod@ndata>.1)>dim(groupedSCS_patient1mod@ndata)[2]*.05
selS = rowSums(current_matrix>0)>dim(groupedSCS_patient1mod@ndata)[2]*.05

gene_names_Joep = strip__chrXX(rownames(groupedSCS_patient1mod@ndata[selJ,]))
gene_names_Seurat = rownames(current_matrix[selS,])

venn_simple_plot_mw(venn_list = list(Joep=gene_names_Joep, Seurat=gene_names_Seurat))
    # slightly different 

cell_names_Joep = colnames(groupedSCS_patient1mod@ndata)
cell_names_Seurat = colnames(current_matrix)

# Have to convert between wellnames to compare .. 
# Not sure I want to ..

# is it easy to translate this?
# I guess so actually
wellplate=data.frame(colnr = rep(1:24, times=16), 
                     rowsym=rep(c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'), each=24),
                     wellcount=1:384)
# now convert
cell_names_Joep_inSeuratstyle =
    sapply(cell_names_Joep, function(x) {
        s=strsplit(x,split = '_')[[1]]; 
        paste0(switch(s[1],JE1='JE01',JE2='JE02',JE3='JE03',JE4='JE04'),'_',
              sprintf("X%03d", wellplate[wellplate$colnr==s[3]&wellplate$rowsym==s[2],]$wellcount))
        })

# just double-check conversion
ggplot(wellplate)+
    geom_text(aes(x=colnr,y=rowsym, label=wellcount))+theme_minimal()+scale_y_discrete(limits=rev)#scale_y_reverse()

venn_simple_plot_mw(list(Seurat=cell_names_Seurat, Joep=cell_names_Joep_inSeuratstyle))












