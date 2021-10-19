# Transform the tables 

library(ggplot2)

dataTable1=
    read.table('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_UmiTools_CountTables_Wang_HCMSCS/ROOIJ/JE2_AHY3WGBGX3_S1_cat_pT.nonRibo_E99_Aligned.out.counts.tsv', header = 1)
dataTable1=dataTable1[order(dataTable1$cell),]
# View(dataTable1)

# just to check all cells are here:
unique(dataTable1$cell)

#library(tidyr)
# dataTable1_df = 
#     pivot_wider(dataTable1, names_from = cell, values_from = count)
# dataTable1_df[is.na(dataTable1_df)]=0

library(reshape2)
dataTable1_df=
    acast(dataTable1, gene~cell, value.var="count")
dataTable1_df[is.na(dataTable1_df)]=0

# Check out some stats
gene_totals = rowSums(dataTable1_df)
cell_totals = colSums(dataTable1_df)
# View(gene_totals)

ggplot(data.frame(celltotals=colSums(dataTable1_df)))+
    geom_histogram(aes(x=log10(.1+celltotals)))+theme_bw()


load('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Resources/Anna_name_conversion.Rdata')

################################################################################

# Obtain gene symbols with ensemble IDs
library(biomaRt)
        
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
all_gene_symbols_ensids <- getBM(
  attributes=c("hgnc_symbol","ensembl_gene_id"),
  mart = mart)
# head(all_gene_symbols_ensids)

# Now conversion table
ens_to_sym_conv_table = all_gene_symbols_ensids$hgnc_symbol
names(ens_to_sym_conv_table) = all_gene_symbols_ensids$ensembl_gene_id

rownames(dataTable1_df)
ens_to_sym_conv_table[rownames(dataTable1_df)]
# show how many gene symbols are redundant
length(unique(ens_to_sym_conv_table[unique(dataTable1$gene)]))/length(unique(dataTable1$gene))
    
rownames(dataTable1_df)=paste0(rownames(dataTable1_df),'_',ens_to_sym_conv_table[rownames(dataTable1_df)])
colnames(dataTable1_df)=paste0('X',sprintf("%03d", as.numeric(colnames(dataTable1_df))))
dataTable1_df=as.data.frame(dataTable1_df)
################################################################################

# Now inspect relationship between mito and pseudogenes over cells

ggplot(data.frame(MT.ATP6=dataTable1_df['ENSG00000198899_MT-ATP6',],
           MTATP6P1=dataTable1_df['ENSG00000248527_MTATP6P1',]),
        aes(x=MT.ATP6,y=MTATP6P1))+
    geom_point()+theme_bw()+geom_abline(slope = 1, intercept = 0)

ggplot(data.frame(MT.ND1=dataTable1_df['ENSG00000198888_MT-ND1',],
                MTND1P23=dataTable1_df['ENSG00000225972_MTND1P23',]),
        aes(x=MT.ND1,y=MTND1P23))+
    geom_point()+theme_bw()+geom_abline(slope = 1, intercept = 0)

#####

dataTable1_anna_multi = 
    #read.table('/Volumes/fastq_m.wehrens/Mapping/HCM_SCS/mapping.93.may25/counttables/JE2_AHY3WGBGX3_S1_cat_nc_total.UFICounts.tsv', header=1, row.names = 1)
    read.table('/Volumes/fastq_m.wehrens/Mapping/HCM_SCS/mapping.93.may25/counttables/JE2_AHY3WGBGX3_S1_cat_pT_uniaggGenes_spliced.UFICounts.tsv', header=1, row.names = 1)
dataTable1_anna = dataTable1_anna_multi[!grepl('-',rownames(dataTable1_anna_multi)),]

# fix names
library(stringr)
rownames(dataTable1_anna)=sapply(str_split(rownames(dataTable1_anna),pattern ='_'), function(x) {paste0(x[[1]],'_',x[[2]])} )
rownames(dataTable1_anna)[grepl('ENSG00000002016_RAD52',rownames(dataTable1_anna))]
rownames(dataTable1_anna)[grepl('ENSG00000006607_FARP2',rownames(dataTable1_anna))]
rowSums(dataTable1_anna[grepl('ENSG00000002016_RAD52',rownames(dataTable1_anna)),])
dataTable1_anna[grepl('ENSG00000006607_FARP2',rownames(dataTable1_anna)),]


# make pseudogene comparison as before
rownames(dataTable1_anna)[grepl('_MT.ND1',rownames(dataTable1_anna))]
ggplot(data.frame(MT.ND1  =as.double(dataTable1_anna['ENSG00000198888_MT.ND1',]),
                     MTND1P23=as.double(dataTable1_anna['ENSG00000225972_MTND1P23',])),
        aes(x=MT.ND1,y=MTND1P23))+
    geom_point()+theme_bw()+geom_abline(slope = 1, intercept = 0)

ggplot(data.frame(MT.ATP6  =as.double(dataTable1_anna['ENSG00000198899_MT-ATP6',]),
                     MTATP6P1=as.double(dataTable1_anna['ENSG00000248527_MTATP6P1',])),
        aes(x=MT.ATP6,y=MTATP6P1))+
    geom_point()+theme_bw()+geom_abline(slope = 1, intercept = 0)
    # Note that this probably doesn't work because these genes were merged ..


# Now compare some cells
gene_totals_anna=rowSums(dataTable1_anna)
gene_totals_umit=rowSums(dataTable1_df)
names(gene_totals_umit)=gsub(pattern = '-', replacement = '\\.', x = names(gene_totals_umit))
gene_names_both = unique(c(names(gene_totals_anna), names(gene_totals_umit)))
compare_df = data.frame(expr_anna = gene_totals_anna[gene_names_both], 
                        expr_umit = gene_totals_umit[gene_names_both],
                        gene=gene_names_both)
compare_df[is.na(compare_df)]=-1
ggplot(compare_df,
            aes(x=log10(2+expr_anna),y=log10(2+expr_umit)))+
    geom_point()+theme_bw()+geom_abline(slope = 1, intercept = 0)
    #geom_vline(xintercept = log10(2+4^6))+geom_hline(yintercept = log10(2+4^6)) # doesn't make sense for rowSums

compare_df$distance=(log10(2+compare_df$expr_anna)-log10(2+compare_df$expr_umit))^2
compare_df$short_name=sapply(str_split(compare_df$gene,pattern ='_'), function(x) {paste0(x[[2]])} )

library(ggrepel)
ggplot(compare_df,
            aes(x=log10(2+expr_anna),y=log10(2+expr_umit)))+
    geom_point()+theme_bw()+geom_abline(slope = 1, intercept = 0)+
    # geom_vline(xintercept = log10(2+4^6))+geom_hline(yintercept = log10(2+4^6))+ # doesn't make sense for rowSums
    geom_point(data=compare_df[compare_df$expr_anna>0 & compare_df$expr_umit>0 &  compare_df$distance>1,], color='red')+
    geom_text_repel(data=compare_df[compare_df$expr_anna>0 & compare_df$expr_umit>0 &  compare_df$distance>1,], aes(label=short_name), color='red')+
    geom_text_repel(data=compare_df[(log10(2+compare_df$expr_anna)>2.5 | log10(2+compare_df$expr_umit)>2.5) & (compare_df$expr_anna<=0 | compare_df$expr_umit<=0),], aes(label=short_name), color='purple')
    #geom_text_repel(data=compare_df[compare_df$distance>3,], aes(label=gene), color='red')

ggplot(compare_df,
            aes(x=log10(2+expr_anna),y=log10(2+expr_umit)))+
    geom_point()+theme_bw()+geom_abline(slope = 1, intercept = 0)+
    # geom_vline(xintercept = log10(2+4^6))+geom_hline(yintercept = log10(2+4^6))+ # doesn't make sense for rowSums
    geom_point(data=compare_df[compare_df$expr_anna>0 & compare_df$expr_umit>0 &  compare_df$distance>1,], color='red')+
    geom_text_repel(data=compare_df[compare_df$expr_anna>0 & compare_df$expr_umit>0 &  compare_df$distance>.5,], aes(label=short_name), color='red', size=2)+
    geom_text_repel(data=compare_df[(log10(2+compare_df$expr_anna)>2.5 | log10(2+compare_df$expr_umit)>2.5) & (compare_df$expr_anna<=0 | compare_df$expr_umit<=0),], aes(label=short_name), color='purple')
    #geom_text_repel(data=compare_df[compare_df$distance>3,], aes(label=gene), color='red')






