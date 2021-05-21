
library(ggplot2)
library(umap)
library(Rtsne)

###############################################################

shorthand_loadRawData = function(filepath) {
    
    countTable=read.table(filepath, row.names = 1, header=1)
    #View(countTable[1:100,1:100])
    
    # some minimal processing
    readsPerWell      = apply(countTable,2,sum)
    median_readsPerWell = median(readsPerWell)
    GenesHowManyCells = apply(countTable>0,1,sum)
    gene_old_selection = sum(GenesHowManyCells>3)
    
    countTable_scaled = 
        sapply(1:dim(countTable)[2], function(X) {
            countTable[,X]/readsPerWell[X]*median_readsPerWell})
    rownames(countTable_scaled) = rownames(countTable)
    colnames(countTable_scaled) = colnames(countTable)
    
    apply(countTable_scaled,2,sum)
    
    #View(countTable_scaled)
    
    # now make sure gene names are also compatible
    # For old data
    # rownames(countTable_scaled)
    #rownames(countTable_scaled) = sapply(strsplit(rownames(countTable_scaled),'_'), function(X){X[1]})
    rownames(countTable_scaled) = sapply(strsplit(rownames(countTable_scaled),'_'), function(X){X[2]})
    #rownames(countTable) = sapply(strsplit(rownames(countTable),'_'), function(X){X[1]})
    
    return(list(countTable=countTable,countTable_scaled=countTable_scaled,readsPerWell=readsPerWell,GenesHowManyCells=GenesHowManyCells))
}

###############################################################
###############################################################
# Comparing JE HCM SCS old/new
###############################################################
###############################################################

# Load new data
# all
JE7_new_data = shorthand_loadRawData('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_rebuttal_analysis/new_count_tables/run2-GRCh38.81/final_tsv/JE7_AHFL7NBGX5_S16_cat_pT_total.TranscriptCounts.tsv')

# load old data 
JE7_old_data = shorthand_loadRawData('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_countdata/files_sent_by_anne_2020-08/Datafiles/JE7_TranscriptCounts.tsv')

JE7_old_CountTable=read.table(oldfilepath, row.names = 1, header=1)
View(JE7_old_CountTable[1:100,1:100])

# some minimal processing
JE7_old_readsPerWell      = apply(JE7_old_CountTable,2,sum)
JE7_old_median_readsPerWell = median(JE7_old_readsPerWell)
JE7_old_GenesHowManyCells = apply(JE7_old_CountTable>0,1,sum)
gene_old_selection = sum(JE7_old_GenesHowManyCells>3)

JE7_old_CountTable_scaled = 
    sapply(1:dim(JE7_old_CountTable)[2], function(X) {
        JE7_old_CountTable[,X]/JE7_old_readsPerWell[X]*JE7_old_median_readsPerWell})
rownames(JE7_old_CountTable_scaled) = rownames(JE7_old_CountTable)
colnames(JE7_old_CountTable_scaled) = colnames(JE7_old_CountTable)

apply(JE7_old_CountTable_scaled,2,sum)

View(JE7_old_CountTable_scaled)

# now make sure gene names are also compatible
# For old data
# rownames(JE7_old_CountTable_scaled)
rownames(JE7_old_CountTable_scaled) = sapply(strsplit(rownames(JE7_old_CountTable_scaled),'_'), function(X){X[1]})




################################################################################

# now get some data out to compare

# simple stats
sum(JE7_old_CountTable)/sum(JE7_new_CountTable)

# first test, TTN gene
cnames= intersect(colnames(JE7_old_CountTable_scaled), colnames(JE7_new_CountTable_scaled))
ggplot(data.frame(old=JE7_old_CountTable_scaled['TTN',cnames], new=JE7_new_CountTable_scaled['TTN',cnames],
                totalReadsPerWell = JE7_readsPerWell[cnames]))+
    geom_point(aes(x=old, y=new, color=totalReadsPerWell>1000))+theme_bw()+ggtitle('TTN gene counts')+
    coord_fixed(ratio = 1)+geom_abline(intercept = 0, slope = 1)

GENE='MYH7'
ggplot(data.frame(old=JE7_old_CountTable_scaled[GENE,cnames], new=JE7_new_CountTable_scaled[GENE,cnames],
                totalReadsPerWell = JE7_readsPerWell[cnames]))+
    geom_point(aes(x=old, y=new, color=totalReadsPerWell>1000))+theme_bw()+ggtitle(paste0(GENE,' gene counts'))+
    coord_fixed(ratio = 1)+geom_abline(intercept = 0, slope = 1)


################################################################################

# combine old and new

# find shared gene names
geneSeenOld = apply(JE7_old_CountTable_scaled>0,1,sum)
geneSeenNew = apply(JE7_new_CountTable_scaled>0,1,sum)
shared_gene_names = intersect(rownames(JE7_old_CountTable_scaled)[geneSeenOld>2], rownames(JE7_new_CountTable_scaled)[geneSeenNew>2])

newcolnames = c(paste0('o.',cnames), paste0('n.',cnames))
combined_df = 
    cbind(JE7_old_CountTable_scaled[shared_gene_names,cnames], JE7_new_CountTable_scaled[shared_gene_names,cnames])
colnames(combined_df) = newcolnames

library(umap)
library(Rtsne)
#library(tsne)

umap_out = t(umap(t(combined_df)))
Rtsne_out = t(Rtsne(t(combined_df), check_duplicates = F))

source=c( rep('old', length(cnames)), rep('new', length(cnames)) )
pairs=rep(cnames,2)

ggplot(data.frame(u1=umap_out[[1]][,1], u2=umap_out[[1]][,2], source=annotation))+
    geom_point(aes(x=u1, y=u2, col=source))


plot_df = data.frame(u1=umap_out[[1]][,1], u2=umap_out[[1]][,2], 
                    tsne1=Rtsne_out[[2]][,1], tsne2=Rtsne_out[[2]][,2],
                     pairs=as.factor(pairs), source=as.factor(source),
                     totalReads=rep(JE7_old_readsPerWell[cnames],2))
ggplot(plot_df[plot_df$totalReads>1000,])+
    geom_point(aes(x=u1, y=u2, color=pairs))+
    geom_line(aes(x=u1, y=u2, color=pairs))+theme_bw()+theme(legend.position='none')
ggplot(plot_df[plot_df$totalReads>1000,])+
    geom_line(aes(x=tsne1, y=tsne2, color=pairs))+
    geom_point(aes(x=tsne1, y=tsne2, color=pairs))+
    theme_bw()+theme(legend.position='none')


#ggplot(data.frame(x=c(10,2,3,4,100,200),y=c(3,7,5,1.3,100,100), cat=as.factor(c(1,1,2,2,3,3))))+
#    geom_line(aes(x=x, y=y, color=cat))+theme_bw()


plot_df$MALAT1 = c(JE7_old_CountTable_scaled['MALAT1',cnames], JE7_new_CountTable_scaled['MALAT1',cnames])
ggplot(plot_df)+
    geom_point(aes(x=tsne1, y=tsne2, color=MALAT1))
    #geom_line(aes(x=tsne1, y=tsne2, color=pairs))+theme_bw()+theme(legend.position='none')


################################################################################

# Let's compare X160, a cell which has many reads

df_compare_HCM_singleCell = data.frame(combined_df[,c('n.X160', 'o.X160')])

ggplot(df_compare_HCM_singleCell, aes(x=n.X160, y=o.X160))+
    geom_point()+ggtitle('Comparing single cell')+theme_bw()

###############################################################
###############################################################
# Comparing Wang old/new
###############################################################


###############################################################

# Cells we've mapped again:
# - SC_92563_0_12  SRR6641031

CURRENT_CELL_NAME = 'SC_92563_0_12'

# Their mapped and processed data
Wang_Data1 = read.table('/Volumes/workdrive_m.wehrens_hubrecht/data/2020_04_Wang-heart/GSE109816_normal_heart_umi_matrix.csv',row.names = 1, header=1, sep=',')

data_WangTest =
    shorthand_loadRawData(
    #        '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_rebuttal_analysis/wang_count_tables/test_1SRA/SRR6640430_nc_total.UFICounts.tsv')
        # only single mappers
             '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_rebuttal_analysis/wang_count_tables/test_1SRA/SRR6640920_nc_uniaggGenes_total.UFICounts.tsv')
             #'/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_rebuttal_analysis/wang_count_tables/test_1SRA/SRR6641031_nc_total.UFICounts.tsv')

#View(data_WangTest$countTable_scaled)
View(data_WangTest$countTable)
sum(data_WangTest$countTable)

# Compare totals
sum(Wang_Data1[,CURRENT_CELL_NAME,drop=F])
sum(data_WangTest$countTable)
# SC_92563_0_12 SRR6641031: 52103 & 60124
# SC_92563_0_23

# previous stuff
#data_WangTest =
#    shorthand_loadRawData(
#        #'/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_rebuttal_analysis/wang_count_tables/test_1SRA/SRR6640920_nc_total.UFICounts.tsv')
#        '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_rebuttal_analysis/wang_count_tables/test_1SRA/SRR6640430_nc_uniaggGenes_total.UFICounts.tsv')

#data_WangTest =
#    shorthand_loadRawData(
#        # unstranded
#        #'/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_rebuttal_analysis/wang_count_tables/test_1SRA/mapping_unstranded/SRR6640920_nc_total.UFICounts.tsv'
#        # assuming stranded
#        '/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_rebuttal_analysis/wang_count_tables/test_1SRA/SRR6640920_nc_total.UFICounts.tsv')

# Compare with original processed table (see other script for loading)
Current_Cell_Original = Wang_Data1[,CURRENT_CELL_NAME,drop=F]
# Make rownames consistent with Anna
rownames(Current_Cell_Original) = gsub(pattern = '-', replacement = '.', x = rownames(Current_Cell_Original))

# Initialize new cell
Current_Cell_New = data_WangTest$countTable

# determine initial new names, but these will have doubles
simpler_GeneNames_New =
    sapply(strsplit(rownames(Current_Cell_New),'_'), function(X){X[2]})
simpler_GeneNames_New_Cnts = table(simpler_GeneNames_New)
double_collection = names(simpler_GeneNames_New_Cnts[simpler_GeneNames_New_Cnts>1])    

# calculate totals for the doubles
thecolname = colnames(Current_Cell_New)
Current_Cell_New_Update = data.frame(double())
colnames(Current_Cell_New_Update) = colnames(Current_Cell_New)
for (double_name in double_collection) {
    
    # calculating sum and just merging them
    current_totalCount = data.frame(
        row.names=double_name,
        sum(Current_Cell_New[simpler_GeneNames_New==double_name,]))
    colnames(current_totalCount) = thecolname
    Current_Cell_New_Update = rbind(Current_Cell_New_Update, current_totalCount)
    
# renaming   
#    simpler_GeneNames_New[simpler_GeneNames_New==double_name] = 
#        paste0(simpler_GeneNames_New[simpler_GeneNames_New==double_name],'__nr',1:length(simpler_GeneNames_New[simpler_GeneNames_New==double_name]))
}

# remove doubules
Current_Cell_New=Current_Cell_New[!(simpler_GeneNames_New %in% double_collection),,drop=F]
# now update rownames
rownames(Current_Cell_New) = sapply(strsplit(rownames(Current_Cell_New),'_'), function(X){X[2]})

# add merged entries
Current_Cell_New=rbind(Current_Cell_New,Current_Cell_New_Update)

#rownames(Current_Cell_New) = simpler_GeneNames_New
View(Current_Cell_New)

# by intersecting genes
matching_genes = 
    intersect(rownames(Current_Cell_New), rownames(Current_Cell_Original))

df_toPlotWangComp = data.frame(n=1:length(matching_genes), gene=matching_genes, ct.wang.byme=Current_Cell_New[matching_genes,1],ct.wang.bywang=Current_Cell_Original[matching_genes,1])
df_toPlotWangComp$mismatch = df_toPlotWangComp$ct.wang.byme-df_toPlotWangComp$ct.wang.bywang

# just merging any positive count
newNames = rownames(Current_Cell_New[Current_Cell_New[,1]>0,,drop=F])
oriNames = rownames(Current_Cell_Original[Current_Cell_Original[,1]>0,,drop=F])
sel_genes = unique(c(newNames,oriNames))
df_toPlotWangComp_all = data.frame(n=1:length(sel_genes), gene=sel_genes, ct.wang.byme=Current_Cell_New[sel_genes,1],ct.wang.bywang=Current_Cell_Original[sel_genes,1])
df_toPlotWangComp_all[is.na(df_toPlotWangComp_all)] = 0
df_toPlotWangComp_all$mismatch = df_toPlotWangComp_all$ct.wang.byme-df_toPlotWangComp_all$ct.wang.bywang
df_toPlotWangComp_all$max=sapply(1:dim(df_toPlotWangComp_all)[1], function(X) {max(df_toPlotWangComp_all$ct.wang.byme[X], df_toPlotWangComp_all$ct.wang.bywang[X])})
df_toPlotWangComp_all$mismatch_rel = (df_toPlotWangComp_all$ct.wang.byme-df_toPlotWangComp_all$ct.wang.bywang)/df_toPlotWangComp_all$max

View(df_toPlotWangComp_all)

ggplot(df_toPlotWangComp, aes(x=ct.wang.byme,y=ct.wang.bywang))+
    geom_point()

# for intersect
library(ggplot2)
library(ggrepel)
ggplot(df_toPlotWangComp, aes(x=ct.wang.byme,y=ct.wang.bywang))+
    geom_point()+theme_bw()+
    geom_point(data=df_toPlotWangComp[df_toPlotWangComp$mismatch>20,], color='red')+
    geom_text_repel(data=df_toPlotWangComp[df_toPlotWangComp$mismatch>20,], aes(label=gene), color='red')+
    geom_text_repel(data=df_toPlotWangComp[df_toPlotWangComp$ct.new>60&!df_toPlotWangComp$mismatch>20,], aes(label=gene), color='blue')+
    geom_abline(slope = 1, intercept = 0)

# for all
library(ggplot2)
library(ggrepel)
weird_genes = df_toPlotWangComp_all$max>100&abs(df_toPlotWangComp_all$mismatch_rel)>.95
ggplot(df_toPlotWangComp_all, aes(x=ct.wang.byme+1,y=ct.wang.bywang+1))+
    geom_point()+theme_bw()+
    geom_point(data=df_toPlotWangComp_all[weird_genes,], color='purple')+
    #geom_text_repel(data=df_toPlotWangComp_all[df_toPlotWangComp$mismatch>20,], aes(label=gene), color='red')+
    geom_text_repel(data=df_toPlotWangComp_all[df_toPlotWangComp_all$max>100 & !weird_genes,], aes(label=gene), color='blue')+
    geom_text_repel(data=df_toPlotWangComp_all[weird_genes,], aes(label=gene,x=ct.wang.byme+1,y=ct.wang.bywang+1), color='purple')+
    geom_abline(slope = 1, intercept = 0)+
    scale_x_log10()+scale_y_log10() 


# conclusion:
# for some genes it has worked, while for other genes it has not..

ggplot(df_toPlotWangComp, aes(x=n, y=mismatch))+
    geom_point()+theme_bw()


View(Current_Cell_New)

# amount of genes that have more than 10 reads
sum(df_toPlotWangComp_all$max>20)
# calculate how many genes that have >10 reads have a >.95 rel. mismatch 
sum(df_toPlotWangComp_all$max>20&abs(df_toPlotWangComp_all$mismatch_rel)>.95)

# Buuuuut .. the ratio of total reads between old and new is not so far off ..
sum(Current_Cell_New)/sum(Current_Cell_Original)

