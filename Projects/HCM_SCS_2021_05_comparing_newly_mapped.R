
library(ggplot2)


JE8_test_CountTable=read.table('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_rebuttal_analysis/new_count_tables/run1/JE8_AHFL7NBGX5_S17_cat_pT_total.TranscriptCounts.tsv')
View(JE8_test_CountTable[1:100,1:100])

# some minimal processing
JE8_readsPerWell      = apply(JE8_test_CountTable,2,sum)
JE8_median_readsPerWell = median(JE8_readsPerWell)
JE8_GenesHowManyCells = apply(JE8_test_CountTable>0,1,sum)
gene_selection = sum(JE8_GenesHowManyCells>3)

JE8_test_CountTable_scaled = 
    sapply(1:dim(JE8_test_CountTable)[2], function(X) {
        JE8_test_CountTable[,X]/JE8_readsPerWell[X]*JE8_median_readsPerWell})
rownames(JE8_test_CountTable_scaled) = rownames(JE8_test_CountTable)
colnames(JE8_test_CountTable_scaled) = colnames(JE8_test_CountTable)


#JE8_test_CountTable_scaled = scale(JE8_test_CountTable, center = F, scale = JE8_readsPerWell)


apply(JE8_test_CountTable_scaled,2,sum)


# TAKE_JE8_DATA_FROM_HERE: groupedData$patient2?

##########

# Old count data
JE8_old_CountTable=read.table('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_countdata/files_sent_by_anne_2020-08/Datafiles/JE8_TranscriptCounts.tsv', row.names = 1, header=1)

JE8_old_readsPerWell      = apply(JE8_old_CountTable,2,sum)
JE8_old_median_readsPerWell = median(JE8_old_readsPerWell)
JE8_old_GenesHowManyCells = apply(JE8_old_CountTable>0,1,sum)
gene_selection_old = sum(JE8_old_GenesHowManyCells>3)

JE8_old_CountTable_scaled = 
    sapply(1:dim(JE8_old_CountTable)[2], function(X) {
        JE8_old_CountTable[,X]/JE8_old_readsPerWell[X]*JE8_old_median_readsPerWell})
rownames(JE8_old_CountTable_scaled) = rownames(JE8_old_CountTable)
colnames(JE8_old_CountTable_scaled) = colnames(JE8_old_CountTable)

apply(JE8_old_CountTable_scaled,2,sum)

##########

# now get some data out to compare

# colnames are already compatible
colnames(JE8_test_CountTable_scaled)
colnames(JE8_old_CountTable_scaled)

# now make sure gene names are also compatible
# For new data
rownames(JE8_test_CountTable_scaled)
rownames(JE8_test_CountTable_scaled) = sapply(strsplit(rownames(JE8_test_CountTable_scaled),'_'), function(X){X[2]})
# For old data
rownames(JE8_old_CountTable_scaled)
rownames(JE8_old_CountTable_scaled) = sapply(strsplit(rownames(JE8_old_CountTable_scaled),'_'), function(X){X[1]})

# first test, TTN gene
cnames=colnames(JE8_old_CountTable_scaled)
ggplot(data.frame(old=JE8_old_CountTable_scaled['TTN',cnames], new=JE8_test_CountTable_scaled['TTN',cnames],
                totalReadsPerWell = JE8_readsPerWell[cnames]))+
    geom_point(aes(x=old, y=new, color=totalReadsPerWell>1000))+theme_bw()+ggtitle('TTN gene counts')+
    coord_fixed(ratio = 1)+geom_abline(intercept = 0, slope = 1)

GENE='MYH7'
ggplot(data.frame(old=JE8_old_CountTable_scaled[GENE,cnames], new=JE8_test_CountTable_scaled[GENE,cnames],
                totalReadsPerWell = JE8_readsPerWell[cnames]))+
    geom_point(aes(x=old, y=new, color=totalReadsPerWell>1000))+theme_bw()+ggtitle(paste0(GENE,' gene counts'))+
    coord_fixed(ratio = 1)+geom_abline(intercept = 0, slope = 1)


