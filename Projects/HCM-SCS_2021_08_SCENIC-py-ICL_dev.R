
# Gather some info about the reliability of the link between gene sets and their TF
add.info.meta = read.csv('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/SCENIC/R.P1/reg.csv', head=F, nrows=3)
add.info = read.csv('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/SCENIC/R.P1/reg.csv', head=F, skip = 3)
colnames(add.info)=c(add.info.meta[3, 1:2], add.info.meta[2, 3:length(add.info)])
    # View(add.info)

# Now summarize this to have NES scores available per regulon
myRegScores = aggregate(list(NES=add.info$NES), by = list(TF=add.info$TF), FUN=median)

# Now also get gene weights
adj.info = read.csv('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/SCENIC/R.P1/adj.csv', head=1)
adj.info[adj.info$TF=='MEF2A',]

# Importance: higher = more important
regulon_gene_importance = 
    lapply(unique(adj.info$TF), function(tf) {adj.info[adj.info$TF==tf,c('target', 'importance')]})
names(regulon_gene_importance) = unique(adj.info$TF)


adj_YY1=adj.info[adj.info$TF=='YY1',]
adj_YY1[order(adj_YY1$importance, decreasing = T),][1:10,]

adj_MEF2A=adj.info[adj.info$TF=='MEF2A',]
adj_MEF2A[order(adj_MEF2A$importance, decreasing = T),][1:10,]

adj_ZEB1=adj.info[adj.info$TF=='ZEB1',]
adj_ZEB1[order(adj_ZEB1$importance, decreasing = T),][1:10,]



############################################################

# 
View(add.info[add.info$TF=='ZEB1',])

ggplot(data.frame(NES=add.info[add.info$TF=='MEF2A',]$NES))+
    geom_histogram(aes(x=NES))
