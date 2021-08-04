


raceid2_directory = '/Users/m.wehrens/Documents/git_repos/SCS_Joep/'
source(paste0(raceid2_directory,"/Functions/RaceID2_StemID_class.R"))




# Let's do a quick RaceID2 analysis for comparison
current_analysis$ROOIJonly_RID2l_RACEID2 = SCseq(expdata = as.data.frame(current_analysis$ROOIJonly_RID2l@assays$RNA@counts))

# Normalization and downsampling
# note:
# expdata = raw data
# ndata = normalized data
# fdata = downsampled data
current_analysis$ROOIJonly_RID2l_RACEID2 =
    filterdata(
          current_analysis$ROOIJonly_RID2l_RACEID2, 
          mintotal = 1000, 
          minexpr = 3, 
          minnumber = 1, 
          maxexpr = Inf, 
          downsample = T, 
          dsn = 1, 
          rseed = 17000
    )

# Clustering based on correlation
# Note: with "pearson", distances are calculated by
# as.dist( 1 - cor(t(x), method='pearson'), to which 
# clustering is then applied
current_analysis$ROOIJonly_RID2l_RACEID2 =
    clustexp(
      current_analysis$ROOIJonly_RID2l_RACEID2,
      clustnr=30,
      bootnr=50,
      metric="pearson",
      do.gap=FALSE,
      sat=TRUE,
      SE.method="Tibs2001SEmax",
      SE.factor=.25,
      B.gap=50,
      cln=0,
      rseed=17000,
      FUNcluster="kmedoids"
    )

# tSNE
# Note: this uses the "tsne" function from the "tsne" package.
# (Which is a rather slow implementation.)
current_analysis$ROOIJonly_RID2l_RACEID2 <- comptsne(current_analysis$ROOIJonly_RID2l_RACEID2, rseed=15555)    
    
# Creates a 2nd reconsidered cluster grouping 
# (Stored in current_analysis$ROOIJonly_RID2l_RACEID2@cluster$kpart)
# (Original is in current_analysis$ROOIJonly_RID2l_RACEID2@cpart)
current_analysis$ROOIJonly_RID2l_RACEID2 <- findoutliers(
      current_analysis$ROOIJonly_RID2l_RACEID2, 
      outminc=5,
      outlg=4,
      probthr=1e-7,
      thr=2**-(1:40),
      outdistquant=.95
    )

# Visualze RaceID2 tsne
ggplot(data.frame(tsne1=current_analysis$ROOIJonly_RID2l_RACEID2@tsne[,1], tsne2=current_analysis$ROOIJonly_RID2l_RACEID2@tsne[,2],cl=as.factor(current_analysis$ROOIJonly_RID2l_RACEID2@cluster$kpart)),
        aes(x=tsne1, y=tsne2, color=cl))+
    geom_point()+theme_bw()

# How many cells / cluster
table(current_analysis$ROOIJonly_RID2l_RACEID2@cluster$kpart)

# Can we create a UMAP for consistency?
library(umap)
testUMAP = umap(current_analysis$ROOIJonly_RID2l_RACEID2@distances)

# Plot UMAP
ggplot(data.frame(umap1=testUMAP$layout[,1], umap2=testUMAP$layout[,2],cl=as.factor(current_analysis$ROOIJonly_RID2l_RACEID2@cluster$kpart)),
        aes(x=umap1, y=umap2, color=cl))+
    geom_point()+theme_bw()

View()

################################################################################
# Project the RaceID2 clustering on the Seurat UMAP

current_analysis$ROOIJonly_RID2l$Race_clusters=
    current_analysis$ROOIJonly_RID2l_RACEID2@cluster$kpart[colnames(current_analysis$ROOIJonly_RID2l)]

# Cl 3
current_analysis$ROOIJonly_RID2l$Race_cluster2=
    current_analysis$ROOIJonly_RID2l_RACEID2@cluster$kpart[colnames(current_analysis$ROOIJonly_RID2l)]==2


DimPlot(current_analysis$ROOIJonly_RID2l, group.by = 'Race_clusters')
DimPlot(current_analysis$ROOIJonly_RID2l, group.by = 'Race_cluster2')

Markers_Cl2 =
    FindMarkers(current_analysis$ROOIJonly_RID2l, ident.1 = colnames(current_analysis$ROOIJonly_RID2l)[current_analysis$ROOIJonly_RID2l$Race_cluster2])


current_analysis$ROOIJonly_RID2l = RunTSNE(object = current_analysis$ROOIJonly_RID2l, tsne.method = 'FIt-SNE') # tsne is not supported any more
current_analysis$ROOIJonly_RID2l = RunUMAP(object = current_analysis$ROOIJonly_RID2l, tsne.method = 'FIt-SNE') # tsne is not supported any more




