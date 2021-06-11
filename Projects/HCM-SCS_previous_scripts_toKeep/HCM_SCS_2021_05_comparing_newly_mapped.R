library(dplyr)
library(umap)
library(ggplot2)
library(Rtsne)

data_dir1='/Volumes/fastq_m.wehrens/Mapping/HCM_SCS/mapping.93.may25/counttables/'
data_dir2='/Volumes/fastq_m.wehrens/Mapping/WANG2/counttables/'
dataset_list_paths=list('AL1'=paste0(data_dir1, 'HUB-AL-s001_HG25TBGXF_S5_cat_pT_total.TranscriptCounts.tsv'),
                        'AL2'=paste0(data_dir1, 'HUB-AL-s002_HG25TBGXF_S6_cat_pT_total.TranscriptCounts.tsv'),
                        'WANG13'=paste0(data_dir2, 'GSM3449619_N13_cat_nc_total.TranscriptCounts.tsv'))
# unique and spliced reads only:
dataset_list_paths=list('AL1'=paste0(data_dir1, 'HUB-AL-s001_HG25TBGXF_S5_cat_pT_uniaggGenes_spliced.TranscriptCounts.tsv'),
                        'JE5'=paste0(data_dir1, 'JE5_AHFL77BGX5_S6_cat_pT_uniaggGenes_spliced.TranscriptCounts.tsv'),
                        'WANG13'=paste0(data_dir2, 'GSM3449619_N13_cat_nc_uniaggGenes_spliced.TranscriptCounts.tsv'),
                        'WANG14'=paste0(data_dir2, 'GSM3449620_N14_cat_nc_uniaggGenes_spliced.TranscriptCounts.tsv'))


# let's take those for now
groupedData_mw = 
    loadData_MW(dataset_list_paths)

lapply(groupedData_mw, dim)

rn1 = rownames(groupedData_mw$AL1)
rn2 = rownames(groupedData_mw$AL2)
rn3 = rownames(groupedData_mw$WANG13)

rn123 = unique(c(rn1, rn2, rn3))

test = 
    groupedData_mw$AL1[rn123,]
rownames(test) = rn12

pooled_df =
    do.call(cbind, lapply(groupedData_mw, function(df) {df[rn123,]}))

testdata = list(s1=data.frame(s1.c1=c(1,2,3),s1.c2=c(5,2,1), row.names = c('gene1','gene2','gene3')),
                s2=data.frame(s2.c1=c(1,2,2),s2.c2=c(4,2,4), row.names = c('gene2','gene3','gene4')),
                s3=data.frame(s3.c1=c(2,2,3),s3.c2=c(5,2,12), row.names = c('gene3','gene4','gene5')))

merge(testdata[[1]], testdata[[2]], all.x = T, )
merge(testdata[[1]], testdata[[2]], by.x=0, by.y=0, all=TRUE)

#####

test_pool=
    merge(groupedData_mw$AL1, groupedData_mw$AL2, by.x=0, by.y=0, all=TRUE)

test_pool2 = Reduce(
      function(x,y) {
        x <- merge(x,y,by.x=0,by.y=0,all=TRUE)
        rownames(x) <- x[,1]
        x[,"Row.names"] <- NULL
        x[is.na(x)] = 0
        return(x)
      },
      list(groupedData_mw$AL1, groupedData_mw$JE5, groupedData_mw$WANG13, groupedData_mw$WANG14)
    )

test_pool2_scaled_out = 
    manual_scale_table(test_pool2)

# let's look at potential to throw some stuff out
sum(test_pool2_scaled_out$GenesHowManyCells>4)
sum(test_pool2_scaled_out$GenesHowManyCells<5)

# create annotations for cells
annotations = rep(NA, dim(test_pool2_scaled_out$countTable_scaled)[2])
annotations[grepl('AL1.',colnames(test_pool2_scaled_out$countTable_scaled))] = 'AL1'
annotations[grepl('JE5.',colnames(test_pool2_scaled_out$countTable_scaled))] = 'JE5'
annotations[grepl('WANG13.',colnames(test_pool2_scaled_out$countTable_scaled))] = 'WANG13'
annotations[grepl('WANG14.',colnames(test_pool2_scaled_out$countTable_scaled))] = 'WANG14'

# let's plot some stats
ggplot(data.frame(counts=test_pool2_scaled_out$readsPerWell))+
  geom_freqpoly(aes(x=log10(.1+counts), color=annotations))+theme_bw()+
  geom_vline(xintercept = log10(3000+.1))
  
# make a selection
test_pool2_scaled_out$countTable_scaled_sel = 
  test_pool2_scaled_out$countTable_scaled[test_pool2_scaled_out$GenesHowManyCells>4,test_pool2_scaled_out$readsPerWell>3000]
annotation_sel = annotations[test_pool2_scaled_out$readsPerWell>3000]

# now create umap of the selection
umap_out = t(umap(t(test_pool2_scaled_out$countTable_scaled_sel)))
# plot
ggplot(data.frame(u1=umap_out[[1]][,1], u2=umap_out[[1]][,2], source=annotation_sel))+
    geom_point(aes(x=u1, y=u2, col=source))+theme_bw()

# also create a tsne 
Rtsne_out = t(Rtsne(t(test_pool2_scaled_out$countTable_scaled_sel), check_duplicates = F))
# plot
p1=ggplot(data.frame(tsne1=Rtsne_out[[2]][,1], tsne2=Rtsne_out[[2]][,2], source=annotation_sel))+
    geom_point(aes(x=tsne1, y=tsne2, col=source))+theme_bw()
print(p1)

library("RColorBrewer")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))


shorthand_expression_plot = function(which_gene) {
  expr=test_pool2_scaled_out$countTable_scaled_sel[which_gene,]
  
  p=ggplot(data.frame(tsne1=Rtsne_out[[2]][,1], tsne2=Rtsne_out[[2]][,2], source=annotation_sel, 
              expr=expr))+
    geom_point(aes(x=tsne1, y=tsne2, col=expr))+theme_bw()+
    scale_colour_gradientn(colours = myPalette(100), limits=c(0,max(expr)))+
    ggtitle(which_gene)
    
  return(p)

}

rownames(test_pool2_scaled_out$countTable_scaled_sel)[grepl('NPPB',rownames(test_pool2_scaled_out$countTable_scaled_sel))]

library(patchwork)
p1+shorthand_expression_plot('ENSG00000155657_TTN_ProteinCoding')
p1+shorthand_expression_plot('ENSG00000092054_MYH7_ProteinCoding')
p1+shorthand_expression_plot('ENSG00000175206_NPPA_ProteinCoding')
p1+shorthand_expression_plot('ENSG00000120937_NPPB_ProteinCoding')

####
# Also without mitochondrial genes
rownames(test_pool2)[grepl('_MT\\.',rownames(test_pool2))]
mito_remove_sel = !grepl('_MT\\.',rownames(test_pool2))
test_pool2_selM=test_pool2[mito_remove_sel,]

test_pool2_selM_scaled_out = 
    manual_scale_table(test_pool2_selM)

# make a selection again on the mito-removed table
test_pool2_selM_scaled_out$countTable_scaled_sel = 
  test_pool2_selM_scaled_out$countTable_scaled[test_pool2_selM_scaled_out$GenesHowManyCells>4,test_pool2_selM_scaled_out$readsPerWell>3000]
annotation_selM = annotations[test_pool2_selM_scaled_out$readsPerWell>3000]

# Let's also determine mitochondrial contents
mitochondrial_Percentage = apply(test_pool2[!mito_remove_sel,],2,sum)/apply(test_pool2, 2, sum)
  # note this is an estimate, since we're also taking into account thrown out genes for sum
  # though probably pretty accurate, because unlikely we threw out mito genes

Rtsne_out_selM = t(Rtsne(t(test_pool2_selM_scaled_out$countTable_scaled), check_duplicates = F))
p1m=ggplot(data.frame(tsne1=Rtsne_out_selM[[2]][,1], tsne2=Rtsne_out_selM[[2]][,2], source=annotation_sel))+
    geom_point(aes(x=tsne1, y=tsne2, col=source))+theme_bw()
print(p1m)

# now create umap of the selection
umap_out_selM = t(umap(t(test_pool2_selM_scaled_out$countTable_scaled)))
# plot
ggplot(data.frame(u1=umap_out_selM[[1]][,1], u2=umap_out_selM[[1]][,2], source=annotation_sel))+
    geom_point(aes(x=u1, y=u2, col=source))+theme_bw()

# plot of mitochondrial percentage
ggplot(data.frame(u1=umap_out_selM[[1]][,1], u2=umap_out_selM[[1]][,2], source=annotation_sel, mito_perc=mitochondrial_Percentage))+
  geom_point(aes(x=u1, y=u2, col=mito_perc))+theme_bw()+
  scale_colour_gradientn(colours = myPalette(100), limits=c(0,1))

ggplot(data.frame(tsne1=Rtsne_out_selM[[2]][,1], tsne2=Rtsne_out_selM[[2]][,2], source=annotation_sel, mito_perc=mitochondrial_Percentage))+
  geom_point(aes(x=tsne1, y=tsne2, col=mito_perc))+theme_bw()+
  scale_colour_gradientn(colours = myPalette(100), limits=c(0,1))



#####
rnames_all = unique(unlist(lapply(groupedData_mw, rownames)))
ldf=lapply(groupedData_mw, function(df) {df[rnames_all,]})
#####

#####



groupedData_mw_pooled = 
    pool_df_mw(groupedData_mw)

#####



