
# most of this should be already implemented elsewhere



RidgePlot(object = pbmc_small, features = 'PC_1')


marker = 'NPPA'
pRidge = RidgePlot(object = mySeuratObject, features = marker, group.by = 'annotation_paper_str') #, group.by = 'from_paper')
ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_4_Ridge_markers_',marker,'.png'), plot = pRidge, height=15, width=15, units='cm')

expr=
mySeuratObject@assays$RNA@scale.data[marker,]
paper=
mySeuratObject$annotation_paper_fct

pFreq=ggplot(data.frame(expression=mySeuratObject@assays$RNA@scale.data[marker,], 
                        source=mySeuratObject$annotation_paper_fct))+
    geom_freqpoly(aes(x=expression, color=source, after_stat(density)))
ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_4_freqpoly_markers_',marker,'.png'), plot = pFreq, height=15, width=15, units='cm')

# Also create histograms with 2% limits
currentlims=list()
for (source in levels(mySeuratObject$annotation_paper_fct)) {
    #currentlims[[source]]=calc_limits(mySeuratObject@assays$RNA@scale.data[marker,mySeuratObject$annotation_paper_fct==source], percentile = .02)
    pHist=ggplot(data.frame(expression=mySeuratObject@assays$RNA@scale.data[marker,mySeuratObject$annotation_paper_fct==source]),
                    aes(ymax=max(..count..)))+
        #geom_histogram(aes(x=expression))+#, after_stat(density)))+
        stat_bin(aes(x=expression))+
        give_better_textsize_plot(10)
        #xlim(c(0,currentlims[[source]][2]))
    ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_4_histogram98Lim_markers_',marker,'_',source,'.png'), plot = pHist, height=10, width=7.5, units='cm')
}

# More complicated histogram, split, x-lims @98% of data, such that shape of curve is visible in case of outliers
currentlims=list()
for (source in levels(mySeuratObject$annotation_paper_fct)) {
    currentlims[[source]]=calc_limits(mySeuratObject@assays$RNA@scale.data[marker,mySeuratObject$annotation_paper_fct==source], percentile = .02)
}
# then calculate breaks
max_val = max(sapply(currentlims, function(x) {x[2]}))
currentbreaks = seq(from=-(max_val/29)/2,by=(max_val/29),to=max_val+(max_val/29)/2) 
# And make histogram, use density as y-value
pHist=ggplot(data.frame(expression=mySeuratObject@assays$RNA@scale.data[marker,],
                   source=mySeuratObject$annotation_paper_fct),
             )+
    geom_histogram(aes(x=expression, y=..density.., fill=source), breaks=currentbreaks)+
    give_better_textsize_plot(10)+
    facet_grid(rows='source')+theme_bw()+theme(legend.position = 'none')
ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_4_histogram98Lim_markers_',marker,'.png'), plot = pHist, height=10, width=7.5, units='cm')


# Custom Violin that takes 
pHist=ggplot(data.frame(expression=mySeuratObject@assays$RNA@scale.data[marker,], 
                        source=mySeuratObject$annotation_paper_fct))+
    geom_histogram(aes(x=expression, fill=source, after_stat(density)))+
    facet_grid(rows='source')+theme_bw()+theme(legend.position = 'none', )+give_better_textsize_plot(10)
ggsave(filename = paste0(base_dir,'Rplots/',run_name,'_4_histogram_markers_',marker,'.png'), plot = pHist, height=10, width=7.5, units='cm')





