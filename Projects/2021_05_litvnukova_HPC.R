# Commands that I ran on the HPC

# After installation, I followed:
# - https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
# - https://mojaveazure.github.io/seurat-disk/reference/Convert.html
# Convert(source, dest, assay, overwrite = FALSE, verbose = TRUE, ...)

# Goal: load the h5ad file into R for further processing

################################################################################
# libs

library(Seurat)
library(SeuratDisk)
library(SeuratData)

library(ggplot2)
library(umap)
library(Rtsne)

library(RColorBrewer)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

################################################################################
# commands performed at the HPC

# First convert to h5seurat 
mysource='/hpc/hub_oudenaarden/mwehrens/data/Teichmann/hca_heart_ventricular_CM_raw.h5ad'
Convert(
  mysource,
  dest = "h5seurat",
  assay = "RNA",
  overwrite = FALSE,
  verbose = TRUE
)

# Then can be further converted
# First load it:
# https://mojaveazure.github.io/seurat-disk/articles/h5Seurat-load.html

# Can be viewed using Connect without loading into memory
Teichmann = SeuratDisk::Connect('/hpc/hub_oudenaarden/mwehrens/data/Teichmann/hca_heart_ventricular_CM_raw.h5seurat')
Teichmann
Teichmann$index()
Teichmann$close_all()

Teichmann <- LoadH5Seurat('/hpc/hub_oudenaarden/mwehrens/data/Teichmann/hca_heart_ventricular_CM_raw.h5seurat')
Teichmann

# > object_size(Teichmann)
# 5.02 GB

# To access the raw data:
Teichmann@assays$RNA[1:100,1:100]
any(grepl('septum',rownames(Teichmann@assays$RNA)))

# There appear to be 15K septal cells in this data
> sum(grepl('septum',colnames(Teichmann@assays$RNA)))
# [1] 15710

> dim(Teichmann@assays$RNA)
#[1]  33538 125289

# Let's see if we can create a somewhat smaller sample of the cells
# to play around with later
# sample(x=dim(Teichmann@assays$RNA)[2],size=4000)
some_sampled_colnames=colnames(Teichmann@assays$RNA)[sample(x=dim(Teichmann@assays$RNA)[2],size=4000)]
sum(grepl('septum',some_sampled_colnames))
# 532

# rsync -avh --progress gw2hpct01:/hpc/hub_oudenaarden/mwehrens/data/Teichmann/Teichmann_subset.Rdata /Volumes/workdrive_m.wehrens_hubrecht/data/2020_12_Litvinukova2020/HPC_sync

################################################################################
# commands that we perform lcoally

# Now load the thing here..
load('/Volumes/workdrive_m.wehrens_hubrecht/data/2020_12_Litvinukova2020/HPC_sync/Teichmann_subset.Rdata')
Teichmann_well_totals = apply(Teichmann_Sampled, 2, sum)
Teichmann_gene_total_counts = apply(Teichmann_Sampled, 1, sum)

Teichmann_Sampled_scaledMan_out = 
    manual_scale_table(Teichmann_Sampled)

ggplot(data.frame(read_totals=Teichmann_Sampled_scaledMan_out$readsPerWell))+
    geom_freqpoly(aes(x=read_totals))+theme_bw()+geom_vline(xintercept = 3000)

# actually quite some read counts are low, so let's see what remains after filtering
sum(Teichmann_Sampled_scaledMan_out$readsPerWell>3000)
sum(grepl('septum',colnames(Teichmann_Sampled_scaledMan_out$countTable_scaled[,Teichmann_Sampled_scaledMan_out$readsPerWell>3000])))
    # still 275 septal cells, OK

# Look at data

Teichmann_well_selection = Teichmann_Sampled_scaledMan_out$readsPerWell>3000
Teichmann_gene_selection = Teichmann_Sampled_scaledMan_out$GenesHowManyCells>4

Teichmann_Sampled_scaledMan_out$countTable_scaled_sel = 
    Teichmann_Sampled_scaledMan_out$countTable_scaled[Teichmann_gene_selection,Teichmann_well_selection]

# Create umap

Teichmann_Sampled_scaledMan_out$umap = 
    t(umap(t(Teichmann_Sampled_scaledMan_out$countTable_scaled_sel)))

annotation=unlist(lapply(  strsplit(colnames(Teichmann_Sampled_scaledMan_out$countTable_scaled_sel), split='_'), 
    function(splitted_string) {splitted_string[2]}))

ggplot(data.frame(u1=Teichmann_Sampled_scaledMan_out$umap[[1]][,1], u2=Teichmann_Sampled_scaledMan_out$umap[[1]][,2], source=annotation))+
    geom_point(aes(x=u1, y=u2, color=source))+theme_bw()

# create tsne

# also create a tsne 
Teichmann_Sampled_scaledMan_out$rtsne = t(Rtsne(t(Teichmann_Sampled_scaledMan_out$countTable_scaled_sel), check_duplicates = F))
# plot
p1=ggplot(data.frame(tsne1=Teichmann_Sampled_scaledMan_out$rtsne[[2]][,1], tsne2=Teichmann_Sampled_scaledMan_out$rtsne[[2]][,2], source=annotation))+
    geom_point(aes(x=tsne1, y=tsne2, col=source))+theme_bw()
print(p1)

#########

shorthand_expression_plot_Teich = function(which_gene=NULL, expr=NULL) {
    if (is.null(expr)) {    
        expr=Teichmann_Sampled_scaledMan_out$countTable_scaled_sel[which_gene,]
    } else {
     if (!is.null(which_gene)) { stop('Double input given') }   
    }
    
    
  p=ggplot(data.frame(tsne1=Teichmann_Sampled_scaledMan_out$rtsne[[2]][,1], tsne2=Teichmann_Sampled_scaledMan_out$rtsne[[2]][,2], source=annotation, 
              expr=expr))+
    geom_point(aes(x=tsne1, y=tsne2, col=expr))+theme_bw()+
    scale_colour_gradientn(colours = myPalette(100), limits=c(0,max(expr)))+
    ggtitle(which_gene)
    
  return(p)

}


rownames(Teichmann_Sampled_scaledMan_out$countTable_scaled_sel)[grepl('MALAT',rownames(Teichmann_Sampled_scaledMan_out$countTable_scaled_sel))]


shorthand_expression_plot_Teich('MALAT1')
shorthand_expression_plot_Teich('TTN')
shorthand_expression_plot_Teich('NPPA')
shorthand_expression_plot_Teich('NPPB')

# Looking at mitochondrial reads is a bit useless, since these are nuclei (does confirm we're looking at nuclei)
rownames(Teichmann_Sampled_scaledMan_out$countTable_scaled_sel)[grepl('MT-',rownames(Teichmann_Sampled_scaledMan_out$countTable_scaled_sel))]
expr_MT = apply(Teichmann_Sampled[grepl('MT-',rownames(Teichmann_Sampled)),], 2, sum)

################################################################################                
################################################################################
# Let's now not take a random selection, but instead take all vCM cells

# Well .. actually, the file I'm using (hca_heart_ventricular_CM_raw)
# of course only contains ventricular cardiomyocytes.
# 
# So, it's probably most convenient to later just add the other
# data to this data





