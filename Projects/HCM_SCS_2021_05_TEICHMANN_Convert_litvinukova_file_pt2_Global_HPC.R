
# This is an edited copy of 
# HCM_SCS_2021_05_TEICHMANN_Convert_litvinukova_file_HPC.R
# 
# It deals with converting other files from the Litviňuková data set to 
# data files I can read.
#
# m.wehrens@hubrecht.eu

# At HPC, e.g.:
# srun --nodes=1 -c 1 --mem=50G --time=2:00:00 --pty bash -i
# R
# etc ..

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
mysource='/hpc/hub_oudenaarden/mwehrens/data/Teichmann/Litvinukova_global_raw.h5ad'
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
Teichmann_all = SeuratDisk::Connect('/hpc/hub_oudenaarden/mwehrens/data/Teichmann/Litvinukova_global_raw.h5seurat')
Teichmann_all
Teichmann_all$index()
Teichmann_all$close_all()

################################################################################
# OK we're done.

################################################################################
# We can test it:

Teichmann_all <- LoadH5Seurat('/hpc/hub_oudenaarden/mwehrens/data/Teichmann/Litvinukova_global_raw.h5seurat')
Teichmann_all

# > object_size(Teichmann_all)
# 5.02 GB




################################################################################
# The rest of the code is old code I used to play around with data;
# but not part of final analysis
################################################################################



# Playing around a bit more

# To access the raw data:
Teichmann_all@assays$RNA[1:100,1:100]
any(grepl('septum',rownames(Teichmann_all@assays$RNA)))
sum(grepl('septum',rownames(Teichmann_all@assays$RNA)))

# There appear to be 15K septal cells in this data
> sum(grepl('septum',colnames(Teichmann_all@assays$RNA)))
# [1] 15710

> dim(Teichmann_all@assays$RNA)
#[1]  33538 125289

# Let's see if we can create a somewhat smaller sample of the cells
# to play around with later
# sample(x=dim(Teichmann_all@assays$RNA)[2],size=4000)
some_sampled_colnames=colnames(Teichmann_all@assays$RNA)[sample(x=dim(Teichmann_all@assays$RNA)[2],size=4000)]
sum(grepl('septum',some_sampled_colnames))
# 532

# rsync -avh --progress gw2hpct01:/hpc/hub_oudenaarden/mwehrens/data/Teichmann/Teichmann_subset.Rdata /Volumes/workdrive_m.wehrens_hubrecht/data/2020_12_Litvinukova2020/HPC_sync




