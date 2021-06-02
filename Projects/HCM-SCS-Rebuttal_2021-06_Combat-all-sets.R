
# OK, so Seurat didn't give the best results for this dataset in terms of correction, so 
# let's try the "Combat" package.

# See: 
# Zhang, Y., Parmigiani, G., & Johnson, W. E. (2020). ComBat-seq: batch effect adjustment for RNA-seq count data. NAR Genomics and Bioinformatics, 2(3), 1–10. https://doi.org/10.1093/nargab/lqaa078
# See also e.g.:
# Tutorials:
# https://rnabio.org/module-03-expression/0003/05/01/Batch-Correction/
# https://towardsdatascience.com/how-to-batch-correct-single-cell-7bad210c7ae1
# https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html
# Andrews, T. S., Kiselev, V. Y., McCarthy, D., & Hemberg, M. (2021). Tutorial: guidelines for the computational analysis of single-cell RNA sequencing data. Nature Protocols, 16(1). https://doi.org/10.1038/s41596-020-00409-w

################################################

source('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Functions_2021/Load-Pool-Scale_Simple_MW.R')

# checking script performance/progress
library(pryr)
library(beepr)

################################################

# Let's first load the data here again, now also including the selection of Teichmann cells
data_dir1='/Volumes/fastq_m.wehrens/Mapping/HCM_SCS/mapping.93.may25/counttables/'
data_dir2='/Volumes/fastq_m.wehrens/Mapping/WANG2/counttables/'
dataset_list_paths=list('AL1'=paste0(data_dir1, 'HUB-AL-s001_HG25TBGXF_S5_cat_pT_uniaggGenes_spliced.TranscriptCounts.tsv'),
                        'JE5'=paste0(data_dir1, 'JE5_AHFL77BGX5_S6_cat_pT_uniaggGenes_spliced.TranscriptCounts.tsv'),
                        'WANG13'=paste0(data_dir2, 'GSM3449619_N13_cat_nc_uniaggGenes_spliced.TranscriptCounts.tsv'),
                        'WANG14'=paste0(data_dir2, 'GSM3449620_N14_cat_nc_uniaggGenes_spliced.TranscriptCounts.tsv'))
# Load data
SCS_df_list_data = loadData_MW(dataset_list_paths, toPool = NULL)
# Teichmann is for now in saved workspace (small subselection)
load('/Volumes/workdrive_m.wehrens_hubrecht/data/2020_12_Litvinukova2020/HPC_sync/Teichmann_subset.Rdata')

# Note however, that gene names are now not compatible..
rownames(Teichmann_Sampled)
rownames(SCS_df_list_data$AL1)

# So we'll only take along uniquely mapped genes for now ..
# And we'll convert those names to only the symbol
preprocess_countdf = function(df) {
    # df = SCS_df_list_data$AL1
    # Take only uniquely identified genes
    the_duplicity=(sapply(rownames(df), str_count, pattern='-')+1)
    df=df[the_duplicity==1,]
    # Now take the symbol name
    symbol_names=sapply(strsplit(rownames(df),'_'), function(x) {x[[2]]})
    
    # For now we'll just sum the duplicates that remain
    # (Note: those will have e.g. a different ENSG ID, 
    # probably due to nomenclature issues.)
    # (Note2: presumably, their mapping also was different, so 
    # they will be different molecules, and their UMIs can be
    # added to each other.)
    # (Note3: the issue might also arise due to different mappings
    # in different samples, which is another argument for adding them)
    table_symbol_names_freq = table(symbol_names)
    symbol_names_plus_duplicity = data.frame(name=symbol_names, freq=as.vector(table_symbol_names_freq[symbol_names]))
    for (name_to_correct in unique(symbol_names_plus_duplicity$name[symbol_names_plus_duplicity$freq>1])) {
           summed_data = apply(df[symbol_names_plus_duplicity$name==name_to_correct,,drop=F],2,sum)
           df=df[symbol_names_plus_duplicity$name!=name_to_correct,]
           df[name_to_correct,]=summed_data
           print(paste0('Correct: ',name_to_correct)
    }
    
    
    
}



# Now combine the list
SCS_df_list_data$LITVI_samp = as.data.frame(Teichmann_Sampled)
colnames(SCS_df_list_data$LITVI_samp)=paste0('LITVI.',colnames(SCS_df_list_data$LITVI_samp))



# And now pool those (will take a while)
# This might be something that needs to be done @ the HPC, 
# since Rsession uses 40 GB of mem, operation takes 20 mins,
# SCS_df_pooled_data is 12 GB.
# Note: and this is only a small selection of data!
start_time <- Sys.time()
SCS_df_pooled_data = pool_df_mw(SCS_df_list_data)
end_time <- Sys.time()
print(end_time - start_time)
beepr::beep()
object_size(SCS_df_pooled_data)
object_size(SCS_df_list_data) # much smaller

# End of custom code
########################################################################

########################################################################
# Start of Seurat / Combat analysis

# Now convert this object to Seurat 

# Before starting, let's deduce some annotation
WRL_annotation_samples=unlist(lapply(  strsplit(colnames(SCS_df_pooled_data), split='\\.'), 
    function(splitted_string) {splitted_string[1]}))

# Initialize object
WRL_object = CreateSeuratObject(SCS_df_pooled_data, project = 'WRL')
object_size(WRL_object) 
    # clearly better mem management
    # perhaps eventually I should use Seurat merging to merge stuff?

# Add annotation
WRL_object[["annotation_sample"]]     <- as.factor(WRL_annotation_samples)

# Count mitochondrial reads
rownames(WRL_object)[grepl("-MT\\.",rownames(WRL_object))]
WR_object[["percent.mt"]] <- PercentageFeatureSet(WRL_object, pattern = "-MT\\.")

VlnPlot(WRL_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'annotation_sample')

# Count mito reads
# TODO: check whether this is workign
rownames(WRL_object)[grepl("-MT\\.",rownames(WRL_object))]
WRL_object[["percent.mt"]] <- PercentageFeatureSet(WRL_object, pattern = "-MT\\.")

# Remove some cells
# Note that percent.mt < 80 is not very stringent
WRL_object <- subset(WRL_object, subset = nFeature_RNA > 200 & nCount_RNA > 3000 & percent.mt < 80)
# Normalize data
WRL_object <- NormalizeData(WRL_object, normalization.method = "LogNormalize", scale.factor = 10000) # values are defaults
# Find variable features
WRL_object <- FindVariableFeatures(WRL_object, selection.method = "vst", nfeatures = 2000)

# ### some plots var features
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(WRL_object), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(WRL_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot2)

all.genes <- rownames(WRL_object)
WRL_object <- ScaleData(WRL_object, features = all.genes)

WRL_object <- RunPCA(WRL_object, npcs = 30, verbose = FALSE)
WRL_object <- RunUMAP(WRL_object, reduction = "pca", dims = 1:30)
WRL_object <- FindNeighbors(WRL_object, reduction = "pca", dims = 1:30)
WRL_object <- FindClusters(WRL_object, resolution = 0.5)


# Analysis done, show some plots

print(WRL_object[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(WRL_object, dims = 1:5, reduction = "pca", ncol = 5)

DimPlot(WRL_object, reduction = "pca")

JackStrawPlot(WRL_object, dims = 1:5) # doesn't work for some reason?
ElbowPlot(WRL_object)

DimPlot(WRL_object, reduction = "umap")

# Some plots relevant to our scientific questions

# First sample origin
DimPlot(WRL_object, group.by = 'annotation_sample', cols = col_vector_60,
            label = T, repel = T, label.size = 5)+theme(legend.position = 'none')
#DimPlot(WRL_object, group.by = 'annotation_sample')

# short_gene_name = 'MALAT1'
shorthand_expand_gene_name = function(the_scs_object, short_gene_name, showcounts=NULL) {
    gene_names = 
        rownames(the_scs_object)[grepl(paste0('-',short_gene_name,'-'), rownames(the_scs_object))]
    if(length(gene_names)>1) {
        print('Multiple genes found witht this identifier, unique chosen from:')
        print(paste0(gene_names,';\n'))
        the_duplicity=(sapply(gene_names, str_count, pattern='-')+1)/3
        if (!is.null(showcounts)) {
            p=ggplot(data.frame(expr=apply(the_scs_object@assays$RNA[gene_names,],1,sum),
                              gene=gene_names,
                              duplicity=the_duplicity,
                              duplicity_f=as.factor((sapply(gene_names, str_count, pattern='-')+1)/3),
                              duplicity_xlab=paste0(
                                  sprintf("%03d",       (sapply(gene_names, str_count, pattern='-')+1)/3     ),
                                  '_(',     sprintf("%03d",     1:length(gene_names)     )    ,')'
                                  )
                                ))+
                #geom_bar(aes(x=duplicity_xlab, y=expr), stat='identity')+
                geom_boxplot(aes(x=duplicity_f, y=expr))+
                geom_jitter(aes(x=duplicity_f, y=expr))+theme_bw()+
                    theme(axis.text.x = element_text(angle = 90))
            #+coord_flip()
            print(paste0('Unique one: ',gene_names[the_duplicity==1]))
            print(p)
        } else {
            return(gene_names[the_duplicity==1])    
        }
        return(NULL)
    } else {
        return(gene_names)   
    }
}


shorthand_expand_gene_name(WRL_object,'NPPB')
shorthand_expand_gene_name(WRL_object,'NPPA', showcounts=T)
shorthand_expand_gene_name(WRL_object,'TTN', showcounts=T)

FeaturePlot(WRL_object, features = c("ENSG00000251562-MALAT1-lincRNA"), cols = rainbow_colors)
FeaturePlot(WRL_object, features = c("ENSG00000155657-TTN-ProteinCoding"), cols = rainbow_colors)
FeaturePlot(WRL_object, features = c("ENSG00000092054-MYH7-ProteinCoding"), cols = rainbow_colors)
FeaturePlot(WRL_object, features = c("ENSG00000175206-NPPA-ProteinCoding"), cols = rainbow_colors)
FeaturePlot(WRL_object, features = c("ENSG00000120937-NPPB-ProteinCoding"), cols = rainbow_colors)

FeaturePlot(WRL_object, features = markers_ext$CM, cols = rainbow_colors)
FeaturePlot(WRL_object, features = markers_ext$FB, cols = rainbow_colors)
FeaturePlot(WRL_object, features = markers_ext$SMC, cols = rainbow_colors)
FeaturePlot(WRL_object, features = markers_ext$EC, cols = rainbow_colors)






