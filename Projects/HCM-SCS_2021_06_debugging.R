

# Debugging why mito count failed earlier. Turned out Teichmann didn't have nCount_RNA and nFeature_RNA fields.
sel_teichmann = RHL_SeuratObject_merged$annotation_paper_str=='Teichmann'
mito_genes_idx =
    which(grepl('^MT-',rownames(RHL_SeuratObject_merged_TEICHMANNonly@assays$RNA@counts)))

test = colSums(x = GetAssayData(object = RHL_SeuratObject_merged,
        assay = 'RNA', slot = "counts")[mito_genes_idx, , drop = FALSE])
test[sel_teichmann][1:10]
test2 = RHL_SeuratObject_merged$nCount_RNA * 100
test2[sel_teichmann][1:10]

RHL_SeuratObject_merged_TEICHMANNonly[["percent.mt.test"]] <- PercentageFeatureSet(RHL_SeuratObject_merged_TEICHMANNonly, pattern = "^MT-")
RHL_SeuratObject_merged_TEICHMANNonly$percent.mt.test

RHL_SeuratObject_merged[["percent.mt.test"]] <- PercentageFeatureSet(RHL_SeuratObject_merged, pattern = "^MT-")
RHL_SeuratObject_merged$percent.mt.test[RHL_SeuratObject_merged$annotation_paper_str=='Teichmann'][1:100]


RHL_SeuratObject_merged$percent.mt[RHL_SeuratObject_merged$annotation_paper_str=='Teichmann'][1:100]
