
venn_simple_plot_mw(list(
    TTN=gene_lists_customcorrelated_reorganized$`posCorrWith_ENSG00000155657:TTN`,
    XIRP2=gene_lists_customcorrelated_reorganized$`posCorrWith_ENSG00000163092:XIRP2`))

sum(gene_lists_customcorrelated_reorganized$`posCorrWith_ENSG00000155657:TTN` %in% gene_lists_customcorrelated_reorganized$`posCorrWith_ENSG00000163092:XIRP2`)
# sanity: sum(gene_lists_customcorrelated_reorganized$`posCorrWith_ENSG00000163092:XIRP2` %in% gene_lists_customcorrelated_reorganized$`posCorrWith_ENSG00000155657:TTN`)
length(gene_lists_customcorrelated_reorganized$`posCorrWith_ENSG00000163092:XIRP2`)
