
################################################################################

# Data from 
# Desjardins, C. A., & Naya, F. J. (2017). Antagonistic regulation of cell-cycle and differentiation gene programs in neonatal cardiomyocytes by homologous MEF2 transcription factors. Journal of Biological Chemistry, 292(25), 10613–10629. https://doi.org/10.1074/jbc.M117.776153
#
# Note: 	log2-transformed, RMA-normalized Entrez Gene expression values
#
# See also https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2438843

library(GEOquery) # BiocManager::install("GEOquery")
# library(limma)
# library(umap)

################################################################################

# load series and platform data from GEO

if (F) {
    gset <- getGEO("GSE92861", GSEMatrix =TRUE, getGPL=FALSE)
    if (length(gset) > 1) idx <- grep("GPL17799", attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    
    ex <- exprs(gset)
}

# Sample names
colname_conversion = gset@phenoData@data$title
names(colname_conversion) = gset@phenoData@data$geo_accession
# 
cols_annot = colname_conversion[colnames(ex)]
        
save(list=c('gset','ex'),file=paste0(base_dir, 'Rdata/data_MEF_siRNA_GSE92861.Rdata')) # gset, ex

################################################################################

source('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Functions/Conversion_mouse_human_genes.R')

if (F) {
    
    mart_rat <- useEnsembl(biomart = "ensembl", 
                       dataset = "rnorvegicus_gene_ensembl")#, 
                       #version = paste0(ENSEMBL_VERSION_BIOTYPES))
    
    conversion_table_hcgn_entrez_rat <- getBM(
      attributes=c("ensembl_gene_id", 'entrezgene_id','rgd_symbol'), # "hgnc_symbol",
      mart = mart_rat 
      )
    
    head(conversion_table_hcgn_entrez_rat)
    
    save(list='conversion_table_hcgn_entrez_rat',file=paste0(base_dir, 'Rdata/conversion_table_hcgn_entrez_rat.Rdata')) # conversion_table_hcgn_entrez_rat
}

load(file=paste0(base_dir, 'Rdata/conversion_table_hcgn_entrez_rat.Rdata')) # conversion_table_hcgn_entrez_rat

conversion_table_hcgn_entrez_rat_ENTRtoSYMB = conversion_table_hcgn_entrez_rat$rgd_symbol
names(conversion_table_hcgn_entrez_rat_ENTRtoSYMB) = conversion_table_hcgn_entrez_rat$entrezgene_id

################################################################################
# Convert names to human names

gene_names_orig = gsub('_at$','',rownames(ex))
gene_names_symb = conversion_table_hcgn_entrez_rat_ENTRtoSYMB[gene_names_orig]
    
if (F) {
    
    #View(conversion_table_hcgn_entrez_rat)
    #conversion_table_hcgn_entrez_rat_ENTRtoSYMB
    
    rathuman_conversion_table = convertMouseGeneList(gene_names_symb, the_other_animal = 'rat', table_dump = T)
    
    rathuman_conversion_lookup = rathuman_conversion_table$HGNC.symbol
    names(rathuman_conversion_lookup) = rathuman_conversion_table$RGD.symbol
    
    save(list='rathuman_conversion_lookup',file=paste0(base_dir, 'Rdata/rathuman_conversion_lookup.Rdata')) # rathuman_conversion_lookup
}
# gene_names_symb_humanized = convertMouseGeneList(gene_names_symb, the_other_animal = 'rat', make_unique = F)

load(file=paste0(base_dir, 'Rdata/rathuman_conversion_lookup.Rdata')) # rathuman_conversion_lookup

gene_names_symb_humanized = rathuman_conversion_lookup[gene_names_symb]



################################################################################
# Humanize the expression table

ex_humonly = ex[!is.na(gene_names_symb_humanized),]
rownames(ex_humonly) = gene_names_symb_humanized[!is.na(gene_names_symb_humanized)]

View(ex_humonly)

# check whether log2fc
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
LogC # = F --> log already done

# normalization seems to have been done
colSums(ex_humonly)
colSums(ex)

# So to create log2fc values, need inside vs. outside sets
cntrl_idxs = grepl('LacZ',cols_annot)
mef2a_idxs = grepl('Mef2A',cols_annot)
mef2c_idxs = grepl('Mef2C',cols_annot)
mef2d_idxs = grepl('Mef2D',cols_annot)

log2fc_mef2a = rowMeans(ex_humonly[,mef2a_idxs]) - rowMeans(ex_humonly[,cntrl_idxs])
log2fc_mef2c = rowMeans(ex_humonly[,mef2c_idxs]) - rowMeans(ex_humonly[,cntrl_idxs])
log2fc_mef2d = rowMeans(ex_humonly[,mef2d_idxs]) - rowMeans(ex_humonly[,cntrl_idxs])

log2fc_table = data.frame(log2fc_mef2a, log2fc_mef2c, log2fc_mef2d, gene=rownames(ex_humonly))

################################################################################

top100_mef2a = log2fc_table[order(log2fc_table$log2fc_mef2a, decreasing = T),][1:100,]$gene
top100_mef2c = log2fc_table[order(log2fc_table$log2fc_mef2c, decreasing = T),][1:100,]$gene
top100_mef2d = log2fc_table[order(log2fc_table$log2fc_mef2d, decreasing = T),][1:100,]$gene

lapply(core_regulons_sorted_shortname, function(X) {sum(X %in% top100_mef2a)/length(X)})
lapply(core_regulons_sorted_shortname, function(X) {sum(X %in% top100_mef2c)/length(X)})
lapply(core_regulons_sorted_shortname, function(X) {sum(X %in% top100_mef2d)/length(X)})


for (reg_name in names(core_regulons_sorted_shortname)){

    the_length=length(core_regulons_sorted_shortname[[reg_name]])
    module2_matrix_mef_ = log2fc_table[log2fc_table$gene %in% core_regulons_sorted_shortname[[reg_name]][1:min(the_length, 20)],]
    name_table =  table(  module2_matrix_mef_[,4]  )
    non_duplicates = names(name_table)[ name_table == 1]
    module2_matrix_mef = module2_matrix_mef_[module2_matrix_mef_[,4] %in% non_duplicates ,1:3]
    rownames(module2_matrix_mef) = module2_matrix_mef_[module2_matrix_mef_[,4] %in% non_duplicates , 4]
    
    
    p1=pheatmap((module2_matrix_mef<log2(.8))*1, main = reg_name, cluster_cols = F, cluster_rows = F)
    
    nrow_effective=nrow(module2_matrix_mef)
    
    ggsave(filename = paste0(base_dir,'Rplots/','Desjardins_customModules_style1_',reg_name,'.pdf'), 
                             plot = p1, height=min(172,3*nrow_effective+40), width=1/3*172-4, units = 'mm', device = cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5
    
    
    breaksList = seq(-3, 3, by = .1)

    # Plots the first heatmap
    p2=pheatmap(1*(module2_matrix_mef), # Plots the first 10 genes of the dataset
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
             breaks = breaksList, 
             main = reg_name, cluster_cols = F, cluster_rows = F)

    ggsave(filename = paste0(base_dir,'Rplots/','Desjardins_customModules_style2_',reg_name,'.pdf'), 
                         plot = p2, height=min(172,3*nrow_effective+40), width=1/3*172-4, units = 'mm', device = cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5

    
}

# SCENIC, top20
for (reg_name in names(SCENIC_reg_top_genes_sorted_full)){

    module2_matrix_mef_ = log2fc_table[log2fc_table$gene %in% SCENIC_reg_top_genes_sorted_full[[reg_name]][1:20],]
    name_table =  table(  module2_matrix_mef_[,4]  )
    non_duplicates = names(name_table)[ name_table == 1]
    module2_matrix_mef = module2_matrix_mef_[module2_matrix_mef_[,4] %in% non_duplicates ,1:3]
    rownames(module2_matrix_mef) = module2_matrix_mef_[module2_matrix_mef_[,4] %in% non_duplicates , 4]
    
    
    p1=pheatmap((module2_matrix_mef<log2(.8))*1, main = reg_name, cluster_cols = F, cluster_rows = F)
    
    nrow_effective=nrow(module2_matrix_mef)
    
    ggsave(filename = paste0(base_dir,'Rplots/','Desjardins_style1_',reg_name,'.pdf'), 
                             plot = p1, height=min(172,3*nrow_effective+40), width=1/3*172-4, units = 'mm', device = cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5
            
    
    #if(F){ 
    breaksList = seq(-3, 3, by = .1)

    # Plots the first heatmap
    p2=pheatmap(1*(module2_matrix_mef), # Plots the first 10 genes of the dataset
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
             breaks = breaksList, 
             main = reg_name, cluster_cols = F, cluster_rows = F)
    #}
    ggsave(filename = paste0(base_dir,'Rplots/','Desjardins_style2_',reg_name,'.pdf'), 
                             plot = p2, height=min(172,3*nrow_effective+40), width=1/3*172-4, units = 'mm', device = cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5
    
    
}







