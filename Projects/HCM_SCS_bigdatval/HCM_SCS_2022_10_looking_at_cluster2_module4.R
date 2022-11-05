

########################################################################

# For convenience, on HPC, this script can be sourced by:
# script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')
#
# Local:
# LOCAL=1; script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Projects/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

########################################################################


# Let's retrieve list of cluster 2 and module 4 genes;
load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3b/Rdata/ROOIJonly.sp.bt_RID2l_core_regulons_sorted.Rdata') 
    # core_regulons_sorted
load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3b/Rdata/DE_cluster_ROOIJonly.sp.bt_RID2l.Rdata')
    # DE_cluster$ROOIJonly.sp.bt_RID2l$`2`
    genes_DE_cluster2 = DE_cluster$ROOIJonly.sp.bt_RID2l$`2`
load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3b/Rdata/ROOIJonly.sp.bt_RID2l__SCENIC_reg_top_genes_sorted_full.Rdata')
    
# Module 4, TFs or ligands
core_regulons_sorted$s.R.4
toString(shorthand_cutname( core_regulons_sorted$s.R.4[shorthand_cutname(core_regulons_sorted$s.R.4, PART1OR2 = 1) %in% human_TFs_ens] ))
toString(shorthand_cutname( core_regulons_sorted$s.R.4[shorthand_cutname(core_regulons_sorted$s.R.4, PART1OR2 = 2) %in% ligand_list_unique] ))

# Cluster 2
sum(genes_DE_cluster2$avg_log2FC>0)
genes_DE_cluster2$sym = shorthand_cutname( rownames(genes_DE_cluster2), PART1OR2 = 2)
genes_DE_cluster2$ens = shorthand_cutname( rownames(genes_DE_cluster2), PART1OR2 = 1)
genes_DE_cluster2_pos = genes_DE_cluster2[genes_DE_cluster2$avg_log2FC>0,]

# Which ones are TF or ligands?
toString(genes_DE_cluster2_pos$sym[genes_DE_cluster2_pos$sym %in% ligand_list_unique])
toString(genes_DE_cluster2_pos$sym[genes_DE_cluster2_pos$ens %in% human_TFs_ens])

# Overlap between regulon 4 and cluster 2
venn_simple_plot_mw(list(module4=core_regulons_sorted$s.R.4, cluster2= rownames(genes_DE_cluster2_pos)))

venn_simple_plot_mw(list(module4=core_regulons_sorted$s.R.4[1:10], cluster2= rownames(genes_DE_cluster2_pos)[1:10]))

# top 10 is completely in other set, though nog necessarily in top 10
all(core_regulons_sorted$s.R.4[1:10] %in% rownames(genes_DE_cluster2_pos))
all(rownames(genes_DE_cluster2_pos)[1:10] %in% core_regulons_sorted$s.R.4)

# joint top 10:
toString(shorthand_cutname( unique(c(core_regulons_sorted$s.R.4[1:10], rownames(genes_DE_cluster2_pos)[1:10])) ))

########################################################################
# Comparison with Alejandro's genes

genes_Alejandro_TAB = 
    c('XIRP2', 'XIRP1', 'NRAP', 'NPPB', 'NPPA', 'MYH7', 'LMOD2', 'HSPB7', 'DUSP27', 'DES', 'CSRP3', 'ANKRD1', 'SYNPO2L', 'SORBS2', 'FILIP1', 'CDH2')


genes_module4_sym = shorthand_cutname( core_regulons_sorted$s.R.4 )
genes_cl2_sym     = shorthand_cutname( rownames(genes_DE_cluster2_pos) )
genes_cl2mod4_sym = unique(c(genes_cl2_sym, genes_module4_sym))

venn_simple_plot_mw(list(TAB_16sel=genes_Alejandro_TAB, cl2mod4= unique(c(genes_cl2_sym, genes_module4_sym))  ))

toString(genes_Alejandro_TAB[genes_Alejandro_TAB %in% genes_cl2mod4_sym])
    # XIRP2, NRAP, MYH7, DES, ANKRD1, SYNPO2L



#######

shorthand_whichReg = function(GENE) {
        names(SCENIC_reg_top_genes_sorted_full)[
            sapply(SCENIC_reg_top_genes_sorted_full, function(X) {GENE %in% X})
        ]
    }

toString(shorthand_whichReg('RTN4'))
shorthand_whichReg('ANKRD1')
shorthand_whichReg('CTGF') # = CCN2

shorthand_whichReg('CMYA5')

TOPX=25
membership_mapping_mod4_SCENIC = 
    sapply(shorthand_cutname( core_regulons_sorted$s.R.4 )[1:TOPX], shorthand_whichReg)

module4_members_topX = shorthand_cutname( core_regulons_sorted$s.R.4 )[1:TOPX]
scenic_regs_of_interest = unique(unlist(membership_mapping_mod4_SCENIC))
membership_mod4_SCENIC_heatmap = 
    sapply(module4_members_topX, function(G) {
        
        mymemberlist=rep(NA, length(scenic_regs_of_interest))
        names(mymemberlist) = scenic_regs_of_interest
        mymemberlist[shorthand_whichReg(G)]=1
        return(mymemberlist)
        
    })
membership_mod4_SCENIC_heatmap[is.na(membership_mod4_SCENIC_heatmap)]=0

pheatmap(membership_mod4_SCENIC_heatmap, color = c('#FFFFFF','#000000'))
shorthand_whichReg('SVIL')

# Quickly check how CCN2 is known in my dataset.

if (F){
    
    DATASET_NAME='ROOIJonly.sp.bt_RID2l'
    current_analysis = list()
    current_analysis[[DATASET_NAME]] =
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))
    
    # http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000118523;r=6:131948176-131951372;t=ENST00000367976
    # CCN2 = ENSG00000118523
    shorthand_seurat_fullgenename(seuratObject = current_analysis[[DATASET_NAME]], gene_names = 'ENSG00000118523')
    rownames(current_analysis[[DATASET_NAME]])[grepl('ENSG00000118523',rownames(current_analysis[[DATASET_NAME]]))]

}










