
################################################################################

# Code to be executed at the HPC
# This looks into gene expression patterns of 
# genes that were identified in the tomo-seq project.

# Note to self:
# inter-active session can be started @ HPC using
# srun --nodes=1 -c 1 --mem=64G --time=2:00:00 --pty bash -i

# one-liner that loads the main HCM SCS script w/ relevant functions/libraries
script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

################################################################################

# Genes of interest are from:
    # intersecting_genes_pc1_pc2_perc25 
    # Determined in ~/Documents/git_repos/Tomoseq_mw/projects/project_Anne/__ANNE_analysis__PCA_and-more.R
intersecting_genes_pc1_pc2_perc25_HARDCODED = 
    c('ACTC1', 'ACTN2', 'ADI1', 'ADSSL1', 'AK3', 'ALDOA', 'ANK1', 'ANKRD1', 'ANKRD2', 'APLN', 'ART3', 'ATP5A1', 'ATP5C1', 'ATP5EP2', 'ATP5F1', 
      'ATP5J2', 'ATP5J', 'ATP5L', 'ATP5O', 'ATP5SL', 'BANF1', 'BBS9', 'BCKDK', 'BZW2', 'C11orf84', 'C14orf159', 'C14orf2', 'C1orf43', 'C3orf43', 
      'C7orf73', 'C9orf123', 'C9orf3', 'CACNA1C', 'CADPS', 'CAMK2B', 'CAMK2D', 'CAMTA1', 'CCNB1IP1', 'CDK18', 'CISD1', 'CISD3', 'CKMT2', 'CMYA5', 
      'COL21A1', 'CORO6', 'COX17', 'COX4I1', 'COX7B', 'COX8A', 'CRYAB', 'CRYM', 'CTNND2', 'CYB5R1', 'DDA1', 'DES', 'DLAT', 'DNAJB5', 'DNAJC4', 
      'DPH5', 'DUSP26', 'ECH1', 'ENAM', 'ENO3', 'EXOG', 'FAM189A2', 'FHIT', 'FHL2', 'FKBP3', 'FRMD3', 'GBAS', 'GJA4', 'GLRX2', 'GNG5', 'GYG1', 
      'HIST2H2BE', 'HMGN3', 'HRC', 'HSPB1', 'HSPB3', 'HYAL1', 'IDH2', 'ISCU', 'KCNH2', 'KLHL41', 'LDB3', 'LOC100131138', 'LOC100507537', 'LRRC39', 
      'LTBP1', 'LYPLAL1', 'MAP4', 'MDH1', 'MGST3', 'MLLT11', 'MPC2', 'MRPL15', 'MRPL33', 'MRPL48', 'MYH7B', 'MYL12A', 'MYL2', 'MYL3', 'MYOM1', 
      'MYOT', 'MYOZ2', 'MYPN', 'NCAM1', 'NDUFA12', 'NDUFA4', 'NDUFA6', 'NDUFB1', 'NDUFB2', 'NDUFB3', 'NDUFB9', 'NDUFC1', 'NDUFS7', 'NDUFV2', 
      'NKX2-5', 'NMRK2', 'NPY6R', 'NUDT19', 'NUDT7', 'PACSIN3', 'PCGF5', 'PFKM', 'PHB', 'PKIA', 'PKP2', 'POPDC3', 'PPP1R12B', 'PPP1R3C', 'PPP2CA', 
      'PPP2R1A', 'PPP2R4', 'PRDX2', 'PSMA5', 'PTS', 'RAB3A', 'RBM24', 'RNF115', 'RPL38', 'RYR2', 'SAMD4A', 'SAP18', 'SDHB', 'SDHD', 'SLC25A4', 
      'SLC25A5', 'SLC2A4', 'SLC39A14', 'SLCO3A1', 'SLIRP', 'SMYD2', 'SNTA1', 'SORBS2', 'SPEG', 'SRL', 'SSBP2', 'ST3GAL6', 'SUPT5H', 'SVIL',
      'SYNPO2L', 'TBX5', 'TECRL', 'THBS4', 'TIMM21', 'TMEM67', 'TNNC1', 'TNNT2', 'TOM1L2', 'TPM1', 'TRIM54', 'TTC1', 'TXLNB', 'TYRP1', 'UQCR10',
      'UQCR11', 'UQCRQ', 'USP13', 'UTP11L', 'VDAC1', 'VDAC2', 'VDAC3', 'VPS37A', 'XIRP1', 'XIRP2', 'YAF2', 'ZNF853')

# The 55 targets:
    # COUPTFII_genes
    # Determined in ~/Documents/git_repos/Tomoseq_mw/projects/project_Anne/__ANNE_analysis__HOMER-target-genes.R
COUPTFII_genes_HARDCODED = 
    c('CCNB1IP1', 'LDB3', 'MYPN', 'HSPB1', 'XIRP1', 'SPEG', 'HRC', 'SUPT5H', 'PRDX2', 'SRL', 'ACTC1', 'YAF2', 'ANKRD2', 'CISD1', 'DES', 'TRIM54', 
      'CRYAB', 'LYPLAL1', 'APLN', 'FKBP3', 'VDAC2', 'UQCRQ', 'HYAL1', 'MYL3', 'FHL2', 'LTBP1', 'GLRX2', 'SLC25A4', 'TECRL', 'SLC25A5', 'COX7B', 
      'NDUFS7', 'MYL12A', 'CISD3', 'ENO3', 'IDH2', 'ADSSL1', 'SLIRP', 'TBX5', 'ISCU', 'NDUFA12', 'PTS', 'COX8A', 'NDUFB9', 'BBS9', 'NPY6R', 'PPP2CA',
      'NDUFC1', 'CADPS', 'PPP1R12B', 'HIST2H2BE', 'RNF115', 'PSMA5', 'DPH5', 'SDHB')

# And the 19 genes determined by Kees:
COUPTFII_genes_overlapOtherWork_HARDCODED = 
    c('ACTC1', 'MYL12A', 'LTBP1', 'PSMA5', 'LYPLAL1', 'GLRX2', 'IDH2', 'FHL2', 'DPH5', 'PPP1R12B',  'ANKRD2', 'TECRL', 'NDUFS7', 'HSPB1', 'LDB3', 'MYL3', 'XIRP1', 'CRYAB', 'DES')
    # Determined by Kees


Tomoseq_genes_of_interest = 
    list(genespc1pc2_188=intersecting_genes_pc1_pc2_perc25_HARDCODED,
         genespc1pc2COUPTFII_55=COUPTFII_genes_HARDCODED,
         genespc1pc2COUPTFIIoverlap_19=COUPTFII_genes_overlapOtherWork_HARDCODED)

################################################################################
# Load the SCS HCM data

# Copied piece that sets up some relevant stuff
# --- load / set up data

# CURRENT_RUNNAME='all_RID2l_VAR'
# CURRENT_RUNNAME='ALL.SP_RID2l_clExtended'
CURRENT_RUNNAME='ALL.SP_btypSel_RID2l_clExtended'

current_analysis = list()
current_analysis[[CURRENT_RUNNAME]] = 
    # LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_merged_noMito_sel.',CURRENT_RUNNAME,'.h5seurat'))
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',CURRENT_RUNNAME,'.h5seurat'))

# Add one-letter shorthand for dataset
name_change = c("vRooij"="R", "Hu"="H", "Teichmann"="T")
current_analysis[[CURRENT_RUNNAME]]$annotation_paper_oneletter_fct = factor(
    name_change[current_analysis[[CURRENT_RUNNAME]]$annotation_paper_str], levels=c('R','H','T'))
# Custom labeling
name_change = c("vRooij"="HCM", "Hu"="Ctrl1", "Teichmann"="Ctrl2")
current_analysis[[CURRENT_RUNNAME]]$annotation_paper_beatified = factor(
    name_change[current_analysis[[CURRENT_RUNNAME]]$annotation_paper_str], levels=c('HCM','Ctrl1','Ctrl2'))
# Custom labeling also for patient
current_analysis[[CURRENT_RUNNAME]]$annotation_patient_str_beauty = gsub("^R\\.", 'HCM\\.',current_analysis[[CURRENT_RUNNAME]]$annotation_patient_str)
current_analysis[[CURRENT_RUNNAME]]$annotation_patient_str_beauty = gsub("^H\\.", 'Ctrl1\\.',current_analysis[[CURRENT_RUNNAME]]$annotation_patient_str_beauty)
current_analysis[[CURRENT_RUNNAME]]$annotation_patient_str_beauty = gsub("^T\\.", 'Ctrl2\\.',current_analysis[[CURRENT_RUNNAME]]$annotation_patient_str_beauty)
# Create custom order for annotation
current_analysis[[CURRENT_RUNNAME]]$annotation_paper_fct = 
    factor(current_analysis[[CURRENT_RUNNAME]]$annotation_paper_str, levels=c("vRooij", "Hu", "Teichmann"))

# --- end of load / set up data

################################################################################
# Create the plots
# Modified version of parts of main script

if (!dir.exists(paste0(base_dir,'Rplots/Tomo-seq/'))) { dir.create(paste0(base_dir,'Rplots/Tomo-seq/'), recursive = T) }





# Re-create the reference plot for convenience
custom_colors = c('#bd0020','#9d9d9c','#575756')
p=DimPlot(current_analysis[[CURRENT_RUNNAME]], group.by = 'annotation_paper_beatified', cols = custom_colors, label = T, repel = T, label.size = 7/.pt, pt.size = NULL)+
        theme_void()+ggtitle(element_blank())+theme(legend.position = 'none')
    # test
ggsave(plot = p, filename = paste0(base_dir,'Rplots/Tomo-seq/',CURRENT_RUNNAME,'_2_umapLabeled_by_','annotation_paper_beatified','_style4-customColors.pdf'), height=172/3-4, width=172/3-4, units='mm')
    


# Comparison plot (bars per condition)
p_lists=shorthand_custom_compositeplot(seuratObject_list=current_analysis, 
                                 gene_lists=Tomoseq_genes_of_interest, 
                                 seuratObjectNameToTake=CURRENT_RUNNAME, 
                                 group.by='annotation_paper_beatified', 
                                 group.by2='annotation_patient_fct',
                                 zscore=T, SUBDIR = 'Tomo-seq/') 
# Comparison plot with bars per gene
# bit useless for datastes with many genes, as we can't display all 188 genes ..
# but might be useful for set w/ 19 genes
shorthand_custom_boxplot(seuratObject_list=current_analysis,
                                 gene_lists=Tomoseq_genes_of_interest,
                                 seuratObjectNameToTake=CURRENT_RUNNAME,
                                 group.by='annotation_paper_oneletter_fct',
                                 topX=25, mylimits=.01,
                                 SUBDIR = 'Tomo-seq/')





# Expression projected on UMAP
p_list_modules = lapply(names(Tomoseq_genes_of_interest), function(sublist_name) {
            shorthand_seurat_custom_expr(seuratObject = current_analysis[[CURRENT_RUNNAME]], 
                             gene_of_interest = Tomoseq_genes_of_interest[[sublist_name]],
                             textsize = 6, pointsize = .5, custom_title = sublist_name, mymargin = .1, zscore = T) 
    })
names(p_list_modules) = names(Tomoseq_genes_of_interest)

# Overview plot of all three gene sets
p_modules=wrap_plots(p_list_modules, nrow=1)
ggsave(filename = paste0(base_dir, 'Rplots/Tomo-seq/', CURRENT_RUNNAME, '_9_customUMAPs_pc1pc2_couptfii.pdf'), plot = p_modules,
       height=(PANEL_WIDTH*3-4)/3, width=(PANEL_WIDTH*3-4), units='mm', device=cairo_pdf)    
ggsave(filename = paste0(base_dir, 'Rplots/Tomo-seq/', CURRENT_RUNNAME, '_9_customUMAPs_pc1pc2_couptfii.png'), plot = p_modules,
       height=(PANEL_WIDTH*3-4)/3, width=(PANEL_WIDTH*3-4), units='mm', dpi=1200) 

# Make a plot of the selected UMAP
ggsave(filename = paste0(base_dir, 'Rplots/Tomo-seq/', CURRENT_RUNNAME, '_9_customUMAPs_pc1pc2_couptfii_overlapExternal.png'), 
       plot = p_list_modules$genespc1pc2COUPTFIIoverlap_19,
       height=(PANEL_WIDTH-4), width=(PANEL_WIDTH-4), units='mm', dpi=1200) 
ggsave(filename = paste0(base_dir, 'Rplots/Tomo-seq/', CURRENT_RUNNAME, '_9_customUMAPs_pc1pc2_couptfii_overlapExternal.pdf'), 
       plot = p_list_modules$genespc1pc2COUPTFIIoverlap_19,
       height=(PANEL_WIDTH-4), width=(PANEL_WIDTH-4), units='mm', dpi=1200, device = cairo_pdf) 



# Now the single genes of interest

# enriched genes expr on UMAPs
ZOOM_FACTOR=2
p_list_tomo19 = lapply(Tomoseq_genes_of_interest$genespc1pc2COUPTFIIoverlap_19, function(gene_name) {
            shorthand_seurat_custom_expr(seuratObject = current_analysis[[CURRENT_RUNNAME]], 
                             gene_of_interest = gene_name,
                             textsize = 6*ZOOM_FACTOR, pointsize = .5, custom_title = gene_name, mymargin = .5*ZOOM_FACTOR, zscore = T) 
                                # note: text size twice as large, because i save at zoom 200%, as trick to reduce point size
    })
# p_clusters=wrap_plots(p_list_clusters, nrow=1)
# ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_CL.pdf'), plot = p_clusters,
#        height=(PANEL_WIDTH*3-4)/5, width=(PANEL_WIDTH*3-4), units='mm', device=cairo_pdf)    
# ggsave(filename = paste0(base_dir, 'Rplots/', CURRENT_RUNNAME, '_9_customUMAPs_CL.png'), plot = p_clusters,
#        height=(PANEL_WIDTH*3-4)/5, width=(PANEL_WIDTH*3-4), units='mm', dpi=1200)    
p_tomo19=wrap_plots(p_list_tomo19, nrow=4)
MYWIDTH=(175/2-4)
ggsave(filename = paste0(base_dir, 'Rplots/Tomo-seq/', CURRENT_RUNNAME, '_9_customUMAPs_Tomo19genes.png'), plot = p_tomo19,
       height=MYWIDTH*ZOOM_FACTOR*4/5, width=MYWIDTH*ZOOM_FACTOR, units='mm', dpi=1200)    
ggsave(filename = paste0(base_dir, 'Rplots/Tomo-seq/', CURRENT_RUNNAME, '_9_customUMAPs_Tomo19genes.pdf'), plot = p_tomo19,
       height=MYWIDTH*ZOOM_FACTOR*4/5, width=MYWIDTH*ZOOM_FACTOR, units='mm', device=cairo_pdf)

for (idx in 1:length(Tomoseq_genes_of_interest$genespc1pc2COUPTFIIoverlap_19)) {
    
    gene_name = Tomoseq_genes_of_interest$genespc1pc2COUPTFIIoverlap_19[idx]
    p = p_list_tomo19[[idx]]
    
    print(paste0('Saving plot for ', gene_name))
    
    ggsave(filename = paste0(base_dir, 'Rplots/Tomo-seq/', CURRENT_RUNNAME, '_9_customUMAPs_Tomo19genes_',gene_name,'.png'), plot = p,
       height=(PANEL_WIDTH/3-4)*ZOOM_FACTOR, width=(PANEL_WIDTH/3-4)*ZOOM_FACTOR, units='mm', dpi=1200)    
    ggsave(filename = paste0(base_dir, 'Rplots/Tomo-seq/', CURRENT_RUNNAME, '_9_customUMAPs_Tomo19genes_',gene_name,'.pdf'), plot = p,
       height=(PANEL_WIDTH/3-4)*ZOOM_FACTOR, width=(PANEL_WIDTH/3-4)*ZOOM_FACTOR, units='mm', device=cairo_pdf)

}






