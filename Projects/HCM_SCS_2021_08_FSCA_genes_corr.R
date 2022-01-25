

################################################################################
# Setting up 

library(dplyr)

# DATASET_NAME = 'ROOIJonly_RID2l_clExtended'
DATASET_NAME = 'ROOIJonly.sp.bt_RID2l_clExtended'
DATASET_NAME_FSCA = paste0(DATASET_NAME, '_FSCA')
TRIAGE_DATASET_NAME='ROOIJonly_TRIAGE'

if (!exists('current_analysis')) {current_analysis = list()}
current_analysis[[DATASET_NAME]] =
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))

# Note afterwards -- TRIAGE wasn't intensely used for analysis
# Triage parts can be skipped in this script
current_analysis[[TRIAGE_DATASET_NAME]] =
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',TRIAGE_DATASET_NAME,'.h5seurat'))




# Also required: ---->
load(paste0(base_dir,'Rdata/ROOIJonly.sp.bt_RID2l_clExtended_core_regulons_sorted.Rdata'))
    # core_regulons_sorted
    # core_regulons



# Load index data
plateJE10 = read.csv('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_analyzed_by_Joep/2020_04_AnneJoepData/Copy-of-data-VersionAnne/IndexDatafiles/2018-244 exp cardiologie_plate2.csv', header=T, sep=",", row.names=NULL, skip=13)
plateJE11 = read.csv('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_analyzed_by_Joep/2020_04_AnneJoepData/Copy-of-data-VersionAnne/IndexDatafiles/2018-244 exp cardiologie_plate3merged.csv', header=T ,sep=",", row.names=NULL, skip=13)
plateAL1 = read.csv('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_analyzed_by_Joep/2020_04_AnneJoepData/Copy-of-data-VersionAnne/IndexDatafiles/AL1.csv', header=T ,sep=",", row.names=NULL, skip=15)
plateAL2 = read.csv('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_analyzed_by_Joep/2020_04_AnneJoepData/Copy-of-data-VersionAnne/IndexDatafiles/AL2.csv', header=T ,sep=",", row.names=NULL, skip=15)
# Fix well names
plateJE10$Cell = unlist(lapply(plateJE10$Well, function(x) {return(paste0('JE10_',substr(x,1,1),'_',substr(x,2,3)))}))
plateJE11$Cell = unlist(lapply(plateJE11$Well, function(x) {return(paste0('JE11_',substr(x,1,1),'_',substr(x,2,3)))}))
plateAL1$Cell = unlist(lapply(plateAL1$Well, function(x) {return(paste0('AL01_',substr(x,1,1),'_',substr(x,2,3)))}))
plateAL1 = plateAL1[colnames(plateJE10)] # Throw away ridiculous column overload
plateAL2$Cell = unlist(lapply(plateAL2$Well, function(x) {return(paste0('AL02_',substr(x,1,1),'_',substr(x,2,3)))}))
plateAL2 = plateAL2[colnames(plateJE10)] # Throw away ridiculous column overload

indexData_MW = rbind(plateAL1, plateAL2, plateJE10, plateJE11)
# Fix column names
colnames(indexData_MW) <- c('Well','Events','Parent','FSC_A','FSC_W','FSC_H','SSC_A','SSC_W','SSC_H','BV421_A', 'Cell')
# Convert measurement columns to numeric
for (column_name in c('FSC_A','FSC_W','FSC_H','SSC_A','SSC_W','SSC_H','BV421_A')) {
    indexData_MW[[column_name]] = as.numeric(unlist(lapply(indexData_MW[[column_name]], function(x) {str_replace(x, ',','.')})))    
}


# Now also generate the cell names compatible with Seurat naming
# plateAL1$Well[1:10] # 16[A-P]x24[1-24]
# colnames(current_analysis$ROOIJonly_RID2l)[1:10] 
# First generate conversion table
alph=toupper(c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p'))
wellsnamed    = paste0(rep(alph, each=24), '_', rep(1:24, times=16))
wellsnumbered = paste0('X',sprintf("%03d", 1:384))
conversion_name_to_nr = wellsnumbered
names(conversion_name_to_nr) = wellsnamed
indexData_MW$Cell_newname =
    sapply(str_split(string = indexData_MW$Cell, pattern = '_'), function(x) {
        paste0(x[[1]],'_',conversion_name_to_nr[paste0(x[[2]],'_',x[[3]])])
    }
    )
rownames(indexData_MW) = indexData_MW$Cell_newname

save(list = 'indexData_MW', file = paste0(base_dir,'Rdata/FSCA__indexData_MW.Rdata')) # indexData_MW, colnames match cell names
        # load(file = paste0(base_dir,'Rdata/FSCA__indexData_MW.Rdata')) # indexData_MW
        # e.g. use indexData_MW[colnames(current_analysis[[DATASET_NAME]]),]$FSC_A

################################################################################
# Insert data into clExtended

if (F) {
    
    current_analysis[[DATASET_NAME_FSCA]] = current_analysis[[DATASET_NAME]]
    current_analysis[[DATASET_NAME_FSCA]][['FSCA']]=indexData_MW[colnames(current_analysis[[DATASET_NAME]]),]$FSC_A
    current_analysis[[DATASET_NAME_FSCA]][['DAPI']]=indexData_MW[colnames(current_analysis[[DATASET_NAME]]),]$BV421_A
    
    SaveH5Seurat(object = current_analysis[[DATASET_NAME_FSCA]], filename = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME_FSCA,'.h5seurat'))
        # (Currently, this file is not used by other scripts)
    
} 

# to load
if (F) {
    current_analysis[[DATASET_NAME_FSCA]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME_FSCA,'.h5seurat'))
}

################################################################################
# Project data on UMAP

FeaturePlot(object = current_analysis[[DATASET_NAME_FSCA]], features = 'FSCA', cols = rainbow_colors, pt.size = 2)
DimPlot(object = current_analysis[[DATASET_NAME_FSCA]], group.by='annotation_patient_fct')
VlnPlot(object = current_analysis[[DATASET_NAME_FSCA]], features = 'FSCA')#, cols = rainbow_colors, pt.size = 2)

# # Triage clusters
# current_analysis$ROOIJonly_TRIAGE_clExt[['FSCA']]=indexData_MW[colnames(current_analysis[[DATASET_NAME]]),]$FSC_A
# FeaturePlot(object = current_analysis$ROOIJonly_TRIAGE_clExt, features = 'FSCA', cols = rainbow_colors, pt.size = 2)
# VlnPlot(object = current_analysis$ROOIJonly_TRIAGE_clExt, features = 'FSCA')#, cols = rainbow_colors, pt.size = 2)

################################################################################

# Now perform the correlations

GENE_MIN_PERCENTAGE=.1
GENE_MIN_PERCENTAGE=.33

# First subset on where we have data for
selected_cells = colnames(current_analysis[[DATASET_NAME]])[colnames(current_analysis[[DATASET_NAME]]) %in% indexData_MW$Cell_newname]
current_analysis$ROOIJonly_RID2l_FSCA = subset(current_analysis[[DATASET_NAME]], cells = selected_cells)
# Re-arrange df such that cell order matches (redundant because already in right order, but to make sure..)
rownames(indexData_MW) = indexData_MW$Cell_newname
indexData_MW=indexData_MW[selected_cells,]
# Sanity check on order of cells
all(selected_cells==colnames(current_analysis$ROOIJonly_RID2l_FSCA))
all(selected_cells==indexData_MW[selected_cells,]$Cell_newname)

# Now select genes
# Let's do this per patient 
the_patients = unique(current_analysis$ROOIJonly_RID2l_FSCA$annotation_patient_str)
correlations_FSCA_per_patient_list = 
    lapply(the_patients, function(current_patient) {
        
        print(paste0('Analyzing patient ',current_patient))
        
        cell_sel2 = current_analysis$ROOIJonly_RID2l_FSCA$annotation_patient_str == current_patient
        genes = rownames(current_analysis$ROOIJonly_RID2l_FSCA@assays$RNA@data)[rowMeans(current_analysis$ROOIJonly_RID2l_FSCA@assays$RNA@data[,cell_sel2]>0)>GENE_MIN_PERCENTAGE]
        
        # Perform test
        correlations_FSCA = as.data.frame(t(data.frame(lapply(genes, function(g) {
                cor_out = cor.test(current_analysis$ROOIJonly_RID2l_FSCA@assays$RNA@data[g,cell_sel2], indexData_MW$FSC_A[cell_sel2])        
                return(c(cor_out$estimate, p=cor_out$p.value))
            }))))
        rownames(correlations_FSCA) = genes
        correlations_FSCA$fdr=p.adjust(correlations_FSCA$p, method='fdr')
        correlations_FSCA$p.adjust=p.adjust(correlations_FSCA$p, method='BH')
        
        return(correlations_FSCA)
        
    })
names(correlations_FSCA_per_patient_list) = the_patients

View(correlations_FSCA_per_patient_list$R.P5)

# Combine the two frames
#correlations_FSCA_per_patient_list$R.P4$patient = 'R.P4'
#correlations_FSCA_per_patient_list$R.P5$patient = 'R.P5'
genes_overlap = intersect(rownames(correlations_FSCA_per_patient_list$R.P4), rownames(correlations_FSCA_per_patient_list$R.P5))
colnames(correlations_FSCA_per_patient_list$R.P4) = paste0('P4_',colnames(correlations_FSCA_per_patient_list$R.P4))
colnames(correlations_FSCA_per_patient_list$R.P5) = paste0('P5_',colnames(correlations_FSCA_per_patient_list$R.P5))
correlations_FSCA_per_patient_combined = cbind(correlations_FSCA_per_patient_list$R.P4[genes_overlap,], correlations_FSCA_per_patient_list$R.P5[genes_overlap,])
correlations_FSCA_per_patient_combined$gene_symbol = shorthand_cutname( rownames(correlations_FSCA_per_patient_combined) )
correlations_FSCA_per_patient_combined$gene = rownames(correlations_FSCA_per_patient_combined) 

tresholds_p4=calc_limits(correlations_FSCA_per_patient_combined$P4_cor, 0.1)[2]
tresholds_p4
tresholds_p5=calc_limits(correlations_FSCA_per_patient_combined$P5_cor, 0.1)[2]
tresholds_p5

#selected_genes=correlations_FSCA_per_patient_combined$P4_cor>.2&correlations_FSCA_per_patient_combined$P5_cor>.1
selected_genes=correlations_FSCA_per_patient_combined$P4_cor>tresholds_p4 & correlations_FSCA_per_patient_combined$P5_cor>tresholds_p5
shorthand_cutname(rownames(correlations_FSCA_per_patient_combined[selected_genes,]))

p=ggplot(correlations_FSCA_per_patient_combined, aes(x=P4_cor, y=P5_cor))+
    geom_point(size=.5, shape=1, color='grey')+
    geom_point(size=.5, shape=1, data=correlations_FSCA_per_patient_combined[selected_genes,], color='red')+
    geom_text_repel(data=correlations_FSCA_per_patient_combined[selected_genes,], color='red',
                    mapping=aes(label=gene_symbol), 
                    max.overlaps = Inf, min.segment.length = 0, force = 4, size=6/.pt, segment.size=.25)+
    theme_bw()+give_better_textsize_plot(8)+
    xlim(c(min(correlations_FSCA_per_patient_combined$P4_cor)*1.1,max(correlations_FSCA_per_patient_combined$P4_cor)*1.3))+
    ylim(c(min(correlations_FSCA_per_patient_combined$P5_cor)*1.1,max(correlations_FSCA_per_patient_combined$P5_cor)*1.3))+
    xlab('Gene-cell size correlation patient 4')+ylab('Gene-cell size correlation patient 5')
    #geom_smooth(method = 'lm')+ylim(c(-.5,.5))+xlim(c(-.5,.5))
p
ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',DATASET_NAME,'_8_FSCA_corrs_Both-10pct.pdf'), 
       width = 172/3-4, height= 172/3-4, units='mm', device = cairo_pdf)
ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',DATASET_NAME,'_8_FSCA_corrs_Both-10pct_L.pdf'), 
       width = 172/2-4, height= 172/2-4, units='mm', device = cairo_pdf)

# With fitted line
correlations_FSCA_per_patient_combined$gene_symbol_italic = paste0('italic(',correlations_FSCA_per_patient_combined$gene_symbol,')')
p=ggplot(correlations_FSCA_per_patient_combined, aes(x=P4_cor, y=P5_cor))+
    geom_point(size=.5, shape=1, color='grey')+
    geom_point(size=.5, shape=1, data=correlations_FSCA_per_patient_combined[selected_genes,], color='red')+
    geom_text_repel(data=correlations_FSCA_per_patient_combined[selected_genes,], color='red',
                    mapping=aes(label=gene_symbol_italic), 
                    max.overlaps = Inf, min.segment.length = 0, force = 4, size=6/.pt, segment.size=.25, parse = T)+
    theme_bw()+give_better_textsize_plot(8)+
    xlim(c(min(correlations_FSCA_per_patient_combined$P4_cor)*1.01,max(correlations_FSCA_per_patient_combined$P4_cor)*1.01))+
    ylim(c(min(correlations_FSCA_per_patient_combined$P5_cor)*1.01,max(correlations_FSCA_per_patient_combined$P5_cor)*1.01))+
    xlab('Gene-cell size correlation patient 4')+ylab('Gene-cell size correlation patient 5')+
    geom_smooth(method = 'lm', size=.25)+
    ggtitle(' ') # This is a little "hack" to make sure this panel has equal size as other titled panels
    #+ylim(c(-.5,.5))+xlim(c(-.5,.5))
p
ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',DATASET_NAME,'_8_FSCA_corrs_Both-10pct_withFit.pdf'), 
       width = 172/3-4, height= 172/3-4, units='mm', device = cairo_pdf)
ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',DATASET_NAME,'_8_FSCA_corrs_Both-10pct_L_withFit.pdf'), 
       width = 172/2-4, height= 172/2-4, units='mm', device = cairo_pdf)


# Save data for later use
save(list = 'correlations_FSCA_per_patient_combined', file = paste0(base_dir,'Rdata/FSCA__correlations_FSCA_per_patient_combined.Rdata')) # correlations_FSCA_per_patient_combined
save(list = 'correlations_FSCA_per_patient_list', file = paste0(base_dir,'Rdata/FSCA__correlations_FSCA_per_patient_list.Rdata')) # correlations_FSCA_per_patient_list

# Export data for supp database
openxlsx::write.xlsx(x = correlations_FSCA_per_patient_combined, file = paste0(base_dir, 'Rplots/ROOIJonly_RID2l_FSCA_correlations_FSCA_per_patient_combined.xlsx'))


# Export list for analyis
gene_list_FSCA = shorthand_cutname( rownames(correlations_FSCA_per_patient_combined[selected_genes,]) )
# Export extended list for analysis
tresholds_p4_extended=calc_limits(correlations_FSCA_per_patient_combined$P4_cor[correlations_FSCA_per_patient_combined$P4_cor>0], 0.5)[2]
tresholds_p4_extended
tresholds_p5_extended=calc_limits(correlations_FSCA_per_patient_combined$P5_cor[correlations_FSCA_per_patient_combined$P5_cor>0], 0.5)[2]
tresholds_p5_extended
selected_genes_extended = correlations_FSCA_per_patient_combined$P4_cor>tresholds_p4_extended&correlations_FSCA_per_patient_combined$P5_cor>tresholds_p5_extended

################################################################################
# Now show some regulons
load(file=paste0(base_dir, 'Rdata/SCENIC_regulons_core_genes_sel.Rdata'))
load(file=paste0(base_dir, 'Rplots/ROOIJonly.sp.bt_RID2l_core_regulons_sorted_shortname.Rdata'))
# If already generated, load data generated above
if (F) {
    load(file = paste0(base_dir,'Rdata/FSCA__correlations_FSCA_per_patient_combined.Rdata')) # correlations_FSCA_per_patient_combined
    #load(file = paste0(base_dir,'Rdata/FSCA__correlations_FSCA_per_patient_list.Rdata')) # correlations_FSCA_per_patient_list   
}


REGULON_SETS = list(custom=core_regulons_sorted_shortname, SCENIC=SCENIC_regulons_core_genes_sel)

for (reg_set_name in names(REGULON_SETS)) {

    # reg_set_name='custom'
    
    current_regulon_set = REGULON_SETS[[reg_set_name]]
    
    list_p=list()
    for (reg_name in names(current_regulon_set)) {
        
        # reg_name = 'ATF4(+)'
        # reg_name = "s.R.2"
        
        current_highlight_genes = current_regulon_set[[reg_name]]
        
        correlations_FSCA_per_patient_combined$regulon_membership = 'no'
        correlations_FSCA_per_patient_combined$regulon_membership[correlations_FSCA_per_patient_combined$gene_symbol %in% current_highlight_genes] = 'yes'
        correlations_FSCA_per_patient_combined$regulon_membership = as.factor(correlations_FSCA_per_patient_combined$regulon_membership)
        # list_p[[reg_name]]=ggplot(correlations_FSCA_per_patient_combined %>% arrange(regulon_membership), aes(x=P4_cor, y=P5_cor))+
        #     geom_point(aes(color=regulon_membership))+
        #     geom_density_2d(aes(color=regulon_membership))+theme_bw()+give_better_textsize_plot(6)+ggtitle(gsub('\\(\\+\\)','',reg_name))+
        #     scale_color_manual(values=c('grey','black'))+theme(legend.position = 'none')
        
        lims_p4=c(min(c(-.21, correlations_FSCA_per_patient_combined$P4_cor)), max(c(.21, correlations_FSCA_per_patient_combined$P4_cor)))
        lims_p5=c(min(c(-.21, correlations_FSCA_per_patient_combined$P5_cor)), max(c(.21, correlations_FSCA_per_patient_combined$P5_cor)))
        p=ggplot(correlations_FSCA_per_patient_combined %>% arrange(regulon_membership), aes(x=P4_cor, y=P5_cor))+
            #geom_point(shape=1,alpha=.1,aes(color=regulon_membership))+
            geom_point(data=correlations_FSCA_per_patient_combined %>% subset(regulon_membership=='no'), color='grey', size=.25)+
            geom_density_2d(data=correlations_FSCA_per_patient_combined %>% subset(regulon_membership=='no'), color='dodgerblue2', size=.1)+
            geom_point(data=correlations_FSCA_per_patient_combined %>% subset(regulon_membership=='yes'), color='black', size=.25)+
            theme_bw()+give_better_textsize_plot(6)+ggtitle(gsub('\\(\\+\\)','',reg_name))+
            theme(legend.position = 'none')+
            scale_x_continuous(breaks = seq(-1, 1, by = .2), limits = lims_p4)+
            scale_y_continuous(breaks = seq(-1, 1, by = .2), limits = lims_p5)
        if (sum(correlations_FSCA_per_patient_combined$gene_symbol %in% current_highlight_genes)>1) {
            p=p+geom_density_2d(data=correlations_FSCA_per_patient_combined %>% subset(regulon_membership=='yes'), color='black', size=.1)
        }
        #p
        list_p[[reg_name]]=p
        
    }
    
    if (reg_set_name=='custom') {
        p=wrap_plots(list_p, ncol=1)&theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"))
        #p
        if (!exists('nosave')) {
            ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',DATASET_NAME,'_8_FSCA_corrs_',reg_set_name,'_RegulonsHighlighted.pdf'), 
                   width = (172/6)-4, height= (172*5/6)-4, units='mm', device = cairo_pdf)
        }
    }
    if (reg_set_name=='SCENIC') {
        p=wrap_plots(list_p, ncol=5)&theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"))
        #p
        if (!exists('nosave')) {
            ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',DATASET_NAME,'_8_FSCA_corrs_',reg_set_name,'_RegulonsHighlighted.pdf'), 
                   width = (172*5/6)-4, height= (172*5/6)-4, units='mm', device = cairo_pdf)
        }
    } 
    
}

# Show selected regulons
REG_NAMES=c('NFE2L1(+)', 's.R.2')
REG_SETS =c('SCENIC', 'custom')
wil_out=list(); tt_out = list()
pval_overview_df = data.frame(patient=character(), test=character(), regulon=character(), pval=numeric())
for (idx in 1:2) {
    
    # idx=2
    
    REG_NAME = gsub('\\(\\+\\)','',REG_NAMES[idx])
    REG_SET = REG_SETS[idx]
        
    current_highlight_genes = REGULON_SETS[[REG_SET]][[REG_NAME]]
    correlations_FSCA_per_patient_combined$regulon_membership = 'no'
    correlations_FSCA_per_patient_combined$regulon_membership[correlations_FSCA_per_patient_combined$gene_symbol %in% current_highlight_genes] = 'yes'
    correlations_FSCA_per_patient_combined$regulon_membership = as.factor(correlations_FSCA_per_patient_combined$regulon_membership)
    
    p=ggplot(correlations_FSCA_per_patient_combined, aes(x=P4_cor, y=P5_cor))+
        geom_point(size=.5, shape=1, color='grey')+
        geom_density_2d(data=correlations_FSCA_per_patient_combined %>% subset(regulon_membership=='no'), color='dodgerblue2', size=.1, bins=3)+
        geom_density_2d(data=correlations_FSCA_per_patient_combined %>% subset(regulon_membership=='yes'), color='red', size=.1, bins=3)+
        geom_point(size=.5, shape=1, data=correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$gene_symbol %in% current_highlight_genes,], color='red')+
        geom_text_repel(data=correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$gene_symbol %in% current_highlight_genes[1:min(20,length(current_highlight_genes))],], color='black',
                        mapping=aes(label=gene_symbol_italic), 
                        max.overlaps = Inf, min.segment.length = 0, force = 6, size=4/.pt, segment.size=.25, parse=T)+
        theme_bw()+give_better_textsize_plot(8)+
        xlim(c(min(correlations_FSCA_per_patient_combined$P4_cor)*1.01,max(correlations_FSCA_per_patient_combined$P4_cor)*1.01))+
        ylim(c(min(correlations_FSCA_per_patient_combined$P5_cor)*1.01,max(correlations_FSCA_per_patient_combined$P5_cor)*1.01))+
        xlab('Gene-cell size correlation patient 4')+ylab('Gene-cell size correlation patient 5')+
        ggtitle(paste0(
            gsub( 's\\.R\\.','Module ', gsub('\\(\\+\\)','',REG_NAME) )
            ,' genes highlighted'))
        #geom_smooth(method = 'lm', size=.25)#+ylim(c(-.5,.5))+xlim(c(-.5,.5))
    p
    
    if (!exists('nosave')) {
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',DATASET_NAME,'_8_FSCA_corrs_Reg',REG_NAME,'_withFit.pdf'), 
               width = 172/3-4, height= 172/3-4, units='mm', device = cairo_pdf)
        ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',DATASET_NAME,'_8_FSCA_corrs_Reg',REG_NAME,'_L_withFit.pdf'), 
               width = 172/2-4, height= 172/2-4, units='mm', device = cairo_pdf)
    }
 
    print(paste0('Wilcox test P4 --', REG_NAME))
    wil_out[[REG_NAME]][['P4']] = wilcox.test(correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$regulon_membership=='no',]$P4_cor, correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$regulon_membership=='yes',]$P4_cor)
    tt_out[[REG_NAME]][['P4']] = t.test(correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$regulon_membership=='no',]$P4_cor, correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$regulon_membership=='yes',]$P4_cor)
    pval_overview_df = rbind(pval_overview_df, data.frame(patient='P4', test='Wilcoxon', regulon=REG_NAME, pval=wil_out[[REG_NAME]][['P4']]$p.value))
    pval_overview_df = rbind(pval_overview_df, data.frame(patient='P4', test='t-test', regulon=REG_NAME, pval=tt_out[[REG_NAME]][['P4']]$p.value))
    
    print(paste0('Wilcox test P5 --', REG_NAME))
    wil_out[[REG_NAME]][['P5']] = wilcox.test(correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$regulon_membership=='no',]$P5_cor, correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$regulon_membership=='yes',]$P5_cor)
    tt_out[[REG_NAME]][['P5']] = t.test(correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$regulon_membership=='no',]$P5_cor, correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$regulon_membership=='yes',]$P5_cor)
    pval_overview_df = rbind(pval_overview_df, data.frame(patient='P5', test='Wilcoxon', regulon=REG_NAME, pval=wil_out[[REG_NAME]][['P5']]$p.value))
    pval_overview_df = rbind(pval_overview_df, data.frame(patient='P5', test='t-test', regulon=REG_NAME, pval=tt_out[[REG_NAME]][['P5']]$p.value))
    
    
    
}

openxlsx::write.xlsx(x = pval_overview_df, file = paste0(base_dir, 'Rplots/ROOIJonly_RID2l_FSCA_Regulons_pvals.xlsx'))

# Now do a statistical test on this
ggplot(correlations_FSCA_per_patient_combined)+
    geom_freqpoly(aes(x=P4_cor, color=regulon_membership))+theme_bw()
ggplot(correlations_FSCA_per_patient_combined)+
    geom_freqpoly(aes(x=P5_cor, color=regulon_membership))+theme_bw()

# Also export lists
gene_list_FSCA_extended = shorthand_cutname( rownames(correlations_FSCA_per_patient_combined[selected_genes_extended,]) )
write.table(x = gene_list_FSCA_extended, file = paste0(base_dir,'Rdata/',DATASET_NAME,'__gene_list_FSCA_extended.txt'), quote = F, row.names = F, col.names = F)


################################################################################
# Old version for custom regulon analysis only


regulon_membership = 
        Reduce(f=rbind, x=
            lapply(names(core_regulons_sorted[sapply(core_regulons_sorted, length)>0]), function(n) {data.frame(gene=core_regulons_sorted[[n]], regulon=n)}))
    membership_mapping = regulon_membership$regulon
    names(membership_mapping) = regulon_membership$gene
    correlations_FSCA_per_patient_combined$regulon_membership = membership_mapping[correlations_FSCA_per_patient_combined$gene]
    ggplot(correlations_FSCA_per_patient_combined, aes(x=P4_cor, y=P5_cor))+
        geom_point(shape=1,alpha=.1)+
        geom_density_2d(color='black')+
        geom_point(data=correlations_FSCA_per_patient_combined[!is.na(correlations_FSCA_per_patient_combined$regulon_membership),], aes(color=regulon_membership))+theme_bw()

ggplot(correlations_FSCA_per_patient_combined, aes(x=P4_cor, y=P5_cor))+
    geom_point(shape=1,alpha=.1)+
    geom_density_2d(aes(color=regulon_membership), bins=3)+
    geom_point(data=correlations_FSCA_per_patient_combined[!is.na(correlations_FSCA_per_patient_combined$regulon_membership),], aes(color=regulon_membership))+theme_bw()

TOPX_REG=15
correlations_FSCA_per_patient_combined_reg2 = correlations_FSCA_per_patient_combined[!is.na(correlations_FSCA_per_patient_combined$regulon_membership)&correlations_FSCA_per_patient_combined$regulon_membership=='s.R.2',]
p=ggplot(correlations_FSCA_per_patient_combined, aes(x=P4_cor, y=P5_cor))+
    geom_point(shape=1,color='lightgrey')+
    geom_density_2d(bins=4,color='grey')+
    geom_point(shape=1,color='orangered1',data=correlations_FSCA_per_patient_combined_reg2)+
    geom_density_2d(bins=4,color='orangered1',data=correlations_FSCA_per_patient_combined_reg2)+
    geom_text_repel(data=correlations_FSCA_per_patient_combined_reg2[correlations_FSCA_per_patient_combined_reg2$gene %in% core_regulons_sorted$s.R.2[1:TOPX_REG],],aes(label=gene_symbol), 
                    max.overlaps = Inf, min.segment.length = 0, force = 4, size=6/.pt, segment.size=.25)+
    # geom_text_repel(data=correlations_FSCA_per_patient_combined_reg2,aes(label=gene_symbol), 
    #                 max.overlaps = Inf, min.segment.length = 0, force = 4, size=4/.pt, segment.size=.1)+
    theme_bw()+give_better_textsize_plot(8)+xlab('Gene-cell size correlation patient 4')+ylab('Gene-cell size correlation patient 5')
p
ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',DATASET_NAME,'_8_FSCA_corrs_reg2.pdf'), 
       width = 172/3-4, height= 172/3-4, units='mm')
ggsave(plot = p,filename = paste0(base_dir, 'Rplots/',DATASET_NAME,'_8_FSCA_corrs_reg2_L.pdf'), 
       width = 172/2-4, height= 172/2-4, units='mm')
    


    #geom_density_2d(aes(color=regulon_membership), bins=3)+
    #geom_point(data=correlations_FSCA_per_patient_combined[!is.na(correlations_FSCA_per_patient_combined$regulon_membership),], aes(color=regulon_membership))+theme_bw()










