################################################################################

# Processing qPCR data to check correlations

################################################################################
library(pheatmap)

filename_data = 'qPCR Raw_shared regulon 2_topgenes_MWe_vClean.xlsx' 
  # will be made available as supp dataset.

################################################################################

qPCR_Data_Maya =
    openxlsx::read.xlsx(paste0('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3b/qPCR/',filename_data), sheet = 'relevant_data_CT')

# Little preprocessing
qPCR_Data_Maya$metadata_type = NA
qPCR_Data_Maya$metadata_type[grepl('Control', qPCR_Data_Maya$metadata_annotation)] = 'Ctrl'
qPCR_Data_Maya$metadata_type[!grepl('Control', qPCR_Data_Maya$metadata_annotation)] = 'HCM'

gene_names_noMYL2 = colnames(qPCR_Data_Maya)[!grepl('^metadata|^MYL2$',colnames(qPCR_Data_Maya))]
gene_names_all = colnames(qPCR_Data_Maya)[!grepl('^metadata',colnames(qPCR_Data_Maya))]
print(paste0('For reference, genes: ',toString(sort(gene_names_all))))

# additional table with 2^-dCT
qPCR_Data_Maya_2dCT = qPCR_Data_Maya
qPCR_Data_Maya_2dCT[,gene_names_all]=2^-qPCR_Data_Maya[,gene_names_all]

# Additional table with -dCT
qPCR_Data_Maya_minusdCT = qPCR_Data_Maya
qPCR_Data_Maya_minusdCT[,gene_names_all]=-qPCR_Data_Maya_minusdCT[,gene_names_all]

# Organize data
qPCR_dataFrame = list(qPCR_Data_Maya=qPCR_Data_Maya, qPCR_Data_Maya_2dCT=qPCR_Data_Maya_2dCT,qPCR_Data_Maya_minusdCT=qPCR_Data_Maya_minusdCT)

################################################################################

QPCRSET = 'qPCR_Data_Maya'
QPCRSET = 'qPCR_Data_Maya_2dCT'
QPCRSET = 'qPCR_Data_Maya_minusdCT'

qPCR_dataFrame[['qPCR_Data_Maya']]$metadata_annotation[qPCR_dataFrame[['qPCR_Data_Maya']]$MYL2 < -15 & !is.na(qPCR_dataFrame[['qPCR_Data_Maya']]$MYL2)]
blacklist='HCM106'
    # HCM106 is blacklisted below

qPCR_Data_Maya_HCM = qPCR_dataFrame[[QPCRSET]][qPCR_dataFrame[[QPCRSET]]$metadata_type=='HCM',]
qPCR_Data_Maya_Ctrl = qPCR_dataFrame[[QPCRSET]][qPCR_dataFrame[[QPCRSET]]$metadata_type=='Ctrl',]

my_unit = list('qPCR_Data_Maya'='(ΔCT)', 'qPCR_Data_Maya_2dCT'='(2^-ΔCT)', 'qPCR_Data_Maya_minusdCT'='(-ΔCT)')[[QPCRSET]]
p_list =
    lapply(sort(gene_names_noMYL2), function(current_gene) {
        #current_gene='ACTC1'
        #current_gene='CRYAB'
        # current_gene='ACTA1'
        qPCR_Data_Maya_sel =     qPCR_dataFrame[[QPCRSET]][(!is.na(qPCR_dataFrame[[QPCRSET]][['MYL2']])) & (!is.na(qPCR_dataFrame[[QPCRSET]][[current_gene]])) &
                                                    (!(qPCR_dataFrame[[QPCRSET]]$metadata_annotation %in% blacklist)),] # blacklist
        qPCR_Data_Maya_HCM_sel = qPCR_Data_Maya_HCM[(!is.na(qPCR_Data_Maya_HCM[['MYL2']])) & (!is.na(qPCR_Data_Maya_HCM[[current_gene]])) &
                                                        (!(qPCR_Data_Maya_HCM$metadata_annotation %in% blacklist)),] # blacklist
        
        lm_out = lm(formula = as.formula(paste0(current_gene, ' ~ MYL2')), data = qPCR_Data_Maya_HCM_sel) 
        cor.test_out = cor.test(qPCR_Data_Maya_HCM_sel$MYL2, qPCR_Data_Maya_HCM_sel[[current_gene]])
        R=cor.test_out$estimate; p=cor.test_out$p.value
        
        p_string = if (p >= .001) {paste0('p=',round(p,3))} else {paste0('p=10^',round(log10(p),1))}
        
        ggplot(qPCR_Data_Maya_sel, aes_string(x='MYL2', y=current_gene))+#, color='metadata_type'))+
            geom_point(data=qPCR_Data_Maya_sel[qPCR_Data_Maya_sel$metadata_type=='HCM',],size=.05, color='#bd0020')+
            geom_point(data=qPCR_Data_Maya_sel[qPCR_Data_Maya_sel$metadata_type=='Ctrl',],size=.05, color='black')+
            theme_bw()+give_better_textsize_plot(6)+theme(legend.position = 'none')+
            geom_abline(slope = lm_out$coefficients[2], intercept = lm_out$coefficients[1], size=.25)+
            #scale_color_manual(values = c('black','#bd0020'), breaks = c('Ctrl','HCM'))+
            ggtitle(paste0('R=',round(R,2), '\n', p_string))+
            xlab(paste0('MYL2 (-ΔCT)'))+ylab(paste0(current_gene,' (-ΔCT)'))
            #xlab(expression(paste(italic('MYL2'),' (ΔCT)')))+ylab(expression(paste(italic(current_gene), '(ΔCT)')))
            #annotate(geom='text', -Inf, -Inf, label=paste0('R=',round(R,2), '\n', p_string), hjust = -0.05, vjust = -.3)
        })

p=wrap_plots(p_list, nrow = 2)
p

if (!exists('nosave')) {
  #ggsave(filename = paste0(base_dir,'Rplots/qPCR_scatters_MYL2.pdf'), 
  #        plot = p, width=172-4, height=PANEL_HEIGHT*2, units='mm', device = cairo_pdf) # 184.6/3*2-4
  ggsave(filename = paste0(base_dir,'Rplots/qPCR_scatters_MYL2_',gsub('Δ','d',my_unit),'.pdf'), 
          plot = p, width=172-4, height=76, units='mm', device = cairo_pdf) # 184.6/3*2-4
  
  ggsave(filename = paste0(base_dir,'Rplots/qPCR_scatters_MYL2_',gsub('Δ','d',my_unit),'_style2.pdf'), 
          plot = p, width=172-4, height=(172-4)/8*2+15, units='mm', device = cairo_pdf) # 184.6/3*2-4
}

# Show sample sizes
sample_sizes = sapply(sort(gene_names_noMYL2), function(current_gene) {
    # current_gene = 'ACTA1'
    c(all=sum((!is.na(qPCR_dataFrame[[QPCRSET]][['MYL2']])) & (!is.na(qPCR_dataFrame[[QPCRSET]][[current_gene]])) &
                                                    (!(qPCR_dataFrame[[QPCRSET]]$metadata_annotation %in% blacklist))),
        HCM=sum((!is.na(qPCR_Data_Maya_HCM[['MYL2']])) & (!is.na(qPCR_Data_Maya_HCM[[current_gene]])) &
                                                    (!(qPCR_Data_Maya_HCM$metadata_annotation %in% blacklist))),
        Ctrl=sum((!is.na(qPCR_Data_Maya_Ctrl[['MYL2']])) & (!is.na(qPCR_Data_Maya_Ctrl[[current_gene]])) &
                                                    (!(qPCR_Data_Maya_Ctrl$metadata_annotation %in% blacklist)))
      
      ) # blacklist)
        
        
})
sample_sizes
View(data.frame(colnames(sample_sizes), sample_sizes[1,]))
# n_ACTA1 = 42, n_ACTC1 = 74, n_CKM = 71, n_COX6A2 = 77, n_CRYAB = 80, n_CSRP3 = 73, n_GAPDH = 71, n_HSPB1 = 70, n_MB = 75, n_MYL3 = 68, n_MYL9 = 83, n_SLC25A3 = 68, n_SLC25A4 = 74, n_TNNC1 = 75, n_TPM1 = 78, n_UBC = 61
min(sample_sizes[1,])
max(sample_sizes[1,])

paste0( paste0(colnames(sample_sizes), ', ', sample_sizes[2,],'; ') , collapse = '')

# Test code

lm_out = lm(formula = MYL2 ~ ACTA1, data = qPCR_Data_Maya) 


# Misc plot
ggplot(qPCR_dataFrame[['qPCR_Data_Maya_minusdCT']], aes(x=metadata_type,y=MYL2, color=metadata_type))+
    geom_jitter()+theme_bw()+theme(legend.position='none')

################################################################################
# Creating a full heatmap

qPCR_dataFrame[[QPCRSET]][,gene_names_all]
2:dim(qPCR_dataFrame[[QPCRSET]][,gene_names_all])[2]

gene_idxs = 1:dim(qPCR_dataFrame[[QPCRSET]][,gene_names_all])[2]
corr_matrix = 
    as.data.frame(sapply(gene_idxs, function(idx1) {
        sapply(gene_idxs, function(idx2) {
             # idx1=1; idx2=2
              x=qPCR_dataFrame[[QPCRSET]][qPCR_dataFrame[[QPCRSET]]$metadata_type=='HCM'&(!(qPCR_dataFrame[[QPCRSET]]$metadata_annotation %in% blacklist)),gene_names_all][idx1]
              y=qPCR_dataFrame[[QPCRSET]][qPCR_dataFrame[[QPCRSET]]$metadata_type=='HCM'&(!(qPCR_dataFrame[[QPCRSET]]$metadata_annotation %in% blacklist)),gene_names_all][idx2]
              the_cor = cor(x[(!is.na(x))&(!is.na(y))], y[(!is.na(x))&(!is.na(y))])
              return(the_cor)
              })}))
p_matrix = 
    as.data.frame(sapply(gene_idxs, function(idx1) {
        sapply(gene_idxs, function(idx2) {
             # idx1=1; idx2=2
              x=qPCR_dataFrame[[QPCRSET]][qPCR_dataFrame[[QPCRSET]]$metadata_type=='HCM'&(!(qPCR_dataFrame[[QPCRSET]]$metadata_annotation %in% blacklist)),gene_names_all][idx1]
              y=qPCR_dataFrame[[QPCRSET]][qPCR_dataFrame[[QPCRSET]]$metadata_type=='HCM'&(!(qPCR_dataFrame[[QPCRSET]]$metadata_annotation %in% blacklist)),gene_names_all][idx2]
              cor.test_out = cor.test(x[(!is.na(x))&(!is.na(y))], y[(!is.na(x))&(!is.na(y))])
              return(cor.test_out$p.value)
              })}))
colnames(corr_matrix) = gene_names_all
rownames(corr_matrix) = gene_names_all

# ====
# Some statistics
# Note that there are less pairs than the matrix has squares
# Exclude self-pairing and double pairing
# Not necessary to exclude double pairing when using percentages

p_values_adjusted = p.adjust(as.vector(as.matrix(p_matrix[corr_matrix<.9999999])))
fdr_values        = p.adjust(as.vector(as.matrix(p_matrix[corr_matrix<.9999999])), method='fdr')
sum(fdr_values<.05)
sum(fdr_values<.05)/length(fdr_values)
paste0('Percentage fdr<.05: ',round(sum(fdr_values<.05)/length(fdr_values)*100,0),'%')

fdr_values_Rpos = p.adjust(as.vector(as.matrix(p_matrix[(corr_matrix<.9999999)&(corr_matrix>0)])), method='fdr')
paste0('Of R>0, percentage fdr<.05: ',round(sum(fdr_values_Rpos<.05)/length(fdr_values_Rpos)*100,0),'%')

sum(corr_matrix[corr_matrix<.9999999]>0)/length(corr_matrix[corr_matrix<.9999999])
paste0('Percentage R>0: ',round(sum(corr_matrix[corr_matrix<.9999999]>0)/length(corr_matrix[corr_matrix<.9999999])*100,0),'%')

# Maybe it makes more sense to check what the odds are these values are >0

# Sanity check with scatters
MYL_corrs_check = corr_matrix$MYL2
names(MYL_corrs_check) = colnames(corr_matrix) 
MYL_corrs_check
#ggplot(data.frame(x=qPCR_dataFrame[[QPCRSET]]$MYL2, y=qPCR_dataFrame[[QPCRSET]]$TPM1), aes(x=x, y=y))+
#    geom_point()

# ====
# Heatmaps

my_color_gradient = colorRampPalette(rainbow_colors)(201)

p = pheatmap(corr_matrix, treeheight_row = 0, treeheight_col = 0, fontsize = 8, 
             cellwidth = 10, cellheight = 10, border_color=NA, color=my_color_gradient, breaks=seq(-1,1,.01))

if (!exists('nosave')) {
  ggsave(filename = paste0(base_dir,'Rplots/qPCR_heatmap_all17genes_HCMsamples.pdf'), 
          plot = p, width=172/3*2-4, height=172/3*2-4, units='mm', device = cairo_pdf) # 184.6/3*2-4
        # 20*10/.pt
}

# smaller size
myBreakList=seq(-1,1,.01)
myColors=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(myBreakList))
myColors=colorRampPalette(rainbow_colors)(length(myBreakList))
rainbow_colors
p = pheatmap(corr_matrix, treeheight_row = 0, treeheight_col = 0, fontsize = 5, 
             cellwidth = 5, cellheight = 5, border_color=NA, legend_breaks = seq(-1,1,.5), breaks=myBreakList, color=myColors)
             # 7/.pt # cell height in mm
if (!exists('nosave')) {
  ggsave(filename = paste0(base_dir,'Rplots/qPCR_heatmap_all17genes_HCMsamples-smallCustomCols.pdf'), 
          plot = p, width=172/3-4, height=172/3-4, units='mm', device = cairo_pdf) # 184.6/3*2-4
}

ggplot(data.frame(correlation=corr_matrix[corr_matrix<1]),aes(x=correlation))+
  geom_histogram()+theme_bw()


################################################################################
# Just to remind ourselves, which genes are we checking out

ANALYSIS_NAME = "ROOIJonly_RID2l"
load(paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_sorted_shortname.Rdata')) # core_regulons_sorted_shortname

inclusion_list = core_regulons_sorted_shortname$s.R.2 %in% gene_names_all*1
names(inclusion_list) = core_regulons_sorted_shortname$s.R.2
View(data.frame(gene=names(inclusion_list), included=inclusion_list))

core_regulons_sorted_shortname$s.R.2[1:20] %in% gene_names_all*1

################################################################################
# Looking at raw data


# Processing qPCR data to check correlations
qPCR_Data_Maya_rawCT =
    openxlsx::read.xlsx(paste0('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3b/qPCR/',filename_data), sheet = 'raw_data')

# Little preprocessing
qPCR_Data_Maya_rawCT$metadata_type = NA
qPCR_Data_Maya_rawCT$metadata_type[grepl('Control', qPCR_Data_Maya_rawCT$metadata_annotation)] = 'Ctrl'
qPCR_Data_Maya_rawCT$metadata_type[!grepl('Control', qPCR_Data_Maya_rawCT$metadata_annotation)] = 'HCM'

rawCT_gene_names_noMYL2 = colnames(qPCR_Data_Maya_rawCT)[!grepl('^metadata|^MYL2$',colnames(qPCR_Data_Maya_rawCT))]
rawCT_gene_names_all = colnames(qPCR_Data_Maya_rawCT)[!grepl('^metadata',colnames(qPCR_Data_Maya_rawCT))]



# Raw data
ggplot(qPCR_Data_Maya_rawCT, aes(x=metadata_type,y=-RPL32, color=metadata_type))+
    geom_jitter()+theme_bw()+theme(legend.position='none')+xlab(element_blank())+give_better_textsize_plot(12)+ylab('RPL32 (-CT)')+
    geom_boxplot(alpha=.2,color='black')+
    ylim(c(-max(qPCR_Data_Maya_rawCT$RPL32, na.rm=T)*1.1,0))

ggplot(qPCR_Data_Maya_rawCT, aes(x=metadata_type,y=GAPDH, color=metadata_type))+
    geom_jitter()+theme_bw()+theme(legend.position='none')+xlab(element_blank())+give_better_textsize_plot(12)+ylab('GAPDH (CT)')

p=ggplot(qPCR_Data_Maya_rawCT, aes(x=GAPDH, y=RPL32, color=metadata_type))+
    geom_point()+theme_bw()+theme(legend.position='none')+give_better_textsize_plot(12)+xlab('GAPDH (CT)')+ylab('RPL32 (CT)')
p

if (!exists('nosave')) {
  ggsave(filename = paste0(base_dir,'Rplots/qPCR_misc_GAPDH-vs-RPL32.pdf'), 
          plot = p, width=172/3-4, height=172/3-4, units='mm', device = cairo_pdf) # 184.6/3*2-4
}

p1=ggplot(qPCR_Data_Maya_rawCT, aes(x=metadata_type,y=MYL2, color=metadata_type))+
    geom_jitter()+theme_bw()+theme(legend.position='none')+xlab(element_blank())+give_better_textsize_plot(12)+ylab('MYL2 (CT)')
p2=ggplot(qPCR_Data_Maya_rawCT, aes(x=metadata_type,y=ACTA1, color=metadata_type))+
    geom_jitter()+theme_bw()+theme(legend.position='none')+xlab(element_blank())+give_better_textsize_plot(12)+ylab('ACTA1 (CT)')
p3=ggplot(qPCR_Data_Maya_rawCT, aes(x=metadata_type,y=CRYAB, color=metadata_type))+
    geom_jitter()+theme_bw()+theme(legend.position='none')+xlab(element_blank())+give_better_textsize_plot(12)+ylab('CRYAB (CT)')

p=p1+p2+p3
p

if (!exists('nosave')) {
  ggsave(filename = paste0(base_dir,'Rplots/qPCR_misc_raw-MYL2-ACTA1-CRYAB.pdf'), 
          plot = p, width=172*2/3-4, height=172/3-4, units='mm', device = cairo_pdf) # 184.6/3*2-4
}


