
################################################################################

# Processing qPCR data to check correlations
qPCR_Data_Maya =
    openxlsx::read.xlsx('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/qPCR/qPCR Raw_shared regulon 2_topgenes_MWe.xlsx', sheet = 'relevant_data_CT')

# Little preprocessing
qPCR_Data_Maya$metadata_type = NA
qPCR_Data_Maya$metadata_type[grepl('Control', qPCR_Data_Maya$metadata_annotation)] = 'Ctrl'
qPCR_Data_Maya$metadata_type[!grepl('Control', qPCR_Data_Maya$metadata_annotation)] = 'HCM'

gene_names_noMYL2 = colnames(qPCR_Data_Maya)[!grepl('^metadata|^MYL2$',colnames(qPCR_Data_Maya))]
gene_names_all = colnames(qPCR_Data_Maya)[!grepl('^metadata',colnames(qPCR_Data_Maya))]

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

#ggsave(filename = paste0(base_dir,'Rplots/qPCR_scatters_MYL2.pdf'), 
#        plot = p, width=172-4, height=PANEL_HEIGHT*2, units='mm', device = cairo_pdf) # 184.6/3*2-4
ggsave(filename = paste0(base_dir,'Rplots/qPCR_scatters_MYL2_',gsub('Δ','d',my_unit),'.pdf'), 
        plot = p, width=172-4, height=76, units='mm', device = cairo_pdf) # 184.6/3*2-4

ggsave(filename = paste0(base_dir,'Rplots/qPCR_scatters_MYL2_',gsub('Δ','d',my_unit),'_style2.pdf'), 
        plot = p, width=172-4, height=(172-4)/8*2+15, units='mm', device = cairo_pdf) # 184.6/3*2-4


# Test code

lm_out = lm(formula = MYL2 ~ ACTA1, data = qPCR_Data_Maya) 


# Misc plot
ggplot(qPCR_dataFrame[['qPCR_Data_Maya_minusdCT']], aes(x=metadata_type,y=MYL2, color=metadata_type))+
    geom_jitter()+theme_bw()+theme(legend.position='none')




################################################################################
# Looking at raw data


# Processing qPCR data to check correlations
qPCR_Data_Maya_rawCT =
    openxlsx::read.xlsx('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/qPCR/qPCR Raw_shared regulon 2_topgenes_MWe.xlsx', sheet = 'raw_data')

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

ggsave(filename = paste0(base_dir,'Rplots/qPCR_misc_GAPDH-vs-RPL32.pdf'), 
        plot = p, width=172/3-4, height=172/3-4, units='mm', device = cairo_pdf) # 184.6/3*2-4


p1=ggplot(qPCR_Data_Maya_rawCT, aes(x=metadata_type,y=MYL2, color=metadata_type))+
    geom_jitter()+theme_bw()+theme(legend.position='none')+xlab(element_blank())+give_better_textsize_plot(12)+ylab('MYL2 (CT)')
p2=ggplot(qPCR_Data_Maya_rawCT, aes(x=metadata_type,y=ACTA1, color=metadata_type))+
    geom_jitter()+theme_bw()+theme(legend.position='none')+xlab(element_blank())+give_better_textsize_plot(12)+ylab('ACTA1 (CT)')
p3=ggplot(qPCR_Data_Maya_rawCT, aes(x=metadata_type,y=CRYAB, color=metadata_type))+
    geom_jitter()+theme_bw()+theme(legend.position='none')+xlab(element_blank())+give_better_textsize_plot(12)+ylab('CRYAB (CT)')

p=p1+p2+p3
p

ggsave(filename = paste0(base_dir,'Rplots/qPCR_misc_raw-MYL2-ACTA1-CRYAB.pdf'), 
        plot = p, width=172*2/3-4, height=172/3-4, units='mm', device = cairo_pdf) # 184.6/3*2-4



