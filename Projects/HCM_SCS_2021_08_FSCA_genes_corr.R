

# Also required:
core_regulons_sorted
core_regulons



# Load index data
plateJE10 = read.csv('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_analyzed_by_Joep/2020_04_AnneJoepData/Copy-of-data-VersionAnne/IndexDatafiles/2018-244 exp cardiologie_plate2.csv', header=T, sep=",", row.names=NULL, skip=13)
plateJE11 = read.csv('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_analyzed_by_Joep/2020_04_AnneJoepData/Copy-of-data-VersionAnne/IndexDatafiles/2018-244 exp cardiologie_plate3merged.csv', header=T ,sep=",", row.names=NULL, skip=13)
plateAL1 = read.csv('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_analyzed_by_Joep/2020_04_AnneJoepData/Copy-of-data-VersionAnne/IndexDatafiles/AL1.csv', header=T ,sep=",", row.names=NULL, skip=15)
plateAL2 = read.csv('/Users/m.wehrens/Data/_2019_02_HCM_SCS/_analyzed_by_Joep/2020_04_AnneJoepData/Copy-of-data-VersionAnne/IndexDatafiles/AL2.csv', header=T ,sep=",", row.names=NULL, skip=15)
# Fix well names
plateJE10$Well = unlist(lapply(plateJE10$Well, function(x) {return(paste0('JE10_',substr(x,1,1),'_',substr(x,2,3)))}))
plateJE11$Well = unlist(lapply(plateJE11$Well, function(x) {return(paste0('JE11_',substr(x,1,1),'_',substr(x,2,3)))}))
plateAL1$Well = unlist(lapply(plateAL1$Well, function(x) {return(paste0('AL01_',substr(x,1,1),'_',substr(x,2,3)))}))
plateAL1 = plateAL1[colnames(plateJE10)] # Throw away ridiculous column overload
plateAL2$Well = unlist(lapply(plateAL2$Well, function(x) {return(paste0('AL02_',substr(x,1,1),'_',substr(x,2,3)))}))
plateAL2 = plateAL2[colnames(plateJE10)] # Throw away ridiculous column overload

indexData_MW = rbind(plateAL1, plateAL2, plateJE10, plateJE11)
# Fix column names
colnames(indexData_MW) <- c('Cell','Events','Parent','FSC_A','FSC_W','FSC_H','SSC_A','SSC_W','SSC_H','BV421_A')
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


################################################################################
# Project data on UMAP

current_analysis$ROOIJonly_RID2l[['FSCA']]=indexData_MW[colnames(current_analysis$ROOIJonly_RID2l),]$FSC_A

FeaturePlot(object = current_analysis$ROOIJonly_RID2l, features = 'FSCA', cols = rainbow_colors, pt.size = 2)

################################################################################

# Now perform the correlations

GENE_MIN_PERCENTAGE=.1
GENE_MIN_PERCENTAGE=.33

# First subset on where we have data for
selected_cells = colnames(current_analysis$ROOIJonly_RID2l)[colnames(current_analysis$ROOIJonly_RID2l) %in% indexData_MW$Cell_newname]
current_analysis$ROOIJonly_RID2l_FSCA = subset(current_analysis$ROOIJonly_RID2l, cells = selected_cells)
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

p=ggplot(correlations_FSCA_per_patient_combined, aes(x=P4_cor, y=P5_cor))+
    geom_point(size=.5, shape=1, color='grey')+
    geom_point(size=.5, shape=1, data=correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$P4_cor>.2&correlations_FSCA_per_patient_combined$P5_cor>.1,], color='red')+
    geom_text_repel(data=correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$P4_cor>.2&correlations_FSCA_per_patient_combined$P5_cor>.1,], color='red',
                    mapping=aes(label=gene_symbol), 
                    max.overlaps = Inf, min.segment.length = 0, force = 4, size=6/.pt, segment.size=.25)+
    theme_bw()+give_better_textsize_plot(8)+
    xlim(c(min(correlations_FSCA_per_patient_combined$P4_cor)*1.1,max(correlations_FSCA_per_patient_combined$P4_cor)*1.3))+
    ylim(c(min(correlations_FSCA_per_patient_combined$P5_cor)*1.1,max(correlations_FSCA_per_patient_combined$P5_cor)*1.3))+
    xlab('Gene-cell size correlation patient 4')+ylab('Gene-cell size correlation patient 5')
    #geom_smooth(method = 'lm')+ylim(c(-.5,.5))+xlim(c(-.5,.5))
p
ggsave(plot = p,filename = paste0(base_dir, 'Rplots/ROOIJonly_RID2l_8_FSCA_corrs_Above02.pdf'), 
       width = 184.6/3-4, height= 184.6/3-4, units='mm')
ggsave(plot = p,filename = paste0(base_dir, 'Rplots/ROOIJonly_RID2l_8_FSCA_corrs_Above02_L.pdf'), 
       width = 184.6/2-4, height= 184.6/2-4, units='mm')

# Now with linear fit
# (Important genes (ie. regulons) are consistently high in both patients, but in general 
# patient consistency is not reflected in linear fit)
p=ggplot(correlations_FSCA_per_patient_combined, aes(x=P4_cor, y=P5_cor))+
    geom_point(size=.5, shape=1, color='grey')+
    geom_point(size=.5, shape=1, data=correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$P4_cor>.2&correlations_FSCA_per_patient_combined$P5_cor>.1,], color='red')+
    geom_text_repel(data=correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$P4_cor>.2&correlations_FSCA_per_patient_combined$P5_cor>.1,], color='red',
                    mapping=aes(label=gene_symbol), 
                    max.overlaps = Inf, min.segment.length = 0, force = 4, size=6/.pt, segment.size=.25)+
    theme_bw()+give_better_textsize_plot(8)+
    xlim(c(min(correlations_FSCA_per_patient_combined$P4_cor)*1.1,max(correlations_FSCA_per_patient_combined$P4_cor)*1.3))+
    ylim(c(min(correlations_FSCA_per_patient_combined$P5_cor)*1.1,max(correlations_FSCA_per_patient_combined$P5_cor)*1.3))+
    xlab('Gene-cell size correlation patient 4')+ylab('Gene-cell size correlation patient 5')+
    geom_smooth(method = 'lm')+ylim(c(-.5,.5))+xlim(c(-.5,.5))
p
ggsave(plot = p,filename = paste0(base_dir, 'Rplots/ROOIJonly_RID2l_8_FSCA_corrs_Above02_withfit.pdf'), 
       width = 184.6/3-4, height= 184.6/3-4, units='mm')
ggsave(plot = p,filename = paste0(base_dir, 'Rplots/ROOIJonly_RID2l_8_FSCA_corrs_Above02_L_withfit.pdf'), 
       width = 184.6/2-4, height= 184.6/2-4, units='mm')

# Now show some regulons
regulon_membership = 
    Reduce(f=rbind, x=
        lapply(names(core_regulons[sapply(core_regulons, length)>0]), function(n) {data.frame(gene=core_regulons[[n]], regulon=n)}))
membership_mapping = regulon_membership$regulon
names(membership_mapping) = regulon_membership$gene
correlations_FSCA_per_patient_combined$regulon_membership = membership_mapping[correlations_FSCA_per_patient_combined$gene]
ggplot(correlations_FSCA_per_patient_combined, aes(x=P4_cor, y=P5_cor))+
    geom_point(shape=1,alpha=.1)+
    geom_density_2d(color='black')+
    geom_point(data=correlations_FSCA_per_patient_combined[!is.na(correlations_FSCA_per_patient_combined$regulon_membership),], aes(color=regulon_membership))+theme_bw()
    #geom_point(data=correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$P4_cor>.1&correlations_FSCA_per_patient_combined$P5_cor>.1,], )
    #geom_text_repel(data=correlations_FSCA_per_patient_combined[correlations_FSCA_per_patient_combined$P4_cor>.1&correlations_FSCA_per_patient_combined$P5_cor>.1,], color='red',
    #                mapping=aes(label=gene_symbol))+theme_bw()

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
ggsave(plot = p,filename = paste0(base_dir, 'Rplots/ROOIJonly_RID2l_8_FSCA_corrs_reg2.pdf'), 
       width = 184.6/3-4, height= 184.6/3-4, units='mm')
ggsave(plot = p,filename = paste0(base_dir, 'Rplots/ROOIJonly_RID2l_8_FSCA_corrs_reg2_L.pdf'), 
       width = 184.6/2-4, height= 184.6/2-4, units='mm')
    


    #geom_density_2d(aes(color=regulon_membership), bins=3)+
    #geom_point(data=correlations_FSCA_per_patient_combined[!is.na(correlations_FSCA_per_patient_combined$regulon_membership),], aes(color=regulon_membership))+theme_bw()










