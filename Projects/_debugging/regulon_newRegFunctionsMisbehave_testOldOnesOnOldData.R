

library(heatmap3)
library(openxlsx)
library(randomcoloR)
library(scales)

# Using the old code
source('/Users/m.wehrens/Documents/git_repos/SCS_Joep/Functions/generateCorrelationMatrix.R')
source('/Users/m.wehrens/Documents/git_repos/SCS_Joep/Functions/analyzeRegulons.R')
source('/Users/m.wehrens/Documents/git_repos/SCS_Joep/Functions/RaceID2_StemID_class.R')
source('/Users/m.wehrens/Documents/git_repos/SCS_Joep/Functions/functions_mw_copy/MW_analysis_functions.R')


# Load the old patients as they were stored separately
# groupedSCS_SinglePatients
load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/Rdata_previous/groupedSCS_SinglePatients.Rdata')
# Load the old patients in a pooled fashion
load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/Rdata_previous/groupedSCS_PatientAllMod.Rdata')

# # For annotation purposes
# samples_patient1_J = c('JE1', 'JE2', 'JE3', 'JE4')
# samples_patient2_J = c('JE5', 'JE6', 'JE7', 'JE8')
# samples_patient3_J = c('MW5', 'MW6', 'MW7', 'MW8')
# samples_patient4_J = c('JE10', 'JE11')
# samples_patient5_J = c('AL1', 'AL2')
# samples_Rooij_J = list(samples_patient1_J,samples_patient2_J,samples_patient3_J,samples_patient4_J,samples_patient5_J)
# 
# current_patient=1
# selection_pt = sapply(str_split(colnames(hcm_scs@ndata),'_'),function(x){x[[1]]}) %in% samples_Rooij_J[[current_patient]]

# hcm_scs_test=hcm_scs
# hcm_scs_test@ndata=hcm_scs@ndata[,selection_pt]
# groupedSCS=list(test.p1=hcm_scs_test)

# Tets with original stuff
outputDir=paste0(base_dir,'testing3/')
correlationMatrices = generateCorrelationMatrix(config, groupedSCS_SinglePatients, groupNames=c('patient1Mod'), 
          excludeOutlierCells=T, minCellFraction=0.05, 
          minCellExpression=0.1, desiredPValue=0.00001, adjustP=T, saveMatrices=F, 
          overrideDir=F, outputMode='pdf')
# Determine the regulons    
regulons <- analyzeRegulons(config, correlationMatrices, minCorrelations=40, 
    clusteringMethod='ward.D2', overrideClusterNum=F, useSquish=0.01, 
    overrideDir='Fig4/A', outputMode='pdf', K.max = 50)
# Throw another "analyzeRegulons over it?

save(list = 'correlationMatrices', file =paste0(base_dir,'Rdata/ancient_correlationMatrices.Rdata'))

########################################################################
# So I think the discrepency with earlier results lies in an additional
# filter I applied

outputDir=paste0(base_dir,'testing4_extrafilter/')

# First a little complication, determine genes that are expressed in all 
# patients ..
MINPERCENTAGCELLSEEXPRESSED=.2
all_genes = unique(unlist(sapply(groupedSCS_SinglePatients, function(x) {rownames(x@ndata)} )))
genes_detected = strip__chrXX(rownames(groupedSCS_PatientAllMod$patientAllMod@ndata[apply(groupedSCS_PatientAllMod$patientAllMod@ndata>.1, 1, mean)>MINPERCENTAGCELLSEEXPRESSED,]))
genes_detected_all_patients_mincelltreshold = genes_detected

# Create a manual matrix
#all_medium_expressed_genes = unique(c(genes_pt1, genes_pt2, genes_pt3, genes_pt4, genes_pt5))
p1_yes = as.double(genes_detected %in% strip__chrXX(rownames(groupedSCS_SinglePatients$patient1Mod@ndata)))
p2_yes = as.double(genes_detected %in% strip__chrXX(rownames(groupedSCS_SinglePatients$patient2Mod@ndata)))
p3_yes = as.double(genes_detected %in% strip__chrXX(rownames(groupedSCS_SinglePatients$patient3Mod@ndata)))
p4_yes = as.double(genes_detected %in% strip__chrXX(rownames(groupedSCS_SinglePatients$patient4Mod@ndata)))
p5_yes = as.double(genes_detected %in% strip__chrXX(rownames(groupedSCS_SinglePatients$patient5Mod@ndata)))

manual_upset_matrix = data.frame(p1=p1_yes, p2=p2_yes, p3=p3_yes, p4=p4_yes, p5=p5_yes)
rownames(manual_upset_matrix) = genes_detected # all_medium_expressed_genes
manual_upset_matrix_genesexpressedpatients = manual_upset_matrix

pheatmap(1*manual_upset_matrix)

# Now take those genes (this line just renames them w. the chromosome suffix)
all_genes_patients = rownames(manual_upset_matrix_genesexpressedpatients)[apply(manual_upset_matrix_genesexpressedpatients,1,all)]
all_genes_patients_chrXX = find__chrXX(all_genes_patients, rownames(groupedSCS_PatientAllMod$patientAll@ndata))

# Calculate corr matrix based on lower p-val and only genes expressed in all patients
correlationMatrices_selectedGenes = generateCorrelationMatrix(config, groupedSCS_SinglePatients, 
            which_genes_to_select = all_genes_patients_chrXX, 
          groupNames=c('patient1Mod'), 
          #groupNames=c('patient1Mod', 'patient2Mod', 'patient3Mod', 'patient4Mod', 'patient5Mod'), 
          excludeOutlierCells=T, minCellFraction=0.05, 
          minCellExpression=0.1, desiredPValue=0.001, adjustP=T, 
          saveMatrices=F, overrideDir=F, outputMode='pdf', 
          filename_precursor='genes_allPs')
regulons_selectedGenes <- analyzeRegulons(config, correlationMatrices_selectedGenes, 
    minCorrelations=10, 
    clusteringMethod='ward.D2', overrideClusterNum=F, useSquish=0.01, 
    overrideDir='Fig4/A', outputMode='pdf', 
    K.max = 50, filename_precursor='genes_allPs')

########################################################################

# Tets with original stuff, now with adjusted parameter
outputDir=paste0(base_dir,'testing2_adjPARAM/')
correlationMatrices_adjPARAM = generateCorrelationMatrix(config, groupedSCS_SinglePatients, groupNames=c('patient1Mod'), 
          excludeOutlierCells=T, minCellFraction=0.05, 
          minCellExpression=0.1, desiredPValue=1e-7, adjustP=T, saveMatrices=F, 
          overrideDir=F, outputMode='pdf')
# Determine the regulons    
regulons_adjPARAM <- analyzeRegulons(config, correlationMatrices_adjPARAM, minCorrelations=40, 
    clusteringMethod='ward.D2', overrideClusterNum=F, useSquish=0.01, 
    overrideDir='Fig4/A', outputMode='pdf', K.max = 50)

save(list = 'correlationMatrices', file =paste0(base_dir,'Rdata/ancient_correlationMatrices.Rdata'))

########################################################################

# Tets with original stuff, CHANGING THE CONNECTIVITY
outputDir=paste0(base_dir,'testing2_adjPARAM_CONN/')
correlationMatrices_adjPARAM_CONN = generateCorrelationMatrix(config, groupedSCS_SinglePatients, groupNames=c('patient1Mod'), 
          excludeOutlierCells=T, minCellFraction=0.05, 
          minCellExpression=0.1, desiredPValue=1e-5, adjustP=T, saveMatrices=F, 
          overrideDir=F, outputMode='pdf')
# Determine the regulons    
regulons_adjPARAM_CONN <- analyzeRegulons(config, correlationMatrices_adjPARAM_CONN, 
    minCorrelations=60, 
    clusteringMethod='ward.D2', overrideClusterNum=F, useSquish=0.01, 
    overrideDir='Fig4/A', outputMode='pdf', K.max = 50)

save(list = 'correlationMatrices', file =paste0(base_dir,'Rdata/ancient_correlationMatrices.Rdata'))

dim(regulons_adjPARAM_CONN$patient1Mod$subCorrelationMatrix)
pheatmap(regulons_adjPARAM_CONN$patient1Mod$subCorrelationMatrix, clustering_method = 'ward.D2')

########################################################################

# let's see whether connectivity is involved here
# First we need clust
genes_subM = rownames(regulons$patient1Mod$subCorrelationMatrix)
hclust_out = hclust(as.dist(1-regulons$patient1Mod$subCorrelationMatrix), method='ward.D2')
p=pheatmap(regulons$patient1Mod$subCorrelationMatrix, fontsize = 1, cluster_cols = hclust_out, cluster_rows = hclust_out)
p
ggsave(filename = paste0(base_dir,'testing2/patient1matrix.png'),plot = p, width=40, height=40, units='cm', dpi = 600)
ggsave(filename = paste0(base_dir,'testing2/patient1matrix.pdf'),plot = p)
plot_df = data.frame(count=correlationMatrices$patient1Mod$correlationsPerGene[
                    colnames(regulons$patient1Mod$subCorrelationMatrix)[hclust_out$order]],
                    x=hclust_out$order, 
                    reg=as.factor(regulons$patient1Mod$clustering[hclust_out$order]))
ggplot(plot_df)+
    geom_bar(aes(x=x,y=count,fill=reg), stat='identity')+theme_bw()+
    scale_fill_manual(values = col_vector_60)

# Let's also check out p-vals
# Note: this of course scales with the strength of the correlations, so it's not so relevant ofr
# current question
rownames(correlationMatrices$patient1Mod$pValueMatrix) = rownames(correlationMatrices$patient1Mod$correlationMatrix)
colnames(correlationMatrices$patient1Mod$pValueMatrix) = colnames(correlationMatrices$patient1Mod$correlationMatrix)
rownames(correlationMatrices$patient1Mod$filteredPValueMatrix) = rownames(correlationMatrices$patient1Mod$correlationMatrix)
colnames(correlationMatrices$patient1Mod$filteredPValueMatrix) = colnames(correlationMatrices$patient1Mod$correlationMatrix)

p2=pheatmap(correlationMatrices$patient1Mod$pValueMatrix[genes_subM,genes_subM], fontsize = 1, cluster_cols = hclust_out, cluster_rows = hclust_out)
p2

ggplot(data.frame(R=as.vector(regulons$patient1Mod$subCorrelationMatrix)))+
    geom_freqpoly(aes(x=R))+theme_bw()

pheatmap((regulons$patient1Mod$subCorrelationMatrix>.35)*1, fontsize = 1, cluster_cols = hclust_out, cluster_rows = hclust_out)
pheatmap((regulons$patient1Mod$subCorrelationMatrix>.5)*1, fontsize = 1, cluster_cols = hclust_out, cluster_rows = hclust_out)
# with pvals
pheatmap((correlationMatrices$patient1Mod$filteredPValueMatrix[genes_subM,genes_subM]<10^-6)*1, fontsize = 1, cluster_cols = hclust_out, cluster_rows = hclust_out)


ggplot(data.frame(connectivity=rowSums((regulons$patient1Mod$subCorrelationMatrix>.35))))+
    geom_histogram(aes(x=connectivity))+theme_bw()

# Look at effect of p-val on remaining genes
# Note: the amount of genes is determined by using the originla matrix, of course
subPvals=correlationMatrices$patient1Mod$pValueMatrix#[genes_subM,genes_subM]
resultingGenesLeftForPCutoffs = sapply(0:10, function(x) {sum(rowSums(subPvals<10^-x)>40)})
ggplot(data.frame(pcut=0:10, genes_left=resultingGenesLeftForPCutoffs))+
    geom_line(aes(x=pcut, y=genes_left))+theme_bw()+
    ylim(  c(0,dim(correlationMatrices$patient1Mod$pValueMatrix)[2]+1)  )+
    geom_hline(yintercept = 200)+
    geom_vline(xintercept = c(5,6,7), color='red')+
    scale_y_continuous(trans = 'log10')

ggplot(data.frame(pval=as.vector(correlationMatrices$patient1Mod$filteredPValueMatrix[genes_subM,genes_subM])))+
    geom_histogram(aes(x=pval), bins=100)+theme_bw()#+
    #scale_x_continuous(trans = 'log10')
ggplot(data.frame(pval=as.vector(correlationMatrices$patient1Mod$pValueMatrix[genes_subM,genes_subM])))+
    geom_freqpoly(aes(x=-log10(pval+.1)))+theme_bw()+
    geom_vline(xintercept = 5)+
    geom_vline(xintercept = 6, color='red')

# I guess the question is: how many genes are still in the matrix when 
# for different cutoffs?
# (Note that this was also analyzed in one of the decision plots)






