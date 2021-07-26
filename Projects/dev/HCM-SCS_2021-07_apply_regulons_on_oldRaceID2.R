
source('/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/HCM_SCS_2021_06_SeuratRevised_Regulon_v3.R')
library(Matrix)

load_dir='/hpc/hub_oudenaarden/mwehrens/data/Wang/Rdata/'
load(paste0(load_dir, 'export_groupedSCS_patientAllMod_only.Rdata'))
colnames(hcm_scs@ndata)

# For annotation purposes
samples_patient1_J = c('JE1', 'JE2', 'JE3', 'JE4')
samples_patient2_J = c('JE5', 'JE6', 'JE7', 'JE8')
samples_patient3_J = c('MW5', 'MW6', 'MW7', 'MW8')
samples_patient4_J = c('JE10', 'JE11')
samples_patient5_J = c('AL1', 'AL2')
samples_Rooij_J = list(samples_patient1_J,samples_patient2_J,samples_patient3_J,samples_patient4_J,samples_patient5_J)

# Just checking
sum(sapply(str_split(colnames(hcm_scs@ndata),'_'),function(x){x[[1]]}) %in% samples_Rooij_J[[1]])
sum(sapply(str_split(colnames(hcm_scs@ndata),'_'),function(x){x[[1]]}) %in% samples_Rooij_J[[2]])
sum(sapply(str_split(colnames(hcm_scs@ndata),'_'),function(x){x[[1]]}) %in% samples_Rooij_J[[3]])
sum(sapply(str_split(colnames(hcm_scs@ndata),'_'),function(x){x[[1]]}) %in% samples_Rooij_J[[4]])
sum(sapply(str_split(colnames(hcm_scs@ndata),'_'),function(x){x[[1]]}) %in% samples_Rooij_J[[5]])

ANALYSIS_NAME='original_HCM_SCS_data_sel'

current_analysis=list()

# Convert to usual matrix
whole_matrix =  Matrix(as.matrix(hcm_scs@ndata)-.1, sparse = TRUE)

# First execute the initial filter
gene_selection=rowMeans(whole_matrix>0)>.2
whole_matrix_filtered=whole_matrix[gene_selection,]

# Version that is parallelized
collected_regulon_objects = list()
collected_regulon_objects[[ANALYSIS_NAME]] =
    mclapply(X = 1:5, FUN = function(current_patient) {
            
            print(paste0('Starting patient ',current_patient))
            
            selection_pt =
                sapply(str_split(colnames(whole_matrix_filtered),'_'),function(x){x[[1]]}) %in% samples_Rooij_J[[current_patient]]
            
            # Convert to sparse matrix before starting .. 
            current_matrix = 
                whole_matrix_filtered[,selection_pt]
                #Matrix(as.matrix(whole_matrix_filtered[,selection_pt])-.1, sparse = TRUE)
                #hcm_scs@ndata[,selection_pt]
            
            return( 
                giveMeRegulons_SeuratVersion(run_name=paste0(ANALYSIS_NAME,'_p.',current_patient),
                    base_dir=base_dir,current_matrix=current_matrix, strip_object = T, 
                    MAX_GENES = MAX_GENES)
            )
    
        }, mc.cores = MYMCCORES)
names(collected_regulon_objects[[ANALYSIS_NAME]]) = paste0('p.',1:5)

# Save analysis outcome
save(list='collected_regulon_objects', file = paste0(base_dir,'Rdata/',ANALYSIS_NAME,'_regulons_per_patient.Rdata'))

