

# Generating correlations to genes of interest for each of the datasets
#
#

print('Starting correlation script ..')

library(Matrix) # required to transpose the sparse matrix, dgCMatrix
library(pheatmap)

source(paste0(script_dir,'Functions/MW_general_functions.R'))
# file.edit('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Functions/MW_general_functions.R')
    # Volcano plot functions (other copies of this function exist in other repos)

# For the integrated dataset
# currentSeuratObject_recombined3@assays$integrated@scale.data

################################################################################
# Load the objects

#OBJECTS_TO_ANALYZE = c('ROOIJonly_RID2l', 'ROOIJonly_default_Int1c','HUonly_RID2l', 'HUonly_default_Int1c')
#INTEGRATED_OR_NOT = c('no'              , 'yes'                     ,'no'                 , 'yes')
#names(INTEGRATED_OR_NOT)=OBJECTS_TO_ANALYZE

OBJECTS_TO_ANALYZE = c('ROOIJonly_RID2l','HUonly_RID2l', 'TEICHMANNonly_RID2l')
INTEGRATED_OR_NOT = c('no'              , 'no'                     ,'no')
names(INTEGRATED_OR_NOT)=OBJECTS_TO_ANALYZE

current_analysis=list()
for (analysis_name in OBJECTS_TO_ANALYZE) {

    # analysis_name = 'HUonly_RID2l'
    
    print(paste0('Loading ', analysis_name))
    
    current_analysis[[analysis_name]] = 
        LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',analysis_name,'.h5seurat'))

}

################################################################################
# Some TTN correlations
# To do: automate this such that it just loops over multiple datasets, patients, and the pooled/integrated data

OBJECTS_TO_ANALYZE=OBJECTS_TO_ANALYZE # defined above; this is just a reminder
GENES_OF_INTEREST=c('TTN','NPPA')

Volcano_df_collection=list()
for (analysis_name in OBJECTS_TO_ANALYZE) {
    
    for (gene_name in GENES_OF_INTEREST) {
        
        # analysis_name='ROOIJonly_RID2l'; gene_name=shorthand_seurat_fullgenename('TTN')
        # analysis_name='ROOIJonly_RID2l'; gene_name=shorthand_seurat_fullgenename('NPPA')
        # analysis_name='ROOIJonly_RID2l'; gene_name=shorthand_seurat_fullgenename('CMYA5')
        # analysis_name='TEICHMANNonly_RID2l'; gene_name='TTN'
        
        # Depending on whether we want to look at integrated data or not, Seurat has stored the
        # expression matrix of interest in a different spot
        assayType = switch(INTEGRATED_OR_NOT[[analysis_name]], 'yes'='integrated', 'no'='RNA')
        
        # Now go over patients
        patients_current_dataset = unique(current_analysis[[analysis_name]]$annotation_patient_str)
        for (current_patient_idx in 0:length(patients_current_dataset)) {
            
            if (current_patient_idx>0) {
                current_patient=patients_current_dataset[current_patient_idx]
                cells_from_patient_selection = (current_analysis[[analysis_name]]$annotation_patient_str == current_patient)
            }else{
                current_patient=paste0(substring(analysis_name, 1,1), 'pooled')
                cells_from_patient_selection=rep(T, length(current_analysis[[analysis_name]]$annotation_patient_str))
            }
            
            # Display progress
            print(paste0('Running analysis for: ',analysis_name,', ',gene_name,', ',current_patient,'..'))
            
            # Get correlation df
            # Note we can take "data" here, since correlation is scaled by definition
            Volcano_df_collection[[analysis_name]][[gene_name]][[current_patient]] = 
                get_volcano_df3(expression_matrix=current_analysis[[analysis_name]]@assays[[assayType]]@data[,cells_from_patient_selection],
                                my_gene=gene_name,calc_qvals=T,no_expression_val=0,min_cell_expressed=.1,manual_expression=NULL)
                
            # Make a plot
            p=plot_volcano3(my_corrs_df_current = Volcano_df_collection[[analysis_name]][[gene_name]][[current_patient]],mycex=3,
                            NRLABELED=20,mypvaltreshold=0.01,manual_gene_name = gene_name,
                              mypointsize=.1, mylinesize=.25, mytextsize=10)+
                ggtitle(paste0('Correlations with ',gene_name,'\n(',analysis_name,'); ',current_patient,''))
            
            # Save & export it
            ggsave(filename = paste0(base_dir,'Rplots/',analysis_name,'_6_Volcano_',gene_name,'_',current_patient,'.png'), 
                plot = p, height=7.5, width=7.5, units = 'cm')
            openxlsx::write.xlsx(x = Volcano_df_collection[[analysis_name]][[gene_name]], file = paste0(base_dir,'Rplots/',analysis_name,'_6_Volcano_',gene_name,'.xlsx'))
            
        }    
    }
    
}#; beepr::beep()

save(list = c('Volcano_df_collection'), file = paste0(base_dir,'Rdata/Volcano_df_collection__for_all.Rdata'))


################################################################################

# custom analysis code to compare patients and datasets

# for (gene_name in GENES_OF_INTEREST) {

if (F) {
    
    load('/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis/Rdata/Volcano_df_collection__for_all.Rdata')
    
    current_gene='TTN'

    collected_genes = lapply(OBJECTS_TO_ANALYZE,
        function(current_dataset) {
            patients_current_dataset = names(Volcano_df_collection[[current_dataset]][[current_gene]])
            lapply(patients_current_dataset, function(current_patient) {
                Volcano_df_collection[[current_dataset]][[current_gene]][[current_patient]]$gene_name[order(Volcano_df_collection[[current_dataset]][[current_gene]][[current_patient]]$corr, decreasing = T)[1:20]]
            })
        })
    collected_genes=unique(unlist(collected_genes))
    
    all_correlations_df=data.frame(gene=numeric(),corr=numeric(),patient=character(),dataset=character())
    for (current_dataset in OBJECTS_TO_ANALYZE) {
    
        # current_dataset='ROOIJonly_RID2l'; gene_name='TTN'; current_patient = 'R.P1'
        
        patients_current_dataset = names(Volcano_df_collection[[current_dataset]][[current_gene]])
        for (current_patient in patients_current_dataset) {

            
            current_data=
                data.frame(gene=collected_genes, 
                           corr=Volcano_df_collection[[current_dataset]][[current_gene]][[current_patient]][collected_genes,]$corr,
                           patient=current_patient, dataset=current_dataset)
                
            all_correlations_df=rbind(all_correlations_df,current_data)
            
        }
        
    }
    
    expand.grid()
    
    ggplot(all_correlations_df, aes(x=patient,y=gene,fill=corr))+
        geom_tile()

    mtx <- matrix(NA, nrow=length(collected_genes), ncol=length(unique(all_correlations_df$patient)) )
    dimnames(mtx) <- list( collected_genes, sort(unique(all_correlations_df$patient) ) )
    mtx[cbind(all_correlations_df$gene, all_correlations_df$patient)] <- all_correlations_df$corr
    
    mtx[is.na(mtx)]=0
    p=pheatmap(mtx, fontsize_row = 3)
    p
    ggsave(filename = paste0(base_dir,'Rplots/combined_TTN_corrs.pdf'), plot = p, height=30, width=15, units='cm')
}
################################################################################

# Old code/for manual running

if (F) {
    
    # NPPA Venn
    plot_Venn_MW_2lists_v3(Volcano_df_collection[['ROOIJonly_RID2l']][['NPPA']]$gene_name_short[2:31],
                            Volcano_df_collection[['ROOIJonly_default_Int1c']][['NPPA']]$gene_name_short[2:31],
                            name1='Rooij, RID2',name2='Rooij, Int')
    
    plot_Venn_MW_2lists_v3(Volcano_df_collection[['ROOIJonly_default_Int1c']][['NPPA']]$gene_name_short[2:31],
                            Volcano_df_collection[['HUonly_default_Int1c']][['NPPA']]$gene_name_short[2:31],
                            name1='Rooij, Int',name2='Hu, Int')
    
    Volcano_df_collection[['ROOIJonly_default_Int1c']][['NPPA']]$gene_name_short[2:31][Volcano_df_collection[['ROOIJonly_default_Int1c']][['NPPA']]$gene_name_short[1:30] %in% Volcano_df_collection[['HUonly_default_Int1c']][['NPPA']]$gene_name_short[2:31]]
    
    # TTN Venn
    plot_Venn_MW_2lists_v3(Volcano_df_collection[['ROOIJonly_RID2l']][['TTN']]$gene_name_short[2:31],
                            Volcano_df_collection[['ROOIJonly_default_Int1c']][['TTN']]$gene_name_short[2:31],
                            name1='Rooij, RID2',name2='Rooij, Int')
    
    plot_Venn_MW_2lists_v3(Volcano_df_collection[['ROOIJonly_default_Int1c']][['TTN']]$gene_name_short[2:31],
                            Volcano_df_collection[['HUonly_default_Int1c']][['TTN']]$gene_name_short[2:31],
                            name1='Rooij, Int',name2='Hu, Int')
}

################################################################################






