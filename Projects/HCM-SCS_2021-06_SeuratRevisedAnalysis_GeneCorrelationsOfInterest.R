

# Generating correlations to genes of interest for each of the datasets
#
#

source('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Functions/MW_general_functions.R')
# file.edit('/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Functions/MW_general_functions.R')
    # Volcano plot functions (other copies of this function exist in other repos)

# For the integrated dataset
currentSeuratObject_recombined3@assays$integrated@scale.data

################################################################################
# Load the objects

OBJECTS_TO_ANALYZE = c('ROOIJonly_RID2l', 'ROOIJonly_default_Int1c','HUonly_RID2l', 'HUonly_default_Int1c')
INTEGRATED_OR_NOT = c('no'              , 'yes'                     ,'no'                 , 'yes')
names(INTEGRATED_OR_NOT)=OBJECTS_TO_ANALYZE

current_analysis=list()
for (analysis_name in OBJECTS_TO_ANALYZE) {

    # analysis_name = 'HUonly_RID2l'
    
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
        
        # Depending on whether we want to look at integrated data or not, Seurat has stored the
        # expression matrix of interest in a different spot
        assayType = switch(INTEGRATED_OR_NOT[[analysis_name]], 'yes'='integrated', 'no'='RNA')
        
        # Get correlation df
        # Note we can take "data" here, since correlation is scaled by definition
        Volcano_df_collection[[analysis_name]][[gene_name]] = 
            get_volcano_df3(expression_matrix=current_analysis[[analysis_name]]@assays[[assayType]]@data,
                            my_gene=gene_name,calc_qvals=T,no_expression_val=0,min_cell_expressed=.1,manual_expression=NULL)
            
        # Make a plot
        p=plot_volcano3(Volcano_df_collection[[analysis_name]][[gene_name]],mytextsize=15,mycex=3,manual_gene_name=NULL,
                NRLABELED=20,mypvaltreshold=0.01,custom_highlight_group=NULL,NRLABELED_hlgroup=NULL,
                          mypointsize=.1, mylinesize=.25)+
            ggtitle(paste0('Correlations with ',gene_name,'\n(',analysis_name,')'))+
            give_better_textsize_plot(10)
        
        # Save & export it
        ggsave(filename = paste0(base_dir,'Rplots/',analysis_name,'_6_Volcano_',gene_name,'.png'), 
            plot = p, height=7.5, width=7.5, units = 'cm')
        openxlsx::write.xlsx(x = Volcano_df_collection[[analysis_name]][[gene_name]], file = paste0(base_dir,'Rplots/',analysis_name,'_6_Volcano_',gene_name,'.xlsx'))
        
    }
    
}; beepr::beep()

################################################################################

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

################################################################################






