

# Project Rooij-specific clusters on all data

# Run this to activate main script dependencies:
# script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

if (!exists('current_analysis')) { current_analysis = list() }

current_analysis[['ALL.SP_RID2l']] = 
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_','ALL.SP_RID2l','.h5seurat'))
current_analysis[['ROOIJonly_RID2l_clExtended']] = 
    LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_','ROOIJonly_RID2l_clExtended','.h5seurat'))

clustering_rooij = Idents(current_analysis$ROOIJonly_RID2l_clExtended)
clustering_rooij_projected = clustering_rooij[colnames(current_analysis$ALL.SP_RID2l)]

current_analysis$ALL.SP_RID2l$clustering_rooij_only_projected = clustering_rooij_projected

# p = DimPlot(current_analysis$ALL.SP_RID2l, group.by = 'clustering_rooij_only_projected')+give_better_textsize_plot(18)+theme_void()

p=DimPlot(current_analysis$ALL.SP_RID2l, group.by = 'clustering_rooij_only_projected', cols = rep(col_vector_60,2), 
          label = F, repel = T, label.size = 7/.pt, pt.size = .1)+
            theme_void()+ggtitle(element_blank())+theme(legend.position = 'bottom', legend.text = element_text(family="Arial", size=7), legend.key.size = unit(8,"pt"))
        
ggsave(plot = p,filename = paste0(base_dir, 'Rplots/','ALL.SP_RID2l','_9_custom_RooijCl_on_ALL.SP.pdf'), 
    #   width = 172/3-4, height= 172/3-4, units='mm', device=cairo_pdf)
    width = 172/3-4, height= 172/3+8/.pt, units='mm', device=cairo_pdf)



