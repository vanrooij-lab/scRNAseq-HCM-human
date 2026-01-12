

ANALYSIS_NAME = 'ROOIJonly.sp.bt_RID2l_clExtended'

SETNAME='Module4'
MOD4_genes = sort(c('CMYA5', 'XIRP2', 'ZNF106',
                    'MAP4','NRAP','ANKRD1',
                    'TTN','HIPK2','DES'))
p_combined=wrap_plots(lapply(MOD4_genes, 
                        function(x) {
                                p=shorthand_seurat_custom_expr(current_analysis[[ANALYSIS_NAME]], x, textsize=8, pointsize=.3)
                                return(p)
                            }
                        ), 
                 nrow=3)*theme(legend.position='none')
p_combined
ggsave(plot=p_combined, filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_7_GenesOfInterest_',SETNAME,'.pdf'), 
                height=50, width=50, units='mm', device = cairo_pdf)
