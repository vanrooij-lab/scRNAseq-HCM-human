##########
            # Now same with only Rooij data displayed (could have been made earlier, but OK)
            # Slightly adjusted heatmap in case of more information
            df_melted_sel2$donor_short = gsub('^.\\.','',df_melted_sel2$donor)
            CUSTOM_PATIENT_ORDER_short = gsub('^.\\.','',CUSTOM_PATIENT_ORDER)
            p=ggplot(df_melted_sel2[df_melted_sel2$paper=='R',], aes(x=factor(donor_short, levels=CUSTOM_PATIENT_ORDER_short), y=factor(gene_name_short, levels=rev(gene_order_short)), fill=corr)) +
                geom_tile(size=.5, width=1, height=1)+
                geom_point(data=df_melted_sel2[df_melted_sel2$pval.adj.sign&df_melted_sel2$paper=='R', ])+
                scale_color_manual(values=c('white','black'))+
                scale_fill_gradientn(colours = rainbow_colors, breaks=seq(-1,1,.5), limits=c(-1,1))+
                #geom_text(aes(label=round(corr,2)), color='#666666', size=6/.pt)+
                theme_bw()+ylab(element_blank())+give_better_textsize_plot(8)+
                ggtitle(paste0('Correlations with ',shorthand_cutname(CURRENT_GENE)))+
                theme(legend.position = 'none', legend.key.height = unit(2,"mm"),
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
                xlab('donor')
            # p
            ggsave(filename = paste0(base_dir,'Rplots/customALL',SP_SWITCH,'_6_TableSignCorrelations-style2_',CURRENT_GENE,'_Corr',negpos,'_.pdf'), 
                             plot = p, height=min(172,3*nrow_effective+20), width=2/3*172-4, units = 'mm', device=cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5
            # Now with legend
            p=p+theme(legend.position='bottom') # p
            ggsave(filename = paste0(base_dir,'Rplots/customALL',SP_SWITCH,'_6_TableSignCorrelations-style2_LEGEND_',CURRENT_GENE,'_.pdf'), 
                             plot = p, height=min(172,3.5*nrow_effective+5), width=2/3*172-4, units = 'mm', device=cairo_pdf) # 185/4 || 46*2 || ncol_effective*7.5
            
            