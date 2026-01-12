
################################################################################
# ESSRA

# Huss, J. M., Torra, I. P., Staels, B., Giguère, V., & Kelly, D. P. (2004). Estrogen-Related Receptor α Directs Peroxisome Proliferator-Activated Receptor α Signaling in the Transcriptional Control of Energy Metabolism in Cardiac and Skeletal Muscle. Molecular and Cellular Biology, 24(20), 9079–9091. https://doi.org/10.1128/mcb.24.20.9079-9091.2004
ESSRA_regulated1 = c('Cd36', 'Facl2', 'Acadm', 'Acadvl', 'Hadha', 'Hadha', 'Fabp3', 'Acox1', 'Hsd17b4', 'Scd2', 'Cyb5', 'Aldh1a4', 'Cycs', 'Cox8h', 'ACP', 'Etfdh', 'Mtch2', 'Hk1', 'Pppr2b2', 'Pp2c2', 'Hmgcs1', 'Hmgcr', 'Idi1', 'Lip1', 'Pdp1', 'Ckb', 'Alas1', 'Chk', 'Cyp51', 'Asah', 'B4galt6', 'Got1', 'Csrp2', 'Hmgb2', 'Smarca4', 'Sfrs10', 'Sfrs10', 'Pln', 'Myh6', 'Myh7')

ESSRA_regulated1_humanized = convertMouseGeneList(x = ESSRA_regulated1, the_other_animal = 'rnorvegicus_gene_ensembl')

    # only few get converted

sum(ESSRA_regulated1_humanized %in% SCENIC_reg_top_genes_sorted_full$ESRRA)

    # only one hit 
    #
    # bit of a reach anyhow


