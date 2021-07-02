
# Some plots showing regulon 6 expression on the UMAP
DimPlot(current_analysis$ROOIJonly_RID2l)
FeaturePlot(current_analysis$ROOIJonly_RID2l, features = c('CRYAB', 'MYL2', 'TNNC1', 'TNNI3', 'TPM1', 'ACTA1', 'ACTC1', 'MYL3', 'MB'),
    cols = rainbow_colors)
VlnPlot(current_analysis$ROOIJonly_RID2l, features = 'MYL2')    

FeaturePlot(current_analysis$ROOIJonly_RID2l, features = c('RYR2'))

####

# composite expression of regulon 1
s.Reg.1.Names = c('CFLAR', 'MYH7B', 'HOOK2', 'DDX17', 'PABPN1', 'TNNT1', 'LUC7L3', 'WSB1', 'SLC25A36', 'NKTR', 'PNISR', 'SLC2A11', 'CACNA1C', 'ASAP1', 'EXOG', 'RNF207', 'PCSK7', 'ACOT11', 'DPY19L2', 'PDE7A', 'N4BP2L2', 'NEAT1', 'MALAT1', 'AC010680.1', 'POLR2J3')
current_analysis$ROOIJonly_RID2l$s.Reg.1 = 
    colSums(
        current_analysis$ROOIJonly_RID2l@assays$RNA@data[s.Reg.1.Names,]/
            rowSums(current_analysis$ROOIJonly_RID2l@assays$RNA@data[s.Reg.1.Names,])
        )
FeaturePlot(current_analysis$ROOIJonly_RID2l, features = 's.Reg.1', cols = rainbow_colors)



# s.reg 3+4
s.Reg.34.Names = c('RPS20', 'ACTB', 'RPS12', 'B2M', 'RPL13A', 'HLA-C', 'EEF1A1', 'TMSB4X', 'PTMA', 'HLA-B', 'RPL14', 'RPL41')
current_analysis$ROOIJonly_RID2l$s.Reg.34 = 
    colSums(
        current_analysis$ROOIJonly_RID2l@assays$RNA@data[s.Reg.34.Names,]/
            rowSums(current_analysis$ROOIJonly_RID2l@assays$RNA@data[s.Reg.34.Names,])
        )
FeaturePlot(current_analysis$ROOIJonly_RID2l, features = 's.Reg.34', cols = rainbow_colors)
#FeaturePlot(current_analysis$ROOIJonly_RID2l, features = s.Reg.34.Names)




# composite expression of regulon 5
s.Reg.5.Names = c('MAP4', 'TNS1', 'NFE2L1', 'MYH7', 'SORBS1', 'MYOM1', 'ZNF106', 'ANKRD1', 'TTN', 'XIRP2', 'CMYA5', 'TCAP', 'DES', 'SYNM', 'CSDE1', 'HIPK2', 'KIF1C', 'NES', 'SYNPO', 'SVIL', 'TECRL')
current_analysis$ROOIJonly_RID2l$S.Reg.5 = 
    colSums(current_analysis$ROOIJonly_RID2l@assays$RNA@data[s.Reg.5.Names,])
FeaturePlot(current_analysis$ROOIJonly_RID2l, features = 'S.Reg.5', cols = rainbow_colors)
FeaturePlot(current_analysis$ROOIJonly_RID2l, features = s.Reg.5.Names, cols = rainbow_colors)




# composite expression of regulon 6
s.Reg.6.names =c('CRYAB', 'MYL2', 'TNNC1', 'TNNI3', 'TPM1', 'ACTA1', 'ACTC1', 'MYL3', 'MB')    
#s.Reg.6.names =c('MYL2', 'TNNI3', 'ACTC1', 'SLC25A4', 'CRYAB', 'MYL3', 'TNNC1', 'ACTA1', 'MYL12A', 'COX6A2', 'CSRP3', 'HSPB1', 'CKM', 'ATP5J', 'PLN', 'MB', 'GAPDH', 'TPM1', 'NDUFA4', 'MDH1')
s.Reg.6.names = s.Reg.6.names[s.Reg.6.names %in% rownames(current_analysis$ROOIJonly_RID2l@assays$RNA@data)]
current_analysis$ROOIJonly_RID2l$S.Reg.6 = 
    colSums(
        current_analysis$ROOIJonly_RID2l@assays$RNA@data[s.Reg.6.names,]/
            rowSums(current_analysis$ROOIJonly_RID2l@assays$RNA@data[s.Reg.6.names,])
        )
#FeaturePlot(current_analysis$ROOIJonly_RID2l, features = 'S.Reg.6')
FeaturePlot(current_analysis$ROOIJonly_RID2l, features = 'S.Reg.6', cols = rainbow_colors)





# Old regulon 4
s.Reg.old4.Names = c('GJA1', 'LPL', 'CASQ2', 'HPSA5', 'ASAH1', 'HHATL', 'BSG', 'ATP1B1', 'GPNMB', 'CD36')
s.Reg.old4.Names = s.Reg.old4.Names[s.Reg.old4.Names %in% rownames(current_analysis$ROOIJonly_RID2l)]
current_analysis$ROOIJonly_RID2l$s.Reg.old4 = 
    colSums(
        current_analysis$ROOIJonly_RID2l@assays$RNA@data[s.Reg.old4.Names,]/
            rowSums(current_analysis$ROOIJonly_RID2l@assays$RNA@data[s.Reg.old4.Names,])
        )
FeaturePlot(current_analysis$ROOIJonly_RID2l, features = 's.Reg.old4', cols = rainbow_colors)







FeaturePlot(current_analysis$ROOIJonly_RID2l, features = c('MALAT1','KCNQ1OT1'), cols = rainbow_colors)


rownames(current_analysis$ROOIJonly_RID2l)[grepl('KCNQ',rownames(current_analysis$ROOIJonly_RID2l))]





