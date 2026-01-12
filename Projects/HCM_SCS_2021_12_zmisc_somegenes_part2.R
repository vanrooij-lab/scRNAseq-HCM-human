FOXN3, MAFK en NFE2l1 

SCENIC_reg_top_genes_sorted_full$FOXN3
SCENIC_reg_top_genes_sorted_full$MAFK
SCENIC_reg_top_genes_sorted_full$NFE2L1

Reduce(intersect, list(SCENIC_reg_top_genes_sorted_full$FOXN3, SCENIC_reg_top_genes_sorted_full$MAFK, SCENIC_reg_top_genes_sorted_full$NFE2L1))

SCENIC_reg_top_genes_sorted_full$FOXN3[SCENIC_reg_top_genes_sorted_full$FOXN3 %in% SCENIC_reg_top_genes_sorted_full$MAFK]
SCENIC_reg_top_genes_sorted_full$MAFK[SCENIC_reg_top_genes_sorted_full$MAFK %in% SCENIC_reg_top_genes_sorted_full$NFE2L1]


intersect(SCENIC_reg_top_genes_sorted_full$FOXN3, SCENIC_reg_top_genes_sorted_full$MAFK)
intersect(SCENIC_reg_top_genes_sorted_full$FOXN3, SCENIC_reg_top_genes_sorted_full$NFE2L1)
intersect(SCENIC_reg_top_genes_sorted_full$MAFK, SCENIC_reg_top_genes_sorted_full$NFE2L1)



# MAFK / NFE2L1
intersect(SCENIC_reg_top_genes_sorted_full$MAFK, SCENIC_reg_top_genes_sorted_full$NFE2L1)
# [1] "MAP4"     "CCND1"    "MYH14"    "EIF3A"    "CCND2"    "TRDN"     "CTNNA3"   "TP53INP2" "ATRX"     "PPP1R12B" "SUPT6H"  

# ZEB1 / NFE2L1
intersect(SCENIC_reg_top_genes_sorted_full$ZEB1, SCENIC_reg_top_genes_sorted_full$NFE2L1)
# [1] "MAP4"   "CCND1"  "CCND2"  "ARMCX3" "ZEB1"  

# MAFK / ZEB1
intersect(SCENIC_reg_top_genes_sorted_full$MAFK, SCENIC_reg_top_genes_sorted_full$ZEB1)
# [1] "MAP4"  "PALLD" "CCND1" "AHNAK" "MFN2"  "CCND2"

# All three
Reduce(intersect, list(SCENIC_reg_top_genes_sorted_full$ZEB1, SCENIC_reg_top_genes_sorted_full$MAFK, SCENIC_reg_top_genes_sorted_full$NFE2L1))
# "MAP4"  "CCND1" "CCND2"