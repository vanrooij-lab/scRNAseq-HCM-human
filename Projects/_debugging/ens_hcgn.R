

sapply(str_split(rownames(current_analysis[[CURRENT_RUNNAME]])),function(x){x[[2]]})

bla = shorthand_cutname( rownames(current_analysis[[CURRENT_RUNNAME]]) , PART1OR2 = 2)

length(bla)
length(unique(bla))

table_bla = table(bla)

table_bla[order(table_bla, decreasing = T)][1:20]

affected_genes = names(table_bla[table_bla>1])

names_full = rownames(current_analysis[[CURRENT_RUNNAME]])
names_p1 = shorthand_cutname( names_full , PART1OR2 = 1)
affected_genes_fullname = sort(names_full[names_p1 %in% affected_genes])