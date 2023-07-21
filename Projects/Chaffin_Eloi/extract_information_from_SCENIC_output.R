

# Note that the analysis initially came up with empty output files, 
# (don't remember precisely but i think because library installation issues -- check this)




# Initialize a parameter for later
scenic_regulons_collected_all_patients_list=list()
scenic_regulons_collected_all_patients_list_regnames=list()

# This can be made into a loop, for now I do one patient
CURRENT_PATIENT = 'P1622'

# Connect to the data file for this patient
myLoomConnection = loomR::connect(filename=paste0('/Users/m.wehrens/Desktop/SCENIC_out_',CURRENT_PATIENT,'.loom'), mode = 'r')

# Retrieve the data
df_regulons_compare = as.data.frame(myLoomConnection$row.attrs$Regulons[])
row.names(df_regulons_compare) = myLoomConnection$row.attrs$Gene[]

# Re-organize the data in a more convenient form
scenic_regulons_collected_all_patients_list[[CURRENT_PATIENT]] = 
    lapply(1:ncol(df_regulons_compare), function(reg_idx) {rownames(df_regulons_compare)[df_regulons_compare[,reg_idx]==1]})
names(scenic_regulons_collected_all_patients_list[[CURRENT_PATIENT]]) = colnames(df_regulons_compare) # paste0(CURRENT_PATIENT,'.',colnames(df_regulons_compare))

# Also make a separate list of just the transcription factors that were identified
scenic_regulons_collected_all_patients_list_regnames[[CURRENT_PATIENT]] = colnames(df_regulons_compare)
