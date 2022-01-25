
################################################################################

# Scenic script version 2
# ===
# Using Python/command line on HPC
# This script just prepares loom files with expression matrices of the respective patients
# Note that matrices will be filtered, but not log-transformed 
# I'm following the SCENIC example/walkthrough PBMC10k_SCENIC-protocol-CLI here.
#
# See also https://github.com/aertslab/SCENICprotocol/tree/master/notebooks
# (HTML versions available)
# Locally: /Users/m.wehrens/Documents/Naslag/documentation_code/SCENIC/PBMC10k_SCENIC-protocol-CLI.ipynb
# (Opened with anaconda // Jupyter Notebook)

# NOTE: @HPC, this needs to be run in the base_plusloom conda environment ..
# conda activate base_plusloom

################################################################################
# Load my own main script, to be able to use my Seurat data sets
#
# Note: this also sets "base_dir", I will output to 
# base_dir/SCENIC
#
# Or use following one-liner @HPC
# script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'; desired_command='dummy'; source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R')); rm('desired_command')

# set script dir (note: overwritten by loading SeuratRevisedAnalysis below)
if (exists('LOCAL')) {
    script_dir = '/Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Projects/'
} else {
    script_dir = '/hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/'
}

desired_command='dummy'
source(paste0(script_dir, 'HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R'))

# Also get some fns from my own regulons script to do post-processing
desired_command_regulon='dummy'; source(paste0(script_dir, if (exists('LOCAL')) { 'Projects/' } else { '' }, 'HCM_SCS_2021_06_SeuratRevised_Regulon_v3.R'))


################################################################################
# Libs/more config

library(pheatmap)
library(loomR)

# library(rowr) # has cbind.fill

# DATASET_NAME='ROOIJonly_RID2l'
# DATASET_NAME='TEICHMANN.SP.only_RID2l'
# DATASET_NAME='HUonly_RID2l'

ANALYSIS_NAME = 'ROOIJonly.sp.bt_RID2l'

DATASET_NAMES = c(ANALYSIS_NAME, 
                  'TEICHMANNonly.sp.bt_RID2l', 
                  'HUonly.sp.bt_RID2l')
PATIENT_SETS = list(paste0("R.P", 1:5), 
                    c("T.H5", "T.H6", "T.H3", "T.H2", "T.H7", "T.H4", "T.D2", "T.D3", "T.D4", "T.D5", "T.D6", "T.D7", "T.D11"), # note T.D1 disappeared in sep only
                    c('H.N1', 'H.N2', 'H.N3', 'H.N4', 'H.N5', 'H.N13', 'H.N14'))

################################################################################
# << PRE PROCESSING >>

# NOTE: @HPC, this needs to be run in the base_plusloom conda environment ..
# conda activate base_plusloom
# (Probably It's safe to install loom also in the base environment)
print("NOTE: @HPC, this needs to be run in the base_plusloom conda environment ..")

# Export patients for analysis

if (!dir.exists(paste0(base_dir, 'Ldata/'))) { dir.create(paste0(base_dir, 'Ldata/')) }

for (idx in 1:3) {

    DATASET_NAME = DATASET_NAMES[idx]
    ALL_PATIENTS = PATIENT_SETS[[idx]]
        
    # Rooij patients
    # ALL_PATIENTS = paste0("R.P", 1:5)
    # Teichmann healthy donors
    # ALL_PATIENTS = c('T.H5', 'T.H6', 'T.H3', 'T.H2', 'T.H7', 'T.H4', 'T.D1', 'T.D2', 'T.D3', 'T.D4', 'T.D5', 'T.D6', 'T.D7', 'T.D11')
    # Hu healthy donors
    # ALL_PATIENTS = c('H.N1', 'H.N2', 'H.N3', 'H.N4', 'H.N5', 'H.N13', 'H.N14')
    
    for (CURRENT_PATIENT in ALL_PATIENTS) {
        
        # CURRENT_PATIENT = 'T.H5'
        
        # Load respective dataset
        if (!exists('current_analysis')) {current_analysis = list()}
        current_set = paste0(CURRENT_PATIENT, 'RID2l')
        current_analysis[[current_set]] = LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',current_set,'.h5seurat'))
        
        print(paste0('Converting ',current_set,' ..'))
        
        # Retrieve data from Seurat object 
        # ===
        exprMat_ <- current_analysis[[current_set]]@assays$RNA@counts
        # Simply throw out genes that have duplicate or unknown hgnc names (NA probably not interesting, duplicate symbols are very few, seem uninteresting)
        rownames_exprMat_ = shorthand_cutname(rownames(exprMat_), PART1OR2 = 2)
        rownames(exprMat_) = rownames_exprMat_
        duplicate_genes = rownames_exprMat_[(duplicated(rownames_exprMat_))]
        gene_selection_beforehand = rownames_exprMat_[!(rownames_exprMat_ %in% duplicate_genes)]
        exprMat = exprMat_[gene_selection_beforehand,]
            dim(exprMat)
        rm('rownames_exprMat_') # remove temporary matrix
            
        # Quality control is in principle already done, but many genes might be present because
        # they're expressed in one of the other datasets/patients.
        # So let's here again create a gene filter
        genesStatPositiveCells = rowMeans(exprMat>0)
        exprMat_filtered = exprMat[genesStatPositiveCells>.1,]
            dim(exprMat_filtered) # 5833 x 306
            
        # Creating a loom file
        print(paste0('Writing loom for ', current_set))
        lfile <- loomR::create(filename = paste0(base_dir, 'Ldata/', current_set, '.loom'), 
                               data = exprMat_filtered, overwrite = T) 
                               #feature.attrs = list(gene_info=), 
                               #cell.attrs = list(cell=colnames(current_analysis$R.P1RID2l@assays$RNA@data)))
        
    }
}

# check
for (idx in 1:3) {

    DATASET_NAME = DATASET_NAMES[idx]
    ALL_PATIENTS = PATIENT_SETS[[idx]]
    
    for (CURRENT_PATIENT in ALL_PATIENTS) {

        current_set = paste0(CURRENT_PATIENT, 'RID2l')
        
        if (file.exists(paste0(base_dir, 'Ldata/', current_set, '.loom'))) {
            print(paste0(current_set, '.loom', ' PRESENT'))
        } else {print(paste0(current_set, '.loom', ' NOT present'))}
               
    }
}

################################################################################
# << POST PROCESSING >>
# Now collect all the result matrices per patient
# Also construct regulon lists using this info

DATASET_NAME = DATASET_NAMES[1] # = "ROOIJonly.sp.bt_RID2l"

ALL_PATIENTS = paste0("R.P", 1:5)

ROOIJ_PATIENTS = paste0("R.P", 1:5)
HU_PATIENTS    = c('H.N1', 'H.N2', 'H.N3', 'H.N4', 'H.N5', 'H.N13')
TEICHMANN_PATIENTS = c('T.D2', 'T.D3', 'T.D4', 'T.D5', 'T.D6', 'T.D7', 'T.D11', 'T.H2', 'T.H3', 'T.H4', 'T.H5', 'T.H6', 'T.H7') # 'T.D1', 

ALL_PATIENTS = c(ROOIJ_PATIENTS, HU_PATIENTS, TEICHMANN_PATIENTS)
# 'H.N14' removed because the analysis failed
# 'T.D1' removed because it disappeared after removing non-septal cells

loom_file_list = list()
scenic_regulons_collected_all_patients_list=list()
scenic_regulons_collected_all_patients_list_regnames=list()
for (CURRENT_PATIENT in ALL_PATIENTS) {

    print(paste0('Loading ', CURRENT_PATIENT))
    
    # collect in list
    loom_file_list[[CURRENT_PATIENT]] <- loomR::connect(filename=paste0(base_dir, '/SCENIC/',CURRENT_PATIENT,'/output/SCENIC_out_',CURRENT_PATIENT,'RID2l.loom'), 
                        mode = 'r')

    df_regulons_compare = as.data.frame(loom_file_list[[CURRENT_PATIENT]]$row.attrs$Regulons[])
    row.names(df_regulons_compare)=loom_file_list[[CURRENT_PATIENT]]$row.attrs$Gene[]
    
    scenic_regulons_collected_all_patients_list[[CURRENT_PATIENT]] = 
        lapply(1:ncol(df_regulons_compare), function(reg_idx) {rownames(df_regulons_compare)[df_regulons_compare[,reg_idx]==1]})
    names(scenic_regulons_collected_all_patients_list[[CURRENT_PATIENT]]) = colnames(df_regulons_compare) # paste0(CURRENT_PATIENT,'.',colnames(df_regulons_compare))
    scenic_regulons_collected_all_patients_list_regnames[[CURRENT_PATIENT]] = colnames(df_regulons_compare)
}

# scenic_regulons_collected_all_patients_list -- holds a list per patient
# scenic_regulons_collected_all_patients_list_regnames -- holds the regulon names per patient

# now create also flattened list
scenic_regulons_collected_all_patients = unlist(scenic_regulons_collected_all_patients_list, recursive = F)
scenic_regulons_collected_all_patients_regnames = 
    unlist(sapply(1:length(scenic_regulons_collected_all_patients_list), function(i) {names(scenic_regulons_collected_all_patients_list[[i]])}), recursive = F)
names(scenic_regulons_collected_all_patients_regnames) = names(scenic_regulons_collected_all_patients)

################################################################################
# Now some analysis of this data

length(scenic_regulons_collected_all_patients)
    
regulon_overlap_heatmap(pooled_regulons = scenic_regulons_collected_all_patients, base_dir = base_dir, run_name = 'SCENIC_regulons_collected_all_patients', MYTREECUTTINGHEIGHT = 5.14, myfontsize = 2)

regulon_overlap_heatmap(pooled_regulons = scenic_regulons_collected_all_patients_list_regnames, base_dir = base_dir, run_name = 'SCENIC_regulons_overlap_regnames', MYTREECUTTINGHEIGHT = .1, myfontsize = 2)

# all_regulons = Reduce(x = scenic_regulons_collected_all_patients_list_regnames, f = unique)

all_regulons = unique(unlist(scenic_regulons_collected_all_patients_list_regnames))

################################################################################
# Upset like heatmap
# ===

upset_df = data.frame(lapply(scenic_regulons_collected_all_patients_list_regnames, function(tf_list) {all_regulons %in% tf_list}), row.names=all_regulons)

# Upset plot not so clear
# library(UpSetR)
# upsetR_df = UpSetR::fromList(scenic_regulons_collected_all_patients_list_regnames)
# UpSetR::upset(upsetR_list)

# Heatmap of upset df

# Make the plot for all patients
cellsize=200/nrow(upset_df)
p=pheatmap(1*upset_df, fontsize = cellsize*.pt, cluster_cols = F, cluster_rows = T, color = c('white','black'), legend = F, 
           treeheight_row = 0, cellwidth = cellsize*.pt, cellheight = cellsize*.pt, border_color = NA) # .pt = points/mm
p
if (!exists('NOSAVE')) { ggsave(filename = paste0(base_dir,'Rplots/ALL.SPcustom_7_RegulonsSCENIC_overlapHeatmap.pdf'), 
        plot = p, width=15+ncol(upset_df)*cellsize, height=15+nrow(upset_df)*cellsize, units='mm') # 184.6/3*2-4
}

# Now Rooij-specific 
upset_df_Rooij = upset_df[,c(paste0('R.P',1:5))]
upset_df_Rooij = upset_df_Rooij[apply(upset_df_Rooij, 1, any),]
cellsize=200/nrow(upset_df_Rooij)
p=pheatmap(1*upset_df_Rooij, fontsize = cellsize*.pt, cluster_cols = F, cluster_rows = T, color = c('white','black'), legend = F, 
           treeheight_row = 0, cellwidth = cellsize*.pt, cellheight = cellsize*.pt, border_color = NA) # .pt = points/mm
p
if (!exists('NOSAVE')) { ggsave(filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_RegulonsSCENIC_overlapHeatmap.pdf'), 
        plot = p, width=15+ncol(upset_df_Rooij)*cellsize, height=15+nrow(upset_df_Rooij)*cellsize, units='mm') # 184.6/3*2-4
}

# Specific reviewer's question: are there regulons specific to MYBPC3 patients?
# At least no P2 & P4 specific ones (undisclosed/untested patients might also carry MYBPC3 mutations)
upset_df_Rooij[upset_df_Rooij$R.P2&upset_df_Rooij$R.P4&(!upset_df_Rooij$R.P1)&(!upset_df_Rooij$R.P3)&(!upset_df_Rooij$R.P5),]


colnames(upset_df)

# pt_count_prelim = apply(upset_df_Rooij[,1:5], 1, sum)
upset_df_rooij_expanded = upset_df[rownames(upset_df_Rooij),]
pt_count = apply(upset_df_rooij_expanded[,1:5], 1, sum)


pheatmap(as.matrix(1*upset_df), cluster_cols = F)
pheatmap(as.matrix(1*upset_df[rownames(upset_df_Rooij),]), cluster_cols = F)
pheatmap(as.matrix(1*upset_df_rooij_expanded[order(pt_count, decreasing = T),]), cluster_cols = F, cluster_rows = F)

upset_df_rooij_expanded_sel = upset_df_rooij_expanded[pt_count>3,]
pt_count_sel    = apply(upset_df_rooij_expanded_sel[,1:5], 1, sum)
donor_count = apply(upset_df_rooij_expanded_sel[,6:(dim(upset_df_rooij_expanded_sel)[2])], 1, sum)
pheatmap(as.matrix(1*upset_df_rooij_expanded_sel[order(pt_count_sel, decreasing = T),]), cluster_cols = F, cluster_rows = F)

min(donor_count)

pheatmap(as.matrix(1*upset_df_rooij_expanded_sel[donor_count<3,]), cluster_cols = F, cluster_rows = F)
dim(upset_df_rooij_expanded_sel)[2]-5

################################################################################
# Now retreive the NES scores
#

# Flag to export stuff later on
lets_export=1
    # Note that all the way done in this script, this flag needs to be 1 such that the supp dataset will be printed.

if (exists('lets_export')) {
    all_collected_add.info = list()
}

myRegScores=list()
for (CURRENT_PATIENT in ALL_PATIENTS) {
    
    # CURRENT_PATIENT = 'R.P1'
 
    # Gather some info about the reliability of the link between gene sets and their TF
    add.info.meta = read.csv(paste0(base_dir,'/SCENIC/',CURRENT_PATIENT,'/reg.csv'), head=F, nrows=3)
    add.info = read.csv(paste0(base_dir,'/SCENIC/',CURRENT_PATIENT,'/reg.csv'), head=F, skip = 3)
    colnames(add.info)=c(add.info.meta[3, 1:2], add.info.meta[2, 3:length(add.info)])
        # View(add.info)
    
    # Now summarize this to have NES scores available per regulon
    myRegScores[[CURRENT_PATIENT]] = aggregate(list(NES=add.info$NES), by = list(TF=add.info$TF), FUN=median)
    # myRegScores[[CURRENT_PATIENT]] = aggregate(list(NES=add.info$NES), by = list(TF=add.info$TF), FUN=max)
    rownames(myRegScores[[CURRENT_PATIENT]]) = myRegScores[[CURRENT_PATIENT]]$TF
       
    # If desired, also gather export data
    if (exists('lets_export')) {
        all_collected_add.info[[paste0(CURRENT_PATIENT, '_add.info')]] = add.info
    }
    
    print(paste0(CURRENT_PATIENT,' done ..'))
    
}

# Slightly more sophisticated "upset_df"
upset_df_NES = data.frame(lapply(names(scenic_regulons_collected_all_patients_list_regnames), 
     function(current_patient)  {
        tf_list = scenic_regulons_collected_all_patients_list_regnames[[current_patient]]
        regulons_of_interest = all_regulons[all_regulons %in% tf_list]
        regulons_of_interest_  = gsub(pattern = '\\([+-]\\)',replacement = '', regulons_of_interest)
        values = myRegScores[[current_patient]][regulons_of_interest_,]$NES
        column_w_score = rep(NA,length(all_regulons))
        column_w_score[all_regulons %in% tf_list] = values
        return(column_w_score)
    }), row.names=all_regulons)
colnames(upset_df_NES) = names(scenic_regulons_collected_all_patients_list_regnames)
rownames(upset_df_NES)=gsub('\\(\\+\\)','',rownames(upset_df_NES))


# Make heatmap again
cellsize=200/nrow(upset_df_NES)
upset_df_NES_ = upset_df_NES; upset_df_NES_[is.na(upset_df_NES)]=0
upset_df_NES_ = Reduce(rbind, rev(lapply(split(upset_df_NES_, apply(upset_df[rownames(upset_df_NES_),], 1, sum)), function(X) {
        # order method of the separate parts 
        # X[hclust(dist(X))$order,] # clustering
        X[order(apply(X, 1, mean), decreasing = T),]
    })))
# beautify colnames
theDonorNames = colnames(upset_df_NES_)
theDonorNames=gsub('^R\\.','HCM.',theDonorNames);theDonorNames=gsub('^H\\.','Ctrl1.',theDonorNames);theDonorNames=gsub('^T\\.','Ctrl2.',theDonorNames)
upset_df_NES_btNames_ = upset_df_NES_
colnames(upset_df_NES_btNames_)=theDonorNames
custom_breaks = seq(0, max(upset_df_NES_btNames_), max(upset_df_NES_btNames_)/99)
p=pheatmap(upset_df_NES_btNames_, fontsize = cellsize*.pt, cluster_cols = F, cluster_rows = F, legend = F, 
           treeheight_row = 0, cellwidth = cellsize*.pt, cellheight = cellsize*.pt, border_color = NA, breaks = custom_breaks) # .pt = points/mm
#p
if (!exists('NOSAVE')) { ggsave(filename = paste0(base_dir,'Rplots/','ALL.SPcustom','_7_RegulonsSCENIC_overlapHeatmap_NES.pdf'), 
        plot = p, width=5+ncol(upset_df_NES_btNames_)*cellsize, height=15+nrow(upset_df_NES_btNames_)*cellsize, units='mm') # 184.6/3*2-4
}

# Rather arbitrarily, show 1st X rows, just as a zoom
FIRSTX=70
MARGIN=10
cellsize=(100-MARGIN)/FIRSTX
p=pheatmap(upset_df_NES_btNames_[1:FIRSTX,], fontsize = cellsize*.pt*.9, cluster_cols = F, cluster_rows = F, legend = F, 
           treeheight_row = 0, cellwidth = cellsize*.pt, cellheight = cellsize*.pt, border_color = NA, breaks = custom_breaks) # .pt = points/mm
#p
if (!exists('NOSAVE')) { ggsave(filename = paste0(base_dir,'Rplots/','ALL.SPcustom','_7_RegulonsSCENIC_overlapHeatmap_NES_first',FIRSTX,'.pdf'), 
        plot = p, width=MARGIN+ncol(upset_df_NES_btNames_)*cellsize, height=MARGIN+FIRSTX*cellsize, units='mm') # 184.6/3*2-4
}
# Same but rotate 90 deg
p=pheatmap(t(upset_df_NES_btNames_[1:FIRSTX,]), fontsize = cellsize*.pt*.9, cluster_cols = F, cluster_rows = F, legend = F, 
           treeheight_row = 0, cellwidth = cellsize*.pt, cellheight = cellsize*.pt, border_color = NA, angle_col = '90', breaks = custom_breaks) # .pt = points/mm
p
if (!exists('NOSAVE')) { ggsave(filename = paste0(base_dir,'Rplots/','ALL.SPcustom','_7_RegulonsSCENIC_overlapHeatmap_NES_first',FIRSTX,'-rotated.pdf'), 
        plot = p, height=MARGIN+ncol(upset_df_NES_btNames_)*cellsize, width=MARGIN+FIRSTX*cellsize, units='mm') # 184.6/3*2-4
}



# Now also Rooij-specific again
upset_df_NES_Rooij_ = upset_df_NES_[,c(paste0('R.P',1:5))]
colnames(upset_df_NES_Rooij_) = paste0('P',1:5)
upset_df_NES_Rooij_ = upset_df_NES_Rooij_[apply(upset_df_NES_Rooij_>0, 1, any),]

# Create index to order the regulons (note that it's already sorted, but for all patients, including other data)
nr_patients_w_R = apply(upset_df_NES_Rooij_>0,1,sum)
mean_NES_score  = apply(upset_df_NES_Rooij_,1,mean)
custom_order = unlist(lapply(5:1, function(x) {theorder=order(mean_NES_score, decreasing = T); theorder[theorder %in% which(nr_patients_w_R==x)]}))
# alternative 2
custom_order = order(mean_NES_score, decreasing = T)
# perform sort
upset_df_NES_Rooij_ordered_ = upset_df_NES_Rooij_[custom_order,]

# Heatmap
p=pheatmap(upset_df_NES_Rooij_ordered_, fontsize = cellsize*.pt, cluster_cols = F, cluster_rows = F, legend = F, 
           treeheight_row = 0, cellwidth = cellsize*.pt, cellheight = cellsize*.pt, border_color = NA) # .pt = points/mm
p
if (!exists('NOSAVE')) { ggsave(filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_RegulonsSCENIC_overlapHeatmap_NES.pdf'), 
        plot = p, width=15+ncol(upset_df_NES_Rooij_)*cellsize, height=15+nrow(upset_df_NES_Rooij_)*cellsize, units='mm', device = cairo_pdf) # 184.6/3*2-4
}
# ggsave(filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_RegulonsSCENIC_overlapHeatmap_NES_testMaxNES.pdf'), plot = p, width=15+ncol(upset_df_NES_Rooij_)*cellsize, height=15+nrow(upset_df_NES_Rooij_)*cellsize, units='mm', device = cairo_pdf)

# Now again Rooij-specific, but only regulons in >1 patient, and beautified names
eff_rows = nrow(upset_df_NES_Rooij_ordered_[current_selection,]); eff_cols=ncol(upset_df_NES_Rooij_ordered_[current_selection,])
cellsize=200/eff_rows
current_selection = apply(upset_df_NES_Rooij_ordered_>0,1,sum)>1
psel=pheatmap(upset_df_NES_Rooij_ordered_[current_selection,], fontsize = cellsize*.pt, cluster_cols = F, cluster_rows = F, legend = F, 
           treeheight_row = 0, cellwidth = cellsize*.pt, cellheight = cellsize*.pt, border_color = NA, device = cairo_pdf) # .pt = points/mm
psel
if (!exists('NOSAVE')) { ggsave(filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_RegulonsSCENIC_overlapHeatmap_NES_selMoreThan1Pat.pdf'), 
        plot = psel, width=15+eff_cols*cellsize, height=15+eff_rows*cellsize, units='mm', device = cairo_pdf) # 184.6/3*2-4
}

################################################################################
# Now let's export this data for reference

View(upset_df_NES_)

################################################################################
# Let's quickly check presence of regulons in how # patients per dataset

upset_df_sel = upset_df[apply(upset_df,1,sum)>3,]

patient_sum_df = 
    as.data.frame(
    cbind(apply(upset_df[,ROOIJ_PATIENTS], 1, sum), 
          apply(upset_df[,HU_PATIENTS], 1, sum), 
          apply(upset_df[,TEICHMANN_PATIENTS], 1, sum)))
colnames(patient_sum_df) = c('R','H','T')
patient_sum_df$TF=gsub(pattern = '\\(\\+\\)',replacement = '',x = rownames(patient_sum_df))
patient_sum_df_melted = melt(patient_sum_df, id.vars = 'TF', variable.name = 'dataset')
# patient_sum_df_melted = melt(patient_sum_df, varnames = c('TF','dataset'))

# ggplot(patient_sum_df_melted, aes(x=TF, fill=dataset, y=value))+
#     geom_bar(stat='identity', position='dodge')+
#     facet_wrap(facets = vars(TF))

patient_sum_df_extended=patient_sum_df
patient_sum_df_extended$TF=gsub(pattern = '\\(\\+\\)',replacement = '',x = rownames(patient_sum_df))
p=ggplot(patient_sum_df_extended,aes(x=R, y=T))+
    geom_density_2d()+
    geom_density_2d(data=patient_sum_df_extended[patient_sum_df_extended$R>=1,], color='green')+
    geom_point(size=.5)+
    geom_text_repel(data=patient_sum_df_extended[patient_sum_df_extended$R>=3&patient_sum_df_extended$T<=3, ], 
                    aes(label=TF), max.overlaps = Inf, size=6/.pt, segment.size=.25, color='red')+
    theme_bw()+give_better_textsize_plot(8)
# p
if (!exists('NOSAVE')) { ggsave(filename = paste0(base_dir,'Rplots/','ALL.SPcustom','_7_RegulonsSCENIC_scatterPatCount.pdf'), 
        plot = p, width=172/3-4, height=172/3-4, units='mm', device = cairo_pdf) # 184.6/3*2-4
}

p=ggplot(patient_sum_df_extended,aes(x=R, y=T))+
    geom_density_2d_filled(data=patient_sum_df_extended[patient_sum_df_extended$R>1,])+theme_minimal()+
    geom_density_2d(color='white')+
    geom_point(data=patient_sum_df_extended[patient_sum_df_extended$R>=3&patient_sum_df_extended$T<=3, ],
               size=.5,color='red')+
    geom_text_repel(data=patient_sum_df_extended[patient_sum_df_extended$R>=3&patient_sum_df_extended$T<=3, ], 
                    aes(label=TF), max.overlaps = Inf, size=6/.pt, segment.size=.25, color='red')+
    theme_bw()+give_better_textsize_plot(8)+theme(legend.position = 'none')
p
if (!exists('NOSAVE')) { ggsave(filename = paste0(base_dir,'Rplots/','ALL.SPcustom','_7_RegulonsSCENIC_scatterPatCount_v2.pdf'), 
        plot = p, width=172/3-4, height=172/3-4, units='mm', device = cairo_pdf) # 184.6/3*2-4
}

# Here are all the HCM-specific specific thingies
patient_sum_df_extended[patient_sum_df_extended$R>=3&patient_sum_df_extended$T<=2, ]
patient_sum_df_extended[patient_sum_df_extended$R>=3&patient_sum_df_extended$H<=2, ]
patient_sum_df_extended[patient_sum_df_extended$R>=3&(patient_sum_df_extended$T<=2|patient_sum_df_extended$H<=2), ]
patient_sum_df_extended[patient_sum_df_extended$R>=3&(patient_sum_df_extended$T<=2&patient_sum_df_extended$H<=2), ]

SCENIC_patient_counts_overview = list(patient_totals=patient_sum_df_extended, selected_regulons=
                                          patient_sum_df_extended[patient_sum_df_extended$R>=3&(patient_sum_df_extended$T<=2|patient_sum_df_extended$H<=2), ])

if (!exists('NOSAVE')) { 
    openxlsx::write.xlsx(x= SCENIC_patient_counts_overview, file = paste0(base_dir,'Rplots/','ALL.SPcustom','_SCENIC_patientsStats_regulons.xlsx'), overwrite = T)
}

################################################################################
# Now further delve into v. Rooij Regulons, select first by requiring minimally present in 4 patients
NR_PATIENTS_SCEN_FILTER = 3

# Create some Rooij specific info
# --> Use "upset_df_NES_Rooij_ordered_" for this.

current_selection_regulons_SCENIC_Rooij = rownames(upset_df_NES_Rooij_ordered_[apply(upset_df_NES_Rooij_ordered_>0, 1, sum)>=NR_PATIENTS_SCEN_FILTER,])
# potential other selections:
rownames(upset_df_NES_Rooij_ordered_[apply(upset_df_NES_Rooij_ordered_>0, 1, sum)>=5,])
rownames(upset_df_NES_Rooij_ordered_[apply(upset_df_NES_Rooij_ordered_>0, 1, sum)>=3,])
    
# Export this selection, with patient count overview
SCENIC_regulons_Rooij_selection_overview = cbind(upset_df_NES_Rooij_ordered_, patient_sum_df[paste0(rownames(upset_df_NES_Rooij_ordered_),'(+)'),])
openxlsx::write.xlsx(x= SCENIC_regulons_Rooij_selection_overview, file = paste0(base_dir,'Rplots/','ALL.SPcustom','_SCENIC_overview_Rooij_regulons.xlsx'), overwrite = T)

# Find out gene overlap between those ones
scenic_regulons_collected_all_patients_overlapping_ones_Rooij = 
    scenic_regulons_collected_all_patients[scenic_regulons_collected_all_patients_regnames %in% current_selection_regulons_SCENIC_Rooij]

# regulon_overlap_heatmap(pooled_regulons = scenic_regulons_collected_all_patients_overlapping_ones, base_dir = base_dir, run_name = 'SCENIC_regulons_collected_all_patients_overlapping', MYTREECUTTINGHEIGHT = .1, myfontsize = 2)

# Better calculate the overlap for each respective regulon between the patients
# First create Rooij-specific selection table
scenic_regulons_collected_all_patients_regnames_ROOIJonly = scenic_regulons_collected_all_patients_regnames[
    grepl('R\\.P',names(scenic_regulons_collected_all_patients_regnames))]
scenic_regulons_collected_all_patients_ROOIJonly = scenic_regulons_collected_all_patients[
    grepl('R\\.P',names(scenic_regulons_collected_all_patients))]
overlapping_scores_regulons_Rooij = 
    sapply(current_selection_regulons_SCENIC_Rooij, function(tf) {
        # tf = 'STAT6(+)' 
        # tf = 'MEF2A(+)' # current_selection_regulons_SCENIC_Rooij[1]
        # tf = 'MLX(+)'
        # tf = 'ZNF540(+)'
        current_reg_multi_pt = scenic_regulons_collected_all_patients_ROOIJonly[scenic_regulons_collected_all_patients_regnames_ROOIJonly == paste0(tf,'(+)')]
        all_targets = unique(Reduce(x = current_reg_multi_pt, f = c)) # sort(Reduce(x = current_reg_multi_pt, f = unique))
        upset_df_tf = data.frame(lapply(current_reg_multi_pt, function(target_list) {all_targets %in% target_list}), row.names=all_targets)
        # current_overlap = mean(apply(upset_df_tf, 1, sum)>=3)
        rsize_per_patient=apply(upset_df_tf, 2, sum)
        effective_max_hits = min(sort(rsize_per_patient[rsize_per_patient>0])[1:3]) 
            # OVERLAP IS DETERMINED IN SLIGHTLY SOPHISTICATED WAY
            # rationale: if we were to have a regulon present in 4 patients, and one of those patient-regulons is very small, but the others are big,
            # i'd like to only compare overlap for the three biggest ones
            # for now, something that happens in 3 patients is considered above the magic treshold
        current_overlap = min(1, sum(apply(upset_df_tf, 1, sum)>=3)/effective_max_hits)
        
    })
# rownames(upset_df_tf[apply(upset_df_tf, 1, sum)>=2,])
overlapping_scores_regulons_df_Rooij = data.frame(overlap=overlapping_scores_regulons_Rooij, gene=names(overlapping_scores_regulons_Rooij))
overlapping_scores_regulons_df_Rooij$gene = factor(overlapping_scores_regulons_df_Rooij$gene, levels=overlapping_scores_regulons_df_Rooij$gene[order(overlapping_scores_regulons_df_Rooij$overlap)])
# overlapping_scores_regulons_df_Rooij$score = overlapping_scores_regulons_df_Rooij$gene
# Determine NES scores per TF
#overlapping_scores_regulons_df_Rooij$NES = apply(upset_df_NES_ROOIJonly, 1, function(x) {median(x[x>0&!is.na(x)])})[overlapping_scores_regulons_df_Rooij$gene]
overlapping_scores_regulons_df_Rooij$NES = apply(upset_df_NES_Rooij_ordered_, 1, function(x) {median(x[x>0&!is.na(x)])})[overlapping_scores_regulons_df_Rooij$gene]
    

# Now show gene overlap between those regulons that were identified in the separate patients 
overlapping_scores_regulons_df_Rooij$gene_ = factor(gsub('\\(\\+\\)','',overlapping_scores_regulons_df_Rooij$gene), levels=(gsub('\\(\\+\\)','',levels(overlapping_scores_regulons_df_Rooij$gene))))
p=ggplot(overlapping_scores_regulons_df_Rooij, aes(x=gene_, y=overlap))+
    geom_bar(stat='identity')+coord_flip()+theme_bw()+give_better_textsize_plot(8)+ylim(c(0,1))+ylab('Fraction identical genes\nbetween patients')+xlab('Regulon')
# p
if (!exists('NOSAVE')) { ggsave(filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_RegulonsSCENIC_MemberGeneConsistencyPatients.pdf'), 
        plot = p, width=50, height=8/.pt*nrow(overlapping_scores_regulons_df_Rooij), units='mm', device=cairo_pdf) # 184.6/3*2-4
}

# Now including NES color
p=ggplot(overlapping_scores_regulons_df_Rooij, aes(x=gene_, y=overlap, fill=NES))+
    geom_bar(stat='identity')+coord_flip()+theme_bw()+give_better_textsize_plot(7)+
    ylab('Fraction identical\ngenes between\npatients')+xlab('Regulon')+
    # scale_fill_gradientn(colours = rainbow_colors)+
    scale_fill_gradientn(colours = c('orange','red'))+theme(legend.position='none')+
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .5))
# p
if (!exists('NOSAVE')) { ggsave(filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_RegulonsSCENIC_MemberGeneConsistencyPatients_NES.pdf'), 
        plot = p, width=39, height=8/.pt*nrow(overlapping_scores_regulons_df_Rooij), units='mm', device = cairo_pdf) # 184.6/3*2-4
}


################################################################################
# Now a secondary thing, let's calculate potential overlap with each other
# First determine core regulons, genes present in three regulons, 
# finally only take regulons >10 genes
# Earlier I required also 50% overlap minimum, but I set that to 0 since it has
# no effect on final selection

PERCENTAGE_TRESHOLD=0
SCENIC_regulons_core_genes = 
    sapply(current_selection_regulons_SCENIC_Rooij[overlapping_scores_regulons_Rooij>PERCENTAGE_TRESHOLD], function(tf) { # [overlapping_scores_regulons>.1]
        # tf = current_selection_regulons_SCENIC_Rooij[1]
        current_reg_multi_pt = scenic_regulons_collected_all_patients_ROOIJonly[scenic_regulons_collected_all_patients_regnames_ROOIJonly == paste0(tf,'(+)')]
        all_targets = unique(Reduce(x = current_reg_multi_pt, f = c))
        upset_df_tf = data.frame(lapply(current_reg_multi_pt, function(target_list) {all_targets %in% target_list}), row.names=all_targets)
        current_overlap = rownames(upset_df_tf)[apply(upset_df_tf, 1, sum)>=3]
    })
# Take any with >=10 genes from these
SCENIC_regulons_core_genes_sel = SCENIC_regulons_core_genes[sapply(SCENIC_regulons_core_genes, length)>=10]
    # TFOI_list = c('CREB3L1', 'HLF', 'IRF6', 'KLF15', 'NFYB', 'NPAS2', 'NR1H2', 'SMAD4', 'SRF', 'STAT1', 'USF2', 'ZBED1', 'ZBTB7A', 'GATA4')
    # TFOI_list[paste0(TFOI_list,'(+)') %in% names(SCENIC_regulons_core_genes_sel)]

length(SCENIC_regulons_core_genes_sel)

save(list='SCENIC_regulons_core_genes_sel', file=paste0(base_dir, 'Rdata/SCENIC_regulons_core_genes_sel.Rdata'))
    # load(file=paste0(base_dir, 'Rdata/SCENIC_regulons_core_genes_sel.Rdata')) # SCENIC_regulons_core_genes_sel


################################################################################
# Additionally, also show overlap inbetween regulons

SCENIC_regulons_core_genes_sel_ = SCENIC_regulons_core_genes_sel
names(SCENIC_regulons_core_genes_sel_) = gsub('\\(\\+\\)','',names(SCENIC_regulons_core_genes_sel_))
regulon_overlap_heatmap(pooled_regulons = SCENIC_regulons_core_genes_sel_, base_dir = base_dir, run_name = 'SCENIC_regulons_overlap_core_sel', 
                        MYTREECUTTINGHEIGHT = 1.8, myfontsize = 8, makeallheatmaps=T)
    

################################################################################
# Create lists of top genes for the selected core regulons

# First load and organize weights of the genes

# Now also get gene weights
RELEVANT_REGULONS = gsub(pattern = '\\([+-]\\)', replacement = '', x = names(SCENIC_regulons_core_genes_sel))
regulon_gene_importance = list()
for (CURRENT_PATIENT in ROOIJ_PATIENTS) {
    print(paste0('Looking up selected genes for ',CURRENT_PATIENT))
    
    adj.info = read.csv(paste0(base_dir, '/SCENIC/',CURRENT_PATIENT,'/adj.csv'), head=1)
    adj.info = adj.info[adj.info$TF %in% RELEVANT_REGULONS,]
    
    # Importance: higher = more important
    regulon_gene_importance[[CURRENT_PATIENT]] = 
        lapply(unique(adj.info$TF), function(tf) {adj.info[adj.info$TF==tf,c('target', 'importance')]})
    names(regulon_gene_importance[[CURRENT_PATIENT]]) = unique(adj.info$TF)
    
}
rm('adj.info')


# Now first create a summary parameter for each of the lists 
TOPX=15
SCENIC_reg_top_genes = 
    lapply(names(SCENIC_regulons_core_genes_sel), function(tf) {
    
        # tf = 'ZEB1(+)'
        # tf = 'MEF2A(+)'
        # tf = 'USF2(+)'
        current_member_list = SCENIC_regulons_core_genes_sel[[tf]]
    
        tf_ = gsub('\\([+-]\\)', '', tf)
        
        relevant_weight_values_per_pt = lapply(regulon_gene_importance, function(X) {X[[tf_]]})
        
        importance_df = data.frame(lapply(relevant_weight_values_per_pt, function(X) {if(is.null(X)) {rep(NA, length(current_member_list))} else {rownames(X) = X$target; X[current_member_list,]$importance}}), row.names=current_member_list)
        #importance_df = data.frame(lapply(relevant_weight_values_per_pt, function(X) {rownames(X) = X$target; X[current_member_list,]$importance}), row.names=current_member_list)
        importance_df$median_nona = apply(importance_df, 1, function(x) {median(x[!is.na(x)])})
        top_genes = rownames(importance_df[order(importance_df$median_nona, decreasing = T),])[1:TOPX]
    })
names(SCENIC_reg_top_genes) = names(SCENIC_regulons_core_genes_sel)
SCENIC_reg_top_genes_df = data.frame(SCENIC_reg_top_genes)

openxlsx::write.xlsx(x= SCENIC_reg_top_genes_df[,sort(colnames(SCENIC_reg_top_genes_df))], file = paste0(base_dir,'Rplots/',DATASET_NAME,'_SCENIC_top',TOPX,'_regulons.xlsx'), overwrite = T)
    # file.exists(        paste0(base_dir,'Rplots/',DATASET_NAME,'_SCENIC_top',TOPX,'_regulons.xlsx'))
# NOTE: ALSO GIVE LENGTHS!!
SCENIC_regulon_lengths = 
    data.frame(lapply(SCENIC_regulons_core_genes_sel[sort(names(SCENIC_regulons_core_genes_sel))], length))
colnames(SCENIC_regulon_lengths)=gsub('\\([+-]\\)','',sort(names(SCENIC_regulons_core_genes_sel)))
SCENIC_regulon_lengths[2,] = paste0(colnames(SCENIC_regulon_lengths), ' (', SCENIC_regulon_lengths[1,],')')
openxlsx::write.xlsx(x= SCENIC_regulon_lengths, file = paste0(base_dir,'Rplots/',DATASET_NAME,'_SCENIC_regulon_lengths.xlsx'), overwrite = T)        

regulon_gene_importance$R.P1$MEF2A

SCENIC_regulons_core_genes_sel


# Better calculate the overlap for each respective regulon between the patietns
lists_of_regulon_genes = 
    sapply(current_selection_regulons_SCENIC_Rooij, function(tf) {
        
        # tf="MEF2A(+)"
        
        current_reg_multi_pt = scenic_regulons_collected_all_patients[scenic_regulons_collected_all_patients_regnames == paste0(tf,'(+)')]
        
        
        
    })


################################################################################
# Create new regulon overview where all genes are sorted

# First whole df
SCENIC_reg_top_genes_sorted_full_dfs = 
    lapply(names(SCENIC_regulons_core_genes_sel), function(tf) {
    
        # tf = 'ZEB1(+)'
        # tf = 'MEF2A(+)'
        # tf = 'USF2(+)'
        current_member_list = SCENIC_regulons_core_genes_sel[[tf]]
    
        tf_ = gsub('\\([+-]\\)', '', tf)
        
        relevant_weight_values_per_pt = lapply(regulon_gene_importance, function(X) {X[[tf_]]})
        
        importance_df = data.frame(lapply(relevant_weight_values_per_pt, 
                            function(X) {if(is.null(X)) {rep(NA, length(current_member_list))} else {rownames(X) = X$target; X[current_member_list,]$importance}}), row.names=current_member_list)
        #importance_df = data.frame(lapply(relevant_weight_values_per_pt, function(X) {rownames(X) = X$target; X[current_member_list,]$importance}), row.names=current_member_list)
        importance_df$median_nona = apply(importance_df, 1, function(x) {median(x[!is.na(x)])})
        
        # top_genes = rownames(importance_df[order(importance_df$median_nona, decreasing = T),])[1:100]
        importance_df[order(importance_df$median_nona, decreasing = T),]
    })
names(SCENIC_reg_top_genes_sorted_full_dfs) = gsub('\\([+-]\\)', '', names(SCENIC_regulons_core_genes_sel))

# Now just lists of names
SCENIC_reg_top_genes_sorted_full =
    lapply(names(SCENIC_regulons_core_genes_sel), function(tf) {
            tf_ = gsub('\\([+-]\\)', '', tf)
            rownames(SCENIC_reg_top_genes_sorted_full_dfs[[tf_]])
        })
names(SCENIC_reg_top_genes_sorted_full) = gsub('\\([+-]\\)', '', names(SCENIC_regulons_core_genes_sel))


# Save the list of SCENIC regulon genes
# TO DO: Not sure why I chose the name "full" here, I use selection anyways
save(x='SCENIC_reg_top_genes_sorted_full', file=paste0(base_dir,'Rdata/',DATASET_NAME,'__SCENIC_reg_top_genes_sorted_full.Rdata'))
    # load(file=paste0(base_dir,'Rdata/',DATASET_NAME,'__SCENIC_reg_top_genes_sorted_full.Rdata')) # SCENIC_reg_top_genes_sorted_full
    # file.exists(paste0(base_dir,'Rdata/',DATASET_NAME,'__SCENIC_reg_top_genes_sorted_full.Rdata'))

save(list='SCENIC_reg_top_genes_sorted_full', file=paste0(base_dir, 'Rdata/SCENIC_reg_top_genes_sorted_full.Rdata'))
# load(file=paste0(base_dir, 'Rdata/SCENIC_reg_top_genes_sorted_full.Rdata')) # SCENIC_reg_top_genes_sorted_full

openxlsx::write.xlsx(x= SCENIC_reg_top_genes_sorted_full, file = paste0(base_dir, 'Rplots/SCENIC_reg_top_genes_sorted_full.xlsx'))


if (exists('lets_export')) {
    big_data_export_structure = c(SCENIC_reg_top_genes_sorted_full, scenic_regulons_collected_all_patients, all_collected_add.info)
    openxlsx::write.xlsx(x= big_data_export_structure, file = paste0(base_dir,'Rplots/','ALL.SPcustom','_SCENIC_all_regulons_with_data.xlsx'), overwrite = T)
        # print(paste0('Exported to: ',base_dir,'Rplots/'))
    rm('big_data_export_structure')
    rm('all_collected_add.info')
}

# additional overview
regulon_overview_df = as.data.frame(Reduce(cbind,lapply(SCENIC_reg_top_genes_sorted_full, function(X){c( X[1:(min(100, length(X)))],
                                                                          rep('', (100-min(100, length(X)))))})))
colnames(regulon_overview_df) = names(SCENIC_reg_top_genes_sorted_full)
regulon_overview_df=regulon_overview_df[,sort(colnames(regulon_overview_df))]
openxlsx::write.xlsx(x = regulon_overview_df, file = paste0(base_dir,'Rplots/','ALL.SPcustom','_SCENIC_all_regulons_overview.xlsx'), overwrite = T)

#####

lapply(SCENIC_reg_top_genes_sorted_full, function(X) {which(X %in% 'MYL2')})
lapply(SCENIC_reg_top_genes_sorted_full, length)


# lapply(SCENIC_reg_top_genes_top100, function(X) {which(X %in% 'MYL2')})

# Sanity check:
GENE='SRF'
target_genes= shorthand_seurat_fullgenename_faster(current_analysis[[DATASET_NAME]], rownames(SCENIC_reg_top_genes_sorted_full_dfs[[GENE]]))
gene_of_interest=shorthand_seurat_fullgenename_faster(current_analysis[[DATASET_NAME]], GENE)
correlations=
    sapply(target_genes[target_genes!=gene_of_interest], function(t) {
        #out=cor.test(current_analysis[[DATASET_NAME]]@assays$RNA@data[gene_of_interest,], current_analysis[[DATASET_NAME]]@assays$RNA@data[t,])
        #out$estimate
        cor(current_analysis[[DATASET_NAME]]@assays$RNA@data[gene_of_interest,], current_analysis[[DATASET_NAME]]@assays$RNA@data[t,])
        })

df_toplot=data.frame(corr=correlations, importance=SCENIC_reg_top_genes_sorted_full_dfs[[GENE]][names(target_genes)[target_genes!=gene_of_interest],]$median_nona, name=names(correlations))
high_corr=rownames(df_toplot[order(df_toplot$corr, decreasing = T),][1:15,])
high_importance=rownames(df_toplot[1:15,])
ggplot(df_toplot,aes(x=corr, y=importance))+
    geom_point()+
    geom_text_repel(data=df_toplot[unique(c(high_corr,high_importance)),], aes(label=name), color='red')+
    geom_smooth(method='lm')+theme_bw()+ggtitle(GENE)

################################################################################
# Now determine the relatedness of their regulons to my regulons

# Load my regulons
load(paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_core_regulons_sorted_shortname.Rdata'))
    # core_regulons_sorted_shortname

# Now compare matrix
SCENIC_MW_regulon_comparison = 
    sapply(core_regulons_sorted_shortname, function(core_mw) { # [overlapping_scores_regulons>.1]
        sapply(SCENIC_regulons_core_genes, function(core_scenic) { # [overlapping_scores_regulons>.1]
            overlap = sum(core_scenic %in% core_mw)/min(length(core_scenic), length(core_mw))
        })
    })

pheatmap(SCENIC_MW_regulon_comparison)

# Heatmap again
myfontsize=8
p=pheatmap(SCENIC_MW_regulon_comparison, cluster_rows = T,cluster_cols = F, 
            fontsize = myfontsize, fontsize_col = myfontsize, fontsize_row = myfontsize, treeheight_row = 0, treeheight_col = 0, legend = F, limits=c(0,1), border_color = NA)
            #annotation_colors = list(colors=col_Dark2[1:max(cutree_out)])))
print(p)
if (!exists('NOSAVE')) { ggsave(filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_regulon_SCENIC_MW_gene_overlap-v1b-L.pdf'), 
    plot = p, width=ncol(SCENIC_MW_regulon_comparison)*myfontsize/.pt*1.1+20, height=nrow(SCENIC_MW_regulon_comparison)*myfontsize/.pt*1.1+15, units='mm', limitsize = F)
}

# Now combined comparision heatmap
core_regulons_sorted_shortname_ = core_regulons_sorted_shortname
names(core_regulons_sorted_shortname_) = gsub('s\\.R\\.', 'Module', names(core_regulons_sorted_shortname_))
combined_regulon_scenic_and_mw = c(core_regulons_sorted_shortname_, SCENIC_reg_top_genes_sorted_full)
regulon_overlap_heatmap(pooled_regulons = combined_regulon_scenic_and_mw, base_dir = base_dir, run_name = 'regulons_overlap_both_scenic_and_mw', 
                        MYTREECUTTINGHEIGHT = 1.8, myfontsize = 8, makeallheatmaps=T)

################################################################################
# And additionally we can project their regulons on our UMAP

# DATASET_NAME = "ROOIJonly_RID2l"
# DATASET_NAME = "ROOIJonly.sp.bt_RID2l"

if (F) {
    
    # Load analysis if necessary
    if (!exists('current_analysis')) {current_analysis = list()}
    if (!(DATASET_NAME %in% names(current_analysis))) {
        current_analysis[[DATASET_NAME]] =
            LoadH5Seurat(file = paste0(base_dir,'Rdata/H5_RHL_SeuratObject_nM_sel_',DATASET_NAME,'.h5seurat'))
    }
    
    
    # X=shorthand_cutname(core_regulons$s.R.4)
    
    plotlist=lapply(# 1:length(SCENIC_reg_top_genes_sorted_full),
        order(names(SCENIC_reg_top_genes_sorted_full)), # this makes it alphabetic
        function(X) {
            gname = gsub(pattern = '\\(\\+\\)',replacement = '',x = names(SCENIC_reg_top_genes_sorted_full)[[X]])
            shorthand_seurat_custom_expr(seuratObject = current_analysis[[DATASET_NAME]], 
                                                    gene_of_interest = SCENIC_reg_top_genes_sorted_full[[X]], textsize=6, pointsize=.5, 
                                                    custom_title = paste0('italic(',gname,')'), mymargin = .5, zscore = T) # label:italic (note to self)
            }) 
    #plots_per_row=round((172/2-4)/25)
    #n_rows=ceiling(length(plotlist)/plots_per_row)
    p=wrap_plots(plotlist, nrow = 5)
    
    #p
    if (!exists('NOSAVE')) { 
        ggsave(filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_RegulonsSCENIC_UMAP_compositeExpr.pdf'), 
            plot = p, width=90, height=90, units='mm', device=cairo_pdf) # 184.6/3*2-4
        ggsave(filename = paste0(base_dir,'Rplots/',DATASET_NAME,'_7_RegulonsSCENIC_UMAP_compositeExpr-v2.pdf'), 
            plot = p, width=172/2-4, height=20*n_rows, units='mm', device=cairo_pdf) # 184.6/3*2-4
    }
    
}

################################################################################
# And also the underlying TFs

# TO DO: update this list
regulon_controllers = sort(c('CEBPB', 'ESRRA', 'ETV1', 'FOXN3', 'JUND', 'MAFK', 'MEF2A', 'MEF2D', 'MXI1', 'NR3C1', 'SREBF2', 'SRF', 'USF2', 'YY1', 'ZEB1', 'ZNF91'))
regulon_controllers = sort(c('YY1', 'SRF', 'NR3C1', 'NFE2L1', 'MXI1', 'MITF', 'MEF2D', 'MEF2A', 'MAFK', 'JUND', 'JUN', 'FOXP1', 'FOXN3', 'ETV1', 'ATF6'))
regulon_controllers = sort(gsub('\\(\\+\\)','',names(SCENIC_regulons_core_genes)))

p_list=lapply(regulon_controllers, function(gene) {
    shorthand_seurat_custom_expr(seuratObject = current_analysis[[ANALYSIS_NAME]], 
                                                    gene_of_interest = gene, textsize=4, pointsize=.25, 
                                                    custom_title = paste0(gene), mymargin = .5, zscore = T)})
p=wrap_plots(p_list, ncol=5)#ceiling(sqrt(length(p_list))))
if (!exists('NOSAVE')) { 
    ggsave(filename = paste0(base_dir,'Rplots/',ANALYSIS_NAME,'_9_custom_SCENIC_tfs.pdf'), 
           plot = p, width=20*5, height=5*20, units='mm', device = cairo_pdf) # 184.6/3*2-4
}
################################################################################
# Now also sort their regulons' genes, just using the same method as for my
# own regulon analysis
# 
# df_core_regulon_overview_list=list()
# 
# function(core_regulon_list, current_analysis, DATASET_NAME) {
#     
#     for (core_reg_name in names(core_regulon_list)) {
#         
#         if (length(core_regulon_list[[core_reg_name]])<1) { break }
#         
#         # core_reg_name='s.R.4'
#         
#         all_patients = unique(current_analysis[[DATASET_NAME]]$annotation_patient_str)
#         #all_patients = names(collected_regulon_objects$TEICHMANNonly_RID2l)
#         df_reg_genescores = data.frame(gene=character(), corr=numeric())
#         for (current_patient in all_patients) {
#         
#             print(paste0('Calculating ',core_reg_name, ' - ',current_patient))
#             
#             patient_selection = current_analysis[[DATASET_NAME]]$annotation_patient_str==current_patient
#             
#             expr_mat = as.matrix(current_analysis[[DATASET_NAME]]@assays$RNA@data)
#             rownames(expr_mat) = shorthand_cutname(rownames(expr_mat))
#             
#             # current_analysis[[DATASET_NAME]]@assays$RNA@counts[core_regulon_list[[core_reg_name]], ]
#             # Let's do this per patient again
#             cor_out = cor(t(expr_mat[core_regulon_list[[core_reg_name]], patient_selection]))
#             cor_reg_avg = colMeans(cor_out)
#             #View(cor_reg_avg)
#             
#             df_reg_genescores=rbind(df_reg_genescores, data.frame(gene=names(cor_reg_avg), corr=cor_reg_avg))
#             
#         }
#         df_reg_genescores_sum       = aggregate(x = list(median_corr=df_reg_genescores$corr), by = list(gene=df_reg_genescores$gene), median)
#         
#         # !!!! CONTINUE EDITING HERE!!!! 
#         XXXXX
#         df_reg_genescores_sum$nr_patients = XXXX core_regulons_gene_info_table[[core_reg_name]][df_reg_genescores_sum$gene,]$total_occurence XXXX
#         XXXXX
#         # !!!! CONTINUE EDITING HERE!!!!
#         rownames(df_reg_genescores_sum) = df_reg_genescores_sum$gene
#         
#         # Sort first by patients, then by correlation (median)
#         df_reg_genescores_sum$nr_patients <- ordered(df_reg_genescores_sum$nr_patients, levels = sort(unique(df_reg_genescores_sum$nr_patients), decreasing = T))
#         df_core_regulon_overview_list[[core_reg_name]]=
#             do.call(rbind, lapply(split(df_reg_genescores_sum, df_reg_genescores_sum$nr_patients), function(x) x[order(x$median_corr,decreasing = T),]))
#         
#     }
#     
#     # Create output
#     core_regulons_sorted = lapply(df_core_regulon_overview_list, function(x) {x$gene})
#     matrix_core_regulons_padded = t(plyr::ldply(core_regulons_sorted, rbind))
#     df_core_regulons_padded = data.frame(matrix_core_regulons_padded[-1,])
#     colnames(df_core_regulons_padded)=matrix_core_regulons_padded[1,]
#     openxlsx::write.xlsx(x= df_core_regulons_padded, file = paste0(base_dir,'Rplots/',DATASET_NAME,'_core_regulons.xlsx'), overwrite = T)
#     
#     df_core_regulons_padded_shortname=as.data.frame(shorthand_cutname_table(df_core_regulons_padded))
#     openxlsx::write.xlsx(x= df_core_regulons_padded_shortname, file = paste0(base_dir,'Rplots/',DATASET_NAME,'_core_regulons_shortname.xlsx'), overwrite = T)
#     
# }
# 
# if (F) {
# 
#     core_regulon_list=SCENIC_regulons_core_genes
#     
#     
#     
# 
# }    
    
################################################################################

    
    
    
    
################################################################################
# Code graveyard
    
# # Now look at results
# lfile <- loomR::connect(filename='/Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/SCENIC/R.P1/output/SCENIC_out_R.P1RID2l.loom', 
#                         mode = 'r')
# 
# lfile[["matrix"]][1:5, 1:5]
# 
# full.matrix <- lfile$matrix[, ]
# dim(x = full.matrix) 
#     # 306 x 5833
#     # note: these are the cells x genes again
# 
# 
# # So, e.g. 
# lfile$row.attrs$Gene[1:3]
#     # This shows genes 1-3 of the original table
# dim(lfile$row.attrs$Regulons[1:3])
#     # This shows for each of the 127 regulons, whether the gene is a member or not (I suppose)
# 
# # So let's create a heatmap of this
# df_regulons_compare = as.data.frame(lfile$row.attrs$Regulons[])
# row.names(df_regulons_compare)=lfile$row.attrs$Gene[]
# 
# # Let's also gather the regulon members
# scenic_regulons_collected = lapply(1:ncol(df_regulons_compare), function(reg_idx) {rownames(df_regulons_compare)[df_regulons_compare[,reg_idx]==1]})
# names(scenic_regulons_collected) = colnames(df_regulons_compare)
# 
# # p=pheatmap(as.matrix(df_regulons_compare))
# 
# 
# regulon_overlap_heatmap(pooled_regulons = scenic_regulons_collected, base_dir = base_dir, run_name = 'SCENIC_test', MYTREECUTTINGHEIGHT = 1, myfontsize = 2)
# 
#     
#     
# 
# # Conversely, we also have regulon "activity" over the cells
# # (±composite expression)
# # lfile$col.attrs$
# 
# 
# 
# # And 
# lfile[["row_attrs/Regulons"]][1:3]
# 
# rownames(full.matrix) = lfile$row.attrs$Regulons[[]]
# lfile[["row_attrs/Regulons"]][1]
# lfile[["row_attrs/Regulons"]][2]
# lfile[["row_attrs/Regulons"]][2]
# 
# lfile[["row_attrs/Genes"]]
# 
# 
# 
# lfile$col.attrs
# 
# df_1 = 
#     data.frame(auc=lfile$col.attrs$RegulonsAUC, 
#             cellid=lfile$col.attrs$CellID)
# 
# regulons = lfile$row.attrs$Regulons
# 
# lfile$row.attrs$Gene
# 
# 
# #####
# 
# attrs <- c("nUMI", "nGene", "orig.ident")
# attr.df <- lfile$get.attribute.df(MARGIN = 2, attribute.names = attrs)
# head(x = attr.df)



