

# Readme

These scripts contain Seurat analysis for the publication "Single-cell transcriptomics provides insights into hypertrophic cardiomyopathy" by Wehrens et al. published in Cell Reports (2022), DOI: [https://doi.org/10.1016/j.celrep.2022.110809](https://doi.org/10.1016/j.celrep.2022.110809).

# Rdata & Seurat data files 

There are Rdata and h5seurat files that contain the backbone processed data of this analysis.
I'll be looking to put those online somewhere too, but currently haven't done that yet.
Contact me (m.wehrens@hubrecht.eu) if you would like to take a look at the data using these files.

## Easy use of this data for Van Rooij lab members

See the script 
`./Projects/howtousedata_HCMSCS/example-load-plot-data.R` and the folder `srv-lnx-varo1:/opt/backup_wehrens/data/Wehrens2022/Rdata-important` to make some plots using this data.

# About the scripts

The most central file in this analysis is `Projects/HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R`, you can start exploring these scripts there.

*More explanation will be added.*

Send an e-mail to m.wehrens@hubrecht.eu if you have any questions.

# Notes about libraries etc. 

The scripts involving the Litvinukova data were executed on a high performance computing cluster.

The necessary software and libraries were installed using Conda.

Four separate environments were created, "base" (default), "MARGE", "LISA", "pyscenic". The base environment was
used for all analysis, except for the MARGE, Lisa and SCENIC analyses. 
A list of installed libraries can be found in the text file "libraries_environments_information.txt", 
this also includes a separate list of all packages that were installed "inside" R.
















