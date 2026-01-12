

# Readme

These scripts contain Seurat analysis for the publication "Single-cell transcriptomics provides insights into hypertrophic cardiomyopathy" by Wehrens et al. published in Cell Reports (2022), DOI: [https://doi.org/10.1016/j.celrep.2022.110809](https://doi.org/10.1016/j.celrep.2022.110809).

# Data repositories

- The raw FASTQ files can be found on GEO, using accession ID [GSE138262](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138262).
- The count tables, Seurat objects on which plots are based, and FACS metadata, can be found via:
[https://doi.org/10.5281/zenodo.18206553](https://doi.org/10.5281/zenodo.18206553).

## Information about the pre-processed files

*This information is also listed at the Zenodo repository.*

#### Description

The files at [https://doi.org/10.5281/zenodo.18206553](https://doi.org/10.5281/zenodo.18206553) should allow you to reproduce analyses and figures from the paper, without re-mapping the data from FASTQ files.  
 
#### Count tables Van Rooij
The file `counttables-ROOIJ.zip` contains the following files:

| GEO | Patient | Plate | Filename |
|-----|---------|--------|----------|
| GSM4103885 | Patient 1 | Plate 1 | rJE1_AHFL7NBGX5_S3_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM4103886 | Patient 1 | Plate 2 | JE2_AHY3WGBGX3_S1_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM4103887 | Patient 1 | Plate 3 | JE3_AHY3WGBGX3_S2_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM4103888 | Patient 1 | Plate 4 | JE4_AHFL7NBGX5_S4_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM4103889 | Patient 2 | Plate 1 | JE5_AHFL77BGX5_S6_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM4103890 | Patient 2 | Plate 2 | JE6_AHFL77BGX5_S7_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM4103891 | Patient 2 | Plate 3 | JE7_AHFL7NBGX5_S16_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM4103892 | Patient 2 | Plate 4 | JE8_AHFL7NBGX5_S17_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM4103893 | Patient 3 | Plate 1 | HUB-MW-005_AH32W2BGX9_S5_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM4103894 | Patient 3 | Plate 2 | HUB-MW-006_AH32W2BGX9_S6_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM4103895 | Patient 3 | Plate 3 | HUB-MW-007_HC3GFBGX9_S6_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM4103896 | Patient 3 | Plate 4 | HUB-MW-008_HC3GFBGX9_S7_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM4103897 | Patient 4 | Plate 1 | HUB-JE-010_HGVN3BGX9_S1_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM4103898 | Patient 4 | Plate 2 | HUB-JE-011_HGVN3BGX9_S2_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM4762107 | Patient 5 | Plate 1 | HUB-AL-s001_HG25TBGXF_S5_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM4762108 | Patient 5 | Plate 2 | HUB-AL-s002_HG25TBGXF_S6_cat_pT.nonRibo_E99_Aligned.out.counts.table.tsv |
 
 
#### FACS index data
The file `FACS_INDEX_data.zip` contains additional recorded data from FACS sorting, related to single cell sizes.
 
Specifically, the file `FSCA__indexData_MW.Rdata` in this archive contains an R object, which contains FACS information (such as FSC_A, FSC_W, FSC_H, SSC_A, SSC_W, SSC_H), and an entry `Cell_newname`, which relates the observations to cellnames in the Seurat objects.
 
- Additional files, where filenames start with a prefix indicating the GEO id:
    - GSM4103897_indexdata_patient4-plate1__2018-244_exp_cardiologie_plate2.csv
   - GSM4103898_indexdata_patient4-plate2__2018-244_exp_cardiologie_plate3merged.csv
   - GSM4762107_indexdata_patient5-plate1_AL1.csv
   - GSM4762108_indexdata_patient5-plate2_AL2.csv
 
#### Count tables Hu et al.
The file `counttables-HU.zip` contains the following files:

| GEO id | Filename |
|-------------|--------------------------------------------------------------------------|
| GSM2970361 | p.N1.plate.97474.part.1_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM2970361 | p.N1.plate.97474.part.2_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM2970358 | p.N2.plate.97452.part.1_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM2970358 | p.N2.plate.97452.part.2_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM2970358 | p.N2.plate.97493.part.1_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM2970358 | p.N2.plate.97493.part.2_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM2970362 | p.N3.plate.97438.part.1_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM2970362 | p.N3.plate.97438.part.2_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM2970366 | p.N4.plate.97461.part.1_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM2970366 | p.N4.plate.97461.part.2_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM2970360 | p.N5.plate.97458.part.1_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM2970360 | p.N5.plate.97458.part.2_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM3449619 | p.N13.plate.100355.part.1_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM3449619 | p.N13.plate.100355.part.2_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM3449620 | p.N14.plate.104720.part.1_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |
| GSM3449620 | p.N14.plate.104720.part.2_cat_nc.nonRibo_E99_Aligned.out.counts.table.tsv |

This is data from Wang et al (re-calculated count tables determined by Wehrens et al. paper from selected data as described in the Wehrens et al. methods section).
 
#### Seurat object (Van Rooij only)
 
The file: `H5_RHL_SeuratObject_nM_sel_ROOIJonly.sp.bt_RID2l_clExtended.Rds` contains a Seurat object with only the processed data from Wehrens et al samples:
```
GSM4103885, GSM4103886, GSM4103887, GSM4103888, GSM4103889, GSM4103890, GSM4103891, GSM4103892,
GSM4103893, GSM4103894, GSM4103895, GSM4103896, GSM4103897, GSM4103898, GSM4762107, GSM4762108.
```
 
This Seurat object was used to generate the figures in the paper.
 
#### Seurat object (all)
 
The file `H5_RHL_SeuratObject_nM_sel_ALL.SP_btypSel_RID2l_clExtended.Rds` contains data processed using Seurat as described in Wehrens et al., Wang et al., and Litvinukova et al. Data processing is described in the methods section.
 
This file contains data from:
```
GSM4103885, GSM4103886, GSM4103887, GSM4103888, GSM4103889, GSM4103890, GSM4103891, GSM4103892, GSM4103893, GSM4103894, GSM4103895, GSM4103896, GSM4103897, GSM4103898, GSM4762107, GSM4762108.
GSM2970361, GSM2970361, GSM2970358, GSM2970358, GSM2970358, GSM2970358, GSM2970362, GSM2970362, GSM2970366, GSM2970366, GSM2970360, GSM2970360, GSM3449619, GSM3449619, GSM3449620, GSM3449620
ERP123138
```
 
This Seurat object was used to generate the figures in the paper.

## Easy use of this data for Van Rooij lab members

See the script 
`./Projects/howtousedata_HCMSCS/example-load-plot-data.R` and the folder `srv-lnx-varo1:/opt/backup_wehrens/data/Wehrens2022/Rdata-important` to make some plots using this data.

## About the scripts

The most central file in this analysis is `Projects/HCM_SCS_2021_06_SeuratRevisedAnalysis_v2_UmiTools.R`, you can start exploring these scripts there.

*More explanation will be added.*

Send an e-mail to m.wehrens@hubrecht.eu if you have any questions.

## Notes about libraries etc. 

The scripts involving the Litvinukova data were executed on a high performance computing cluster.

The necessary software and libraries were installed using Conda.

Four separate environments were created, "base" (default), "MARGE", "LISA", "pyscenic". The base environment was
used for all analysis, except for the MARGE, Lisa and SCENIC analyses. 
A list of installed libraries can be found in the text file "libraries_environments_information.txt", 
this also includes a separate list of all packages that were installed "inside" R.


## Contact information

Either contact Prof. Van Rooij (e.vanrooij@hubrecht.eu) or Martijn Wehrens (m.wehrens@uva.nl) in case you have questions.















