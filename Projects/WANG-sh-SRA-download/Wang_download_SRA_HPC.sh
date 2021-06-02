#!/bin/bash
#

# This script will use the SRA toolkit to download data of interest.
#
# It will require that you make directory with the name <IDENTIFIER>,
# which contains an accession list containing the SRA IDs (obtain
# them from [https://www.ncbi.nlm.nih.gov/Traces/study/?]), accession list should
# be named like <IDENTIFIER>_SRR_Acc_List.txt.
#
# @HPC, use e.g. to execute:
# sbatch Wang_download_SRA.sh
#
# Locally, use 
# sh /Users/m.wehrens/Documents/git_repos/SCS_More_analyses/Projects/Wang_download_SRA.sh
#
# Edit basepath appropriately below

scriptpath='./'

datapath='/hpc/hub_oudenaarden/mwehrens/fastq/WANG/'
#datapath=/Volumes/workdrive_m.wehrens_hubrecht/Fastq__raw_files/Wang/

for IDENTIFIER in GSM2970361_N1_LV GSM2970358_N2_LV GSM2970362_N3_LV GSM2970366_N4_LV GSM2970360_N5_LV GSM2970359_N6_LA GSM2970369_N7_LA GSM2970368_N8_LA GSM2970365_N9_LA GSM2970363_N10_LA GSM2970364_N11_LA GSM2970367_N12_LA GSM3449619_N13 GSM3449620_N14
do

    sbatch --time=48:00:00 --mem=10G  --job-name=${IDENTIFIER} --export=ALL,datapath="${datapath}",IDENTIFIER="${IDENTIFIER}" Wang_download_SRA_onesubset.sh 

done




