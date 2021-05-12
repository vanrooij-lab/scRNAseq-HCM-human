#!/bin/bash

cd /Volumes/workdrive_m.wehrens_hubrecht/Fastq__raw_files/Wang/

filename=GSM3449619_SRR_Acc_List.txt

linecount=$(cat $filename | wc -l)

n=1
while read line; do
  
  # reading each line
  echo "Line No. $n/$linecount"
  n=$((n+1))
  
  # executing command
  fasterq-dump $line
  
done < $filename