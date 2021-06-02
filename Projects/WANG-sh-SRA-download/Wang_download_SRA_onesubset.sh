#!/bin/bash


if [[ datapath == "" ]]; then
  datapath=$1
  IDENTIFIER=$2
fi

cd $datapath/${IDENTIFIER}/

filename=${IDENTIFIER}_SRR_Acc_List.txt

linecount=$(cat $filename | wc -l)

echo "Starting job ${IDENTIFIER} .." > log_${IDENTIFIER}.txt

n=1
while read line; do
  
  # reading each line
  echo "Line No. $n/$linecount" >> log_${IDENTIFIER}.txt
  n=$((n+1))
  
  # executing command
  fasterq-dump $line 
  
  #fastq-dump $line 
  # NOTE: as opposed to fasterq-dump, fastq-dump doesn't automatically split
  # files, you need to add some options for that; anyways, I re-installed the
  # newer version of SRA toolkit on the server now, which has fasterq-dump also,
  # instead of only fastq-dump (which was the case for the conda install).
  
done < $filename

echo "DONE. ${IDENTIFIER} should be fully downloaded.." >> log_${IDENTIFIER}.txt