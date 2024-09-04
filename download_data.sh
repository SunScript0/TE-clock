#!/bin/bash

mkdir -p data
mkdir -p data/annotation
mkdir -p data/infinium
mkdir -p data/processed

wget https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv -O data/annotation/full_probe_info_450.csv
tail -n +8 data/annotation/full_probe_info_450.csv | head -n -851 > data/annotation/short_probe_info_450.csv

# Read the TSV file line by line
while IFS=$'\t' read -r col1 col2 col3
do
    # Create the folder structure
    mkdir -p "data/$col1/$col2"
    # wget and gunzip files
    wget -P "data/$col1/$col2" "$col3"
	gunzip "data/$col1/$col2/*.gz"
done < "data_sources.tsv"