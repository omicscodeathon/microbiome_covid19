#!/usr/bin/env bash

# This script grabs the first 105 accession ids from the Metadata file arranged by country and writes it to a new .txt file.

input_file="../accessions/microbiome_data.csv"
echo "Input File: $input_file"

output_file="../accessions/accessions_usa_16s.txt"
echo "Input File: $output_file"

num_accessions=105
echo "Number of Accessions to Grab from SRA Run Table: $num_accessions"

accessions_plus_header=$(($num_accessions + 1))
echo "Accessions including Header line from SRA Run Table: $accessions_plus_header"

head -n $accessions_plus_header $input_file | tail -n $num_accessions | cut -f2 -d',' > $output_file
