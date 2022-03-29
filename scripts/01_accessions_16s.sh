#!/usr/bin/env bash

# This script grabs the first 105 accession ids from the Metadata file arranged by country and writes it to a new .txt file.
# Countries in the Metadata: usa_wc, italy_milan, france, germany, india, northamerica_ny, japan_tokyo, usa_nashville, china


country_metagenome=$1

input_file="../accessions/microbiome_data.csv"
echo "Input File: $input_file"

echo
echo "Processing Country: $country_metagenome"
if [ country_metagenome == "usa_wc" ]
then
{
  output_file="../accessions/accessions_usa_wc_16s.txt";
  field_num=2
  echo "Output Filename: $output_file"
  echo "Field Number: $field_num"
}
elif [ country_metagenome == "italy_milan" ]
then
{
  output_file="../accessions/accessions_italy_milan_16s.txt"
  field_num=3
}
elif [ country_metagenome == "france" ]
then
{
  output_file="../accessions/accessions_france_16s.txt"
  field_num=4
}
elif [ country_metagenome == "germany" ]
then
{
  output_file="../accessions/accessions_germany_16s.txt"
  field_num=5
}
elif [ country_metagenome == "india" ]
then
{
  output_file="../accessions/accessions_india_16s.txt"
  field_num=6
}
elif [ country_metagenome == "northamerica_ny" ]
then
{
  output_file="../accessions/accessions_northamerica_ny_16s.txt"
  field_num=7
}
elif [ country_metagenome == "japan_tokyo" ]
then
{
  output_file="../accessions/accessions_japan_tokyo_16s.txt"
  field_num=8
}
elif [ country_metagenome == "usa_nashville" ]
then
{
  output_file="../accessions/accessions_usa_nashville_16s.txt"
  field_num=9
}
elif [ country_metagenome == "china" ]
then
{
  output_file="../accessions/accessions_china_16s.txt"
  field_num=10
}

fi

echo "Output Filename: $output_file"
echo "Field Number: $field_num"

num_accessions=105
echo "Number of Accessions to Grab from SRA Run Table: $num_accessions"

accessions_plus_header=$(($num_accessions + 1))
echo "Accessions including Header line from SRA Run Table: $accessions_plus_header"

#head -n $accessions_plus_header $input_file | tail -n $num_accessions | cut -f$field_num -d','
