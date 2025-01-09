#!/bin/bash

input_file=$1
output_file=$2

awk -F'\t' '{if($2 > $3) {t=$2; $2=$3; $3=t} print}' "${input_file}" | tr ' ' '\t' > "${output_file}_temp1"
bedtools sort -i "${output_file}_temp1" > "${output_file}_temp2"
bedtools merge -i "${output_file}_temp2" > "${output_file}_temp3"
mv "${output_file}_temp3" "${output_file}"
rm "${output_file}_temp1" "${output_file}_temp2"

