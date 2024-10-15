#!/bin/bash

#species=( athaliana celegans dmelanogaster hsapiens pfalciparum scerevisiae )
properties=( MGW Shear Stretch Stagger Buckle ProT Opening Shift Slide Rise Tilt Roll HelT )

input_dir="/local/home/quee4387/raw_promoter_sequences/"
output_dir="/local/home/quee4387/dna_shape/"
spec=$1
echo "$spec"

for prop in "${properties[@]}"
    do
        echo "Predicting $prop in $spec"
        python /local/home/quee4387/.conda/envs/myenv/bin/deepDNAshape \
        --file ${input_dir}${spec}_200.fa \
        --feature $prop \
        --output $output_dir${spec}_${prop}_200.txt
    done


