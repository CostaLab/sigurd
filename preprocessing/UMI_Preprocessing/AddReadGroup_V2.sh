#!/bin/bash

echo "Running script:"
echo "AddReadGroup.sh"
echo "Starting time:"
date

# Variables
# Where do you want the output to be?
output="/Path/To/Your/Output/Folder/"
# Where is the central input file?
molecule_input_path="/Path/To/Your/molecule_info.h5"
# The sample to be analysed.
sample_use="YourSample"
# How many reads should a UMI have in a specific cell?
# UMIs in different cells are not relevant.
min_reads_per_umi_use=10
# What is the length of your UMIs?
# 10 for 3'v2 or 5' chemistry and 12 for 3'v3 chemistry.
umi_length_use=12
# The BAM file to be subset.
bam_use="/Path/To/The/possorted_genome_bam.bam"

echo "Using output:             ${output}"
echo "Using Central Input File: ${central_input_file}"
echo "Minimum Reads per UMI:    ${min_reads_per_umi_use}"
echo "UMI length:               ${umi_length_use}"


echo "We subset the UMIs."
Rscript ~/sigurd/preprocessing/UMI_Preprocessing/ReadMoleculeInfo.R \
	--molecule_input_path $molecule_input_path \
	--output $output \
	--sample $sample_use \
	--min_reads_per_umi $min_reads_per_umi_use \
	--umi_length $umi_length_use


python ~/sigurd/preprocessing/UMI_Preprocessing/AddReadGroup.py \
       $bam_use \
       ${output}/${sample_use}/possorted_genome_filtered.bam \
       ${output}/${sample_use}/UGs_Filtered.csv
samtools index ${output}/${sample_use}/possorted_genome_filtered.bam

echo "Finishing time:"
date
