#!/bin/bash

echo "Running script:"
echo "AddReadGroup.sh"
echo "Starting time:"
date

# Variables
# Where do you want the output to be?
output="/Path/To/Your/Output/Folder/"
# What is the sample you are using?
sample_use="SampleID"
# Where is the cell barcodes file?
cells_use="/Path/To/Your/CellRanger/barcodes.tsv"
# Where is the respective molecule_info.h5 file from cellranger?
molecule_info_use="/Path/To/Your/molecule_info.h5"
# Where is the original BAM file?
bams_use="/Path/To/Your/CellRanger/possorted_genome_bam.bam"
# How many reads should a UMI have in a specific cell?
# UMIs in different cells are not relevant.
min_reads_per_umi_use=10
# What is the length of your UMIs?
# 10 for 3'v2 or 5' chemistry and 12 for 3'v3 chemistry.
umi_length_use=12

echo "Using output:          ${output}"
echo "Doing sample:          ${sample_use}"
echo "Doing cells:           ${cells_use}"
echo "Doing molecule info:   ${molecule_info_use}"
echo "Doing BAM:             ${bams_use}"
echo "Minimum Reads per UMI: ${min_reads_per_umi_use}"
echo "UMI length:            ${umi_length_use}"

# Creating the output folder for the sample.
mkdir -p ${output_path}/$sample_use/
cd ${output_path}/${sample_use}/

echo "We subset the UMIs."
Rscript ~/sigurd/preprocessing/UMI_Preprocessing/ReadMoleculeInfo.R \
	--molecule_info $molecule_info_use \
	--barcodes_path $cells_use \
	--output $output_path \
	--sample $sample_use \
	--min_reads_per_umi $min_reads_per_umi_use \
	--umi_length $umi_length_use


python ~/sigurd/preprocessing/UMI_Preprocessing/AddReadGroup.py \
       $bam_use \
       $output_path/${sample_use}/ \
       ${output_path}/$sample_use/UGs_Filtered.csv
samtools index $output_path/${sample_use}/possorted_genome_UG_tagged.bam

echo "Finishing time:"
date
