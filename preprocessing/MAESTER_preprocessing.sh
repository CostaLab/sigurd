#!/bin/bash

echo "Running script:"
echo "MAESTER_preprocessing.sh"
date

# Variables
# How many cores do you want to use?
NCORES=7
# Where do you want the output to be?
OUTPUT="/Path/To/Your/Output/Folder/"
BWA_INDEX="bwa_index/genome_MT/genome_MT.fa"
# What is the read quality threshold?
ALIGNMENT_QUALITY=30
BASE_QUALITY=30
MIN_BARCODE_READS=3

sample_use="SampleID"
cells_use="/Path/To/Your/CellRanger/barcodes.tsv" # Not necessary, but it recommended.
bams_use="/Path/To/Your/CellRanger/possorted_genome_bam.bam"
mkdir -p $OUTPUT/$sample_use
cd ${OUTPUT}${sample_use}

maegatk bcall \
       --ncores $NCORES \
       --snake-stdout \
       --min-barcode-reads $MIN_BARCODE_READS \
       --nsamples 1024 \
       --barcode-tag "CB" \
       --umi-barcode "UB" \
       --input $bams_use \
       --barcodes $cells_use \
       --output $OUTPUT$sample_use \
       --alignment-quality $ALIGNMENT_QUALITY \
       --base-qual $BASE_QUALITY \
       --mito-genome $BWA_INDEX

# This automatically removes the temporary files.
# This is not necessary, but once overthing works reliably it removes clutter.
# rm -r .snakemake
# rm -r logs
# rm -r qc
echo "Finishing time:"
date
