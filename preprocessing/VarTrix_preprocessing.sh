#!/bin/bash

echo "Running script:"
echo "VarTrix_preprocessing.sh"
echo "Starting time:"
date

# Variables
# How many cores do you want to use?
NCORES=7
# Where do you want the output to be?
OUTPUT="/Path/To/Your/Output/Folder/"
# What is the read quality threshold?
MAPQ=30


# Input files
# The reference fasta file.
# Choose one of the following.
REFERENCE_FASTA="~/gene_annotation/refdata-gex-GRCh38-2020-A/fasta/genome_numericMT.fa"
REFERENCE_FASTA="~/bwa_index/genome_MT/genome_MT.fa"

# These are your positions of interest.
# These files are provided as example data.
GenesOfInterest="inst/extdata/ALFA_subset_MAF2_noprefix.vcf"
GenesOfInterest="inst/extdata/ALFA_subset_MAF2_prefix.vcf"
GenesOfInterest="inst/extdata/MT_Input_VCF_NoMAF_Filtering.vcf"
GenesOfInterest="inst/extdata/chrM_Input_VCF_NoMAF_Filtering.vcf"

# We create the necessary folders, if they don't already exist.
mkdir -p $OUTPUT

# We get the loci names.
awk '{print $1,$2,$4,$5}' $GenesOfInterest > $OUTPUT/SNV.loci.txt
sed -i 's/\s/:/g' $OUTPUT/SNV.loci.txt

sample_use="SampleID"
bam_use="/Path/To/Your/CellRanger/possorted_genome_bam.bam"
cells_use="/Path/To/Your/CellRanger/barcodes.tsv"
mkdir -p $OUTPUT/$sample_use

echo "Coverage"
vartrix_linux --threads $NCORES --mapq $MAPQ --bam $bam_use --umi \
              --vcf $GenesOfInterest --scoring-method "coverage" \
              --ref-matrix $OUTPUT/$sample_use/ref_matrix_coverage.mtx \
              --cell-barcodes $cells_use --fasta $REFERENCE_FASTA \
              --out-matrix $OUTPUT/$sample_use/out_matrix_coverage.mtx
echo "Consensus"
vartrix_linux --threads $NCORES --mapq $MAPQ --bam $bam_use --umi \
              --vcf $GenesOfInterest --scoring-method "consensus" \
              --cell-barcodes $cells_use --fasta $REFERENCE_FASTA \
              --out-matrix $OUTPUT/$sample_use/out_matrix_consensus.mt
echo "VarTrix is finished"

echo "Finishing time:"
date

