#!/bin/bash

echo "Running script:"
echo "VarTrix_preprocessing.sh"
echo "Starting time:"
date

# Variables
# How many cores do you want to use?
NCORES=1
# Where do you want the output to be?
OUTPUT="/data/mg000001/sigurd_results/test/"
# What is the read quality threshold?
MAPQ=30
# The sample to be analysed.
SAMPLE="AB58_J"
# The BAM file to be analysed.
BAM="/data/Schneider_lab/biopsy_MPN_CML/MF/data/amplicon_Cellranger_output/Trimming_WithIndex_WithUGSubsetting/AB58_J/outs/possorted_genome_bam.bam"
# How many reads should a UMI have in a specific cell?
# UMIs in different cells are not relevant.
MIN_READS_PER_UMI=10
# What is the length of your UMIs?
# 10 for 3'v2 or 5' chemistry and 12 for 3'v3 chemistry.
UMI_LENGTH=12
# The reference fasta file.
REFERENCE_FASTA="/data/mg000001/gene_annotation/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
# Where is the central input file?
MOLECULE_INPUT_PATH="/data/Schneider_lab/biopsy_MPN_CML/MF/data/amplicon_Cellranger_output/Trimming_WithIndex_WithUGSubsetting/AB58_J/outs/molecule_info.h5"

# These are your positions of interest.
# These files are provided as example data.
GENES_OF_INTEREST="/data/mg000001/MPN/Mutations/Cosmic/results/AdamBenabid_CML_MF_Mutations/20231206_AdamBenabid_VariantsOfInterest_AddedByHand.vcf"




echo "We create the necessary folders, if they don't already exist."
mkdir -p ${OUTPUT}/${SAMPLE}

echo "We get the loci names."
awk '{print $1,$2,$4,$5}' ${GENES_OF_INTEREST} > ${OUTPUT}/SNV.loci.txt
sed -i 's/\s/:/g' ${OUTPUT}/SNV.loci.txt

echo "We subset the UMIs."
Rscript ReadMoleculeInfo.R \
        --molecule_input_path ${MOLECULE_INPUT_PATH} \
        --output ${OUTPUT}/${SAMPLE}/ \
        --sample ${SAMPLE} \
        --min_reads_per_umi ${MIN_READS_PER_UMI} \
        --umi_length ${UMI_LENGTH}

echo "We filter the reads."
python AddReadGroup.py \
       ${BAM} \
       ${OUTPUT}/${SAMPLE}/possorted_genome_filtered.bam \
       ${OUTPUT}/${SAMPLE}/Molecules_${SAMPLE}_Filtered.tsv
samtools index ${OUTPUT}/${SAMPLE}/possorted_genome_filtered.bam

echo "Coverage"
vartrix_linux --threads ${NCORES} --mapq ${MAPQ} --bam ${BAM} --umi \
              --vcf ${GENES_OF_INTEREST} --scoring-method "coverage" \
              --ref-matrix ${OUTPUT}/${SAMPLE}/ref_matrix_coverage.mtx \
              --cell-barcodes ${OUTPUT}/${SAMPLE}"/Molecules_"${SAMPLE}"_Filtered.tsv" \
	      --fasta ${REFERENCE_FASTA} \
              --out-matrix ${OUTPUT}/${SAMPLE}/out_matrix_coverage.mtx
echo "Consensus"
vartrix_linux --threads ${NCORES} --mapq ${MAPQ} --bam ${BAM} --umi \
              --vcf ${GENES_OF_INTEREST} --scoring-method "consensus" \
              --cell-barcodes ${OUTPUT}/${SAMPLE}"/Molecules_"${SAMPLE}"_Filtered.tsv" \
	      --fasta ${REFERENCE_FASTA} \
              --out-matrix ${OUTPUT}/${SAMPLE}/out_matrix_consensus.mt
echo "VarTrix is finished"

echo "Finishing time:"
date
