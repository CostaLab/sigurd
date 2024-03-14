#!/bin/bash

echo "Running script:"
echo "VarTrix_preprocessing.sh"
echo "Starting time:"
date

# Variables
# How many cores do you want to use?
NCORES=1
# Where do you want the output to be?
OUTPUT="/Path/To/Your/Output/Folder/"
# What is the read quality threshold?
MAPQ=30
# The sample to be analysed.
SAMPLE="SampleID"
# The BAM file to be analysed.
BAM="/Path/To/Your/CellRanger/possorted_genome_bam.bam"
# The cell barcodes.
CELLS="/Path/To/Your/CellRanger/barcodes.tsv"
# How many reads should a UMI have in a specific cell?
# UMIs in different cells are not relevant.
MIN_READS_PER_UMI=100
# What is the length of your UMIs?
# 10 for 3'v2 or 5' chemistry and 12 for 3'v3 chemistry.
UMI_LENGTH=12
# The reference fasta file.
REFERENCE_FASTA="/Path/To/Your/gene_annotation/refdata-gex-GRCh38-2020-A/fasta/genome_numericMT.fa"
# Where is the central input file?
MOLECULE_INPUT_PATH="/Path/To/Your/molecule_info.h5"

# These are your positions of interest.
# These files are provided as example data.
GENES_OF_INTEREST="../../inst/extdata/MT_Input_VCF_NoMAF_Filtering.vcf"
GENES_OF_INTEREST="../../inst/extdata/chrM_Input_VCF_NoMAF_Filtering.vcf"
GENES_OF_INTEREST="../../inst/extdata/JAK2V617F.vcf"



echo "We create the necessary folders, if they don't already exist."
mkdir -p ${OUTPUT}/${SAMPLE}

echo "We get the loci names and the chromosomes present in the input VCF file."
awk '{print $1,$2,$4,$5}' ${GENES_OF_INTEREST} > ${OUTPUT}/SNV.loci.txt
sed -i 's/\s/:/g' ${OUTPUT}/SNV.loci.txt
CHROMOSOMES_PRESENT=$(bcftools query -f '%CHROM\n' ${GENES_OF_INTEREST} | sort | uniq)

echo "We subset the UMIs."
Rscript ReadMoleculeInfo.R \
        --molecule_input_path ${MOLECULE_INPUT_PATH} \
        --output ${OUTPUT}/${SAMPLE}/ \
        --sample ${SAMPLE} \
	--cell_barcodes ${CELLS} \
        --min_reads_per_umi ${MIN_READS_PER_UMI} \
        --umi_length ${UMI_LENGTH}

echo "We filter the reads."
# The flag -F indicates that reads with the flag 4 (unmapped) should be removed.
# -q means that all reads with a mapping quality less than MAPQ should be removed.
# -b means that samtools should generate a BAM file as output.
# -o indicates the output file.
samtools view -b -F 4 -q ${MAPQ} -o ${OUTPUT}/${SAMPLE}/possorted_genome_mapped_reads.bam ${BAM} ${CHROMOSOMES_PRESENT}
samtools index -@ ${NCORES} ${OUTPUT}/${SAMPLE}/possorted_genome_mapped_reads.bam
python AddReadGroup.py \
       ${OUTPUT}/${SAMPLE}/possorted_genome_mapped_reads.bam \
       ${OUTPUT}/${SAMPLE}/possorted_genome_filtered.bam \
       ${OUTPUT}/${SAMPLE}/Molecules_${SAMPLE}_Filtered.tsv
samtools index ${OUTPUT}/${SAMPLE}/possorted_genome_filtered.bam
rm ${OUTPUT}/${SAMPLE}/possorted_genome_mapped_reads.bam ${OUTPUT}/${SAMPLE}/possorted_genome_mapped_reads.bam.bai

echo "Coverage"
vartrix_linux --threads ${NCORES} --mapq ${MAPQ} --bam ${OUTPUT}/${SAMPLE}/possorted_genome_filtered.bam --umi \
              --vcf ${GENES_OF_INTEREST} --scoring-method "coverage" \
              --ref-matrix ${OUTPUT}/${SAMPLE}/ref_matrix_coverage.mtx \
              --cell-barcodes ${OUTPUT}/${SAMPLE}"/CBs_"${SAMPLE}"_Filtered.tsv" \
	      --fasta ${REFERENCE_FASTA} \
              --out-matrix ${OUTPUT}/${SAMPLE}/out_matrix_coverage.mtx
echo "Consensus"
vartrix_linux --threads ${NCORES} --mapq ${MAPQ} --bam ${OUTPUT}/${SAMPLE}/possorted_genome_filtered.bam --umi \
              --vcf ${GENES_OF_INTEREST} --scoring-method "consensus" \
              --cell-barcodes ${OUTPUT}/${SAMPLE}"/CBs_"${SAMPLE}"_Filtered.tsv" \
	      --fasta ${REFERENCE_FASTA} \
              --out-matrix ${OUTPUT}/${SAMPLE}/out_matrix_consensus.mt
echo "VarTrix is finished"

echo "Finishing time:"
date
