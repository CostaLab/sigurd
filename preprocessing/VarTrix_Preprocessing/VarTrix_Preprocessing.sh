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
MIN_READS_PER_UMI=2
# The reference fasta file.
REFERENCE_FASTA="/Path/To/Your/gene_annotation/refdata-gex-GRCh38-2020-A/fasta/genome.fa"

# These are your positions of interest.
# These files are provided as example data.
GENES_OF_INTEREST="../../inst/extdata/MT_Input_VCF_NoMAF_Filtering.vcf"
GENES_OF_INTEREST="../../inst/extdata/chrM_Input_VCF_NoMAF_Filtering.vcf"
GENES_OF_INTEREST="../../inst/extdata/JAK2V617F.vcf"


echo "We create the necessary folders, if they don't already exist."
mkdir -p ${OUTPUT}/${SAMPLE}

echo "We generate a BED file to subset the BAM file."
bcftools query -f'%CHROM\t%POS0\t%POS\t%REF,%ALT\n' ${GENES_OF_INTEREST} > ${OUTPUT}/${SAMPLE}/genes_of_interest.bed

echo "We subset the BAM files."
samtools view -@ ${NCORES} -h -b -F 4 -q ${MAPQ} -L ${OUTPUT}/${SAMPLE}/genes_of_interest.bed -o ${OUTPUT}/${SAMPLE}/bam_filtered.bam ${BAM}
samtools index -@ ${NCORES} ${OUTPUT}/${SAMPLE}/bam_filtered.bam

echo "We get the loci names and the chromosomes present in the input VCF file."
awk '{print $1,$2,$4,$5}' ${GENES_OF_INTEREST} > ${OUTPUT}/SNV.loci.txt
sed -i 's/\s/:/g' ${OUTPUT}/SNV.loci.txt

echo "We add the molecule tag and count the number of apperances."
# It consists of the cell barcode (CB) and the UMI (UB).
python Get_Reads_Per_Molecule.py \
       ${OUTPUT}/${SAMPLE}/bam_filtered.bam \
       ${OUTPUT}/${SAMPLE}/ \
       "bam_tagged" \
       ${SAMPLE}
samtools index ${OUTPUT}/${SAMPLE}/bam_tagged.bam
rm ${OUTPUT}/${SAMPLE}/bam_filtered.bam ${OUTPUT}/${SAMPLE}/bam_filtered.bam.bai

echo "We get only molecules with at least ${MIN_READS_PER_UMI} reads."
awk -F, -v threshold="$MIN_READS_PER_UMI" 'NR > 1 && $5 >= threshold {print $1}' "${OUTPUT}/${SAMPLE}/bam_tagged.csv" > "${OUTPUT}/${SAMPLE}/bam_tagged_subset.csv"

echo "We filter the reads."
samtools view -@ ${NCORES} -h -b -F 4 -q ${MAPQ} --tag-file=UG:${OUTPUT}/${SAMPLE}/bam_tagged_subset.csv -o ${OUTPUT}/${SAMPLE}/bam_filtered_tagged.bam ${OUTPUT}/${SAMPLE}/bam_tagged.bam
samtools index -@ ${NCORES} ${OUTPUT}/${SAMPLE}/bam_filtered_tagged.bam
rm ${OUTPUT}/${SAMPLE}/bam_tagged.bam ${OUTPUT}/${SAMPLE}/bam_tagged.bam.bai

echo "Coverage"
vartrix_linux --threads ${NCORES} --mapq ${MAPQ} --bam ${OUTPUT}/${SAMPLE}/bam_filtered_tagged.bam --umi \
              --vcf ${GENES_OF_INTEREST} --scoring-method "coverage" \
              --ref-matrix ${OUTPUT}/${SAMPLE}/ref_matrix_coverage.mtx \
              --cell-barcodes ${CELLS} \
	      --fasta ${REFERENCE_FASTA} \
              --out-matrix ${OUTPUT}/${SAMPLE}/out_matrix_coverage.mtx
echo "Consensus"
vartrix_linux --threads ${NCORES} --mapq ${MAPQ} --bam ${OUTPUT}/${SAMPLE}/bam_filtered_tagged.bam --umi \
              --vcf ${GENES_OF_INTEREST} --scoring-method "consensus" \
              --cell-barcodes ${CELLS} \
	      --fasta ${REFERENCE_FASTA} \
              --out-matrix ${OUTPUT}/${SAMPLE}/out_matrix_consensus.mt
echo "VarTrix is finished"

echo "Finishing time:"
date
