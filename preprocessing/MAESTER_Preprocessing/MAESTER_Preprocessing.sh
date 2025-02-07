#!/bin/bash
echo "Running script:"
echo "MAESTER_preprocessing.sh"
date

# Variables
# Input folder of the FASTQ files.
INPUT_FOLDER_FASTQ="/Path/To/Your/FASTQ/Files/"
# Where should the results be saved?
# Per sample a directory is created inside this folder.
OUTPUT_FOLDER_PATH="/Path/To/Your/Output/Folder/"
# Sample to use.
SAMPLE="SampleID"
# Path to the cell barcodes.
CELL_BARCODES_PATH="/Path/To/Your/CellRanger/barcodes.tsv"
# The length of the cell barcodes supplied in CELL_BARCODES_PATH.
# 10X uses 16.
CB_LENGTH=16
# The length of the UMIs.
# 10X uses 12.
UMI_LENGTH=12
# Where is the reference genome located?
REFERENCE_GENOME_DIRECTORY="/Path/To/The/Reference/Genome/"
# Memory in GB
MEMORY_IN_GB=64
# What is the mitochondrial chromosome name? MT or chrM?
MT_NAME="MT"
BWA_INDEX="bwa_index/genome_MT/genome_MT.fa"
# How many cores do you want to use?
NCORES=1
# How many barcodes should a cell have to be included?
MIN_BARCODE_READS=3
# How many reads should a UMI in a cell have to be included?
MIN_READS=3
# What is the read quality threshold?
BASE_QUALITY=30
# What is the alignment quality threshold?
MAPQ=30
# Do you want to remvoe unneeded files after running the pipeline? Yes or No.
CLEAN_UP="No"



echo "We assemble the FASTQ file for sample ${SAMPLE}."
mkdir -p ${OUTPUT}/${SAMPLE}
cd ${OUTPUT}${SAMPLE}
Rscript "AssembleFASTQ.R" \
  --Input_Folder_Path ${INPUT_FOLDER_FASTQ} \
  --Sample ${SAMPLE} \
  --Cell_Barcodes_Path ${CELL_BARCODES_PATH} \
  --CB_Length ${CB_LENGTH} \
  --UMI_Length ${UMI_LENGTH} \
  --Output_Folder_Path "${OUTPUT_FOLDER_PATH}${SAMPLE}/"


echo "We trim the FASTQ file."
homerTools trim -5 24 "${OUTPUT_FOLDER_PATH}${SAMPLE}/${SAMPLE}.fastq.gz"
# homerTools does not produce a nice output file. We rename it.
mv "${OUTPUT_FOLDER_PATH}${SAMPLE}/${SAMPLE}.fastq.gz.trimmed" "${OUTPUT_FOLDER_PATH}${SAMPLE}/${SAMPLE}_trimmed.fastq"
gzip "${OUTPUT_FOLDER_PATH}${SAMPLE}/${SAMPLE}_trimmed.fastq"


echo "Alignment using STAR."
STAR --genomeDir ${REFERENCE_GENOME_DIRECTORY} \
     --runThreadN ${NCORES} \
     --readFilesIn "${OUTPUT_FOLDER_PATH}${SAMPLE}/${SAMPLE}_trimmed.fastq.gz" \
     --outFileNamePrefix ${OUTPUT_FOLDER_PATH}${SAMPLE}/ \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM \
     --readFilesCommand zcat \
     --limitBAMsortRAM ${MEMORY_IN_GB}000000000
samtools index -@ ${NCORES} "${OUTPUT_FOLDER_PATH}${SAMPLE}/Aligned.sortedByCoord.out.bam"


echo "Converting ${OUTPUT_FOLDER_PATH}${SAMPLE}/Aligned.sortedByCoord.out.bam into ${OUTPUT_FOLDER_PATH}${SAMPLE}/Aligned.sortedByCoord.out.10x.bam"
samtools view -h ${INPUT} | awk 'BEGIN{FS="\t"; OFS="\t"} {
	if (substr($1,1,1) == "@") {
		print $0
	} else {
		split($1, a, "_")
		$1=""
		print a[1]"_"a[2]$0"\tCB:Z:"a[3]"-1\tUB:Z:"a[4]
	} }' | samtools view -bh > "${OUTPUT_FOLDER_PATH}${SAMPLE}/Aligned.sortedByCoord.out.10x.bam"
samtools index -@ ${NCORES} "${OUTPUT_FOLDER_PATH}${SAMPLE}/Aligned.sortedByCoord.out.10x.bam"


echo "We subset the BAM file to only include mitochondrial reads, remove unmapped reads and reads with a mapping quality below ${MAPQ}."
samtools view -@ ${NCORES} -b -F 4 -q ${MAPQ} -o "${OUTPUT_FOLDER_PATH}${SAMPLE}/Aligned.sortedByCoord.out.10x.MT.bam" \
	 "${OUTPUT_FOLDER_PATH}${SAMPLE}/Aligned.sortedByCoord.out.10x.bam" ${MT_NAME}
samtools index -@ 16 "${OUTPUT_FOLDER_PATH}${SAMPLE}/Aligned.sortedByCoord.out.10x.MT.bam"


echo "We start MAEGATK."
maegatk bcall \
       --ncores ${NCORES} \
       --snake-stdout \
       --min-barcode-reads ${MIN_BARCODE_READS} \
       --min-reads ${MIN_READS} \
       --barcodes ${CELL_BARCODES_PATH} \
       --barcode-tag "CB" \
       --umi-barcode "UB" \
       --input "${OUTPUT_FOLDER_PATH}${SAMPLE}/Aligned.sortedByCoord.out.10x.MT.bam" \
       --output ${OUTPUT} \
       --base-qual ${BASE_QUALITY} \
       --alignment-quality ${MAPQ} \
       --mito-genome ${BWA_INDEX}


if [ ${CLEAN_UP} = "Yes" ]; then
  echo "Removing unneeded files."
  rm "${OUTPUT_FOLDER_PATH}${SAMPLE}/${SAMPLE}_trimmed.fastq.gz"
  rm "${OUTPUT_FOLDER_PATH}${SAMPLE}/Aligned.sortedByCoord.out.bam"
  rm "${OUTPUT_FOLDER_PATH}${SAMPLE}/Aligned.sortedByCoord.out.10x.bam"
fi
