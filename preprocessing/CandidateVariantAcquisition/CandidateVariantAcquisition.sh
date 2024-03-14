#!/bin/bash


# Where are the input files?"
INPUT_PATH="/Path/To/Cosmic/Public/Files/CosmicCodingMuts.vcf.gz"

# Where are the results supposed to go?
RESULT_PATH="/Path/To/Output/Folder/"

# In how many samples should a variant have been detected?
# This values is only used to subset variants that are from a gene or a genomic region.
# Variants retrieved by their COSMIC ID are always returned.
SAMPLE_COUNT=10

# The list of genes of interest.
# The GenesOfInterest.txt file should only contain one column with no header. The column should contain the gene names.
PATH_TO_GENES_OF_INTEREST="Example_Genes_Of_Interest.tsv"

# COSMIC Variants of interest.
PATH_TO_COSMIC_IDS="Example_COSMIC_IDs.tsv"

# The regions of interest.
PATH_TO_REGIONS="Example_Regions_Of_Interest.tsv"

echo "The used variables:"
echo "The input files:     ${INPUT_PATH}"
echo "The output folder:   ${RESULT_PATH}"
echo "The candidate genes: ${PATH_TO_GENES_OF_INTEREST}"
echo "The COSMIC IDs:      ${PATH_TO_COSMIC_IDS}"
echo "The regions:         ${PATH_TO_REGIONS}"




echo "We get the header of the COSMIC file."
tabix -p vcf ${INPUT_PATH} # We index the vcf.gz file. -p means that we set the input format as vcf.
tabix -H ${INPUT_PATH} > ${RESULT_PATH}/CosmicSubset_Intermediary.vcf # -H only get the header.

echo "We get variants for the genes of interest."
readarray -t genes_array < ${PATH_TO_GENES_OF_INTEREST}
for gene in ${genes_array[*]}; do
    echo ${gene}
    bcftools view -H -i "CNT>=${SAMPLE_COUNT} && GENE='${gene}'" ${INPUT_PATH} >> ${RESULT_PATH}/CosmicSubset_Intermediary.vcf
done

echo "We get variants by their COSMIC IDs."
readarray -t cosmic_array < ${PATH_TO_COSMIC_IDS}
for cosmic in ${cosmic_array[*]}; do
    echo ${cosmic}
    bcftools view -H -i "ID='${cosmic}'" ${INPUT_PATH} >> ${RESULT_PATH}/CosmicSubset_Intermediary.vcf
done

echo "We get variants in specified regions."
bcftools view -H -R ${PATH_TO_REGIONS} ${INPUT_PATH} >> ${RESULT_PATH}/CosmicSubset_Intermediary.vcf

echo "We remove any duplicates."
bcftools norm --rm-dup all ${RESULT_PATH}/CosmicSubset_Intermediary.vcf > ${RESULT_PATH}/CosmicSubset.vcf
rm ${RESULT_PATH}/CosmicSubset_Intermediary.vcf
