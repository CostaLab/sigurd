#!/bin/bash

echo "Supply a file with the genes, the COSMIC IDs, the regions and the histologies for which variants should be selected."
echo "Simply set the value of a file is NO to skip the part."

echo "Where are the input files?"
echo "The following files from COSMIC are necessary:"
echo "CosmicCodingMuts.vcf.gz"
echo "CosmicGenomeScreensMutantExport.tsv, for the histology selection."
COSMIC_CODING_MUTATIONS="/Path/To/Cosmic/Public/Files/Cosmic_GenomeScreensMutant_Normal_v99_GRCh38.vcf.gz"
COSMIC_CANCER_MUTATION_CENSUS="/Path/To/Cosmic/Public/Files/cmc_export.tsv"

echo "Where are the results supposed to go?"
RESULT_PATH="/Path/To/Output/Folder/"

echo "In how many samples should a variant have been detected?"
echo "This values is only used to subset variants that are from a gene or a genomic region."
echo "Variants retrieved by their COSMIC ID are always returned."
SAMPLE_COUNT=10

echo "The list of genes of interest."
echo "The GenesOfInterest.txt file should only contain one column with no header. The column should contain the gene names."
PATH_TO_GENES_OF_INTEREST="Example_Genes_Of_Interest.tsv"

echo "COSMIC IDs for candidate variants."
PATH_TO_COSMIC_IDS="Example_COSMIC_IDs.tsv"

echo "The regions of interest."
PATH_TO_REGIONS="Example_Regions_Of_Interest.tsv"

echo "The histologies of interest."
PATH_TO_HISTOLOGIES="Example_Histologies_Of_Interest.tsv"

echo "The used variables:"
echo "The coding mutations: ${COSMIC_CODING_MUTATIONS}"
echo "The CMC:              ${COSMIC_CANCER_MUTATION_CENSUS}"
echo "The output folder:    ${RESULT_PATH}"
echo "The sample count:     ${SAMPE_COUNT}"
echo "The candidate genes:  ${PATH_TO_GENES_OF_INTEREST}"
echo "The COSMIC IDs:       ${PATH_TO_COSMIC_IDS}"
echo "The regions:          ${PATH_TO_REGIONS}"
echo "The histologies:      ${PATH_TO_HISTOLOGIES}"




echo "We index the VCF file."
tabix -p vcf ${COSMIC_CODING_MUTATIONS} # We index the vcf.gz file. -p means that we set the input format as vcf.
bcftools view -h ${COSMIC_CODING_MUTATIONS} > ${RESULT_PATH}/CosmicSubset_Intermediary.vcf

if [ "${PATH_TO_GENES_OF_INTEREST}" != "NO" ]; then
  echo "We get variants for the genes of interest."
  readarray -t genes_array < ${PATH_TO_GENES_OF_INTEREST}
  for gene in ${genes_array[*]}; do
      echo ${gene}
      bcftools view -H -i "SAMPLE_COUNT>=${SAMPLE_COUNT} && GENE='${gene}'" ${COSMIC_CODING_MUTATIONS} >> ${RESULT_PATH}/CosmicSubset_Intermediary.vcf
  done
fi

if [ "${PATH_TO_COSMIC_IDS}" != "NO" ]; then
  echo "We get variants by their COSMIC IDs."
  readarray -t cosmic_array < ${PATH_TO_COSMIC_IDS}
  for cosmic in ${cosmic_array[*]}; do
      echo ${cosmic}
      bcftools view -H -i "ID='${cosmic}'" ${COSMIC_CODING_MUTATIONS} >> ${RESULT_PATH}/CosmicSubset_Intermediary.vcf
  done
fi

if [ "${PATH_TO_REGIONS}" != "NO" ]; then
  echo "We get variants in specified regions."
  bcftools view -H -R ${PATH_TO_REGIONS} ${COSMIC_CODING_MUTATIONS} >> ${RESULT_PATH}/CosmicSubset_Intermediary.vcf
fi

if [ "${PATH_TO_HISTOLOGIES}" != "NO" ]; then
  echo "We get variants for specified histologies."
  readarray -t histologies_array < ${PATH_TO_HISTOLOGIES}
  for histology in ${histologies_array[*]}; do
      echo ${histology}
      grep ${histology} ${COSMIC_CANCER_MUTATION_CENSUS} | awk -F '\t' '{print $19}' | xargs -I{} grep {} ${COSMIC_CODING_MUTATIONS} >> ${RESULT_PATH}/CosmicSubset_Intermediary.vcf
  done
fi

echo "We remove any duplicates."
bcftools norm --rm-dup all ${RESULT_PATH}/CosmicSubset_Intermediary.vcf > ${RESULT_PATH}/CosmicSubset.vcf
rm ${RESULT_PATH}/CosmicSubset_Intermediary.vcf
