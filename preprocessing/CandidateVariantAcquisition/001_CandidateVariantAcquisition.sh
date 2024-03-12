#!/bin/bash


echo "Where is our data?"
data_path="/Path/To/Cosmic/Files/"
echo $data_path
echo "Where are the results supposed to go?"
result_path="/Path/To/Result/Folder/"
echo $result_path

echo "The list of genes of interest."
# The GenesOfInterest.txt file should only contain one column with no header. The column should contain the gene names.
path_to_genes_of_interest="/Path/To/GenesOfInterest.tsv"
echo $path_to_genes_of_interest
readarray -t genes_array < $path_to_genes_of_interest


echo "We get the header of the COSMIC file."
tabix -p vcf ${data_path}/CosmicCodingMuts.vcf.gz # We index the vcf.gz file. -p means that we set the input format as vcf.
tabix -H ${data_path}/CosmicCodingMuts.vcf.gz > ${result_path}/CosmicSubset.vcf # -H only get the header.


# We get the lines with a gene from our list.
for gene in ${genes_array[*]}; do
    echo $gene
    grep $gene ${data_path}/CosmicCodingMuts.vcf >> ${result_path}/CosmicSubset.vcf
done


# We load the VCF file into R and subset it.
Rscript SubsetVCF_From_CosmicCodingMuts.R \
	--input_file ${result_path}/CosmicSubset.vcf \
	--output_file ${result_path}/CosmicSubset_filtered.vcf


