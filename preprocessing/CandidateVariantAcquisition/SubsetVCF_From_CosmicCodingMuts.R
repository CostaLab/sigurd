# We load the VCF file and subset it.

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(VariantAnnotation))


option_list <- list(
  optparse::make_option("--input_file",                     type = "character", default = "",   help = "The input tsv file."),
  optparse::make_option("--output_file",                    type = "character", default = "",   help = "The output vcf file."),
  optparse::make_option("--output_file_UCSC",               type = "character", default = "",   help = "The output vcf file using the UCSC format."),
  optparse::make_option("--sample_count",                   type = "numeric",   default = 10,   help = "The minimum sample count."),
  optparse::make_option("--remove_alternative_transcripts", type = "logical",   default = TRUE, help = "Should alternative transcripts be removed?")
)
opt_parser <- optparse::OptionParser(option_list=option_list)
opt <- optparse::parse_args(opt_parser)
input_file <- opt$input_file
output_file <- opt$output_file
output_file_UCSC <- opt$output_file_UCSC
minimum_sample_count <- opt$sample_count
remove_alternative_transcripts <- opt$remove_alternative_transcripts


print("We load the VCF file.")
vcf <- VariantAnnotation::readVcf(input_file)


print("We remove all variants with a sample count below 10.")
vcf <- vcf[info(vcf)$CNT >= minimum_sample_count, ]


if(remove_alternative_transcripts){
  print("We remove the ENST genes.")
  vcf <- vcf[!grepl("_ENST", info(vcf)$GENE),]
}


print("We save the filtered VCF file.")
VariantAnnotation::writeVcf(vcf, output_file)


print("We change the sequence level style to UCSC.")
GenomeInfoDb::seqlevelsStyle(vcf) <- "UCSC"
VariantAnnotation::writeVcf(vcf, output_file_UCSC)
