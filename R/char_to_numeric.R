#'char_to_numeric
#'@description
#'A function to convert the heterozygous/homozygous information from the VCF to the consensus information from VarTrix.
#'It is only used in LoadingVCF_typewise.R.
#'@param char_value What is the genotype encoding you want to convert? 
#'@export
char_to_numeric <- function(char_value) {
  if(char_value == "1/1") return(2)
  if(char_value %in% c("1/0", "0/1")) return(2)
  if(char_value == "0/0") return(1)
  return(0)
}
