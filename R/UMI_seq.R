#'UMI_seq
#'@description
#'This function converts the numeric UMI representation found in a 10X molecule_info.h5 file to a character string.
#'@param number The numeric expression of a UMI.
#'@param umi_length The length of the UMI. Default = 12
#'@export
UMI_seq <- function(number, umi_length = 12){
  if(is.factor(number)) stop("Your input UMI is a factor.")
  if(!is.numeric(number)) stop("Your UMI is not numeric.")
  if(any(number <= 0)) stop("Your input UMI is <= 0.")
  if(!is.numeric(umi_length)) stop("Your umi_length is not numeric.")
  if(umi_length <= 0) stop("Your umi_length is <= 0.")
  string_length <- umi_length * 2
  bases <- c("00" = "A", "01" = "C", "10" = "G", "11" = "T")
  result <- lapply(number, function(x){
    # First, we convert the number to a vector of binaries.
    # Then we convert the binary numbers to numerics.
    # We only take the first umi_length*2, since that is the number of bits we need to encode a UMI of this length.
    # Each position needs two bits since we have 4 options.
    # Then we revert the vector, since the intToBits function returns the least important values first.
    # This basically means, that the output from intToBits is in the wrong order.
    # Then we paste every first with every second element together.
    # Then we replace the new pseudo bits with their according base.
    res <- rev(as.numeric(intToBits(x))[1:string_length])
    res <- paste0(res[seq(1, length(res), by = 2)], res[seq(2, length(res), by = 2)])
    res <- paste0(bases[res], collapse = "")
    return(res)
  })
  result <- unlist(result)
  return(result)
}
