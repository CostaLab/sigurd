getRefMatrix <- function(SE_object, letter, ref_alleles){
  mat <- (assays(SE_object)[[paste0(letter, "_counts_fw")]] + assays(SE_object)[[paste0(letter, "_counts_rev")]])
  rownames(mat) <- paste0(as.character(1:dim(mat)[1]), "_", toupper(ref_alleles), ">", letter)
  return(mat[toupper(ref_alleles) == letter,])
}