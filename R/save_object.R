#'save_object
#'@description
#'Saving function to save the RDS files quicker.
#'Source:https://github.com/CostaLab/CimpleG
#'@import archive
#'@param object The R object to be save.
#'@param file_name The path were the file shall be save.
#'@param file_format The format of the save file. Has to be one of: zstd, lz4, gzip, bzip2, xz, nocomp.
#'@export
save_object <- function(object, file_name, file_format = "zstd"){
  stopifnot(file_format %in% c("zstd", "lz4", "gzip", "bzip2", "xz", "nocomp"))
  stopifnot(length(file_format) == 1)
  if(file_format %in% c("zstd", "lz4")){
    if(requireNamespace("archive", quietly = TRUE)){
      con <- archive::file_write(file = file_name, filter = file_format)
      open(con)
      saveRDS(object = object, file = con)
      close(con)
    }else{
      warning("Package 'archive' needs to be installed to compress files in formats 'zstd' and 'lz4'.\n Saving object with default 'saveRDS()' function instead.")
      saveRDS(object = object, file = file_name, compress = TRUE)
    }
    return(invisible(NULL))
  }
  saveRDS(
    object = object,
    file = file_name,
    compress = ifelse(file_format %in% "nocomp", FALSE, file_format)
  )
}

