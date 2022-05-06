#'We define a loading function to save the RDS files quicker.
#'@import archive
#'@param object The R object to be save.
#'@param file_name The path were the file shall be save.
#'@param file_format The format of the save file. Has to be one of: zstd, lz4, gzip, bzip2, xz, nocomp.
#'@export
save_object <- function(object, file_name, file_format = NULL){
  stopifnot(file_format %in% c("zstd", "lz4", "gzip", "bzip2", "xz", "nocomp"))
  if(file_format %in% "nocomp"){
    saveRDS(object = object, file = file_name, compress = FALSE)
    return(invisible(NULL))
  }
  if(file_format %in% c("zstd", "lz4")){
    con <- archive::file_write(file = file_name, filter = file_format)
    open(con)
    saveRDS(object = object, file = con)
    close(con)
  }else{
    saveRDS(object = object, file = file_name, compress = file_format)
  }
}
