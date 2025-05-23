#'load_object
#'@description
#'A loading function to load the RDS files quicker.
#'Source: https://github.com/CostaLab/CimpleG
#'@import archive
#'@param file_name The path to the file.
#'@export
load_object <- function(file_name){
  if(!file.exists(file_name)) stop(paste0("File '",file_name,"' not found."))
  
  if(requireNamespace("archive", quietly = TRUE)){
    con <- archive::file_read(file = file_name)
    res <- readRDS(file = con)
    close(con)
    return(res)
  }
  res <- readRDS(file = file_name)
}
