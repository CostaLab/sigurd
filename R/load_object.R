#'load_object
#'@description
#'A loading function to load the RDS files quicker.
#'Source: https://github.com/CostaLab/CimpleG
#'@import archive
#'@param file_name The path to the file.
#'@export
load_object <- function(file_name){
  con <- archive::file_read(file = file_name)
  res <- readRDS(file = con)
  close(con)
  return(res)
}
