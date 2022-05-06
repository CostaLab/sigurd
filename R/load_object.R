#'We define a loading function to load the RDS files quicker.
#'@import archive
#'@param file_name
load_object <- function(file_name){
  con <- archive::file_read(file = file_name)
  res <- readRDS(file = con)
  close(con)
  return(res)
}
