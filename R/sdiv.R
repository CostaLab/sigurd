#'Division of sparse matrix.
#'@importFrom Matrix sparseMatrix
#'@param X First sparse matrix.
#'@param Y Second sparse matrix.
#'@param names The dimension names (dimnames(X)).
#'@export
sdiv <- function(X, Y, names = dimnames(X)) {
  sX <- summary(X)
  sY <- summary(Y)
  sRes <- merge(sX, sY, by = c("i", "j"))
  result <- Matrix::sparseMatrix(i = sRes[,1], j = sRes[,2], x = sRes[,3] / sRes[,4], dimnames = names)
  return()
}
