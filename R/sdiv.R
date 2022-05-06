#'Division of sparse matrix.
#'@import Matrix
#'@param X Y names
sdiv <- function(X, Y, names = dimnames(X)) {
  sX <- summary(X)
  sY <- summary(Y)
  sRes <- merge(sX, sY, by = c("i", "j"))
  sparseMatrix(i = sRes[,1], j = sRes[,2], x = sRes[,3] / sRes[,4],
               dimnames = names)
}
