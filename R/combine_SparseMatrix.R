#'combine_sparseMatrix
#'@description
#'We combine two sparse matrices
#'@importFrom Matrix sparseMatrix
#'@param matrix_1 Your first sparse matrix.
#'@param matrix_2 Your second matrix.
#'@export
combine_SparseMatrix <- function(matrix_1, matrix_2){
  # We have two spare matrices. We want to now join them by rows and columns.
  # We get the new row names by getting the unique row names from both matrices.
  variants_unique <- unique(c(rownames(matrix_1), rownames(matrix_2)))
  # We get a new vector with where all the rownames get a row number starting from 0.
  new_rows <- 1:length(variants_unique)
  names(new_rows) <- variants_unique

  cells_unique <- unique(c(colnames(matrix_1), colnames(matrix_2)))
  new_cols <- 1:length(cells_unique)
  names(new_cols) <- cells_unique

  # Now, we get the old rows.
  old_rows_1 <- 1:nrow(matrix_1)
  names(old_rows_1) <- rownames(matrix_1)
  old_cols_1 <- 1:ncol(matrix_1)
  names(old_cols_1) <- colnames(matrix_1)
  old_rows_2 <- 1:nrow(matrix_2)
  names(old_rows_2) <- rownames(matrix_2)
  old_cols_2 <- 1:ncol(matrix_2)
  names(old_cols_2) <- colnames(matrix_2)

  # We get the old positions.
  positions_1 <- summary(matrix_1)
  positions_1[,"i"] <- names(old_rows_1)[positions_1[,"i"]]
  positions_1[,"j"] <- names(old_cols_1)[positions_1[,"j"]]
  positions_2 <- summary(matrix_2)
  positions_2[,"i"] <- names(old_rows_2)[positions_2[,"i"]]
  positions_2[,"j"] <- names(old_cols_2)[positions_2[,"j"]]
  
  positions_1[,"i"] <- new_rows[positions_1[,"i"]]
  positions_1[,"j"] <- new_cols[positions_1[,"j"]]
  positions_2[,"i"] <- new_rows[positions_2[,"i"]]
  positions_2[,"j"] <- new_cols[positions_2[,"j"]]
  
  positions_combined <- rbind(positions_1, positions_2)
  result <- Matrix::sparseMatrix(i = positions_combined[,"i"], j = positions_combined[,"j"], x = positions_combined[,"x"], 
                                 dimnames = list(variants_unique, cells_unique), dims = c(length(variants_unique), length(cells_unique)))
}





