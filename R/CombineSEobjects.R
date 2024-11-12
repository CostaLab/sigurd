#'CombineSEobjects
#'@description
#'We combine two SummarizedExperiment objects. This is originally intended to add the mitochondrial genotyping information to the somatic genotyping information.
#'The function can combine any two objects and it is not requried to be from two different genotyping results.
#'@importFrom SummarizedExperiment assays colData rowData SummarizedExperiment
#'@importFrom S4Vectors DataFrame
#'@param se_1 SummarizedExperiment object 1.
#'@param se_2 SummarizedExperiment object 2.
#'@param suffixes The suffixes you want to add to the meta data.frame.
#'@param patient_check Should the patient column in the data be check, if they are the same?
#'@export
CombineSEobjects <- function(se_1, se_2, suffixes = c("_1", "_2"), patient_check = FALSE){
  # We check if the assays are equally named.
  assay_names_1 <- names(SummarizedExperiment::assays(se_1))
  assay_names_2 <- names(SummarizedExperiment::assays(se_2))
  if(!all(assay_names_1 == assay_names_2)){
    stop("Your assays are not equally named or ordered.")
  }
  patients_1 <- unique(se_1$Patient)
  if(length(patients_1) > 1) stop("The se_1 object has more than 1 patient.")
  patients_2 <- unique(se_2$Patient)
  if(length(patients_2) > 1) stop("The se_2 object has more than 1 patient.")
  if(patient_check){
    if(patients_1 != patients_2) stop("The objects are not from the same patient.")
  }

  features <- sigurd::combine_NAMES(rownames(se_1), rownames(se_2))
  cells    <- sigurd::combine_NAMES(colnames(se_1), colnames(se_2))

  meta_data_1 <- SummarizedExperiment::colData(se_1)
  meta_data_2 <- SummarizedExperiment::colData(se_2)
  meta_data   <- merge(meta_data_1, meta_data_2, by = "Cell", all = TRUE, suffixes = suffixes)
  meta_data   <- meta_data[match(cells, meta_data$Cell),]
  # We check if a patient value (the entry for a specific cells) is NA. If not, we use the patient value from the first column (se_1), otherwise we use the second column (se_2).
  patient_values <- meta_data[, paste0("Patient", suffixes)]
  patient_values <- ifelse(is.na(patient_values[,1]), patient_values[,2], patient_values[,1])
  # We check if a sample value (the entry for a specific cells) is NA. If not, we use the sample value from the first column (se_1), otherwise we use the second column (se_2).
  sample_values <- meta_data[, paste0("Sample", suffixes)]
  sample_values <- ifelse(is.na(sample_values[,1]), sample_values[,2], sample_values[,1])
  # We generate a new DataFrame with a new Cell column, a Patient column showing the patient value for all cells and a Sample column for the sample value.
  # The cells that were only in one SE object now have a Patient value. We also remove the original Pati
  columns_to_keep <- colnames(meta_data)
  columns_to_keep <- columns_to_keep[!columns_to_keep %in% c("Cell", paste0("Patient", suffixes[1]), paste0("Sample", suffixes[2]), paste0("Patient", suffixes[1]), paste0("Sample", suffixes[2]))]
  meta_data <- S4Vectors::DataFrame(Cell = meta_data$Cell, Patient = patient_values, Sample = sample_values, meta_data[,columns_to_keep])

  meta_row_1 <- SummarizedExperiment::rowData(se_1)
  meta_row_2 <- SummarizedExperiment::rowData(se_2)
  if(ncol(meta_row_1) > 0 & ncol(meta_row_2) > 0){
    meta_row <- merge(meta_row_1, meta_row_2, by = "VariantName", all = TRUE, suffixes = suffixes)
    meta_row <- meta_row[match(features, meta_row$VariantName),]
    rownames(meta_row) <- meta_row$VariantName
  } else if(ncol(meta_row_1) == 0 & ncol(meta_row_2) > 0){
    meta_row_1 <- matrix(NA, nrow = nrow(meta_row_1), ncol = ncol(meta_row_2))
    rownames(meta_row_1) <- rownames(se_1)
    colnames(meta_row_1) <- colnames(meta_row_2)
    meta_row_1 <- S4Vectors::DataFrame(meta_row_1)
    meta_row_1$VariantName <- rownames(meta_row_1)
    meta_row <- merge(meta_row_1, meta_row_2, by = "VariantName", all = TRUE, suffixes = suffixes)
    meta_row <- meta_row[match(features, meta_row$VariantName),]
  } else if(ncol(meta_row_1) > 0 & ncol(meta_row_2) > 0){
    meta_row_2 <- matrix(NA, nrow = nrow(meta_row_2), ncol = ncol(meta_row_1))
    rownames(meta_row_2) <- rownames(se_2)
    colnames(meta_row_2) <- colnames(meta_row_1)
    meta_row_2 <- S4Vectors::DataFrame(meta_row_2)
    meta_row_2$VariantName <- rownames(meta_row_2)
    meta_row <- merge(meta_row_1, meta_row_2, by = "VariantName", all = TRUE, suffixes = suffixes)
    meta_row <- meta_row[match(features, meta_row$VariantName),]    
  }

  assays_combined <- lapply(assay_names_1, function(x){
    result <- sigurd::combine_SparseMatrix(matrix_1 = SummarizedExperiment::assays(se_1)[[x]], matrix_2 = SummarizedExperiment::assays(se_2)[[x]])
  })
  names(assays_combined) <- assay_names_1

  se_combined <- SummarizedExperiment::SummarizedExperiment(assays = assays_combined, colData = meta_data, rowData = meta_row)
  return(se_combined)
}

