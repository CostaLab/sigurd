#'CalculateStrandCorrelation
#'@description
#'We calculate the correlation between the amount of forward and reverse reads per variant.
#'@importFrom SummarizedExperiment rowRanges assays
#'@importFrom data.table data.table
#'@importFrom dplyr sym
#'@param SE SummarizedExperiment object.
#'@param chromosome_prefix List of matrices for the alternative reads.
#'@export
CalculateStrandCorrelation <- function(SE, chromosome_prefix = "chrM"){
  ref_allele <- as.character(SummarizedExperiment::rowRanges(SE)$refAllele)
  variants_A <- paste0(chromosome_prefix, "_", 1:length(ref_allele), "_", ref_allele, "_A")
  variants_A <- variants_A[ref_allele != "A"]
  reads_A_fw  <- SummarizedExperiment::assays(SE)[["A_counts_fw"]]
  reads_A_rev <- SummarizedExperiment::assays(SE)[["A_counts_rev"]]
  rownames(reads_A_fw)  <- paste0(chromosome_prefix, "_", 1:nrow(reads_A_fw), "_", ref_allele, "_A")
  rownames(reads_A_rev) <- paste0(chromosome_prefix, "_", 1:nrow(reads_A_rev), "_", ref_allele, "_A")
  reads_A_fw  <- reads_A_fw[ref_allele != "A",]
  reads_A_rev <- reads_A_rev[ref_allele != "A",]
  dt <- merge(data.table::data.table(summary(reads_A_fw)), 
              data.table::data.table(summary(reads_A_rev)), 
              by.x = c("i", "j"), by.y = c("i", "j"), 
              all = TRUE)
  dt <- dt[dt[,x.x] > 0 | dt[,x.y] > 0,]
  dt <- data.table::data.table(variant = variants_A[dt[[1]]],
                               cell_id = dt[[2]],
		                           fw = dt[[3]], rev = dt[[4]])
  cor_result_A <- dt[, .(cor = suppressWarnings(stats::cor(c(fw), c(rev), method = "pearson", use = "pairwise.complete"))), by = list(variant)]

  variants_C <- paste0(chromosome_prefix, "_", 1:length(ref_allele), "_", ref_allele, "_C")
  variants_C <-	variants_C[ref_allele != "C"]
  reads_C_fw  <- SummarizedExperiment::assays(SE)[["C_counts_fw"]]
  reads_C_rev <- SummarizedExperiment::assays(SE)[["C_counts_rev"]]
  rownames(reads_C_fw)  <- paste0(chromosome_prefix, "_", 1:nrow(reads_C_fw), "_", ref_allele, "_C")
  rownames(reads_C_rev) <- paste0(chromosome_prefix, "_", 1:nrow(reads_C_rev), "_", ref_allele, "_C")
  reads_C_fw  <- reads_C_fw[ref_allele != "C",]
  reads_C_rev <- reads_C_rev[ref_allele != "C",]
  dt <- merge(data.table::data.table(summary(reads_C_fw)),
              data.table::data.table(summary(reads_C_rev)),
              by.x = c("i", "j"), by.y = c("i", "j"),
              all = TRUE)
  dt <- dt[dt[,x.x] > 0 | dt[,x.y] > 0,]
  dt <- data.table::data.table(variant = variants_C[dt[[1]]],
                               cell_id = dt[[2]],
                               fw = dt[[3]], rev = dt[[4]])
  cor_result_C <- dt[, .(cor = suppressWarnings(stats::cor(c(fw), c(rev), method = "pearson", use = "pairwise.complete"))), by = list(variant)]

  variants_G <- paste0(chromosome_prefix, "_", 1:length(ref_allele), "_", ref_allele, "_G")
  variants_G <- variants_G[ref_allele != "G"]
  reads_G_fw  <- SummarizedExperiment::assays(SE)[["G_counts_fw"]]
  reads_G_rev <- SummarizedExperiment::assays(SE)[["G_counts_rev"]]
  rownames(reads_G_fw)  <- paste0(chromosome_prefix, "_", 1:nrow(reads_G_fw), "_", ref_allele, "_G")
  rownames(reads_G_rev) <- paste0(chromosome_prefix, "_", 1:nrow(reads_G_rev), "_", ref_allele, "_G")
  reads_G_fw  <- reads_G_fw[ref_allele != "G",]
  reads_G_rev <- reads_G_rev[ref_allele != "G",]
  dt <- merge(data.table::data.table(summary(reads_G_fw)),
              data.table::data.table(summary(reads_G_rev)),
              by.x = c("i", "j"), by.y = c("i", "j"),
              all = TRUE)
  dt <- dt[dt[,x.x] > 0 | dt[,x.y] > 0,]
  dt <- data.table::data.table(variant = variants_G[dt[[1]]],
                                         cell_id = dt[[2]],
                                         fw = dt[[3]], rev = dt[[4]])
  cor_result_G <- dt[, .(cor = suppressWarnings(stats::cor(c(fw), c(rev), method = "pearson", use = "pairwise.complete"))), by = list(variant)]

  variants_T <- paste0(chromosome_prefix, "_", 1:length(ref_allele), "_", ref_allele, "_T")
  variants_T <- variants_T[ref_allele != "T"]
  reads_T_fw  <- SummarizedExperiment::assays(SE)[["T_counts_fw"]]
  reads_T_rev <- SummarizedExperiment::assays(SE)[["T_counts_rev"]]
  rownames(reads_T_fw)  <- paste0(chromosome_prefix, "_", 1:nrow(reads_T_fw), "_", ref_allele, "_T")
  rownames(reads_T_rev) <- paste0(chromosome_prefix, "_", 1:nrow(reads_T_rev), "_", ref_allele, "_T")
  reads_T_fw  <- reads_T_fw[ref_allele != "T",]
  reads_T_rev <- reads_T_rev[ref_allele != "T",]
  dt <- merge(data.table::data.table(summary(reads_T_fw)),
              data.table::data.table(summary(reads_T_rev)),
              by.x = c("i", "j"), by.y = c("i", "j"),
              all = TRUE)
  dt <- dt[dt[,x.x] > 0 | dt[,x.y] > 0,]
  dt <- data.table::data.table(variant = variants_T[dt[[1]]],
                               cell_id = dt[[2]],
                               fw = dt[[3]], rev = dt[[4]])
  cor_result_T <- dt[, .(cor = suppressWarnings(stats::cor(c(fw), c(rev), method = "pearson", use = "pairwise.complete"))), by = list(variant)]

  cor_results <- rbind(cor_result_A, cor_result_C, cor_result_G, cor_result_T)
  
  variants <- c(variants_A, variants_C, variants_G, variants_T)
  result <- rep(NA, length(variants))
  names(result) <- variants
  result[cor_results$variant] <- cor_results$cor
  return(result)
}
