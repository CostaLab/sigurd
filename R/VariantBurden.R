VariantBurden <- function(se){
  suppressPackageStartupMessages(library(Matrix))
  suppressPackageStartupMessages(library(SummarizedExperiment))

  MT_check <- grepl("^chrM_|^chrMT_|^MT_|^chrMt_|^chrM_|^Mt_", rownames(se))
  se_subset <- se[MT_check,]
  burden_MT <- colSums(assays(se_subset)$fraction)
  se_subset <- se[!MT_check,]
  burden_somatic <- colSums(assays(se_subset)$fraction)

  colData(se)[,"Burden_Somatic"] <- burden_somatic
  colData(se)[,"Burden_MT"] <- burden_MT
  return(se)
}
