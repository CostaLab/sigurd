# Analysis of Blastic Plasmacytoid Dendritic Cell Neoplasm (BPDCN)
# Here, we compare the analysis of cells from a patient with BPDCN.
# The original script is available here: https://github.com/petervangalen/MAESTER-2021/blob/main/4_CH_sample/4.2_Variant_Selection.R
# In this setup, the cells are all from a single donor and do not have to be separated into categories. 
# Then small clones that might characterize the clonal lineages in the sample are selected.

# Loading necessary packages.

print("Libraries for SIGURD.")
suppressPackageStartupMessages(library(sigurd))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(ggplot2))

# Loading the data using SIGURD.
# The design matrix contains information for the genotyping data.
genotyping <- LoadingMAEGATK_typewise(patient = "BPDCN712", samples_file = "MAESTER_Reproduction.csv", type_use = "Amplicon_MT", verbose = FALSE)
# Loading the scRNA-seq data. The sample ID has been prepended to the cell names.
scrna <- readRDS("BPDCN712_Seurat_with_TCR_Renamed.rds")


# Generating a block list
# Variants that are also detected in the cell mixture data is treated as possible false positives. They are used as a blacklist and are removed from the results.
genotyping_tenx <- readRDS("TenX_CellMixture_Genotyping.rds")
blocklist <- AllelFrequencyFoldChange(genotyping_tenx, group_of_interest = "CellType", group1 = "K562", group2 = "BT142", maximum_foldchange = 5, minimum_coverage = 5, minimum_allele_freq = 0.001, maximum_allele_freq = 0.999)$Variant
# We add the variant chrM_1583_A_G by hand. It is identified as misleading based on downstream analysis.
blocklist <- c(blocklist, "chrM_1583_A_G")


# Selecting the Variants of Interest
# Now, we select the variants of interest from the loaded data.
voi.ch.sigurd <- VariantSelection_TopCells(genotyping, min_coverage = 5, quantiles = 0.9, thresholds = 0, top_cells = 10, top_VAF = 0.5, min_quality = 30, remove_nocall = FALSE, verbose = FALSE)
voi.ch.sigurd <- voi.ch.sigurd[!voi.ch.sigurd %in% blocklist]
print(voi.ch.sigurd)


# Visualisation using SIGURD
genotyping <- Filtering(SE = genotyping, cells_include = colnames(scrna))
colData(genotyping)$CellType <- scrna$CellType
HeatmapVoi(SE = genotyping, voi = voi.ch.sigurd, annotation_trait = "CellType", sort_cells = TRUE, remove_empty_cells = TRUE, minimum_allele_freq = 0.01)


# Cell Type Enrichment 2593G>A
# We check the cell type enrichment for the variant 2593G>A.
result_enrichment <- Enrichment_FisherTest(se = genotyping, variant = "chrM_2593_G_A", use_nocall = FALSE, trait = "CellType")
ComplexHeatmap::Heatmap(as.matrix(result_enrichment$Heatmap_Values), cluster_rows = FALSE, cluster_columns = FALSE, col= circlize::colorRamp2(c(0, max(result_enrichment$Heatmap_Values, na.rm = TRUE)), c("white", "red")))


# Cell Type Enrichment 6243G>A
result_enrichment <- Enrichment_FisherTest(se = genotyping, variant = "chrM_6243_G_A", use_nocall = FALSE, trait = "CellType")
ComplexHeatmap::Heatmap(as.matrix(result_enrichment$Heatmap_Values), cluster_rows = FALSE, cluster_columns = FALSE, col= circlize::colorRamp2(c(0, max(result_enrichment$Heatmap_Values, na.rm = TRUE)), c("white", "red")))


# Enrichment with somatic variants
# We get the number of 2593_G>A cells that are also mutated for a different somatic variant.
somatic_variants <- c("ASXL1.G642fs.1", "ASXL1.G642fs.2", "TET2.S792X", "TET2.Q1034X", "TET2.R1216X", "TET2.H1380Y")
# Loading the GoT results.
variants <- LoadingRawMatrix_typewise(samples_file = "MAESTER_Reproduction.csv", patient = "BPDCN712", variant_use = somatic_variants, matrix_column_separator = "\t", verbose = FALSE)
# Selecting only our variants of interest.
genotyping <- genotyping[voi.ch.sigurd, ]
# Getting the intersection between the MT genotyped cells and the somatic results.
genotyping <- Filtering(genotyping, cells_include = colnames(variants))
variants <- Filtering(variants, cells_include = colnames(scrna), fraction_threshold = 0.01, reject_value = "Reference")
# We combine both genotyping objects. The object now contains the information for all variants and cells.
genotyping <- CombineSEobjects(se_1 = genotyping, se_2 = variants)
# We define a clone a carrying the 2593_G>A variant, but not the 6243_G>A variant.
genotyping <- ClonalDefinition(se = genotyping, variants_ls = list("chrM_2593_G_A"), grouping = NULL, identities = NULL, explicit = TRUE, explicit_not = list("chrM_6243_G_A"), explicit_min_vaf = 0.01, verbose = TRUE)
# Subsetting to only include the clonal cells.
genotyping <- genotyping[, genotyping$Clones == "C1"]
# Using the Fisher test to check for a significant enrichment. 
res <- VariantWiseFisherTest(variants_list = RowWiseSplit(se = genotyping))
print(res)
