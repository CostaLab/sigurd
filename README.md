# SIngle cell Genotyping Using Rna Data (SIGURD)

Martin Graßhoff<sup>1</sup>
Ivan G. Costa<sup>1</sup>

<sup>1</sup>Institute for Computational Genomics, Faculty of Medicine, RWTH Aachen University, Aachen, 52074 Germany

**Motivation:** With the advent of single RNA seq assays, it became possible to determine the mutational status for each individual cell. Single cell RNA seq data is by its very nature sparse and the probability of hitting a specific variants of interest is therefore very low. While this issue can be overcome using modified amplicon assays, it is also possible to impute the mutational status using the correlation between detected mitochondrial and somatic variants. 

**Results:** Sigurd is an R package for the analysis of single cell data. We determine the overall variant burden per cell and also the number of interesting mitochondrial variants using previously published approaches.\
We employ a imputation approach that utilizes the correlation between mitochondrial variants and somatic variants. Mitochondrial mutations that are significantly associated to somatic mutations are used as stand-ins.

## Installation

You can install sigurd using the following code. The vignette requires data that is currently not published, but is provided as a reference.

```{r}

install.packages("devtools")
devtools::install_github("https://github.com/CostaLab/sigurd.git", build_vignettes = FALSE)
require(sigurd)

```

# SIGURD 

We have provided a small example data set for SIGURD. It consists of chromosome 9 and MT for one MPN sample.

```{r}
# This will be included for published data.
# vignette('sigurd')

```

# Current Features v0.3.14

- Loading data from VarTrix and MAEGATK.
- Transforming the data to be compatible for joint analysis.
- Calculating the variant burden per cell.
- Thresholding variants using the approach described by Miller et al. [2]
- Finding associated variants using correlation or the Fisher Test.

# Sources

This package implements approaches from the following packages and respositories:
- https://github.com/petervangalen/MAESTER-2021 
-- Variant Thresholding and functions for the loading of MAEGATK data.
- https://github.com/CostaLab/CimpleG
-- The loading and saving function.

# Future 
- Memory optimization
- Loading of CB sniffer results
- Providing data for the vignette

# References
[1] VarTrix. [github](https://github.com/10XGenomics/vartrix)

[2] Miller, T.E., et al. Mitochondrial variant enrichment from high-throughput single-cell RNA sequencing resolves clonal populations. Nat Biotechnol (2022). [link](https://doi.org/10.1038/s41587-022-01210-8). See also: [MAEGATK Analysis](https://github.com/petervangalen/MAESTER-2021), [Data](https://vangalenlab.bwh.harvard.edu/resources/maester-2021/)
