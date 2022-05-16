# SIngle cell Genotyping Using Rna Data (SIGURD)

Martin Graßhoff<sup>1</sup>
Ivan G. Costa<sup>1</sup>

<sup>1</sup>Institute for Computational Genomics, Faculty of Medicine, RWTH Aachen University, Aachen, 52074 Germany

**Motivation:** With the advent of single RNA seq assays, it became possible to determine the mutational status for each individual cell. Single cell RNA seq data is by its very nature sparse and the probability of hitting a specific variants of interest is therefore very low. While this issue can be overcome using modified amplicon assays, it is also possible to impute the mutational status using the correlation between detected mitochondrial and somatic variants. 

**Results:** Sigurd is an R package for the analysis of single cell data. We determine the overall variant burden per cell and also the number of interesting mitochondrial variants using previously published approaches. 
We employ a imputation approach that utilizes the correlation between mitochondrial variants and somatic variants. Mitochondrial mutations that are significantly correlated to somatic mutations are used as stand-ins.

## Installation

You can install 
