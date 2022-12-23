# Seeker for homologous genes with highly variable expression
🎯 Keeping the potential homologous genes set with biological meaning is of great importance before comparative analysis between species. To this end, we sought to take advantage of the concept of highly variable genes (HVGs), of which are widely used in single-cell RNA-seq analysis and may related to genuine biological variation. Furthermore, HVGs can be identified in an unsupervised and low calculation cost manner that are applicable to various kind of development system. Here we introduce HomoSeeker that is designed to identify homologous genes set with highly variable expression (Homo-HVGs) for cross-species analysis while keeping species-specific homologous/non-homologous genes for additional purpose

# Installation
To install *HomoSeeker*, please use:
```r
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("Soap4/HomoSeeker")

# or
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("Soap4/HomoSeeker")
```
# Usage
## 1. Data preperation
