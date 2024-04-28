<h1 align="center">Seeker for highly variable homologous genes</h1>
🎯 Keeping the potential homologous genes set with biological meaning is of great importance before comparative analysis between species. To this end, we sought to take advantage of the concept of highly variable genes (HVGs), of which are widely used in single-cell RNA-seq analysis and may related to genuine biological variation. Furthermore, HVGs can be identified in an unsupervised and low calculation cost manner that are applicable to various kind of development system. Here we introduce HomologySeeker that is designed to identify homologous genes set with highly variable expression (Homology-HVGs) for cross-species analysis while keeping species-specific homologous/non-homologous genes for additional purpose

![image](https://github.com/Soap4/HomologySeeker/blob/master/image/HomologySeeker.png)

# Installation
To install *HomologySeeker*, please use:
```r
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("YenLab/HomologySeeker", upgrade = FALSE, dependencies = TRUE)

## or
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("YenLab/HomologySeeker", upgrade = FALSE, dependencies = TRUE)
```
# Usage
## 1. Data preperation ▶️

***HomologySeeker*** requires single-cell matrics (row as genes and column as samples) from two species to seek Homo-HVGs.

In this vignette, we provied pre-built single-cell RNA-seq dataset sampled from mouse and human midbrain by [Manno et al., *Cell*., 2016](https://doi.org/10.1016/j.cell.2016.09.027):

```r
## Midbrain_singlecell_Manno.RData contains single-cell matrix of mouse and human
## Download it with
download.file("https://github.com/Soap4/Data/raw/main/scRNA-seq/Midbrain_singlecell_Manno.RData","Midbrain_singlecell_Manno.RData")
load("Midbrain_singlecell_Manno.RData")

## Or directly load it with
load(url("https://github.com/Soap4/Data/raw/main/scRNA-seq/Midbrain_singlecell_Manno.RData"))
```

## 2. One-step process 🚀

To use ***HomologySeeker***, you can simply run:
```r
library(HomologySeeker)

midbrain <- HomoSeeker(RefSpec = "mouse",   ## Name for species 1
                       QuySpec = "human",   ## Name for species 2
                       RefSpec_mat = mouse,   ## Single-cell matrix for species 1
                       QuySpec_mat = human,   ## Single-cell matrix for species 2
                       HVGs_method = "seurat_vst")   ## HVG selection method                                   
```
The ```HomoSeeker()```function returns ```midbrain```, a ```HomoHVG object``` that contains 6 slots: (access different slot by using ```midbrain@slot name```)

+ ```Species```: Species offical name  
+ ```HomoHVG```: Species whole/overlapped Homo-HVGs list  
+ ```Matrix_orig```: Original single-cell matrix  
+ ```Matrix_homo```: Single-cell matrix with homologous genes as row  
+ ```HVG_feature```: HVG information  
+ ```Table_homo```: Homology information  

## 3. Step-by-step process 📜

To use ***HomologySeeker*** more flexibly, you can use standard ***HomologySeeker*** pipline as well:

### 3.1 Extract homologous gene between two species
First of all, ***HomologySeeker*** extract homologous information between input species by using ```HomoSelector()``` function:  
+ You can check available species name by using ```GetSpecNames()``` function
```r
library(HomologySeeker)

homo_mat <- HomoSelector(RefSpec = "mouse",   ## Name for species 1
                         QuySpec = "human",   ## Name for species 2
                         usedataset = TRUE)   ## Whether use pre-built dataset
```
The ```homo_mat``` contains homologous information between species 1 and 2 and used as input for further HVG selection.  
+ We implemented pre-built homologous table: [Mouse-Human](https://github.com/Soap4/Data/files/10283572/Orthologues_Human_Mouse.csv), [Mouse-Zebrafish](https://github.com/Soap4/Data/files/10283574/Orthologues_Mouse_Zebrafish.csv) and [Human-Zebrafish](https://github.com/Soap4/Data/files/10283573/Orthologues_Human_Zebrafish.csv) in ***HomologySeeker***. Those datasets are available through ```AvilData()``` function.

### 3.2 Identify HVGs

Next, ***HomologySeeker*** identify Homo-HVGs through ```HVGSelector()``` function: 
```r
midbrain <- HVGSelector(RefSpec = "mouse",   ## Name for species 1
                        QuySpec = "human",   ## Name for species 1
                        RefSpec_mat = mouse,   ## Single-cell matrix for species 1
                        QuySpec_mat = human,   ## Single-cell matrix for species 2
                        RefSpec_gene = "Gene_sym",   ## Type of gene name of single-cell matrix of species 1
                        QuySpec_gene = "Gene_sym",   ## Type of gene name of single-cell matrix of species 2
                        homo_mat = homo_mat,   ## Homologous table returned by HomoSelector()
                        HVGs_method = "seurat_vst",   ## HVG selection method
                        verbose = TRUE)
```
```HVGSelector()``` returns a ```HomoHVG object``` as describe above.

# Maintenance
🧐 Any possible questions or improvements for ***HomologySeeker*** are welcome to post on the [issue page](https://github.com/Soap4/Data/issues)

Shaokang Mo🤡  
soap79022@outlook.com

# Citation
Mo, S., Qu, K., Huang, J., Li, Q., Zhang, W., Yen, K. Cross-species transcriptomics reveals bifurcation point during the arterial-to-hemogenic transition. ***Commun Biol***. 6, 827 (2023). doi: 10.1038/s42003-023-05190-6
