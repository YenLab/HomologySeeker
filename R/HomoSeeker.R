#' @title Seeker for highly varible homologous gene sets between species
#'
#' @description Get homologous information and variable gene sets between two species
#'
#' @include GetSpecNames.R
#' @include HomoHVG_Object.R
#' @include HomoSelector.R
#' @include HVGSelector.R
#'
#' @importFrom methods setClass
#' @importFrom magrittr %>% set_names
#' @importFrom dplyr select mutate
#' @importFrom biomaRt useEnsembl getLDS
#' @import Seurat SeuratObject
#'
#'
#' @param species1,species2 Names of species in comparative analysis. Case is ignored
#' @param species1_mat,species2_mat Single cell expression matrix of species 1 or 2 with row as gene symbol/ID and column as sample ID
#' @param species1_gene,species2_gene Type of gene name of Single cell expression matrix.
#' \itemize{
#'  \item{"Gene_sym"} \strong{:} Gene symbel (default).
#'  \item{"Gene_id"} \strong{:} Gene ID.
#'  }
#' @param method HVG selection method to be used. Currently only Seurat is supported.
#' @param HVGs_method if method = "Seurat", available methods: vst(Default), sct. See Seurat::FindVariableFeatures() and Seurat::SCTransform() for further details
#' @param version Ensembl version to be connected. See HomoSelector() and biomaRt::useEnsembl() for detailed information.
#' @param verbose Whether show calculation progress. Default is TRUE.
#' @param usedataset Whether used pre-built homologou information datasets. See AvilData() for details
#'
#' @return Returns a HomoHVG object with slot:
#' #' \itemize{
#'  \item{Species} \strong{:} Query species names.
#'  \item{HomoHVG} \strong{:} Highly variable gene sets for species 1 and 2.
#'  \item{Matrix_orig} \strong{:} Original single cell matrices for species 1 and 2.
#'  \item{Matrix_homo} \strong{:} Single cell matrices with one-to-one homologous genes between species 1 and 2 as row.
#'  \item{HVG_feature} \strong{:} Returned HVG information of all homologous genes for species 1 and 2.
#'  \item{Table_homo} \strong{:} Returned homologous genes table between species 1 and 2.
#'  }
#' @export
#'
HomoSeeker <- function(species1,
                       species2,
                       species1_mat,
                       species2_mat,
                       species1_gene = "Gene_sym",
                       species2_gene = "Gene_sym",
                       method = "seurat",
                       HVGs_method = "vst",
                       version = NULL,
                       usedataset = TRUE,
                       verbose = TRUE){
  if(usedataset){
    if(all(GetSpecNames(c(species1,species2))[,2] %in% c("Mouse","Human"))){

      message(paste0("Loading existing Human2Mouse datasets"))
      homo_mat <- read.csv("https://github.com/Soap4/Data/files/10098277/Orthologues_Human2Mouse_v108.csv",row.names = 1)
      homo_mat <- homo_mat[homo_mat %>% select(contains('type')) == "ortholog_one2one",]
      homo_mat <- homo_mat[homo_mat %>% select(contains('confi'))==1,]

    }else if(all(GetSpecNames(c(species1,species2))[,2] %in% c("Zebrafish","Human"))){

      message(paste0("Loading existing Human2Zebrafish datasets"))
      homo_mat <- read.csv("https://github.com/Soap4/Data/files/10098276/Orthologues_Human2Zebrafish_v108.csv",row.names = 1)
      homo_mat <- homo_mat[homo_mat %>% select(contains('type')) == "ortholog_one2one",]
      homo_mat <- homo_mat[homo_mat %>% select(contains('confi'))==1,]

    }else if(all(GetSpecNames(c(species1,species2))[,2] %in% c("Zebrafish","Mouse"))){

      message(paste0("Loading existing Mouse2Zebrafish datasets"))
      Mouse2Zebrafish <- read.csv("https://github.com/Soap4/Data/files/10098275/Orthologues_Mouse2Zebrafish_v108.csv",row.names = 1)
      homo_mat <- homo_mat[homo_mat %>% select(contains('type')) == "ortholog_one2one",]
      homo_mat <- homo_mat[homo_mat %>% select(contains('confi'))==1,]

    }else{
      if(is.null(version)){
        version = 105
      }else{
        next
        }
      homo_mat <- HomoSelector(species1 = species1,
                               species2 = species2,
                               homotype = "ortholog_one2one",
                               version = version)
  }}
  object <- HVGSelector(species1 = species1,
                        species2 = species2,
                        species1_mat = species1_mat,
                        species2_mat = species2_mat,
                        homo_mat = homo_mat,
                        species1_gene = species1_gene,
                        species2_gene = species2_gene,
                        method = method,
                        HVGs_method = HVGs_method,
                        verbose = verbose)
  return(object)
}
