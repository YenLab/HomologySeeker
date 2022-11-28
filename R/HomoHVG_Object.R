#' @title  Create a HomoHVG object
#'
#' @description HomoHVG object containing species matrics, homologous genes and variable genes information
#'
#' @importFrom methods setClass
#' @importFrom magrittr %>% set_names
#' @importFrom dplyr select mutate
#' @importFrom biomaRt useEnsembl getLDS
#' @import Seurat SeuratObject
#'
#' @param species1,species2 Names of species in comparative analysis. Case is ignored
#' @param species1_mat,species2_mat Single cell expression matrix of species 1 or 2 with row as gene symbol/ID and column as sample ID
#' @param species1_mat_homo,species2_mat_homo Single cell expression matrix of species 1 or 2 with row as homologous gene name/ID and column as sample ID
#' @param spec1_HVG,spec2_HVG Homo-HVG list of species 1 or 2
#' @param Table_homo Table of one to one homologous genes list between species 1 and species 2
#'
#' @return HomoHVG object
#' @export
#'
HomoHVGObject <- function(species1,
                          species2,
                          species1_mat,
                          species2_mat,
                          species1_mat_homo,
                          species2_mat_homo,
                          spec1_HVG,
                          spec2_HVG,
                          Table_homo){
  setClass(Class = "HomoHVG",representation(Species = "list",
                                            HomoHVG = "list",
                                       Matrix_orig = "list",
                                       Matrix_homo = "list",
                                       HVG_feature = "list",
                                       Table_homo = "data.frame"))
  object <- new(
    Class = 'HomoHVG',
    Species = list(species1, species2) %>%
      set_names(c("Species1","Species2")),
    HomoHVG = list(rownames(spec1_HVG)[spec1_HVG$is.HomoHVGs=="TRUE"],
                   rownames(spec2_HVG)[spec2_HVG$is.HomoHVGs=="TRUE"]) %>%
      set_names(paste0(c(species1,species2),"_HomoHVGs")),
    Matrix_orig = list(species1_mat, species2_mat) %>%
      set_names(paste0(c(species1, species2),"_Matrix")),
    Matrix_homo = list(species1_mat_homo, species2_mat_homo) %>%
      set_names(paste0(c(species1, species2),"_HomoMatrix")),
    HVG_feature = list(spec1_HVG, spec2_HVG) %>%
      set_names(paste0(c(species1,species2),"_HVGFeatures")),
    Table_homo = Table_homo
    )
  return(object)
}
