#' @title  Create a HomoHVG object
#'
#' @description HomoHVG object containing species matrics, homologous genes and variable genes information
#'
#' @importFrom methods setClass
#' @importFrom magrittr %>% set_names
#'
#' @param species1,species2 Names of species in comparative analysis. Case is ignored
#' @param species1_mat,species2_mat Single cell expression matrix of species 1 or 2 with row as gene symbol/ID and column as sample ID
#' @param species1_mat_homo,species2_mat_homo Single cell expression matrix of species 1 or 2 with row as homologous gene name/ID and column as sample ID
#' @param spec1_HVG,spec2_HVG Homo-HVG list of species 1 or 2
#' @param Table_homo Table of one to one homologous genes list between species 1 and species 2
#'
#' @slot Species List of species name
#' @slot HomoHVG List of Homo-HVGs
#' @slot Matrix_orig List of original single cell matrix
#' @slot Matrix_homo List of single cell matrix with homologous genes as row
#' @slot HVG_feature List of HVG information table
#' @slot Table_homo Table of homologous information
#'
#' @name HomoHVG-class
#' @concept objects
#' @exportClass HomoHVG
#'
#' @return HomoHVG object
#' @export
#'

setClass(Class = "HomoHVG",
         slots =list(Species = "list",
                     HomoHVG = "list",
                     Matrix_orig = "list",
                     Matrix_homo = "list",
                     HVG_feature = "list",
                     Table_homo = "data.frame"))

setMethod(f = "show",
          signature = "HomoHVG",
          definition = function(object){
            cat("A HomoHVG object containing homologous information between",
                object@Species[[1]],"and",object@Species[[2]],"\n")
            cat("Species:",object@Species[[1]],object@Species[[2]],"\n")
            cat(nrow(object@Matrix_homo[[1]]),"one-to-one homologous genes pairs found in total","\n")
            cat(length(object@HomoHVG[[1]]),"Homo-HVGs in",object@Species[[1]],"\n")
            cat(length(object@HomoHVG[[2]]),"Homo-HVGs in",object@Species[[2]],"\n")
          })

HomoHVGObject <- function(species1,
                          species2,
                          species1_mat,
                          species2_mat,
                          species1_mat_homo,
                          species2_mat_homo,
                          spec1_HVG,
                          spec2_HVG,
                          Table_homo){
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
    Table_homo = Table_homo)
  return(object)
}
