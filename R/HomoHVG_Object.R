#' @title  Create a HomoHVG object
#'
#' @description HomoHVG object containing species matrics, homologous genes and variable genes information
#'
#' @importFrom methods setClass
#' @importFrom magrittr %>% set_names
#'
#' @param RefSpec,QuySpec Names of reference(ref) and query(quy) species in comparative analysis. Case is ignored
#' @param RefSpec_mat,QuySpec_mat Single cell expression matrix of reference or query species with row as gene symbol/ID and column as sample ID
#' @param RefSpec_mat_homo,QuySpec_mat_homo Single cell expression matrix of reference or query species with row as homologous gene name/ID and column as sample ID
#' @param RefSpec_HVGs,QuySpec_HVGs Homo-HVG list of reference or query species
#' @param Table_homo Table of one to one homologous genes list between reference and query species
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
                     Table_homo = "data.frame",
                     HVGs_method = "character"))

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

HomoHVGObject <- function(RefSpec,
                          QuySpec,
                          RefSpec_mat,
                          QuySpec_mat,
                          RefSpec_mat_homo,
                          QuySpec_mat_homo,
                          RefSpec_HVG,
                          QuySpec_HVG,
                          Table_homo,
                          HVGs_method){
  object <- new(
    Class = 'HomoHVG',
    Species = list(RefSpec, QuySpec) %>%
      set_names(c("RefSpec","QuySpec")),
    HomoHVG = list(rownames(RefSpec_HVG)[RefSpec_HVG$is.HomoHVGs=="TRUE"],
                   rownames(QuySpec_HVG)[QuySpec_HVG$is.HomoHVGs=="TRUE"]) %>%
      set_names(paste0(c(RefSpec,QuySpec),"_HomoHVGs")),
    Matrix_orig = list(RefSpec_mat, QuySpec_mat) %>%
      set_names(paste0(c(RefSpec, QuySpec),"_Matrix")),
    Matrix_homo = list(RefSpec_mat_homo, QuySpec_mat_homo) %>%
      set_names(paste0(c(RefSpec, QuySpec),"_HomoMatrix")),
    HVG_feature = list(RefSpec_HVG, QuySpec_HVG) %>%
      set_names(paste0(c(RefSpec, QuySpec),"_HVGFeatures")),
    Table_homo = Table_homo,
    HVGs_method = HVGs_method)
  return(object)
}
