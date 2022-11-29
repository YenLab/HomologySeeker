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
#' @importFrom scuttle logNormCounts
#' @importFrom scran modelGeneVar
#' @import Seurat SeuratObject
#'
#'
#' @param RefSpec,QuySpec Names of reference(ref) and query(quy) species in comparative analysis. Case is ignored
#' @param RefSpec_mat,QuySpec_mat Single cell expression matrix of reference or query species with row as gene symbol/ID and column as sample ID
#' @param RefSpec_gene,QuySpec_gene Type of gene name of Single cell expression matrix.
#' \itemize{
#'  \item{Gene_sym} \strong{:} Gene symbel (default).
#'  \item{Gene_id} \strong{:} Gene ID.
#'  }
#' @param method HVG selection method to be used. Currently Supported methods:
#' \itemize{
#'  \item{seurat}
#'  \item{scran}
#'  }
#' @param HVGs_method if method = "Seurat", available methods:
#' \itemize{
#'  \item{vst (Default)} \strong{:} See Seurat::FindVariableFeatures().
#'  \item{sct} \strong{:} See Seurat::SCTransform().
#'  }
#' @param version Ensembl version to be connected. See HomoSelector() and biomaRt::useEnsembl() for detailed information.
#' @param verbose Whether show calculation progress. Default is TRUE.
#'
#' @return Returns a HomoHVG object with slot:
#' \itemize{
#'  \item{Species} \strong{:} Query species names.
#'  \item{HomoHVG} \strong{:} Highly variable gene sets for species 1 and 2.
#'  \item{Matrix_orig} \strong{:} Original single cell matrices for species 1 and 2.
#'  \item{Matrix_homo} \strong{:} Single cell matrices with one-to-one homologous genes between species 1 and 2 as row.
#'  \item{HVG_feature} \strong{:} Returned HVG information of all homologous genes for species 1 and 2.
#'  \item{Table_homo} \strong{:} Returned homologous genes table between species 1 and 2.
#'  }
#' @export
#'
HomoSeeker <- function(RefSpec,
                       QuySpec,
                       RefSpec_mat,
                       QuySpec_mat,
                       RefSpec_gene = "Gene_sym",
                       QuySpec_gene = "Gene_sym",
                       method = "seurat",
                       HVGs_method = "vst",
                       version = NULL,
                       verbose = TRUE){

      homo_mat <- HomoSelector(RefSpec = RefSpec,
                               QuySpec = QuySpec,
                               homotype = "ortholog_one2one",
                               version = version)
      object <- HVGSelector(RefSpec = RefSpec,
                            QuySpec = QuySpec,
                            RefSpec_mat = RefSpec_mat,
                            QuySpec_mat = QuySpec_mat,
                            homo_mat = homo_mat,
                            RefSpec_gene = RefSpec_gene,
                            QuySpec_gene = QuySpec_gene,
                            method = method,
                            HVGs_method = HVGs_method,
                            verbose = verbose)
  return(object)
}
