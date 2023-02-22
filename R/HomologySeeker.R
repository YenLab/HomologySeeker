#' @title Seeker for highly varible homologous gene sets between species
#'
#' @description Get homologous information and variable gene sets between two species
#'
#' @include GetSpecNames.R
#' @include HomologyHVG_Object.R
#' @include HomologySelector.R
#' @include HVGSelector.R
#'
#' @importFrom methods setClass
#' @importFrom magrittr %>% set_names
#' @importFrom dplyr select mutate
#' @importFrom biomaRt useEnsembl getLDS
#' @import Seurat SeuratObject
#'
#'
#' @param RefSpec,QuySpec Names of reference(ref) and query(quy) species in comparative analysis. Case is ignored.
#' See GetSpecNames() for detailed list.
#' @param RefSpec_mat,QuySpec_mat Single cell expression matrix of reference or query species with row as gene
#' symbol/ID and column as sample ID.
#' @param homo_mat Table of one to one homologous genes table between reference and query species. Default to use homologous
#' gene table collected by HomoSelector(). To use self constructed homologous gene table, please input table with:
#' #' \itemize{
#'  \item{1st Column} \strong{:} Gene symbel/ID of reference species.
#'  \item{2nd Column} \strong{:} Corresponding Gene symbel/ID of query species.
#'  }
#' @param RefSpec_gene,QuySpec_gene Type of gene name of Single cell expression matrix.
#' \itemize{
#'  \item{Gene_sym} \strong{:} Gene symbel (default).
#'  \item{Gene_id} \strong{:} Gene ID.
#'  }
#' @param HVGs_method HVG selection method to be used. Supported methods:
#' \itemize{
#'  \item{seurat_vst(default)}
#'  \item{seurat_sct}
#'  \item{scran}
#'  \item{scmap}
#'  \item{ROGUE}
#'  }
#'  See HVGSelector() for detailed description.
#' @param Integrated_mat Whether input matrics are returned by integration
#' @param Integrated_method If Integrated_mat=TRUE, please specify what integration method was used.
#' currently, only seurat_LogNorm and seurat_SCT intergation are supported
#' @param verbose Whether show calculation progress. Default is TRUE.
#'
#' @return Returns a HomologyHVG object with slot:
#' \itemize{
#'  \item{Species} \strong{:} Query species names.
#'  \item{HomologyHVG} \strong{:} Highly variable gene sets for species 1 and 2.
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
                       homo_mat = NULL,
                       RefSpec_gene = "Gene_sym",
                       QuySpec_gene = "Gene_sym",
                       HVGs_method = "seurat_vst",
                       Integrated_mat = FALSE,
                       Integrated_method = NULL,
                       verbose = TRUE){
  if(is.null(homo_mat)){
    homo_mat <- HomoSelector(RefSpec = RefSpec,
                             QuySpec = QuySpec,
                             homotype = "ortholog_one2one")
  }
    object <- HVGSelector(RefSpec = RefSpec,
                          QuySpec = QuySpec,
                          RefSpec_mat = RefSpec_mat,
                          QuySpec_mat = QuySpec_mat,
                          homo_mat = homo_mat,
                          RefSpec_gene = RefSpec_gene,
                          QuySpec_gene = QuySpec_gene,
                          HVGs_method = HVGs_method,
                          verbose = verbose,
                          Integrated_mat = Integrated_mat,
                          Integrated_method = Integrated_method)
  return(object)
}
