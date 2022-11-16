#' Title
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
#' @param species1
#' @param species2
#' @param species1_mat
#' @param species2_mat
#' @param species1_gene
#' @param species2_gene
#' @param method
#' @param HVGs_method
#' @param verbose
#' @param mirror
#' @param version
#' @param host
#'
#' @return
#' @export
#'
#' @examples
HomoSeeker <- function(species1,
                        species2,
                        species1_mat,
                        species2_mat,
                        species1_gene = "Gene_sym",# Gene_sym/Gene_id
                        species2_gene = "Gene_sym",# Gene_sym/Gene_id
                        method = "seurat",
                        HVGs_method = "vst",
                        verbose = T,
                        mirror = NULL,
                        version = NULL,
                        host = "https://dec2021.archive.ensembl.org/"){
  homo_mat <- HomoSelector(species1 = species1,
                           species2 = species2,
                           homotype = "ortholog_one2one",
                           mirror = mirror,
                           version = version,
                           host = host)
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
