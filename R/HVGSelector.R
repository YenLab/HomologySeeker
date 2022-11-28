#' @title Highly variable homologous genes identifier
#'
#' @description Identify highly variable genes using homologous gene set with HVG Methodologies
#'
#' @include HomoHVG_Object.R
#' @include HomoSelector.R
#'
#' @importFrom methods setClass
#' @importFrom magrittr %>% set_names
#' @importFrom dplyr select mutate
#' @importFrom biomaRt useEnsembl getLDS
#' @importFrom scuttle logNormCounts
#' @importFrom scran modelGeneVar
#'
#' @import Seurat SeuratObject
#'
#' @param species1,species2 Names of species in comparative analysis. Case is ignored.
#' @param species1_mat,species2_mat Single cell expression matrix of species 1 or 2 with row as gene symbol/ID and column as sample ID.
#' @param homo_mat Table of one to one homologous genes table between species 1 and species 2. See HomoSelector() for further details.
#' @param species1_gene,species2_gene Type of gene name of Single cell expression matrix.
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
#'  \item{vst(Default)} \strong{:} See Seurat::FindVariableFeatures().
#'  \item{scran} \strong{:} See Seurat::SCTransform().
#'  }
#' @param verbose Whether show calculation progress. Default is TRUE.
#'
#' @return A HomoHVG object with Highly variable homologous gene sets for both species
#' @export
#'
HVGSelector <- function(species1,
                        species2,
                        species1_mat,
                        species2_mat,
                        homo_mat = NULL,
                        species1_gene = "Gene_sym",
                        species2_gene = "Gene_sym",
                        method = "seurat",
                        HVGs_method = "vst",
                        verbose = TRUE){
  if(is.null(homo_mat)){
    homo_mat <- HomoSelector(species1 = species1,species2 = species2)
    mat <- homo_mat[,c("Gene.name","Gene.stable.ID",
                       "Gene.name.1","Gene.stable.ID.1")]
    }else if(all(c("Gene.name","Gene.stable.ID",
                   "Gene.name.1","Gene.stable.ID.1") %in% colnames(homo_mat))){
      mat <- homo_mat[,c("Gene.name","Gene.stable.ID",
                       "Gene.name.1","Gene.stable.ID.1")]
    }else{
      mat <- homo_mat[,1:4] %>% set_colnames(c("Gene.name","Gene.stable.ID",
                                               "Gene.name.1","Gene.stable.ID.1"))
  }
  species1_index <- ifelse(species1_gene == "Gene_sym","Gene.name","Gene.stable.ID")
  species2_index <- ifelse(species2_gene == "Gene_sym","Gene.name.1","Gene.stable.ID.1")
  mat_fil <- mat[mat[,species1_index] %in% rownames(species1_mat) &
                   mat[,species2_index] %in% rownames(species2_mat),]
  species1_mat_homo <- species1_mat[mat_fil[,species1_index],]
  species2_mat_homo <- species2_mat[mat_fil[,species2_index],]

  if(method == "seurat" & HVGs_method == "vst"){
    if(verbose){message(paste0("Identifying HVGs using ",HVGs_method," from ",method))}
    species1_HVGs <- CreateSeuratObject(species1_mat_homo) %>%
      suppressWarnings() %>%
      FindVariableFeatures(verbose = F) %>%
      {.[["RNA"]]@meta.features} %>%
      dplyr::mutate(HomoGene=rownames(species2_mat_homo)) %>%
      {.[order(.$vst.variance.standardized,decreasing = T),]}

    species2_HVGs <- CreateSeuratObject(species2_mat_homo) %>%
      suppressWarnings() %>%
      FindVariableFeatures(verbose = F) %>%
      {.[["RNA"]]@meta.features} %>%
      dplyr::mutate(HomoGene=rownames(species1_mat_homo)) %>%
      {.[order(.$vst.variance.standardized,decreasing = T),]}
    if(verbose){message("HVGs Screening")}
    species1_HVGs$is.HomoHVGs <- ifelse(species1_HVGs$vst.variance.standardized>=
                                          mean(species1_HVGs$vst.variance.standardized),
                                        "TRUE","FALSE")
    species2_HVGs$is.HomoHVGs <- ifelse(species2_HVGs$vst.variance.standardized>=
                                          mean(species2_HVGs$vst.variance.standardized),
                                        "TRUE","FALSE")
  }
  if(method == "seurat" & HVGs_method == "sct"){
    if(verbose){message(paste0("Identifying HVGs using ",HVGs_method," from ",method))}
    species1_HVGs <- CreateSeuratObject(species1_mat_homo) %>%
      suppressWarnings() %>%
      SCTransform(return.only.var.genes = F,verbose = F) %>%
      {.[['SCT']]@SCTModel.list$model1@feature.attributes} %>%
      mutate(HomoGene=rownames(species2_mat_homo)) %>%
      {.[order(.$residual_variance,decreasing = T),]}

    species2_HVGs <- CreateSeuratObject(species2_mat_homo) %>%
      suppressWarnings() %>%
      FindVariableFeatures(verbose = F) %>%
      {.[['SCT']]@SCTModel.list$model1@feature.attributes} %>%
      mutate(HomoGene=rownames(species1_mat_homo)) %>%
      {.[order(.$residual_variance,decreasing = T),]}
    if(verbose){message("Screening HVGs")}
      species1_HVGs$is.HomoHVGs <- ifelse(species1_HVGs$residual_variance>=
                                            mean(species1_HVGs$residual_variance),
                                          "TRUE","FALSE")
      species2_HVGs$is.HomoHVGs <- ifelse(species2_HVGs$residual_variance>=
                                            mean(species2_HVGs$residual_variance),
                                          "TRUE","FALSE")
      }
  if(method == "scran"){
    if(verbose){message(paste0("Identifying HVGs using ",method))}
    species1_HVGs <- SingleCellExperiment::SingleCellExperiment(list(counts = as.matrix(species1_mat_homo))) %>%
      scuttle::logNormCounts() %>%
      scran::modelGeneVar() %>%
      as.data.frame() %>%
      dplyr::mutate(HomoGene=rownames(species2_mat_homo)) %>%
      {.[order(.$bio,decreasing = T),]}

    species2_HVGs <- SingleCellExperiment::SingleCellExperiment(list(counts = as.matrix(species2_mat_homo))) %>%
      scuttle::logNormCounts() %>%
      scran::modelGeneVar() %>%
      as.data.frame() %>%
      dplyr::mutate(HomoGene=rownames(species1_mat_homo)) %>%
      {.[order(.$bio,decreasing = T),]}

    if(verbose){message("Screening HVGs")}
    species1_HVGs$is.HomoHVGs <- ifelse(species1_HVGs$bio>=mean(species1_HVGs$bio),
                                        "TRUE","FALSE")
    species2_HVGs$is.HomoHVGs <- ifelse(species2_HVGs$bio>=mean(species2_HVGs$bio),
                                        "TRUE","FALSE")
  }
  if(verbose){message("Complete Successfully")}
    message("There are totally: ","\n",
            grep("TRUE",species1_HVGs$is.HomoHVGs) %>% length()," Homo-HVGs in Species1","\n",
            grep("TRUE",species2_HVGs$is.HomoHVGs) %>% length()," Homo-HVGs in Species2","\n",
            intersect(species1_HVGs$HomoGene[species1_HVGs$is.HomoHVGs=="TRUE"],
                      rownames(species2_HVGs)[species2_HVGs$is.HomoHVGs=="TRUE"]) %>% length(),
            " Homo-HVGs overlapped")

  Object <- HomoHVGObject(species1 = species1,
                          species2 = species2,
                          species1_mat = species1_mat,
                          species2_mat = species2_mat,
                          species1_mat_homo = species1_mat_homo,
                          species2_mat_homo = species2_mat_homo,
                          spec1_HVG = species1_HVGs,
                          spec2_HVG = species2_HVGs,
                          Table_homo = homo_mat)
  return(Object)
}

