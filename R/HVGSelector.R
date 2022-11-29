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
#' @importFrom scmap selectFeatures
#'
#' @import Seurat SeuratObject
#'
#' @param RefSpec,QuySpec Names of reference(ref) and query(quy) species in comparative analysis. Case is ignored
#' @param RefSpec_mat,QuySpec_mat Single cell expression matrix of reference or query species with row as gene symbol/ID and column as sample ID
#' @param homo_mat Table of one to one homologous genes table between species 1 and species 2. See HomoSelector() for further details.
#' @param RefSpec_gene,QuySpec_gene Type of gene name of Single cell expression matrix.
#' \itemize{
#'  \item{Gene_sym} \strong{:} Gene symbel (default).
#'  \item{Gene_id} \strong{:} Gene ID.
#'  }
#' @param method HVG selection method to be used. Currently Supported methods:
#' \itemize{
#'  \item{seurat}
#'  \item{scran}
#'  \item{scmap}
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
HVGSelector <- function(RefSpec,
                        QuySpec,
                        RefSpec_mat,
                        QuySpec_mat,
                        homo_mat = NULL,
                        RefSpec_gene = "Gene_sym",
                        QuySpec_gene = "Gene_sym",
                        method = "seurat",
                        HVGs_method = "vst",
                        verbose = TRUE){
  if(is.null(homo_mat)){
    homo_mat <- HomoSelector(RefSpec = RefSpec,QuySpec = QuySpec)
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
  RefSpec_index <- ifelse(RefSpec_gene == "Gene_sym","Gene.name","Gene.stable.ID")
  QuySpec_index <- ifelse(QuySpec_gene == "Gene_sym","Gene.name.1","Gene.stable.ID.1")
  mat_fil <- mat[mat[,RefSpec_index] %in% rownames(RefSpec_mat) &
                   mat[,QuySpec_index] %in% rownames(QuySpec_mat),]
  RefSpec_mat_homo <- RefSpec_mat[mat_fil[,RefSpec_index],]
  QuySpec_mat_homo <- QuySpec_mat[mat_fil[,QuySpec_index],]

  if(method == "seurat" & HVGs_method == "vst"){
    if(verbose){message(paste0("Identifying HVGs using ",HVGs_method," from ",method))}
    RefSpec_HVG <- CreateSeuratObject(RefSpec_mat_homo) %>%
      suppressWarnings() %>%
      FindVariableFeatures(verbose = F) %>%
      {.[["RNA"]]@meta.features} %>%
      dplyr::mutate(HomoGene=rownames(QuySpec_mat_homo)) %>%
      {.[order(.$vst.variance.standardized,decreasing = T),]}

    QuySpec_HVG <- CreateSeuratObject(QuySpec_mat_homo) %>%
      suppressWarnings() %>%
      FindVariableFeatures(verbose = F) %>%
      {.[["RNA"]]@meta.features} %>%
      dplyr::mutate(HomoGene=rownames(RefSpec_mat_homo)) %>%
      {.[order(.$vst.variance.standardized,decreasing = T),]}
    if(verbose){message("HVGs Screening")}
    RefSpec_HVG$is.HomoHVGs <- ifelse(RefSpec_HVG$vst.variance.standardized>=
                                          mean(RefSpec_HVG$vst.variance.standardized),
                                        "TRUE","FALSE")
    QuySpec_HVG$is.HomoHVGs <- ifelse(QuySpec_HVG$vst.variance.standardized>=
                                          mean(QuySpec_HVG$vst.variance.standardized),
                                        "TRUE","FALSE")
  }
  if(method == "seurat" & HVGs_method == "sct"){
    if(verbose){message(paste0("Identifying HVGs using ",HVGs_method," from ",method))}
    RefSpec_HVG <- CreateSeuratObject(RefSpec_mat_homo) %>%
      suppressWarnings() %>%
      SCTransform(return.only.var.genes = F,verbose = F) %>%
      {.[['SCT']]@SCTModel.list$model1@feature.attributes} %>%
      mutate(HomoGene=rownames(QuySpec_mat_homo)) %>%
      {.[order(.$residual_variance,decreasing = T),]}

    QuySpec_HVG <- CreateSeuratObject(QuySpec_mat_homo) %>%
      suppressWarnings() %>%
      FindVariableFeatures(verbose = F) %>%
      {.[['SCT']]@SCTModel.list$model1@feature.attributes} %>%
      mutate(HomoGene=rownames(RefSpec_mat_homo)) %>%
      {.[order(.$residual_variance,decreasing = T),]}
    if(verbose){message("Screening HVGs")}
      RefSpec_HVG$is.HomoHVGs <- ifelse(RefSpec_HVG$residual_variance>=
                                            mean(RefSpec_HVG$residual_variance),
                                          "TRUE","FALSE")
      QuySpec_HVG$is.HomoHVGs <- ifelse(QuySpec_HVG$residual_variance>=
                                            mean(QuySpec_HVG$residual_variance),
                                          "TRUE","FALSE")
      }
  if(method == "scran"){
    if(verbose){message(paste0("Identifying HVGs using ",method))}
    RefSpec_HVG <- SingleCellExperiment::SingleCellExperiment(list(counts = as.matrix(RefSpec_mat_homo))) %>%
      scuttle::logNormCounts() %>%
      scran::modelGeneVar() %>%
      as.data.frame() %>%
      dplyr::mutate(HomoGene=rownames(QuySpec_mat_homo)) %>%
      {.[order(.$bio,decreasing = T),]}

    QuySpec_HVG <- SingleCellExperiment::SingleCellExperiment(list(counts = as.matrix(QuySpec_mat_homo))) %>%
      scuttle::logNormCounts() %>%
      scran::modelGeneVar() %>%
      as.data.frame() %>%
      dplyr::mutate(HomoGene=rownames(RefSpec_mat_homo)) %>%
      {.[order(.$bio,decreasing = T),]}

    if(verbose){message("Screening HVGs")}
    RefSpec_HVG$is.HomoHVGs <- ifelse(RefSpec_HVG$bio>=mean(RefSpec_HVG$bio),
                                        "TRUE","FALSE")
    QuySpec_HVG$is.HomoHVGs <- ifelse(QuySpec_HVG$bio>=mean(QuySpec_HVG$bio),
                                        "TRUE","FALSE")
  }
  if(method == "scmap"){
    RefSpec_HVG <- SingleCellExperiment(assays = list(counts = as.matrix(RefSpec_mat_homo)),
                                        rowData = list(feature_symbol = rownames(RefSpec_mat_homo))) %>%
      scuttle::logNormCounts() %>%
      scmap::selectFeatures(suppress_plot = T) %>%
      rowData() %>%
      as.data.frame() %>%
      mutate(HomoGene=rownames(QuySpec_mat_homo)) %>%
      {.[order(.$scmap_scores,decreasing = T),]}

    QuySpec_HVG <- SingleCellExperiment(assays = list(counts = as.matrix(QuySpec_mat_homo)),
                                        rowData = list(feature_symbol = rownames(QuySpec_mat_homo))) %>%
      scuttle::logNormCounts() %>%
      scmap::selectFeatures(suppress_plot = T) %>%
      rowData() %>%
      as.data.frame() %>%
      mutate(HomoGene=rownames(RefSpec_mat_homo)) %>%
      {.[order(.$scmap_scores,decreasing = T),]}

    if(verbose){message("Screening HVGs")}
    RefSpec_HVG$is.HomoHVGs <- ifelse(RefSpec_HVG$scmap_scores>=
                                        mean(RefSpec_HVG$scmap_scores,na.rm = T),
                                      "TRUE","FALSE")
    RefSpec_HVG$is.HomoHVGs[is.na(RefSpec_HVG$is.HomoHVGs)] <- "FALSE"
    QuySpec_HVG$is.HomoHVGs <- ifelse(QuySpec_HVG$scmap_scores>=
                                        mean(QuySpec_HVG$scmap_scores,na.rm = T),
                                      "TRUE","FALSE")
    QuySpec_HVG$is.HomoHVGs[is.na(RefSpec_HVG$is.HomoHVGs)] <- "FALSE"
  }
  if(verbose){message("Complete Successfully")}
    message("There are totally: ","\n",
            grep("TRUE",RefSpec_HVG$is.HomoHVGs) %>% length()," Homo-HVGs in RefSpec","\n",
            grep("TRUE",QuySpec_HVG$is.HomoHVGs) %>% length()," Homo-HVGs in QuySpec","\n",
            intersect(RefSpec_HVG$HomoGene[RefSpec_HVG$is.HomoHVGs=="TRUE"],
                      rownames(QuySpec_HVG)[QuySpec_HVG$is.HomoHVGs=="TRUE"]) %>% length(),
            " Homo-HVGs overlapped")

  Object <- HomoHVGObject(RefSpec = RefSpec,
                          QuySpec = QuySpec,
                          RefSpec_mat = RefSpec_mat,
                          QuySpec_mat = QuySpec_mat,
                          RefSpec_mat_homo = RefSpec_mat_homo,
                          QuySpec_mat_homo = QuySpec_mat_homo,
                          RefSpec_HVG = RefSpec_HVG,
                          QuySpec_HVG = QuySpec_HVG,
                          Table_homo = homo_mat)
  return(Object)
}

