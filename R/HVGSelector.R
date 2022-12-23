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
#' @importFrom scmap selectFeatures
#' @importFrom scran modelGeneVar
#' @importFrom ROGUE SE_fun
#' @importFrom SummarizedExperiment rowData
#' @import Seurat SeuratObject
#'
#' @param RefSpec,QuySpec Names of reference(ref) and query(quy) species in comparative analysis. Case is ignored.
#' See GetSpecNames() for detailed list.
#' @param RefSpec_mat,QuySpec_mat Single cell expression matrix of reference or query species with row as gene symbol/ID and column as sample ID
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
#' @param HVGs_method HVG selection method to be used. Currently Supported methods:
#' \itemize{
#'  \item{\strong{seurat_vst}} \strong{:}  Identify HVGs through variance of standardized gene expression.
#'  Brifly, vst obtain expected variance of gene expression through relationship between logarithmic
#'  variance and logarithmic mean value of gene expression level (local polynomial regression).
#'  The variance of standardized value were used for selecting HVGs. See Seurat::FindVariableFeatures for further information.
#'  \item{\strong{seurat_sct}} \strong{:}  Identify HVGs through Pearson residuals of gene expression using
#'  the regularized negative binomial regression. See Seurat::SCTransform() and sctransform::vst()
#'  for further information.
#'  \item{\strong{scran}} \strong{:}  Identify HVGs through fitting the relationship between variance of genes
#'  logarithmic expression value and mean expression value of genes.
#'  See scran::modelGeneVar() for further information.
#'  \item{\strong{scmap}} \strong{:}  Identify HVGs through fitting the relationship between genes average
#'  expression and dropout rate of single-cell matrix. See scmap::selectFeatures()
#'  for further information.
#'  \item{\strong{ROGUE}} \strong{:}  Identify HVGs through fitting the relationship gene expression entropy and
#'  gene average expression levels. See ROGUE::SE_fun() for further information.
#'  }
#' @param Integrated_mat Whether input matrics are returned by integration
#' @param Integrated_method If Integrated_mat=TRUE, please specify what integration method was used.
#' currently, only seurat_LogNorm and seurat_SCT intergation are supported.
#' If Integrated_method = seurat_LogNorm, please provides 'data' slot of integrated assay from your
#' seurat object. If Integrated_method = seurat_SCT, please provides 'data' slot of integrated
#' assay from your seurat object.
#' @param verbose Whether show calculation progress. Default is TRUE.
#'
#' @return A HomoHVG object with Highly variable homologous gene sets for both species
#' @export
#' @references{
#'
#' }
#'
HVGSelector <- function(RefSpec,
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

  RefSpec_index <- ifelse(RefSpec_gene == "Gene_sym","Gene.name","Gene.stable.ID")
  QuySpec_index <- ifelse(QuySpec_gene == "Gene_sym","Gene.name.1","Gene.stable.ID.1")

  if(is.null(homo_mat)){
    homo_mat <- HomoSelector(RefSpec = RefSpec,QuySpec = QuySpec)
    mat <- homo_mat[,c("Gene.name","Gene.stable.ID",
                       "Gene.name.1","Gene.stable.ID.1")]
    }else if(all(c("Gene.name","Gene.stable.ID",
                   "Gene.name.1","Gene.stable.ID.1") %in% colnames(homo_mat))){
      mat <- homo_mat[,c("Gene.name","Gene.stable.ID",
                       "Gene.name.1","Gene.stable.ID.1")]
    }else{
      mat <- homo_mat[,1:2] %>% set_colnames(c(RefSpec_index, QuySpec_index))
  }
  mat_fil <- mat[mat[,RefSpec_index] %in% rownames(RefSpec_mat) &
                   mat[,QuySpec_index] %in% rownames(QuySpec_mat),]
  RefSpec_mat_homo <- RefSpec_mat[mat_fil[,RefSpec_index],]
  QuySpec_mat_homo <- QuySpec_mat[mat_fil[,QuySpec_index],]

  if(!Integrated_mat){
    if(HVGs_method == "seurat_vst"){
      if(verbose){message("Identifying HVGs using vst from Seurat")}
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
    if(HVGs_method == "seurat_sct"){
      if(verbose){message("Identifying HVGs using sct from Seurat")}
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
    if(HVGs_method == "scran"){
      if(verbose){message(paste0("Identifying HVGs using ",HVGs_method))}
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
    if(HVGs_method == "scmap"){
      if(verbose){message(paste0("Identifying HVGs using ",HVGs_method))}
      RefSpec_HVG <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(RefSpec_mat_homo)),
                                          rowData = list(feature_symbol = rownames(RefSpec_mat_homo))) %>%
        scuttle::logNormCounts() %>%
        scmap::selectFeatures(suppress_plot = T) %>%
        SummarizedExperiment::rowData() %>%
        as.data.frame() %>%
        mutate(HomoGene=rownames(QuySpec_mat_homo)) %>%
        {.[order(.$scmap_scores,decreasing = T),]}

      QuySpec_HVG <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(QuySpec_mat_homo)),
                                          rowData = list(feature_symbol = rownames(QuySpec_mat_homo))) %>%
        scuttle::logNormCounts() %>%
        scmap::selectFeatures(suppress_plot = T) %>%
        SummarizedExperiment::rowData() %>%
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
    if(HVGs_method == "ROGUE"){
      if(verbose){message(paste0("Identifying HVGs using ",HVGs_method))}
      RefSpec_HVG <- ROGUE::SE_fun(RefSpec_mat_homo) %>% #may filter some genes
        as.data.frame() %>%
        magrittr::set_rownames(.$Gene)
      QuySpec_HVG <- ROGUE::SE_fun(QuySpec_mat_homo) %>%
        as.data.frame() %>%
        magrittr::set_rownames(.$Gene)

      index <- mat_fil
      index <- index[index[,RefSpec_index] %in% RefSpec_HVG$Gene &
                       index[,QuySpec_index] %in% QuySpec_HVG$Gene,]

      RefSpec_HVG <- RefSpec_HVG[index[,RefSpec_index],] %>%
        dplyr::mutate(HomoGene=index[,QuySpec_index]) %>%
        {.[order(.$entropy,decreasing = T),]}

      QuySpec_HVG <- QuySpec_HVG[index[,QuySpec_index],] %>%
        dplyr::mutate(HomoGene=index[,RefSpec_index]) %>%
        {.[order(.$entropy,decreasing = T),]}

      if(verbose){message("Screening HVGs")}
      RefSpec_HVG$is.HomoHVGs <- ifelse(RefSpec_HVG$entropy>=
                                          mean(RefSpec_HVG$entropy,na.rm = T),
                                        "TRUE","FALSE")
      RefSpec_HVG$is.HomoHVGs[is.na(RefSpec_HVG$is.HomoHVGs)] <- "FALSE"
      QuySpec_HVG$is.HomoHVGs <- ifelse(QuySpec_HVG$entropy>=
                                          mean(QuySpec_HVG$entropy,na.rm = T),
                                        "TRUE","FALSE")
      QuySpec_HVG$is.HomoHVGs[is.na(RefSpec_HVG$is.HomoHVGs)] <- "FALSE"
    }
  }else{
    if(Integrated_mat & Integrated_method=="seurat_LogNorm"){
      if(verbose){message("Identifying HVGs using vst from Seurat")}
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
    if(Integrated_mat & Integrated_method=="seurat_SCT"){
      if(verbose){message("Identifying HVGs using residual variance")}
      RefSpec_HVG <- data.frame(res_var = apply(RefSpec_mat_homo,1,var))
      RefSpec_HVG$HomoGene <- rownames(QuySpec_mat_homo)
      QuySpec_HVG <- data.frame(res_var = apply(QuySpec_mat_homo,1,var))
      QuySpec_HVG$HomoGene <- rownames(RefSpec_mat_homo)

      RefSpec_HVG <- RefSpec_HVG %>% .[order(.$res_var,decreasing = T),,drop = F]
      QuySpec_HVG <- QuySpec_HVG %>% .[order(.$res_var,decreasing = T),,drop = F]

      if(verbose){message("HVGs Screening")}
      RefSpec_HVG$is.HomoHVGs <- ifelse(RefSpec_HVG$res_var>=mean(RefSpec_HVG$res_var),
                                        "TRUE","FALSE")
      QuySpec_HVG$is.HomoHVGs <- ifelse(QuySpec_HVG$res_var>=mean(QuySpec_HVG$res_var),
                                        "TRUE","FALSE")
    }
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
                          Table_homo = homo_mat,
                          HVGs_method = HVGs_method)
  return(Object)
}

