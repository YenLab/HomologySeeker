#' Title
#' @include HomoHVG_Object.R
#' @include HomoSelector.R
#'
#' @importFrom methods setClass
#' @importFrom magrittr %>% set_names
#' @importFrom dplyr select mutate
#' @importFrom biomaRt useEnsembl getLDS
#' @import Seurat SeuratObject
#'
#' @param species1
#' @param species2
#' @param species1_mat
#' @param species2_mat
#' @param homo_mat
#' @param species1_gene
#' @param species2_gene
#' @param method
#' @param HVGs_method
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
HVGSelector <- function(species1,
                        species2,
                        species1_mat,
                        species2_mat,
                        homo_mat = NULL, # col1: species1 gene name
                                         # col2: species1 gene id
                                         # col3: species2 gene name
                                         # col4: species2 gene id
                        species1_gene = "Gene_sym",# Gene_sym/Gene_id
                        species2_gene = "Gene_sym",# Gene_sym/Gene_id
                        method = "seurat",# Current: seurat
                        HVGs_method = "vst", # Current: vst sct
                        verbose = T){
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
                                          "TRUE","FALSE")}

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

