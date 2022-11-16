#' Title
#'
#' @include GetSpecNames.R
#'
#' @importFrom methods setClass
#' @importFrom magrittr %>% set_names
#' @importFrom dplyr select mutate
#' @importFrom biomaRt useEnsembl getLDS
#' @import Seurat SeuratObject
#'
#' @param species1
#' @param species2
#' @param homotype
#' @param mirror
#' @param version
#' @param host
#'
#' @return
#' @export
#'
#' @examples
HomoSelector <- function(species1,
                         species2,
                         homotype = "ortholog_one2one",
                         mirror = NULL,
                         version = NULL,
                         host = "https://dec2021.archive.ensembl.org/"){
  if(class(try(GetSpecNames(species1),silent = T))=="try-error"){
    GetSpecNames(species1)
  }else if(class(try(GetSpecNames(species2),silent = T))=="try-error"){
    GetSpecNames(species2)
  }else if(GetSpecNames(species1)$Species_Name==GetSpecNames(species2)$Species_Name){
    stop("\n","It appears that species1 and species2 are the same.","\n",
         "Please check your 'species1' and 'species2' input")
  }else{
    species1_name <- GetSpecNames(species1)
    species2_name <- GetSpecNames(species2)

    message(paste0("Peparing available ",species1_name$Species_Name," datasets"))
    species1_mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                                         # GRCh = GRCh,
                                         mirror = mirror,
                                         version = version,
                                         host = host,
                                         dataset = paste0(species1_name$Scientific_Name,"_gene_ensembl"))
    message(paste0("Peparing available ",species2_name$Species_Name," datasets"))
    species2_mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                                         # GRCh = GRCh,
                                         mirror = mirror,
                                         version = version,
                                         host = host,
                                         dataset = paste0(species2_name$Scientific_Name,"_gene_ensembl"))
    message("Extracting homologous information")
    tmp <- getLDS(mart = species1_mart,
                  attributes = c("external_gene_name","ensembl_gene_id","description",
                                 "chromosome_name","start_position","end_position","strand",
                                 paste0(species2_name$Scientific_Name,"_homolog_orthology_type"),
                                 paste0(species2_name$Scientific_Name,"_homolog_orthology_confidence")),
                  martL = species2_mart,
                  attributesL = c("external_gene_name","ensembl_gene_id","description",
                                  "chromosome_name","start_position","end_position","strand"))
    if(is.null(homotype)){
      tmp <- tmp[tmp %>% select(contains('confi'))==1,]
    }else{
      tmp <- tmp[tmp %>% select(contains('type')) == homotype,]
      tmp <- tmp[tmp %>% select(contains('confi'))==1,]
    }
    message("Extraction Complete! :)")
    message(paste0("Total Homologous Genes Number: ",nrow(tmp)))
  }
  return(tmp)
}
