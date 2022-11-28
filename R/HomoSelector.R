#' @title Homologous information collector between two Species
#'
#' @description Get homologous gene set and information between two species
#'
#' @include GetSpecNames.R
#'
#' @importFrom methods setClass
#' @importFrom magrittr %>% set_names
#' @importFrom dplyr select mutate
#' @importFrom biomaRt useEnsembl getLDS
#' @import Seurat SeuratObject
#'
#' @param species1,species2 Names of species in comparative analysis. Case is ignored
#' @param homotype Homologous gene relationship between two species:
#' \itemize{
#'  \item{"ortholog_one2one"} \strong{:} One-to-one orthologues(Default).
#'  \item{"ortholog_one2many"} \strong{:} One-to-many orthologues.
#'  \item{"ortholog_many2many"} \strong{:} Many-to-many orthologues.
#'  }
#' @param version Ensembl version to be connected. See biomaRt::useEnsembl() for detailed information.
#'
#' @return Table of homologous genes list between species 1 and species 2
#' @export
#'
#' @details Due to the query issue after the latest update to Ensembl version 106 as described in
#' 'https://github.com/grimbough/biomaRt/issues/61', currently the version '105' is the most stable one.
#' If other Ensembl version are needed, the online source of BioMart ('http://www.ensembl.org/biomart/martview')
#' are avilable for querying homologous information among species. We also implement available datasets: Human2Mouse,
#' Human2Zebrafish and Mouse2Zebrafish in HomoSeeker (load with data() function) and deposite in our GitHub
#' repository under Ensembl verison 108:'https://github.com/Soap4/Data'.
#'
#' @examples
#' \dontrun{
#' # Get Homologous gene list between mouse and human
#'
#' Homo_mat <- HomoSelector(species1 = "mouse",species1 = "human")
#'
#' # Or get available datasets in HomoSeeker:
#'
#' Human2Mouse <- read.csv("https://github.com/Soap4/Data/files/10098277/Orthologues_Human2Mouse_v108.csv",row.names = 1)
#' Human2Zebrafish <- read.csv("https://github.com/Soap4/Data/files/10098276/Orthologues_Human2Zebrafish_v108.csv",row.names = 1)
#' Mouse2Zebrafish <- read.csv("https://github.com/Soap4/Data/files/10098275/Orthologues_Mouse2Zebrafish_v108.csv",row.names = 1)
#' }
HomoSelector <- function(species1,
                         species2,
                         homotype = "ortholog_one2one",
                         version = 105){
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
    species1_mart <- useEnsembl(biomart = "ensembl",
                                # GRCh = GRCh,
                                # mirror = mirror,
                                # host = host,
                                version = version,
                                dataset = paste0(species1_name$Scientific_Name,"_gene_ensembl"))
    message(paste0("Peparing available ",species2_name$Species_Name," datasets"))
    species2_mart <- useEnsembl(biomart = "ensembl",
                                # GRCh = GRCh,
                                # mirror = mirror,
                                # host = host,
                                version = version,
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
