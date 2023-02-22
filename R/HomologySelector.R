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
#' @param RefSpec,QuySpec Names of reference(ref) and query(quy) species in comparative analysis. Case is ignored.
#' See GetSpecNames() for detailed list.
#' @param homotype Homologous gene relationship between two species:
#' \itemize{
#'  \item{ortholog_one2one} \strong{:} One-to-one orthologues(Default).
#'  \item{ortholog_one2many} \strong{:} One-to-many orthologues.
#'  \item{ortholog_many2many} \strong{:} Many-to-many orthologues.
#'  }
#' @param version Ensembl version to be connected. See biomaRt::useEnsembl() for detailed information.
#' @param usedataset Whether use pre-built datasets from Ensembl gene version 108
#'
#' @return Table of homologous genes list between species 1 and species 2
#' @export
#'
#' @details Due to the query issue after the latest update to Ensembl version 106 as described in
#' 'https://github.com/grimbough/biomaRt/issues/61', currently the version '105' is the most stable one.
#' If other Ensembl version are needed, the online source of BioMart ('http://www.ensembl.org/biomart/martview')
#' are avilable for querying homologous information among species. We also implement available datasets: Human2Mouse,
#' Human2Zebrafish and Mouse2Zebrafish in HomoSeeker (load with data() function) and deposite in our GitHub
#' repository under Ensembl gene verison 108:'https://github.com/Soap4/Data'.
#'
#' @examples
#' \dontrun{
#' # Get Homologous gene list between mouse and human
#'
#' Homo_mat <- HomoSelector(RefSpec = "mouse",RefSpec = "human")
#'
#' # Or get available datasets in HomoSeeker:
#'
#' Human2Mouse <- read.csv("https://github.com/Soap4/Data/files/10098277/Orthologues_Human2Mouse_v108.csv",row.names = 1)
#' Human2Zebrafish <- read.csv("https://github.com/Soap4/Data/files/10098276/Orthologues_Human2Zebrafish_v108.csv",row.names = 1)
#' Mouse2Zebrafish <- read.csv("https://github.com/Soap4/Data/files/10098275/Orthologues_Mouse2Zebrafish_v108.csv",row.names = 1)
#' }
HomoSelector <- function(RefSpec,
                         QuySpec,
                         homotype = "ortholog_one2one",
                         version = 105,
                         usedataset = TRUE){
  if(class(try(GetSpecNames(RefSpec),silent = T))=="try-error"){
    GetSpecNames(RefSpec)
  }
  if(class(try(GetSpecNames(QuySpec),silent = T))=="try-error"){
    GetSpecNames(QuySpec)
  }
  if(GetSpecNames(RefSpec)$Species_Name==GetSpecNames(QuySpec)$Species_Name){
    stop("\n","It appears that RefSpec and QuySpec are the same.","\n",
         "Please check your 'RefSpec' and 'QuySpec' input")
  }
  if(usedataset & all(GetSpecNames(c(RefSpec,QuySpec))[,2] %in% c("Human","Mouse","Zebrafish"))){
    AvilData <- AvilData()
    RefSpec_name <- GetSpecNames(RefSpec)
    QuySpec_name <- GetSpecNames(QuySpec)
    # url <- AvilData$Source[intersect(grep(RefSpec_name$Species_Name,AvilData$Available_dataset),
    #                                  grep(QuySpec_name$Species_Name,AvilData$Available_dataset))]
    used_dataset <- AvilData$Available_dataset[intersect(grep(RefSpec_name$Species_Name,AvilData$Available_dataset),
                                                         grep(QuySpec_name$Species_Name,AvilData$Available_dataset))]
    message("Loading existing ",used_dataset," datasets")
    if(used_dataset=="Human-Mouse"){
      homo_mat <- get("mouse2human")
    }else if(used_dataset=="Human-Zebrafish"){
      homo_mat <- get("human2zebrafish")
    }else{
      homo_mat <- get("mouse2zebrafish")
    }
    colnames(homo_mat)[match(paste0(c("Gene.name","Gene.stable.ID"),".",
                                    rep(c(RefSpec,QuySpec),each = 2)),
                             colnames(homo_mat))] <- c("Gene.name","Gene.stable.ID",
                                                       "Gene.name.1","Gene.stable.ID.1")

  }else{

    RefSpec_name <- GetSpecNames(RefSpec)
    QuySpec_name <- GetSpecNames(QuySpec)

    message(paste0("Peparing available ",RefSpec_name$Species_Name," datasets"))
    RefSpec_mart <- useEnsembl(biomart = "ensembl",
                                version = version,
                                dataset = paste0(RefSpec_name$Scientific_Name,"_gene_ensembl"))
    message(paste0("Peparing available ",QuySpec_name$Species_Name," datasets"))
    QuySpec_mart <- useEnsembl(biomart = "ensembl",
                                version = version,
                                dataset = paste0(QuySpec_name$Scientific_Name,"_gene_ensembl"))
    message("Extracting homologous information")
    homo_mat <- getLDS(mart = RefSpec_mart,
                       attributes = c("external_gene_name","ensembl_gene_id","description",
                                      "chromosome_name","start_position","end_position","strand",
                                      paste0(QuySpec_name$Scientific_Name,"_homolog_orthology_type"),
                                      paste0(QuySpec_name$Scientific_Name,"_homolog_orthology_confidence")),
                       martL = QuySpec_mart,
                       attributesL = c("external_gene_name","ensembl_gene_id","description",
                                       "chromosome_name","start_position","end_position","strand"))
  }
  if(is.null(homotype)){
    homo_mat <- homo_mat[homo_mat %>% select(contains('confi'))==1,]
  }else{
    homo_mat <- homo_mat[homo_mat %>% select(contains('type')) == homotype,]
    homo_mat <- homo_mat[homo_mat %>% select(contains('confi'))==1,]
  }
  message("Homologous Genes Extraction Complete")
  message(paste0("Extracted ",nrow(homo_mat),
                 " Homologous Genes between ",
                 RefSpec_name$Species_Name," and ",QuySpec_name$Species_Name))
  return(homo_mat)
}
