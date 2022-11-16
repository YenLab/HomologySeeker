#' Title
#'
#' @importFrom methods setClass
#' @importFrom magrittr %>% set_names
#' @importFrom dplyr select mutate
#' @importFrom biomaRt useEnsembl getLDS
#' @import Seurat SeuratObject
#'
#' @param Species
#'
#' @return
#' @export
#'
#' @examples
GetSpecNames <- function(Species = NULL){
  '%ni%' <- Negate("%in%")
  data("SpeciesName")
  tmp <- Dataset
  tmp$Species_Name <- tolower(tmp$Species_Name)
  if(is.null(Species)){
    return(SpeciesName)
  }else if(!all((Species %>% tolower()) %in% tmp$Species_Name)){
    stop(paste0("\nSpecies ",
                paste(Species[which((Species %>% tolower()) %ni% tmp$Species_Name)],collapse = ", "),
                " is not available","\n",
                "Please check the available list using \'SpeciesName(Species = NULL)\'"))
  }else{
    return(SpeciesName[match(Species %>% tolower(),tmp$Species_Name),])
  }
}
