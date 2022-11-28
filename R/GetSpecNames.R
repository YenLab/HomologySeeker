#' @title  Get Species official name
#'
#' @description Check available species official scientific name
#'
#' @importFrom methods setClass
#' @importFrom magrittr %>% set_names
#' @importFrom dplyr select mutate
#' @importFrom biomaRt useEnsembl getLDS
#' @import Seurat SeuratObject
#'
#' @param Species Names of species to be queried. Case is ignored.
#'
#' @return list of species scientific name and species name
#' @export
#' @examples
#' \dontrun{
#' # Get scientific name for single or multiple species:
#' species_name <- GetSpecNames(c("mouse","human"))
#'
#' # Or get a full list:
#' species_name <- GetSpecNames(Species = NULL)
#' }
GetSpecNames <- function(Species = NULL){
  '%ni%' <- Negate("%in%")
  SpeciesName <- get("SpeciesName")
  tmp <- SpeciesName
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
