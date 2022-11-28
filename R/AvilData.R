#' @title  Query and download available datasets
#'
#' @description Query and download available datasets in HomoSeeker
#'
#' @export
#'
AvilData <- function(){
    return(data.frame(Available_dataset = c("Human2Mouse",
                                            "Human2Zebrafish",
                                            "Mouse2Zebrafish"),
                      Ensembl_version = rep(108,3)))
}
