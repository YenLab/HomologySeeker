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
                      Ensembl_version = rep(108,3),
           Source = c("https://github.com/Soap4/Data/files/10098277/Orthologues_Human2Mouse_v108.csv",
                   "https://github.com/Soap4/Data/files/10098276/Orthologues_Human2Zebrafish_v108.csv",
                   "https://github.com/Soap4/Data/files/10098275/Orthologues_Mouse2Zebrafish_v108.csv")))
}
