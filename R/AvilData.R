#' @title  Query and download available datasets
#'
#' @description Query and download available datasets in HomoSeeker
#'
#' @export
#'
AvilData <- function(){
    return(data.frame(Available_dataset = c("Human-Mouse",
                                            "Human-Zebrafish",
                                            "Mouse-Zebrafish"),
                      Ensembl_Gene_Version = rep(108,3),
           Source = c("https://github.com/Soap4/Data/files/10283572/Orthologues_Human_Mouse.csv",
                   "https://github.com/Soap4/Data/files/10283573/Orthologues_Human_Zebrafish.csv",
                   "https://github.com/Soap4/Data/files/10283574/Orthologues_Mouse_Zebrafish.csv")))
}
