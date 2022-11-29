#' @title  Query and download available datasets
#'
#' @description Query and download available datasets in HomoSeeker
#'
#' @export
#'
AvilData <- function(){
    return(data.frame(Available_dataset = c("Human2Mouse",
                                            "Human2Zebrafish",
                                            "Mouse2Human",
                                            "Mouse2Zebrafish",
                                            "Zebrafish2Human",
                                            "Zebrafish2Mouse"),
                      Reference = c("Human","Human","Mouse","Mouse","Zebrafish","Zebrafish"),
                      Query = c("Mouse","Zebrafish","Human","Zebrafish","Human","Mouse"),
                      Ensembl_Gene_Version = rep(108,3),
           Source = c("https://github.com/Soap4/Data/files/10098277/Orthologues_Human2Mouse_v108.csv",
                   "https://github.com/Soap4/Data/files/10098276/Orthologues_Human2Zebrafish_v108.csv",
                   "https://github.com/Soap4/Data/files/10110437/Orthologues_Mouse2Human_v108.csv",
                   "https://github.com/Soap4/Data/files/10098275/Orthologues_Mouse2Zebrafish_v108.csv",
                   "https://github.com/Soap4/Data/files/10110442/Orthologues_Zebrafish2Human_v108.csv",
                   "https://github.com/Soap4/Data/files/10110535/Orthologues_Zebrafish2Mouse_v108.csv")))
}
