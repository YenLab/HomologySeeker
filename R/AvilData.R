#' @title  Query and download available datasets
#'
#' @description Query and download available datasets in HomologySeeker
#'
#' @export
#'
#'#' @examples
#' \dontrun{
#' # Download homologous gene list using download.file() function
#'
#' Human2Mouse <- read.csv(url("https://github.com/Soap4/Data/raw/main/scRNA-seq/Orthologues_Human_Mouse.csv"))
#'
#' # Or you can load it with data() or get() function
#'
#' Human2Mouse <- get("mouse2human")
#' Mouse2Zebrafish <- get("mouse2zebrafish")
#' Human2Zebrafish <- get("human2zebrafish")
#'
#' }
AvilData <- function(){
    return(data.frame(Available_dataset = c("Human-Mouse",
                                            "Human-Zebrafish",
                                            "Mouse-Zebrafish"),
                      Ensembl_Gene_Version = rep(108,3),
           Source = c("https://github.com/Soap4/Data/raw/main/scRNA-seq/Orthologues_Human_Mouse.csv",
                   "https://github.com/Soap4/Data/raw/main/scRNA-seq/Orthologues_Human_Zebrafish.csv",
                   "https://github.com/Soap4/Data/raw/main/scRNA-seq/Orthologues_Mouse_Zebrafish.csv")))
}
