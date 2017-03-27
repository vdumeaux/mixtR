#' Read MSigdb gmt files
#'
#' The gmt files downloaded from the MSigDB website are read and formatted
#'
#' @param inputFile ile character string pointing to the file to read
#'
#' @return list of msigdb gene sets
#'
#' @export

read.gmt.msigdb <- function(inputFile) {

  con  <- file(inputFile, open = "r")

  msgidb <- list()
  current.line = 1
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    line.vector <- strsplit(line, "\t")[[1]]
    name = line.vector[1]
    line <- paste(c(paste(line.vector[3:(length(line.vector)-1)],'\t',sep=""),line.vector[length(line.vector)]),collapse="")
    line.vector <- strsplit(gsub("/"," ",gsub("\t"," ", line,fixed=TRUE),fixed=TRUE)," ")[[1]]
    line.vector = line.vector[line.vector != ""]
    msigdb[[current.line]] <- line.vector
    names(msigdb)[current.line] = name
    current.line = current.line + 1
  }
  close(con)
  return (msigdb)
}

