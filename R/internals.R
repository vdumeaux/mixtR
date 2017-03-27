#' Hypergeometric test
#'
#' computation of hypergeometric test so it matches fisher exact test
#'
#' @keywords internal
#'
hyper.fisher <- function(pop1, pop2, bckg)
{
  s <- length(intersect(pop1, bckg))
  e <- length(intersect(pop2, bckg))
  com <- length(intersect(pop1, pop2))
  ret<- sum(stats::dhyper(com:e,s,length(bckg)-s, e))
  return(ret)
}
