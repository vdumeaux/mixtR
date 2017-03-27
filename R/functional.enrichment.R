#' Functional enrichment analyses for genes partitioned by sig.ranksum()
#'
#' Enrichemnt is estimated by the hypergeometric minimum-likelihood p-values,
#' computed with the function ‘dhyper’ (equivalent to one-sided Fisher
#' exact test). P-values are then adjusted for multiple testing.
#'
#' @param mixt.ranksum Output of sig.ranksum()
#' @param tissue character string that provides the name of tissue of interest
#' @param cohort.name character string that provides the name of the patient
#' cohort to analyze, default is set to 'all'
#' @param functional.groups a list of length 2: each object is a named vector of
#' gene set appartenance for each tissue.
#' @param p.adjust.method correction method. See `?stats::p.adjust`
#' @param mc.cores number of cores to use
#'
#' @return The output is a list containing the following objects
#' \item{results}{data frame of enrichment results for each gene set/module.
#' Enrichment is estimated for all genes in the gene set/module and
#' for genes that are positively (up genes) or negatively (dn genes) correlated
#' with the patient ranksum. Number of genes in common with signature are
#' provided for all, up, and down gene groups.}
#' \item{updn.common}{names of genes in common between gene set/module and
#' functional signature.}
#' \item{up.common}{names of genes in common between up genes in gene set/module
#' and functional signature.}
#' \item{up.common}{names of genes in common between dn genes in gene set/module
#' and functional signature.}
#' @export
#'
functional.enrichment <- function (mixt.ranksum, tissue, cohort.name= "all",
                                   functional.groups, p.adjust.method = "BH",
                                   mc.cores = 2)
{
  bs <- mixt.ranksum[[tissue]]

  all.genes <- unique(unlist(lapply(bs, function(x)
    rownames(x[[cohort.name]]$dat))))

  msigdb <- sapply(functional.groups$all.sig, function(sig)
    intersect(all.genes, sig), simplify = FALSE)
  sig.set <- functional.groups$set[sapply(msigdb, length) >= 5]
  msigdb2 <- msigdb[sapply(msigdb, length) >= 5]
  all.genes <- intersect(all.genes, unlist(msigdb2))

  ## report size and set of each signature
  sig.name<-names(msigdb2)
  sig.size<-sapply(msigdb2, length)

  ret <- parallel::mclapply(names(bs), function(mod) {

      bs.mod <- bs[[mod]][[cohort.name]]

      ## msigdb enrichment for all genes in each module
      updn.genes <- rownames(bs.mod$dat)
      updn.genes <- intersect(updn.genes, all.genes)
      updn.intersections <- lapply(msigdb2, function(sig)
        intersect(updn.genes, sig))
      updn.common<-sapply(updn.intersections, length)
      updn.pval <- stats::p.adjust(lapply(msigdb2, function(msig)
        hyper.fisher(pop1=updn.genes, pop2=msig, bckg=all.genes)),
        method=p.adjust.method)

      ## msigdb enrichment for up genes only in each module
      if (!length(bs.mod$up)==0){
      up.genes <- rownames(bs.mod$dat)[which(bs.mod$up.dn > 0)]
      up.genes <- intersect(up.genes, all.genes)
      up.intersections <- lapply(msigdb2, function(sig)
        intersect(up.genes, sig))
      up.common<-sapply(up.intersections, length)
      up.pval<-stats::p.adjust(lapply(msigdb2, function(msig)
        hyper.fisher(pop1=up.genes, pop2=msig, bckg=all.genes)),
        method=p.adjust.method)} else {
          up.pval<-NA
          up.common<-NA
          up.intersections<-NA
        }

      ## msigdb enrichment for dn genes only in each module
      if (!length(bs.mod$dn)==0){
        dn.genes <- rownames(bs.mod$dat)[which(bs.mod$up.dn < 0)]
        dn.genes <- intersect(dn.genes, all.genes)
        dn.intersections <- lapply(msigdb2, function(sig)
          intersect(dn.genes, sig))
        dn.common<-sapply(dn.intersections, length)
        dn.pval <- stats::p.adjust(lapply(msigdb2, function(msig)
          hyper.fisher(pop1=dn.genes, pop2=msig, bckg=all.genes)),
          method=p.adjust.method)} else {
            dn.pval<-NA
            dn.common<-NA
            dn.intersections<-NA
          }

    ## return pvalues and genes in common between gene list and msigdb2 signature
    return(list(results=as.data.frame(cbind(sig.set, sig.size,
                                            updn.common, updn.pval,
                                            up.common, up.pval,
                                            dn.common, dn.pval),
                                      stringsAsFactors = FALSE,
                                      row.names = sig.name),
                updn.common=updn.intersections,
                up.common=up.intersections,
                dn.common=dn.intersections))
      }, mc.cores=mc.cores)

  names(ret) <- names(bs)
  return(ret)
}
