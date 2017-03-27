#' Compute associations between gene set/module and patient clinical variables
#'
#' @description
#' Using ranksums to capture gene set/module expression, we can ask how gene
#' sets/modules in each tissue are differentially expressed
#' according to patient’s clinicopathological variables.
#' The type of the clinicopathological attribute (categorical or continuous)
#' determines the underlying statistical test.
#' Pearson correlation (Student asymptotic p-value) is used to test association
#' between a given module and continuous patient attributes (eg. age).
#' Analysis of Variance (ANOVA) is used to test association
#' between a given module and categorical patient attributes (eg. ER status).
#' @param mixt.dat Data object from matched tissues.
#' @param mixt.ranksum output of sig.ranksum()
#' @param tissue character string that provides the name of tissue we want to
#' look at
#' @param cohort.name character string that provides the name of the patient
#' cohort to analyze, default is set to 'all'
#'
#' @return matrix of p-values determining significance of association between
#' clinical variables (in row) and gene sets/modules (in column)
#'
#' @export
clinical.ranksum <- function(mixt.dat, mixt.ranksum, tissue,
                             cohort.name = "all") {

  bs <- mixt.ranksum[[tissue]]
  ranksum.df <- data.frame(lapply(bs, function(x) unlist(lapply(x, "[", "ranksum"))))
  ranksum.df <- ranksum.df[grep(cohort.name, rownames(ranksum.df)), ]

  patients <- mixt.dat[[tissue]]$cohorts[[cohort.name]]
  dat.cl <- mixt.dat[[tissue]]$clinical[colnames(mixt.dat[[tissue]]$exprs) %in%
                                          patients, ]
  quant.var <- names(dat.cl)[sapply(dat.cl, class) == "numeric"]
  qual.var <- names(dat.cl)[!names(dat.cl) %in% quant.var]

  tmp <- data.frame()
  ## qualitative variables
  if(length(qual.var) > 0){
    dat.cl.qual <- dat.cl [, qual.var, drop = FALSE]
    sel <- sapply(dat.cl.qual, function(data) length(levels(factor(data))) > 1)
    dat.cl.qual <- dat.cl.qual[ , sel, drop = FALSE]

    if (ncol(dat.cl.qual) > 0){
      tmp <- plyr::laply(dat.cl.qual, function(y) {
        plyr::laply(ranksum.df, function(x) {
          stats::anova(stats::lm(x ~ y))$`Pr(>F)`[1]
        })
      })
      rownames(tmp) <- names(dat.cl.qual)
      colnames(tmp) <- names(ranksum.df)
    }
  }

  # quantitiative variables
  if(length(quant.var) > 0){
    tmp.cor = stats::cor(ranksum.df, dat.cl[, quant.var, drop = FALSE], use = "p")
    tmp.cor.pvalue = WGCNA::corPvalueStudent(tmp.cor, length(patients))
    tmp <- rbind(tmp, t(tmp.cor.pvalue))
  }

  clinicalVars = rownames(tmp)
  cols = colnames(tmp)
  tmp = cbind(clinicalVars, tmp)
  colnames(tmp) = c("Clinical", cols)

  return(tmp)
}

