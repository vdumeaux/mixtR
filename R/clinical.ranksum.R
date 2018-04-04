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
#' @param corType a character string indicating which correlation coefficient
#' is to be computed. Default 'p' for pearson
#' @param nRuns number of permutations
#' @param randomSeed seed number for random number generation.
#' Default set as '1234'
#' @param mc.cores number of cores
#' @param verbose numerical. default > 0 show informational text on progress
#'

#'
#' @return matrix of p-values determining significance of association between
#' clinical variables (in row) and gene sets/modules (in column)
#'
#' @export
clinical.ranksum <- function(mixt.dat, mixt.ranksum, tissue,
                             cohort.name = "all",
                             corType = "p",
                             nRuns=10000, randomSeed = 12345,
                             mc.cores = 2,
                             verbose = 2) {

  if (!tissue %in% names(mixt.ranksum))
    stop ("tissue should be a name of mixt.ranksum")
  if (is.null(mixt.dat[[tissue]]$clinical))
    stop ("clinical is missing")
  if (!cohort.name %in% names(mixt.dat[[tissue]]$cohorts))
    stop ("cohort.name does not match names of cohorts in mixt.dat")

  if (is.null(mixt.dat[[tissue]]$cohorts))
    stop("members of cohorts should be named")
  if (!is.null(mixt.dat[[tissue]]$cohorts) & any(is.na(names(mixt.dat[[tissue]]$cohorts))))
    stop("all members of cohorts be named")
  if (!all(unlist(sapply(mixt.dat[[tissue]]$cohorts, function(x) {
    x %in% colnames(mixt.dat[[tissue]]$exprs)})
  )))
    stop ("cohorts should be a list of character vectors of
             sample names as given in column names of exprs")


  bs <- mixt.ranksum[[tissue]]
  ranksum.df <- data.frame(lapply(bs, function(x) unlist(lapply(x, "[", "ranksum"))))
  ranksum.df <- ranksum.df[grep(cohort.name, rownames(ranksum.df)), ]

  patients <- mixt.dat[[tissue]]$cohorts[[cohort.name]]
  dat.cl <- mixt.dat[[tissue]]$clinical[colnames(mixt.dat[[tissue]]$exprs) %in%
                                          patients, , drop = FALSE]
  quant.var <- names(dat.cl)[sapply(dat.cl, class) == "numeric"]
  qual.var <- names(dat.cl)[!names(dat.cl) %in% quant.var]

  if(length(qual.var) > 0){
    dat.cl.qual <- dat.cl [, qual.var, drop = FALSE]
    sel <- sapply(dat.cl.qual, function(data) length(levels(factor(data))) > 1)
    dat.cl.qual <- dat.cl.qual[ , sel, drop = FALSE]}

    nSamples = length(patients)

    seedSaved = FALSE

    if (!is.null(randomSeed))
    {
      if (exists(".Random.seed"))
      {
        seedSaved = TRUE;
        savedSeed = .Random.seed
      }
      set.seed(randomSeed)
    }

    tmp <- NULL

    tmp <- parallel::mclapply(1:nRuns, function(i){
        set.seed(randomSeed + 2*i + 1)

      if (verbose > 0) print(paste("...working on run", i, ".."))

      if (i > 1)
        {
          useSamples = sample(nSamples)
        } else {
          useSamples = c(1:nSamples)
        }

      if (ncol(dat.cl.qual) > 0){
          cl.anova <- data.frame()
          cl.anova <- plyr::laply(dat.cl.qual, function(y) {
            plyr::laply(ranksum.df, function(x) {
              stats::anova(stats::lm(x ~ y[useSamples]))$`F value`[1]
              })
            })
          rownames(cl.anova) <- names(dat.cl.qual)
          colnames(cl.anova) <- names(ranksum.df)}

      if(length(quant.var) > 0){
        cl.cor <- t(stats::cor(ranksum.df,
                               dat.cl[useSamples, quant.var, drop = FALSE],
                               use = "pairwise.complete.obs", method = corType))
        }

      ret <- rbind(cl.anova, cl.cor)
      return(ret)
      }, mc.cores = mc.cores)

      n.list <- NULL
      for (i in 2:length(tmp)){
        n.list[[i-1]]<- abs(tmp[[i]]) >= abs(tmp[[1]])}

      cl.p <-  Reduce('+', n.list)/9999

      return(cl.p)
}

