#' Compute association between gene sets/modules across tissues
#'
#' Interactions between modules across tissues are identified using a random permutation approach
#' based on the correlation between ranksums of gene expression in gene sets/modules across tissues.
#'
#' @param mixt.ranksum output of sig.ranksum()
#' @param tissue1 name of the first tissue the test is performed for.
#' should be a valid name of mixt.ranksum
#' @param tissue2 name of second tissue the test is performed for.
#' should be a valid name of mixt.ranksum.
#' @param corType a character string indicating which correlation coefficient
#' is to be computed.
#' Default 'p' for pearson
#' @param nRuns number of permutations
#' @param randomSeed seed number for random number generation.
#' Default set as '1234'
#' @param mc.cores number of cores
#' @param verbose numerical. default > 0 show informational text on progress
#'
#' @export
#'
stat.ranksum <- function(mixt.ranksum,
                         tissue1,
                         tissue2,
                      corType = "p",
                      nRuns=10000, randomSeed = 12345,
                      mc.cores = 2,
                      verbose = 2)

{
  dat.ranksum <- lapply(mixt.ranksum, function (bs) {
    data.frame(lapply(bs, function(x) unlist(lapply(x, "[", "ranksum"))))
    })

  cohort.name <- names(mixt.ranksum[[1]][[1]])

  perm.cor.p <- lapply(cohort.name, function(cohort){
    dat.ranksum <- lapply(dat.ranksum, function (x) {x[grep(cohort, rownames(x)), ]})

    result = list()
    tissue1 = names(dat.ranksum)[1]
    tissue2 = names(dat.ranksum)[2]
    nSamples = length(dat.ranksum[[tissue1]][[1]])

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

    mods <- parallel::mclapply(1:nRuns, function(i) {
    set.seed(randomSeed + 2*i + 1);
    if (verbose > 0) print(paste("...working on run", i, ".."));

    if (i > 1)
    {
      useTissue1Samples = sample(nSamples)
      useTissue2Samples = sample(nSamples)

    } else {
      useTissue1Samples = c(1:nSamples)
      useTissue2Samples = c(1:nSamples)}

      mods <- cor(dat.ranksum[[tissue1]][useTissue1Samples, ],
                  dat.ranksum[[tissue2]][useTissue2Samples, ], use=corType)

      rownames(mods) <- names(dat.ranksum[[tissue1]])
      colnames(mods) <- names(dat.ranksum[[tissue2]])

      return(mods)
  }, mc.cores = mc.cores)

    n.list <- NULL
    for (i in 2:length(mods)){
    n.list[[i-1]]<- abs(mods[[i]]) >= abs(mods[[1]])}

    ret <-  Reduce('+', n.list)/9999
    return(ret)})

    names(perm.cor.p) <- cohort.name
    return(perm.cor.p)

}


