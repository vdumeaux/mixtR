#' Patient linear ordering based on gene set expression
#'
#' @description
#' \code{sig.ranksum} linearly orders samples based on expression of a defined
#' set of genes. Patients with the same ranksum are assigned their 'average' rank,
#' i.e., they all receive the same rank value, which might not be an integer.
#'
#' Genes in the set are partitioned into two groups around myeloids (𝑀1 and 𝑀2)
#' using correlation as the distance metric. Each gene in 𝑀1 and 𝑀2 is ordered
#' from high to low and low to high expression, respectively.
#' Expression values of each gene are then replaced by their ranks across patients.
#' The sum of gene ranks (ranksum) are then used to linearly ordered patients.
#' Note that 'up' and 'dn' or 'ns' are used to index into exprdata.
#'
#' @param x.dat a list with at least 3 objects. exprs: gene expression data,
#' must be a matrix. Genes are in rows, patients in columns.
#' NAs in exprdata are not supported. clinical: clinical information
#' data frame. cohorts: list of patient subgroups
#' @param up indices of up-regulated genes
#' @param dn indices of down-regulated genes
#' @param ns indices of genes whose directions have not been specified.
#' @param n number of permutations to perform to compute the region of indendence
#' @param middle.range percentile point of the distribution of random patient
#' ranks to define the region of independence
#' @param seed random seed number
#' @param mc.cores number of cores to use
#'
#' @return
#' The output is a list containing the following objects
#' \item{pat.order}{vector providing the patient linear ordering.
#' data[, pat.order] orders the
#' columns of
#'  expression matrix according to the sum of gene ranks (gene ranksums).}
#' \item{ranksum}{vector of gene ranksums acr}
#' m{gene.order}{gene ordering. exprdata[gene.order, ] sorts rows of exprdata
#' so that genes within each group around the two myeloids are ordered according
#' to the strength of their association with the patient linear ordering.
#' Most correlated genes are at each end of the
#' ordering}
#' \item{dat}{matrix of expression of \code{up and dn} or \code{ns} genes
#' reordered by patient and gene ordering.}
#' \item{up.dn}{vector with nrow(dat) values. -1, indicates down-gene
#'        and 1 indicates up-gene.}
#' \item{up}{vector of indexes pointing at M1 genes in exprdata defined as upregulated
#' i.e positively correlated with patient ordering.}
#' \item{dn}{vector of indexes pointing at M2 genes in exprdata defined as downregulated
#' i.e negatively correlated with patient ordering.}
#'
#' @export
#'
sig.ranksum <- function(x.dat, up = NULL, dn = NULL, ns = NULL, n = 1000,
                        middle.range = 0.95, seed = 123456, mc.cores = 2)
  {
  if (length(up) == 0 && length(dn) == 0 && length(ns) == 0)
        stop("no indices were specified")
    if (length(ns) > 0 && (length(up) > 0 || length(dn) > 0))
        stop("both directional and non-directional indices were specified")
    if (is.logical(up))
        up <- which(up)
    if (is.logical(dn))
        dn <- which(dn)
    if (is.logical(ns))
        ns <- which(ns)

    if (is.null(x.dat$cohorts))
      stop ("cohorts is missing")
    if (is.null(x.dat$exprs))
      stop ("exprs is missing")
    if (!identical(class(x.dat$cohorts), "list"))
      stop("cohorts should be a list")
    if (is.null(x.dat$cohorts))
      stop("members of cohorts should be named")
    if (!is.null(x.dat$cohorts) & any(is.na(names(x.dat$cohorts))))
      stop("all members of cohorts be named")
    if (!all(unlist(sapply(x.dat$cohorts, function(x) {
      x %in% colnames(x.dat$exprs)})
      )))
      stop ("cohorts should be a list of character vectors of
             sample names as given in column names of exprs")


    bresat <- lapply(x.dat$cohorts, function(pat) {

        exprdata <- x.dat$exprs[, colnames(x.dat$exprs) %in% pat]

        if (ncol(exprdata) < 2) {
            ret <- list()
            up <- c(up, ns)
            ret$rank <- 1
            ret$pat.order <- 1
            ret$gene.order <- c(up, dn)
            ret$dat <- exprdata[ret$gene.order, , drop = FALSE]
            ret$up.dn <- c(rep(1, length(up)), rep(-1, length(dn)))
            return(ret)
        }

        if (length(ns) > 0) {
            tmp <- matrixStats::rowSds(exprdata[ns, , drop = FALSE]) == 0
            zero.sd.idx <- ns[tmp]
            ns <- ns[!tmp]

            if (length(ns) == 1) {
                up <- ns
            } else if (length(ns) == 2) {
                if (stats::cor(exprdata[ns[1], ], exprdata[ns[2], ],
                               use = "pairwise") < 0) {
                  up <- ns[1]
                  dn <- ns[2]
                } else {
                  up <- ns
                }
            } else if (length(ns) > 2) {
                diss <- 1 - stats::cor(t(exprdata[ns, , drop = FALSE]),
                                       use = "pairwise")
                diss[which(is.na(diss))] <- 1
                diss[diss >= 1] <- diss[diss >= 1] + 1

                clustering <- cluster::pam(diss, k = 2, diss = TRUE,
                                           cluster.only = TRUE)
                up.cluster <- which.max(table(clustering))
                up <- ns[which(clustering == up.cluster)]
                dn <- ns[which(clustering != up.cluster)]
                up.dn.cor <- stats::cor(t(exprdata[up, , drop = FALSE]),
                                        t(exprdata[dn, , drop = FALSE]),
                                        use = "pairwise")
                if (sum(up.dn.cor < 0, na.rm = TRUE) < length(up) * length(dn)/2) {
                  up <- ns
                  dn <- NULL
                }
            }

            if (length(zero.sd.idx) > 0) {
                up <- c(up, zero.sd.idx)
            }
        }

        ranksum <- double(ncol(exprdata))
        col.counts <- rep(0, ncol(exprdata))
        if (length(up) != 0) {
            dat <- exprdata[up, , drop = FALSE]
            ranksum <- rowSums(apply(dat, 1, function(x) {
                rank(x, length(x), "average")
            }))
            col.counts <- colSums(!is.na(dat))
        }
        if (length(dn) != 0) {
            dat <- exprdata[dn, , drop = FALSE]
            ranksum <- ranksum + rowSums(ncol(exprdata) - apply(dat, 1,
                                                              function(x) {
                rank(x, length(x), "average")
            }) + 1)
            col.counts <- col.counts + colSums(!is.na(dat))
        }

        ranksum <- ranksum/col.counts
        rank <- rank(ranksum, length(ranksum), "average")


       if (length(up) == 0)
            up <- NULL
        if (length(dn) == 0)
            dn <- NULL

        ret <- list()
        ret$pat.order <- order(-rank)
        ret$ranksum <- ranksum

        ## computes correlation of a single row of the expression data with the
        ## the patient ordering. Used to obtain ret$gene.order.
        gene.cor <- function(gene.idx, is.up, exprdata, pat.order) {
            gene.expr <- exprdata[gene.idx, pat.order]
            if (is.up == TRUE)
                gene.expr <- -gene.expr
            stats::cor(rank(gene.expr), 1:ncol(exprdata))
        }

        up.cor <- sapply(up, gene.cor, TRUE, exprdata, ret$pat.order)
        dn.cor <- sapply(dn, gene.cor, FALSE, exprdata, ret$pat.order)
        ret$gene.order <- c(up[rev(order(up.cor))], dn[order(dn.cor)])
        ret$dat <- exprdata[ret$gene.order, ret$pat.order, drop = FALSE]
        ret$up.dn <- c(rep(1, length(up)), rep(-1, length(dn)))
        ret$up <- up
        ret$dn <- dn

        datrank.up <- datrank.dn <- NULL
        if (any(ret$up.dn > 0))
            datrank.up <- t(apply(ret$dat[ret$up.dn > 0, , drop = FALSE], 1,
                                  function(x) {
                rank(x, "average", na.last = "keep")
            }))
        if (any(ret$up.dn < 0))
            datrank.dn <- ncol(ret$dat) - t(apply(ret$dat[ret$up.dn < 0, ,
                                                          drop = FALSE], 1,
                                                  function(x) {
                rank(x, "average", na.last = "keep")
            })) + 1
        datrank <- rbind(datrank.up, datrank.dn)
        nvals.cols <- nrow(datrank) - c(colSums(is.na(datrank)), 0)
        nvals.rows <- ncol(datrank) - rowSums(is.na(datrank))
        set.seed(seed)
        random.cols <- t(sapply(1:nrow(datrank), function(i) {
            stats::runif(n, 1, nvals.rows[i] + 1)
        }))
        rand.dist <- unlist(parallel::mclapply(1:n, function(i) {
            datrank.rand.col <- cbind(datrank, random.cols[, i])
            datrank.rand.col <- t(apply(datrank.rand.col, 1, function(x) {
                rank(x, "average", na.last = "keep")
            }))
            ranksums <- colSums(datrank.rand.col, na.rm = TRUE)/nvals.cols
            utils::tail(rank(ranksums, "average", na.last = TRUE), 1)
        }, mc.cores = mc.cores)) - 1
        roi <- graphics::hist(rand.dist, breaks = 0:(ncol(datrank) + 1),
                              plot = FALSE)$counts

        ret$roi <- roi

        random.ranks.cdf <- cumsum(roi)/sum(roi)
        left <- max(c(0, which(random.ranks.cdf < (1 - middle.range)/2)))
        right <- min(which(random.ranks.cdf > 1 - ((1 - middle.range)/2))) - 1
        roi.cat <- rep(2, length(rank))
        roi.cat[rank < left] <- 1
        roi.cat[rank > right] <- 3

        ret$roi.cat <- roi.cat

        return(ret)
    })

    names(bresat) <- names(x.dat$cohorts)
    return(bresat)
}
