#' Plot reordered heatmap according to sig.ranksum output
#'
#' Generate heatmap where genes and patients are reordered
#' according the sig.ranksum() output for a given gene set expressed in a given tissue
#'
#' @param mixt.dat Data object from matched tissues
#' @param mixt.ranksum Output of sig.ranksum()
#' @param tissue character string that provides the name of tissue of interest
#' @param module character string that provides the name of the module of interest
#' @param cohort.name character string that provides the name of the patient
#' cohort we select patients from, defaults is set to 'all'
#' @param orderByModule patient ordering can be based on another module/gene
#' set expression.
#' @param orderByTissue the tissue where the orderByModule module is found.
#' Default is the same tissue.
#' @param cl.height dimension of plotting area for clinical variable. Default = 6
#'
#' @import grid
#'
#' @return heatmap plot with patients and genes reordered according to
#' sig.ranksum()
#'
#' @export
#'
cohort_heatmap <- function(mixt.dat, mixt.ranksum, tissue, module,
                           cohort.name = "all",
                           orderByModule = NULL, orderByTissue = NULL,
                           cl.height = 6) {

  graphics::plot.new()
  title = ""
  if (is.null(orderByModule)) {
    orderByModule = module
  }
  if (is.null(orderByTissue)) {
    orderByTissue = tissue
  }

  if (orderByModule == module && orderByTissue == tissue) {
    title = paste(module, " module from ", tissue, sep = "")
  } else {
    title = paste(cohort.name, module, tissue, "ordered by", orderByModule, orderByTissue)
  }

  # heatmap variables
  col.clust = FALSE
  layout.m = matrix(c("key", "title", "", "", "", "", "", "", "", "", "ranksum.text", "ranksum.line", "", "", "", "row.labels.rjust", "heatmap", "", "", "", "", "", "", "",
                      "", "ranks.text", "ranks", "", "", "", "", "", "", "", "", "clinical.labels.rjust", "clinical", "", "", "", "", "", "", "", ""), nrow = 9, ncol = 5, byrow = TRUE)
  layout.m.sum = matrix(c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "row.labels.rjust", "heatmap", "", "", "", "",
                          "", "", "", "", "", "", "", "", "", "", "", "", "", ""), nrow = 9, ncol = 5, byrow = TRUE)
  layout.m.updn = matrix(c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "heatmap", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
                           "", "", "", "", "", "", "", "", "", "", ""), nrow = 9, ncol = 5, byrow = TRUE)
  widths = c(2, 5, 0.25, 0.25, 0.25, 0.25)
  heights = c(1, 0.25, 0.5, 3, 0.25, 0.25, 0.25, cl.height, 0.25)

  scale = "none"
  min.val = -5
  max.val = 5
  key.min = -5
  key.max = 5

  ## define reordered variables
  bs.order.by <- mixt.ranksum[[orderByTissue]][[orderByModule]][[cohort.name]]

  order.by <- bs.order.by$pat.order
  roi <- bs.order.by$roi
  roi.cat <- bs.order.by$roi.cat

  ## define patients to include
  patients <- mixt.dat[[tissue]]$cohorts[[cohort.name]]

  ## define reordered clinical data
  mclinical = NULL
  cl = NULL
  cl <- mixt.dat[[tissue]]$clinical.colors[colnames(mixt.dat[[tissue]]$exprs) %in% patients, ]
  mclinical = cl[order.by, ]

  ## define reordered expression
  bs <- mixt.ranksum[[tissue]][[module]][[cohort.name]]
  data = bs$dat[, match(colnames(bs.order.by$dat), colnames(bs$dat))]

  ## plot heatmap
  ddrs = heatmap.simple(data, clinical = mclinical, layout.mat = layout.m, widths = widths, heights = heights, col.clust = FALSE, row.clust = FALSE, title = title,
                        row.labels = rownames(data), col.labels = rep("", ncol(data)))

  ## plot updn top left heatmap
  up.dn = as.vector(array(1, dim = c(1, length(bs$gene.order))))
  names(up.dn) = rownames(bs$dat)
  if (length(bs$dn) > 0) {
    up.dn[names(up.dn) %in% rownames(mixt.dat[[tissue]]$exprs)[bs$dn]] = -1
  }
  to.plot = (as.matrix(up.dn, ncol = 1))
  color.scheme = heatmap.color.scheme(low.breaks = c(-1.5, 0), high.breaks = c(0, 1.5))
  heatmap.simple(to.plot, scale = scale, layout.mat = layout.m.updn, widths = widths, heights = heights, col.clust = FALSE, row.clust = FALSE, color.scheme = color.scheme)

  ## plot ranks for top left heatmap
  the.layout <- grid.layout(nrow(layout.m), ncol(layout.m), widths = widths, heights = heights)
  mid.vp <- viewport(layout = the.layout, name = "heatmap.mid.vp")
  pushViewport(mid.vp)
  elem = "ranks"
  idx <- which(layout.m == elem, arr.ind = TRUE)
  pushViewport(viewport(name = elem, layout.pos.row = unique(idx[, 1]), layout.pos.col = unique(idx[, 2])))

  rank.colors <- rev(colorspace::diverge_hcl(n = ncol(bs$dat)))
  names(rank.colors) <- colnames(bs$dat)
  rank.colors <- rank.colors[match(colnames(bs.order.by$dat), colnames(bs$dat))]
  ranksum = t(as.matrix(rank.colors))[, names(rank.colors) %in% patients, drop = FALSE]
  heatmap.clinical(ranksum)
  upViewport()
  elem = "ranks.text"
  idx <- which(layout.m == elem, arr.ind = TRUE)
  pushViewport(viewport(name = elem, layout.pos.row = unique(idx[, 1]), layout.pos.col = unique(idx[, 2])))
  heatmap.labels("ranks", type = "row.labels", just = "right")
  upViewport()

  ## ranksum bottom left
  the.layout <- grid.layout(nrow(layout.m), ncol(layout.m), widths = widths, heights = heights)
  top.vp <- viewport(layout = the.layout, name = "heatmap.top.vp")
  pushViewport(top.vp)
  ranksum.plot = bs$ranksum[order.by][colnames(data) %in% patients]
  elem = "ranksum.line"
  idx <- which(layout.m == elem, arr.ind = TRUE)
  pushViewport(viewport(name = elem, layout.pos.row = unique(idx[, 1]), layout.pos.col = unique(idx[, 2]), xscale = c(0.5, length(ranksum.plot) + 0.5), yscale = range(ranksum.plot)))

  grid.rect(gp = gpar(lwd = 0.1))
  grid.polyline(rep(c(0, 1), 4), rep(c(0.2, 0.4, 0.6, 0.8), each = 2), id.lengths = rep(2, 4), gp = gpar(lwd = 0.1, col = "grey70"))
  xrange <- range(1:length(ranksum.plot))
  n <- length(ranksum.plot)
  grid.segments(unit(1:length(ranksum.plot), "native"), rep(0, n), unit(1:length(ranksum.plot), "native"), unit(ranksum.plot, "native"))
  upViewport()

  elem = "ranksum.text"
  idx <- which(layout.m == elem, arr.ind = TRUE)
  pushViewport(viewport(name = elem, layout.pos.row = unique(idx[, 1]), layout.pos.col = unique(idx[, 2])))

  heatmap.labels("ranksum", type = "row.labels", just = "right")
  upViewport()

  ## plot roi lines top left heatmap
  first.ind = length(which(roi.cat == 3))
  last.ind = first.ind + length(which(roi.cat == 2))

  res.random.dist.begin = first.ind
  res.random.dist.end = last.ind

  elem = "heatmap"
  idx <- which(layout.m == elem, arr.ind = TRUE)
  pushViewport(viewport(name = elem, layout.pos.row = unique(idx[, 1]), layout.pos.col = unique(idx[, 2])))

  graphics::par(new = TRUE, fig = gridFIG(), mar = c(0, 0, 0, 0))
  grid.lines(x = unit(c(first.ind/length(patients), first.ind/length(patients)), "npc"),
             gp = gpar(col = "yellow", lty = 1,lwd = 2))
  grid.lines(x = unit(c(last.ind/length(patients), last.ind/length(patients)), "npc"),
             gp = gpar(col = "yellow", lty = 1,lwd = 2))
  upViewport()
}
