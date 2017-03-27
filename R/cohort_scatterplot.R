#' Plot scatterplot of ranksums
#'
#' Generate scatterplot of ranksums obtained from two gene sets/modules
#'
#' @param mixt.ranksum  output of sig.ranksum()
#' @param mixt.stat output of stat.ranksum()
#' @param x.tissue character string providing the name of tissue on the x axis
#' @param x.module character string providing the name of module on the x axis
#' @param y.tissue character string providing the name of tissue on the y axis
#' @param y.module character string providing the name of module on the y axis
#' @param cohort.name character string that provides the name of the patient cohort.
#' Defaults is set to 'all'
#'
#' @export

cohort_scatterplot <- function(mixt.ranksum, mixt.stat, x.tissue, x.module, y.tissue, y.module,
                               cohort.name = "all") {


  # we need to swap around these to fix accessing the per.cor.p object
  if (x.tissue != names(mixt.ranksum)[1]) {
    tmp = NULL
    xmodule = NULL
    ymodule = NULL
    tmp = x.module
    xmodule = y.module
    ymodule = tmp
  } else {
    xmodule = x.module
    ymodule = y.module
  }

  ### data for scatterplot
  x.ranksum = mixt.ranksum[[x.tissue]][[x.module]][[cohort.name]]$ranksum
  y.ranksum = mixt.ranksum[[y.tissue]][[y.module]][[cohort.name]]$ranksum
  plot.data <- data.frame(x.ranksum = x.ranksum,
                          y.ranksum = y.ranksum)


  p1 <- ggplot2::ggplot(plot.data, ggplot2::aes(x = x.ranksum, y = y.ranksum)) +
    ggplot2::geom_smooth(method = "lm", colour = "white", alpha = 0.2, size = 0.4) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(y = paste(y.module, y.tissue, "module ranksum"), x = paste(x.module,x.tissue, "module ranksum"),
                  title = paste(cohort.name, " patients", " (cor=", as.character(signif(stats::cor.test(plot.data$x.ranksum, plot.data$y.ranksum)$estimate, digits = 1)), ", p=",
                                signif(mixt.stat[[cohort.name]][xmodule, ymodule], digits = 3), ")", sep = "")) +
    ggplot2::theme(legend.position = "none",
                   panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                   axis.line.x = ggplot2::element_line(colour = "grey60"),
                   axis.line.y = ggplot2::element_line(colour = "grey60"),
                   axis.title = ggplot2::element_text(size = 10),
                   plot.title = ggplot2::element_text(hjust = 0, vjust = 1, size = 12))

  print(p1)
}