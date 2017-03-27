#' PLot bxoplot distribution of ranksum for a defined gene set according to ROI categories
#'
#' @param mixt.ranksum Output of sig.ranksum()
#' @param tissue character string that provides the name of tissue we want to look at
#' @param module character string that provides the name of the module we want to generate a heatmap for
#' @param cohort.name character string that provides the name of the patient cohort we select patients from, defaults is set to 'all'
#' @param orderByModule patient ordering can be based on another module/gene set expression.
#' @param orderByTissue the tissue where the orderByModule module is found. Default is the same tissue.
#'
#' @export
cohort_boxplot <- function(mixt.ranksum, tissue, module, cohort.name = "all",
                           orderByTissue=NULL, orderByModule=NULL) {
  if (is.null(orderByModule)) {
    orderByModule <-  module
  }
  if (is.null(orderByTissue)) {
    orderByTissue <- tissue
  }

  bs.order.by <- mixt.ranksum[[orderByTissue]][[orderByModule]][[cohort.name]]
  roi.cat <- mixt.ranksum[[orderByTissue]][[orderByModule]][[cohort.name]]$roi.cat

  bs <- mixt.ranksum[[tissue]][[module]][[cohort.name]]

  plot.data <- data.frame(ranksum = bs$ranksum, cohort = rep(cohort.name, length(bs$ranksum)), roi.cat = roi.cat[match(names(bs.order.by$ranksum), names(bs$ranksum))])
  cat.label <- stats::setNames(c("low", "mid", "high"), 1:3)
  plot.data$roi.cat <- cat.label[plot.data$roi.cat]
  plot.data$roi.cat [is.na(plot.data$roi.cat)] <- "other"
  plot.data$roi.cat.ordered <- factor(plot.data$roi.cat, levels=c("high", "mid", "low", "other"), ordered = TRUE)

  p <- ggplot2::ggplot(data = plot.data, ggplot2::aes(x = roi.cat.ordered, y = ranksum)) +
    ggplot2::geom_boxplot(colour = "black", outlier.shape = "+", outlier.size = 1, fill = NA) +
    ggplot2::geom_point(shape = "+") +
    ggplot2::labs(x = paste(orderByModule, orderByTissue, "ROI module category"), y = paste(module, tissue, "module ranksum", sep = " "),
                  title = paste(cohort.name, " \n(aov p=", as.character(signif(stats::anova(stats::lm(plot.data$ranksum ~ plot.data$roi.cat))$`Pr(>F)`[1], digits = 1)), ")", sep = "")) +
    ggplot2::theme(legend.position = "none",
                   panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                   axis.line.x = ggplot2::element_line(colour = "grey60"),
                   axis.line.y = ggplot2::element_line(colour = "grey60"),
                   axis.title.x = ggplot2::element_text(hjust = 0.2),
                   axis.title = ggplot2::element_text(size = 10),
                   plot.title = ggplot2::element_text(hjust = 0, vjust = 1, size = 12))
  print(p)
}

