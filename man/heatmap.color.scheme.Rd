% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap.R
\name{heatmap.color.scheme}
\alias{heatmap.color.scheme}
\title{Break points and color vector for heatmaps.}
\usage{
heatmap.color.scheme(low.breaks = seq(-3, 0, length.out = 21),
  high.breaks = seq(0, 3, length.out = 21), masked.color = "black",
  col.palette = c("blue.red"), plot = FALSE, gp = gpar())
}
\arguments{
\item{low.breaks}{breakpoints for the low color}

\item{high.breaks}{breakpoints for the low color.}

\item{masked.color}{color to be used for the gap between the low and high colors.}

\item{col.palette}{either a palette name or index, or a full specification of an hcl palette.
See the 'supported.palettes' definition in the code for examples.}

\item{plot}{if TRUE, will plot the HCL color space for the chosen color palette.}

\item{gp}{graphical parameters passed to grid.points and grid.lines when 'plot' is TRUE.}
}
\value{
a list with three fields:
\item{breaks}{vector with break points, size is one more that 'col'}
\item{col}{vector with colors, size one less than 'breaks'}
\item{masked.idx}{idx into 'col' and 'breaks' pointing to the color and
the left break point of the masked region (if such a region exists).}
}
\description{
\code{heatmap.color.scheme} computes break points and color vector for heatmaps. Can also plot
the HCL space for easier customization of colors.
}
\keyword{internal}
