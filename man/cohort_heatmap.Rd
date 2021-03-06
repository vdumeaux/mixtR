% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cohort_heatmap.R
\name{cohort_heatmap}
\alias{cohort_heatmap}
\title{Plot reordered heatmap according to sig.ranksum output}
\usage{
cohort_heatmap(mixt.dat, mixt.ranksum, tissue, module, cohort.name = "all",
  orderByModule = NULL, orderByTissue = NULL, cl.height = 6)
}
\arguments{
\item{mixt.dat}{Data object from matched tissues}

\item{mixt.ranksum}{Output of sig.ranksum()}

\item{tissue}{character string that provides the name of tissue of interest}

\item{module}{character string that provides the name of the module of interest}

\item{cohort.name}{character string that provides the name of the patient
cohort we select patients from, defaults is set to 'all'}

\item{orderByModule}{patient ordering can be based on another module/gene
set expression.}

\item{orderByTissue}{the tissue where the orderByModule module is found.
Default is the same tissue.}

\item{cl.height}{dimension of plotting area for clinical variable. Default = 6}
}
\value{
heatmap plot with patients and genes reordered according to
sig.ranksum()
}
\description{
Generate heatmap where genes and patients are reordered
according the sig.ranksum() output for a given gene set expressed in a given tissue
}
