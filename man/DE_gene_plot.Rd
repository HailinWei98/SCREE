% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DE_gene_plot.R
\name{DE_gene_plot}
\alias{DE_gene_plot}
\title{Plot DE Gene Numbers for Each Perturbation}
\usage{
DE_gene_plot(
  score,
  pval,
  project = "perturb",
  top = 20,
  select = NULL,
  score_cut = 0.2,
  pval_cut = 0.05,
  y_break = c(50, 200),
  y_height = c(1/5, 4/5),
  title.size = 25,
  legend.text.size = 16,
  legend.title.size = 20,
  x.text.size = 16,
  x.title.size = 20,
  y.text.size = 16,
  y.title.size = 20,
  sort_by = "name",
  plot.save = TRUE,
  prefix = ".",
  label = "",
  width = 7,
  height = 7,
  png_res = 720
)
}
\arguments{
\item{score}{Data frame or directory of score from \code{improved_scmageck_lr}, genes in rows and perturbations in columns.}

\item{pval}{Data frame or directory of p_value from \code{improved_scmageck_lr}, genes in rows and perturbations in columns.}

\item{project}{Title of the barplot. Default is "perturb".}

\item{top}{The top perturbations to plot, with the most differential gene numbers. Default is 20, means plot top 20 perturbations.}

\item{select}{Specific perturbations selected to plot. Default is \code{NULL}, which means to use the top perturbations. If provided a vector of perturbations, \code{top} parameter will be ignored.}

\item{score_cut}{Score cutoff of \code{improved_scmageck_lr} results. Default is 0.2.}

\item{pval_cut}{P-value cutoff of \code{improved_scmageck_lr} results. Default is 0.05.}

\item{y_break}{Break of the barplot, the first item means the maximum y-axis of the bottom ggplot while the second item means the minimum y-axis of the top ggplot. If you don't need a break of y-axis, set the second item to a number bigger than the max number of the up-regulated and down-regulated genes, among all perturbations to plot. Default is c(50, 200).}

\item{y_height}{Heights of the two breaked ggplot to be arranged. Default is c(1/5, 4/5).}

\item{title.size}{Numeric, title size of the barplot. Default is 25.}

\item{legend.text.size}{Numeric, legend text size of the barplot. Default is 16.}

\item{legend.title.size}{Numeric, legend title size of the barplot. Default is 20.}

\item{x.text.size}{Numeric, x-axis text size of the barplot. Default is 16.}

\item{x.title.size}{Numeric, x-axis title size of the barplot. Default is 20.}

\item{y.text.size}{Numeric, y-axis text size of the barplot. Default is 16.}

\item{y.title.size}{Numeric, y-axis title size of the barplot. Default is 20.}

\item{sort_by}{How to order the perturbations labeled in x-axis. Can be one of "name" (name of the perturbations), "select" (the same order as the vector of the \code{select} parameter) and "number" (total differential gene numbers).}

\item{plot.save}{Logical, save plots or not. Default is \code{TRUE}.}

\item{prefix}{Path to save the plots. Default is current directory.}

\item{label}{The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".}

\item{width}{Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.}

\item{height}{Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.}

\item{png_res}{The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.}
}
\description{
Plot potential DE gene numbers for each perturbation based on score and p-value from \code{\link[SCREE]{imporved_scmageck_lr}}.
}
