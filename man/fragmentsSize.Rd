% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fragmentsSize.R
\name{fragmentsSize}
\alias{fragmentsSize}
\title{Plot Distribution of Fragments Size}
\usage{
fragmentsSize(
  mtx,
  fragments,
  CBCindex = 4,
  startIndex = 2,
  endIndex = 3,
  maxSize = 1000,
  title.size = 25,
  x.text.size = 16,
  x.title.size = 20,
  y.text.size = 16,
  y.title.size = 20,
  plot.save = TRUE,
  prefix = ".",
  label = "",
  width = 7,
  height = 7,
  png_res = 720
)
}
\arguments{
\item{mtx}{SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows.}

\item{fragments}{Data frame or directory of fragments file.}

\item{CBCindex}{Cell barcode index, indicates the column index of cell barcode. Default is 4, in line with a standard fragments file.}

\item{startIndex}{Start site index, indicates the column index of fragments start site. Default is 2, in line with a standard fragments file.}

\item{endIndex}{End site index, indicates the column index of fragments end site. Default is 3, in line with a standard fragments file.}

\item{maxSize}{Maximum size of fragments. Default is 1000.}

\item{title.size}{Numeric, title size of the plot. Default is 25.}

\item{x.text.size}{Numeric, x-axis text size of the plot. Default is 16.}

\item{x.title.size}{Numeric, x-axis title size of the plot. Default is 20.}

\item{y.text.size}{Numeric, y-axis text size of the plot. Default is 16.}

\item{y.title.size}{Numeric, y-axis title size of the plot. Default is 20.}

\item{plot.save}{Logical, save plots or not. Default is \code{TRUE}.}

\item{prefix}{Path to save the plots. Default is current directory.}

\item{label}{The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".}

\item{width}{Width of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.}

\item{height}{Height of the graphics region of the pdf file in inches, for both png and pdf format. Default is 7.}

\item{png_res}{The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.}
}
\description{
Calculate fragments size and visualize the distribution.
}
