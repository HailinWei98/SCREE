% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ATAC.R
\name{ATAC_scQC}
\alias{ATAC_scQC}
\title{Quality Control for scATAC-seq Based Input}
\usage{
ATAC_scQC(
  mtx,
  chromation.assay = FALSE,
  peak_frac = 0.01,
  nFeature = c(200, 5e+05),
  nCount = 1000,
  FRiP = 0.1,
  blank_NTC = FALSE,
  title.size = 20,
  x.text.size = 15,
  x.title.size = 15,
  y.text.size = 20,
  pt.size = 0.1,
  plot.show = FALSE,
  plot.save = TRUE,
  raster = FALSE,
  prefix = ".",
  label = "",
  width = 8.3,
  height = 8,
  png_res = 720
)
}
\arguments{
\item{mtx}{SeuratObject or directory to rds file of SeuratObject, with cell in columns and features in rows.}

\item{chromation.assay}{Logical, is the assay in the SeuratObject a ChromatinAssay or not. Default is \code{FALSE}}

\item{peak_frac}{A paramter for filtering low accessibility peaks. By default, only peaks that have counts in at least that fractions of cells are kept. Default is 0.01.}

\item{nFeature}{Limitation of detected feature numbers in each cell, in the format like c(200, 500000). Default is c(200, 500000). If you don't need an upper limit, you can set the upper limit to an extremely large number.}

\item{nCount}{Minimal count numbers in each cell. Default is 1000.}

\item{FRiP}{Minimal FRiP of each cell. Default is 0.1.}

\item{blank_NTC}{Logical, use blank control as negative control or not. Default is \code{FALSE}.}

\item{title.size}{Numeric, title size of the violin plot. Default is 20.}

\item{x.text.size}{Numeric, x-axis text size of the violin plot. Default is 15.}

\item{x.title.size}{Numeric, x-axis title size of the violin plot. Default is 15.}

\item{y.text.size}{Numeric, y-axis text size of the violin plot. Default is 20.}

\item{pt.size}{Point size of the violin plot. Default is 0.1.}

\item{plot.show}{Logical, show the Violin plot or not. Default is \code{TRUE}.}

\item{plot.save}{Logical, save plots or not. Default is \code{TRUE}.}

\item{raster}{Logical, convert points to raster format, will be useful to reduce the storage cost of the output figure or pdf. Default is \code{FALSE}.}

\item{prefix}{Path to save the plots. Default is current directory.}

\item{label}{The prefix label of the output file. Notably, there needs a separator between default file names and the label, so label would be better to be like "label_". Default is "".}

\item{width}{Width of the graphics region of the output file in inches, for both png and pdf format. Default is 8.3.}

\item{height}{Height of the graphics region of the output file in inches, for both png and pdf format. Default is 8.}

\item{png_res}{The nominal resolution in ppi of png file. Higher png_res indicates a bigger and more clear png file. Default is 720.}
}
\description{
Perform quality control for scATAC-seq based input based on peaks fraction, nFeature, nCount, FRiP and sgRNA information.
}
