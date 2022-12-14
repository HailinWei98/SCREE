% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/single_gene_matrix_regression.R
\name{single_gene_matrix_regression}
\alias{single_gene_matrix_regression}
\title{function definitions ##### Return a list with 2 matrix:
index matrix, but values of NegCtrl columns are all 1;
genes * cells matrix, similar to scale_data matrix,
but all values over quantile with
setting of outlier_threshold will be replace with the outlier_threshold quantile}
\usage{
single_gene_matrix_regression(
  targetobj,
  ngctrlgene = c("NonTargetingControlGuideForHuman"),
  indmatrix = NULL,
  selected_genes_list = NULL,
  NTC_baseline = TRUE
)
}
\description{
function definitions ##### Return a list with 2 matrix:
index matrix, but values of NegCtrl columns are all 1;
genes * cells matrix, similar to scale_data matrix,
but all values over quantile with
setting of outlier_threshold will be replace with the outlier_threshold quantile
}
