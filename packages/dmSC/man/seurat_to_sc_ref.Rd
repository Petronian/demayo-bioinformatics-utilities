% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{seurat_to_sc_ref}
\alias{seurat_to_sc_ref}
\title{Obtain a single-cell reference matrix compatible with CIBERSORTx from a Seurat object.}
\usage{
seurat_to_sc_ref(st_obj, cluster_colname = NULL, ...)
}
\arguments{
\item{st_obj}{The Seurat object to analyze.}

\item{cluster_colname}{The name of the column in st_obj[[]] that contains cluster
information. If NULL, information obtained from
Seurat::Idents(st_obj). Default is NULL.}

\item{assay}{The name of the assay to pull from. If NULL, gets the default assay
from the Seurat object. Default is NULL.}
}
\value{
A data.frame representing the single-cell reference matrix. Careful, since
these can get pretty large in size!
}
\description{
Note that this function is basically just calling GetAssayData and then reassigning
the column names of the resultant table.
}
