% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{apply_mapper}
\alias{apply_mapper}
\title{Map using a mapper_table}
\usage{
apply_mapper(to_map, to_ensembl = FALSE, mart = NULL, mapper_table = NULL)
}
\arguments{
\item{to_map}{A vector containing values to map from.}

\item{to_ensembl}{If mapper_table is NULL, whether ENSEMBL IDs should be the keys
(FALSE) or values (TRUE) of the created mapping table. Default
is FALSE. If mapper_table is not NULL, this parameter will be
ignored.}

\item{mart}{If mapper_table is NULL, the mart to be used for fetching ENSEMBL IDs
and MGI symbols. Default is NULL. If mapper_table is not NULL, this
parameter will be ignored. Uses package biomaRt.}

\item{mapper_table}{The mapping table to use. Default is NULL; if a table is not
provided, then one will be created using the
\code{create_default_mapper_table} function.}
}
\value{
A vector with mapped values (or 'NA' in the corresponding slot if the
element in to_map is not in the mapping table).
}
\description{
Given a list of gene symbols or ENSEMBL IDs and a mapping table describing
the appropriate relationships between them, converts the entries of \code{to_map}
into their corresponding values.
}
\details{
Note that \code{to_map} must contain entries that fall within the keys of the
mapping table. These keys will then be turned into their corresponding values
within the mapping table. There may be duplicate keys or values in the mapping
table, but a mapping table containing a one-to-one mapping between keys and
values is recommended for compatibility with other preprocessing functions.

Note that this can be used outside the specific context of ENSEMBL IDs and MGI
symbols.

UPDATE: The GetBM function may already carry out these functions. The
functionality of the method may be replaced in the future.
}
