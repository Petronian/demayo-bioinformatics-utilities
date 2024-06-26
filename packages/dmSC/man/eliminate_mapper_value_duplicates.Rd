% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{eliminate_mapper_value_duplicates}
\alias{eliminate_mapper_value_duplicates}
\title{Eliminate duplicates in mapping table values}
\usage{
eliminate_mapper_value_duplicates(mapper_table)
}
\arguments{
\item{mapper_table}{The mapper_table whose values are to be de-duplicated.}
}
\value{
A list containing the following:
* A deduplicated mapping table based on the description above
('mapper_table'),
* A data.table containing the duplicated values that were removed
('duplicate_values').
}
\description{
Eliminates duplicate values in the 'values' (second column) of a mapper_table,
returning both the de-duplicated mapper_table and a data.frame consisting
of the removed duplicate values (and the number of occurrences for each
duplicate value).
}
\details{
At the end of the day, analogous to
\code{eliminate_mapper_key_duplicates(rev(mapper_table))}.
}
