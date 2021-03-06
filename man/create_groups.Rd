% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_create_groups_seperate_overlap_v2.R
\name{create_groups}
\alias{create_groups}
\title{Random Grouping of Predictors}
\usage{
create_groups(
  npred = 1000,
  nassoc = 10,
  nassoc_percent = 10,
  overlap = FALSE,
  overlap_size = 10
)
}
\arguments{
\item{npred}{p, number of prediction variables.}

\item{nassoc}{p_t, number of variables with true effect size = 0.}

\item{nassoc_percent}{an integer from 0 to 100, indicating percentage of variables in
a cluster with true effect size > 0.}

\item{overlap}{boolean, default: FALSE. It determines if the clusters should be
overlapping.}

\item{overlap_size}{an integer, default: 10. It determines the percentage of
variables to overlap across clusters.}
}
\value{
A list of
\itemize{
 \item \code{csd}: Cluster size data, as a data.table containing cluster,
 cluster size, group weights as in the \code{sqrt(cluster size)}, and
 number of variables with true effect size.
 \item \code{groups}: A named list of integers of length length equal to the
 number \code{npred} containing the cluster and named with the variable name.
}
}
\description{
A function to create a list of random variables clustered based
on a predefined cluster size. Currently, cluster sizes are based on
pathway size curated from KEGG database \code{load(pathway_dr_ds_key_value)},
and stored as vector of integers (outlier excluded
\code{cluster_size <- cluster_size[-1]}). In case of overlap, a variable
can be present in more than one group.
}
