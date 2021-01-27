# ------------------------------ Grouping --------------------------------------
#' @title Random Grouping of Predictors
#'
#' @description A function to create a list of random variables clustered based
#' on a predefined cluster size. Currently, cluster sizes are based on
#' pathway size curated from KEGG database \code{load(pathway_dr_ds_key_value)},
#' and stored as vector of integers (outlier excluded
#' \code{cluster_size <- cluster_size[-1]}). In case of overlap, a variable
#' can be present in more than one group.
#'
#' @param npred p, number of prediction variables.
#' @param nassoc p_t, number of variables with true effect size = 0.
#' @param nassoc_percent an integer from 0 to 100, indicating percentage of variables in
#' a cluster with true effect size > 0.
#' @param overlap boolean, default: FALSE. It determines if the clusters should be
#' overlapping.
#' @param overlap_size an integer, default: 10. It determines the percentage of
#' variables to overlap across clusters.
#
#' @return A list of
#' \itemize{
#'  \item \code{csd}: Cluster size data, as a data.table containing cluster,
#'  cluster size, group weights as in the \code{sqrt(cluster size)}, and
#'  number of variables with true effect size.
#'  \item \code{groups}: A named list of integers of length length equal to the
#'  number \code{npred} containing the cluster and named with the variable name.
#' }
#' @importFrom gtools mixedorder
#' @importFrom utils tail
#' @importFrom data.table data.table
#' @export
create_groups <- function(npred = 1000,
                          nassoc = 10,
                          nassoc_percent = 10,
                          overlap = FALSE,
                          overlap_size = 10) {
  csd <- get_cluster_sizes(target_clusters_size = npred)
  subsetted <- vector(mode = "list", length = length(csd$clname))
  predictor_index <- c(1:npred)
  names(predictor_index) <- paste0("V", 1:npred)
  n <- names(predictor_index)[-(1:nassoc)]
  na <- names(predictor_index)[1:nassoc]
  nassoc_nr <- round((nassoc * nassoc_percent) / 100.0)
  # Most important elements are 1 to nassoc (10 or 50, as supplied by the user)
  # Ensure the percentage of most important elements in a single
  # cluster ( subsetted[[i]] ) by nassoc_percent.
  # e.g., 75% of the truly associated variables are in one group

  # first: assign the nassoc into subsets
  i <- 1
  while(length(na) >= nassoc_nr) {
    if (csd$clsize[i] >= nassoc_nr) {
      subsetted[[i]] <- sample(na, size = nassoc_nr)
      na <- na[!(na %in% subsetted[[i]])]
    }
    i <- i + 1
  }

  # handle leftover na
  n <- c(n, na)
  na <- character(0)

  # assign the non-nassoc into subsets
  for (i in 1:length(csd$clname)) {
    non_nassoc_sample_size <- csd$clsize[i] - length(subsetted[[i]])
    if (non_nassoc_sample_size > 0 & length(n) >= non_nassoc_sample_size) {
      subsetted[[i]] <- c(subsetted[[i]],
                          sample(n, size = non_nassoc_sample_size))
      n <- n[!(n %in% subsetted[[i]])]
    }
  }
  # in case groups should overlap
  if (isTRUE(overlap)) {
    target_clusters_overlap_size <- ((overlap_size/100) * npred)
    csd_overlap <- get_cluster_sizes(target_clusters_size = target_clusters_overlap_size)
    # sample of sample of elements to overlap
    subset_for_overlap <- sample(x = names(predictor_index),
                                 size = target_clusters_overlap_size,
                                 replace = FALSE)
    # add leftover n
    # subset_for_overlap <- c(subset_for_overlap, n)
    # n <- character(0)
    # add overlapped elements
    subsetted_overlap <- vector(mode = "list", length = length(csd_overlap$clname))
    for (i in 1:length(csd_overlap$clname)) {
      ## working: one var cannot be in more than two groups
      # subsetted_overlap[[i]] <- sample(subset_for_overlap,
      #                                  size = csd_overlap$clsize[i],
      #                                  replace = FALSE)
      # subset_for_overlap <- subset_for_overlap[!(subset_for_overlap %in% subsetted_overlap[[i]])]

      ## novel:
      # subsetted_overlap[[i]] <- sample(subset_for_overlap,
      #                                  size = csd_overlap$clsize[i],
      #                                  replace = TRUE)
      # s <- subset_for_overlap[1:csd_overlap$clsize[i]]
      # subsetted_overlap[[i]] <- s
      # subset_for_overlap <- subset_for_overlap[-seq(1:csd_overlap$clsize[i])]

      ## without overlap set, sample directly from the npred_names
      subsetted_overlap[[i]] <- sample(names(predictor_index),
                                       size = csd_overlap$clsize[i],
                                       replace = FALSE)
    }
    subsetted <- c(subsetted, subsetted_overlap)
    csd <- rbind(csd, csd_overlap)
    csd$clname <- as.integer(rownames(csd))
  }
  # return
  # save as list for grpreg
  groups_index <- rep(csd$clname, csd$clsize)
  names(groups_index) <- unlist(subsetted)
  groups_index <- groups_index[gtools::mixedorder(names(groups_index))]

  # return number of assoc. vars in each cluster
  n_nassoc <- as.vector(groups_index[names(groups_index) %in% paste0("V", 1:nassoc)])

  csd$n_nassoc <- sapply(seq(1:dim(csd)[1]),
                         function(x) {
                           RGP:::count(x = n_nassoc, n = x)
                         })

  list(csd = csd,
       groups = subsetted,
       groups_index = groups_index)
}

#' @export
get_cluster_sizes <- function(target_clusters_size = 1000) {
  # Read sizes
  cluster_size <- RGP::cluster_size
  # Remove outlier gp with >900 elements
  cluster_size <- cluster_size[-1]
  # Draw randomly C clusters whose cluster sizes cs sum up to p
  sampled_cluster_sizes <- sample(x = cluster_size) # clusters
  sampled_cluster_sizes <- sampled_cluster_sizes[cumsum(sampled_cluster_sizes)
                                                 <= target_clusters_size]
  if (length(sampled_cluster_sizes) == 0) {
    sampled_cluster_sizes <- sample(cluster_size[cluster_size < target_clusters_size],
                                    size = 1)
  }
  # Draw randomly C clusters whose cluster sizes cs sum up to target_clusters_size
  repeat {
    # last element: tail(cumsum(sampled_cluster_sizes))[6] or
    if (tail(cumsum(sampled_cluster_sizes), n = 1) == target_clusters_size) {
      break
    }
    s <- sample(cluster_size, size = 1)
    if ((tail(cumsum(sampled_cluster_sizes), n = 1) + s) <= target_clusters_size)
      sampled_cluster_sizes <- c(sampled_cluster_sizes, s)
  }

  csd <- data.table("clname" = c(1:length(sampled_cluster_sizes)),
                    "clsize" = sampled_cluster_sizes,
                    "gw" = sqrt(sampled_cluster_sizes))
  return(csd)
}
