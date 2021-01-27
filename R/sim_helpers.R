# ---------------------------- find_intercept_iter -----------------------------
#' @title Find the Intercept for a Dataset
#'
#' @description This function approximates intercept (\eqn{\beta_{0}})
#' iteratively to get a balanced case-control ratio in the simulated data.
#' This function uses \code{simstudy} package.
#'
#' @usage find_intercept_iter(true_beta, dtx_old, ccr = 0.5, ...)
#'
#' @param true_beta A vector of true beta values including both assciated and
#' non-associated variables.
#' @param dtx_old Data.table of old dataset to be appended with outcome.
#' @param ccr Numeric, case-control ratio (Default: 0.5).
#' @param cepts A numeric list of intercepts to use, recommmended to be NULL
#' (Default: NULL).
#' @param iter Integer, number of increments between intercept limits (Default: 500).
#'
#' @return A list of
#' \itemize{
#'  \item \code{data}: The new data.table with outcome column \code{data$out}.
#'  \item \code{intercept}: An estimation of the intercept.
#'  \item \code{difference}: The absolute difference between \code{ccr}
#' and \code{mean(data$out)}.
#'  \item \code{sequence}: Sequence of tested intercepts.
#'  \item \code{ccr}: Case-to-control ratio for each intercept.
#' }
#'
#' @examples
#' dtx <- simstudy::genCorGen(n = 100, nvars = 5, params1 = 0.15,
#'  dist = "binary", wide = TRUE, rho = 0.1, corstr = "ar1")
#' find_intercept_iter(c(-2, -2, 2, 2, 2), dtx)
#'
#' @importFrom utils tail
#' @seealso \code{\link[simstudy:defDataAdd]{defDataAdd}} and
#' \code{\link[simstudy:addColumns]{addColumns}} in \code{simstudy} package.
#'
#' @export
find_intercept_iter <- function(true_beta, dtx_old, ccr = 0.5, cepts = NULL,
                                iter = 500) {
  if (is.null(cepts)) {
    tmp <- - round(sum(true_beta))
    cepts <- c(round(tmp - length(true_beta)/10),
      round(tmp + length(true_beta)/10))
  }
  cepts <- seq(cepts[1], cepts[2], length = iter)
  cacor <- vector()
  daten <- list()

  i <- 0
  for (intercept in cepts) {
    esize <- paste(true_beta, "*", names(dtx_old)[2:(length(true_beta) + 1)],
      sep = "", collapse = " + ")
    formel <- paste0(intercept, " + ", esize)

    def <- simstudy::defDataAdd(
      varname = "out", dist = "binary",
      formula = formel, link = "logit"
    )
    dtx_new <- simstudy::addColumns(def, dtx_old)
    i <- i + 1
    cacor[i] <- mean(dtx_new$out)
    daten[[i]] <- dtx_new
  }

  diff <- abs(cacor - ccr)
  out <- list(
    data = daten[[min(which(diff == min(diff)))]],
    intercept = cepts[which(diff == min(diff))],
    difference = diff[which(diff == min(diff))],
    sequence = cepts,
    ccr = cacor
  )
  out
}

# ------------------------------ Grouping --------------------------------------
#' @title Random Grouping of Predictors
#'
#' @description A function to create a list of random variables clustered based
#' on a predefined cluster size. Currently, cluster sizes are based on
#' pathway size curated from KEGG database and stored as vector of integers
#' (outlier excluded \code{cluster_size <- cluster_size[-1]}).
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
#'  cluster size, and group weights as in the \code{sqrt(cluster size)}.
#'  \item \code{groups}: A named list of integers of length length equal to the
#'  number \code{npred} containing the cluster and named with the variable name.
#' }
#' @importFrom gtools mixedorder
#' @export
create_groups_leg <- function(npred = 1000,
                          nassoc = 10,
                          nassoc_percent = 10,
                          overlap = FALSE,
                          overlap_size = 10) {

  # subset sizes
  # load("rgp/R/cluster_size.rda")
  cluster_size <- RGP::cluster_size
  # Remove outlier gp with >900 elements
  cluster_size <- cluster_size[-1]
  # Draw randomly C clusters whose cluster sizes cs sum up to p
  sampled_cluster_sizes <- sample(x = cluster_size) # clusters
  if (isTRUE(overlap)) {
    target_clusters_size <- npred + ((overlap_size/100) * npred)
  } else {
    target_clusters_size <- npred
  }
  sampled_cluster_sizes <- sampled_cluster_sizes[cumsum(sampled_cluster_sizes)
                                                 <= target_clusters_size]

  # Draw randomly C clusters whose cluster sizes cs sum up to npred
  repeat {
    if (tail(cumsum(sampled_cluster_sizes))[6] == target_clusters_size) {
      break
    }
    s <- sample(cluster_size, size = 1)
    if ((tail(cumsum(sampled_cluster_sizes))[6] + s) <= target_clusters_size)
      sampled_cluster_sizes <- c(sampled_cluster_sizes, s)
  }

  csd <- data.table::data.table("clname" = c(1:length(sampled_cluster_sizes)),
    "clsize" = sampled_cluster_sizes,
    "gw" = sqrt(sampled_cluster_sizes))

  predictor_names <- c(1:npred) # all same group called 1 or all indep. called 1:1000
  names(predictor_names) <- paste0("V", 1:npred)

  # Ensure that the percentage of most important elements in a
  # cluster ( subsetted[[i]] ) does not exceed nassoc_percent.
  # Most important elements are 1 to nassoc (10 or 50, as supplied by the user)
  #

  subsetted <- list()
  n <- names(predictor_names)

  if (isTRUE(overlap)) {
    # a sample of elements to overlap
    subset_for_overlap <- sample(n, size = (target_clusters_size - npred))
    # number of overlap elements per cluster
    gap_in_a_cluster <- ceiling((target_clusters_size - npred)/length(csd$clsize))

    for (i in 1:length(csd$clname)) {
      gap <- subset_for_overlap[1:gap_in_a_cluster]

      if (csd$clsize[i] < length(gap) | any(is.na(gap))) {
        # cluster size is smaller than the gap
        subsetted[[i]] <- sample(n, size = csd$clsize[i])
        n <- n[!(n %in% subsetted[[i]])]
      } else {
        # cluster size is >= gap
        subsetted[[i]] <- sample(n, size = (csd$clsize[i] - length(gap)))
        # remove the already sampled names
        n <- n[!(n %in% subsetted[[i]])]
        # add the gap
        subsetted[[i]] <- c(subsetted[[i]], gap)
        # remove the gap elements from the original list that were added to subsetted
        subset_for_overlap <- subset_for_overlap[-c(1:gap_in_a_cluster)]
      }
      print(subsetted[[i]])
    }
  } else {
    for (i in 1:length(csd$clname)) {
      clsize = csd$clsize[i]

      subsetted[[i]] <- sample(n, size = clsize)

      # new code
      clsize_nassoc_percent <- round((clsize * nassoc_percent) / 100.0)

      nassoc_elements <- subsetted[[i]][subsetted[[i]] %in% paste0("V", 1:nassoc)]

      # print("nassoc elements:")
      # print(nassoc_elements)
      if (length(nassoc_elements) > clsize_nassoc_percent) {
        # print("")
        remove_elements_count = length(nassoc_elements) - clsize_nassoc_percent

        # print(sprintf("remove elements count: %d", remove_elements_count))
        remove_elements <- sample(nassoc_elements, remove_elements_count)

        # print(sprintf("elements to remove: %s", remove_elements))
        # print(sprintf("subsetted before removing: %s",
        #               paste0(subsetted[[i]], collapse = ", ")))
        # remove random elements which are in nassoc_elements
        subsetted[[i]] <- subsetted[[i]][!(subsetted[[i]] %in% remove_elements)]
        # print(sprintf("subsetted after removing: %s",
        #               paste0(subsetted[[i]], collapse = ", ")))

        # concatenate new random elements selected from n
        nn <- n[!(n %in% paste0("V", 1:nassoc))]
        subsetted[[i]] <- c(subsetted[[i]],
                            sample(nn, size = remove_elements_count))

        # print(sprintf("subsetted after adding new elements: %s",
        #               paste0(subsetted[[i]], collapse = ", ")))
      }

      n <- n[!(n %in% subsetted[[i]])]
    }
  }
  # save as list for grpreg
  groups_index <- rep(csd$clname, csd$clsize)
  names(groups_index) <- unlist(subsetted)
  groups_index <- groups_index[mixedorder(names(groups_index))]

  list(csd = csd,
       groups = subsetted,
       groups_index = groups_index)
}
