# ---------------------------- create_cohort -----------------------------------
#' @title Simulate Case-Control Data with Correlated Predictors
#'
#' @description This function uses \code{simstudy} package and creates a dataset of
#' a cohort based on given parameter.
#'
#' @param nobs Integer, number of observations.
#' @param npred Integer, number of predictor variables.
#' @param nassoc Integer, number of predictors that are truely
#'               associated with the outcome.
#' @param marginal_prob A numeric vector of size \code{npred} that determines the
#'                      marginal probabilities of the predictor variables.
#' @param true_beta numeric vector of true effects.
#' @param ccr Numeric, case-control ratio
#' @param corr_param either a correlation structure between predictor variables as numeric.
#'                   matrix or rho as numeric for \code{genCorGen()} (Default: NULL).
#' @param bal if TRUE, the case-to-control ratio is printed (Default: FALSE).
#'
#' @return A data.table object of size \code{nobs} x (\code{npred} + 2).
#'
#' @examples
#' library(simstudy)
#'
#' set.seed(4875)
#' nobs <- 1000
#' npred <- 100
#' ccr <- 0.5
#' nassoc <- 5
#' marginal_prob <- rep(0.15, npred)
#' true_beta <- rep(c(log(5), 0), c(nassoc, npred - nassoc))
#'
#' # -- simulate data
#' cohort <- create_cohort(nobs = nobs, npred = npred, nassoc = nassoc, ccr = ccr,
#'                         marginal_prob = marginal_prob, true_beta = true_beta)
#'
#' # -- test marginals and effects
#' apply(cohort, 2, mean)
#' summary(glm(out ~ ., data = cohort[,-1], family = "binomial", maxit = 100))$coefficients
#'
#'
#' # -- simulate correlated data (positive-definite)
#' corr_param <- outer(1:npred, 1:npred,
#'                      function(x,y) {(-1)^(1+abs(x-y)) * 0.95^abs(x-y)})
#'
#' # -- or provide rho
#' corr_param <- 0.1
#'
#' cohort <- create_cohort(nobs = nobs, npred = npred, nassoc = nassoc,
#'                         ccr = ccr, marginal_prob = marginal_prob,
#'                         true_beta = true_beta, corr_param = corr_param)
#'
#' @importFrom simstudy genCorGen
#' @seealso \code{\link[simstudy]{genCorGen}}.
#' @export
create_cohort <- function(nobs, npred, ccr, nassoc,
                          marginal_prob, true_beta,
                          corr_param = NULL, bal = FALSE) {

  if (!(is.null(corr_param))) {
    if (is.matrix(corr_param)) {
      dpred <- simstudy::genCorGen(n = nobs, nvars = npred, params1 = marginal_prob,
                                  dist = "binary", corMatrix = corr_param, wide = TRUE)
    }
    if (length(corr_param) == 1) {
      dpred <- simstudy::genCorGen(n = nobs, nvars = npred, params1 = marginal_prob,
                                   dist = "binary", rho = corr_param, corstr = "ar1",
                                   wide = TRUE)
    }
  } else {
    # NULL, non-correlation matrix
    dpred <- simstudy::genCorGen(n = nobs, nvars = npred, params1 = marginal_prob,
                                 dist = "binary", rho = 0, corstr = "cs", wide = TRUE)
  }

  # --- assign outcome based on best intercept
  dtx <- find_intercept_iter(true_beta, dtx_old = dpred, ccr = ccr)
  cohort <- dtx["data"][[1]]

  if (isTRUE(bal)) {
    attr(cohort, "bal") <- sprintf("Mean of simulated outcome: %f", mean(cohort$out))
  }
  cohort
}
