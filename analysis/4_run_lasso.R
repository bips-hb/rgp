# to run: qsub -I -l mem=90gb,nodes=2:ppn=12
PRJ_DIR <- 
PKG_DIR <-
DATA_DIR <- "pvm13_gepard_data"
RES_DIR <- "results"
ADE <- c("gastro", "icr")
do.this <- 'lasso'
RV <- R.version$major

pkgs <- c('data.table', 'Matrix', 'doParallel', 'glmnet', 'precrec')
lapply(pkgs, require, character.only = TRUE)

# set seed
set.seed(1)

# rsl <- vector(mode = 'list', length = length(ADE))
all_targets_l <- readRDS(file.path(PRJ_DIR, PKG_DIR, "data", "ttd_clean_data", "ttd_pathways_doac.Rds"))

for (i in seq_along(ADE)) {
  ade <- ADE[i]
  # read dcasted cov data
  subcohort <- readRDS(file.path(PRJ_DIR, DATA_DIR, paste0(ade, "_cc_cov.Rds")))

  data <- prep_for_analytics(cohort = subcohort,
    matrix = TRUE,
    sparse_matrix = TRUE,
    keep_sparse_matrix_intercept = FALSE,
    nchar_gkz5 = 3, dummy_gkz5 = TRUE)

  rs <- lasso_wrapper(data, best_lambda = 'lambda.min', type_measure = 'mse')

  saveRDS(rs$model, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_R%s.rds", do.this, "model", ade, RV)))
  saveRDS(rs$performance, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_R%s.rds", do.this, "performance", ade, RV)))

  # variable importance
  # get ttd pathway vector data
  data$blocks <- create_covars_groups(dt_colnames = colnames(data$x_train),
    ftargetsl = all_targets_l, singletons_as_one_group = TRUE)
  nzvars_blocks <- variable2block(variable_names = rs$nzvars$features,
    data_colnames = colnames(data$x_train),
    blocks = data$blocks)
  rs$nzvars$blocks <- sapply(seq_along(nzvars_blocks), function(x) {
    paste0(nzvars_blocks[[x]], collapse = ', ')
  })

  nzvars <- list(
    variables = rs$nzvars
  )

  saveRDS(nzvars, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_R%s.rds", do.this, "nzvars", ade, RV)))
}
