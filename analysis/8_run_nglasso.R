# to run: qsub -I -l mem=90gb,nodes=2:ppn=12 -N ngl

PRJ_DIR <- 

DATA_DIR <- 
RES_DIR <- "results"
ADE <- c("gastro", "icr")
blocktypes <- c('ft', 'atcicd')
do.this <- 'naivegrplasso'
RV <- R.version$major

pkgs <- c('data.table', 'Matrix', 'parallel', 'glmnet', 'precrec')
lapply(pkgs, require, character.only = TRUE)
# set seed
set.seed(1)

all_targets_l <- readRDS(file.path(PRJ_DIR, PKG_DIR, "data", "ttd_clean_data", "ttd_pathways_doac.Rds"))

for (i in seq_along(ADE)) {
  ade <- ADE[i]

  # read dcasted cov data
  subcohort <- readRDS(file.path(PRJ_DIR, DATA_DIR, paste0(ade, "_cc_cov.Rds")))

  data <- prep_for_analytics(cohort = subcohort,
    matrix = FALSE,
    sparse_matrix = FALSE,
    keep_sparse_matrix_intercept = FALSE,
    nchar_gkz5 = 3, dummy_gkz5 = TRUE)

  demographic_vars <- strsplit(as.character(data$model_formula), split = ' + ', fixed = TRUE)[[3]][-1]

  # FTs
  ft <- create_covars_groups(dt_colnames = colnames(data$x_train), demographic_vars = demographic_vars,
  ftargetsl = all_targets_l, singletons_as_one_group = FALSE, create_groups_for_single_covariates = TRUE)

  # ATC/ICD
  atcicd <- create_covars_groups(dt_colnames = colnames(data$x_train), demographic_vars = demographic_vars,
  ftargetsl = NULL, nchar_atc = 4, nchar_icd = 3, create_groups_for_single_covariates = TRUE, singletons_as_one_group = FALSE)

  for (b in seq_along(blocktypes)) {
    blocktype <- blocktypes[b]

    data$blocks <- get(blocktype)

    rs <- naive_grpl_wrapper(data, method = "sum")

    saveRDS(rs$model, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_%s_R%s.rds", do.this, "model", blocktype, ade, RV)))
    saveRDS(rs$performance, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_%s_R%s.rds", do.this, "performance", blocktype, ade, RV)))

    nzvars <- list(
      variables = rs$nzvars
    )

    saveRDS(nzvars, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_%s_R%s.rds", do.this, "nzvars", blocktype, ade, RV)))
  }
}
