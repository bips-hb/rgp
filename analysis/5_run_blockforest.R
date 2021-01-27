# to run:
# $ qsub -I -l mem=250gb,nodes=1:ppn=12 -N bfftgastro
# $ R
# > ade <- "gastro" or "icr"
# > source('/home/Gepard/pv_monitor/task1_3/task13_realdata/src/5_run_blockforest.R')

PRJ_DIR <- 

blocktypes <- c('ft', 'atcicd')
DATA_DIR <- 
RES_DIR <- "results"
do.this <- 'blockforest'
RV <- R.version$major

# devtools::install_local(file.path(PRJ_DIR, 'blockForest'))
pkgs <- c('data.table', 'batchtools', 'blockForest', 'precrec')
lapply(pkgs, require, character.only = TRUE)

# set seed
set.seed(1)

# get ttd pathway vector data
all_targets_l <- readRDS(file.path(PRJ_DIR, PKG_DIR, "data", "ttd_clean_data", "ttd_pathways_doac.Rds"))

#### read and prep data ####
# read dcasted cov data
subcohort <- readRDS(file.path(PRJ_DIR, DATA_DIR, paste0(ade, "_cc_cov.Rds")))

data <- prep_for_analytics(cohort = subcohort, keep_sparse_matrix_intercept = TRUE,
  nchar_gkz5 = 3, dummy_gkz5 = FALSE, sparse_matrix = FALSE, matrix = FALSE)

demographic_vars <- strsplit(as.character(data$model_formula), split = ' + ', fixed = TRUE)[[3]][-1]

# FTs
ft <- create_covars_groups(dt_colnames = colnames(data$x_train), demographic_vars = demographic_vars,
ftargetsl = all_targets_l, singletons_as_one_group = TRUE, create_groups_for_single_covariates = TRUE)

# ATC/ICD
atcicd <- create_covars_groups(dt_colnames = colnames(data$x_train), demographic_vars = demographic_vars,
ftargetsl = NULL, nchar_atc = 4, nchar_icd = 3, create_groups_for_single_covariates = TRUE, singletons_as_one_group = TRUE)

# recode for tuning
dat <- cbind(data$x_train, "CASE" = data$y_train)
dat[, CASE := as.factor(CASE)]
dat[, GKZ5 := as.factor(GKZ5)]

for (b in seq_along(blocktypes)) {
  # b <- 1 # TODO test
  blocktype <- blocktypes[b]
  data$blocks <- get(blocktype)

  # bf parameters
  mtry <- sapply(data$blocks, function(x) sqrt(length(x)))
  nsets <- 50  # 50 n tuning sets was 100x30minutes, stopped after 14 after 1week
  num.trees.pre <- 50 # 50 tuning
  block.method <- "BlockForest"
  splitrule <- "gini" # for binary is best
  save.memory <- FALSE # saves memory, no sorting
  probability <- TRUE # estimate probability instead of hard classification
  num.threads <- 12  # get entire memory of node qsub -I -l mem=90gb,nodes=1:ppn=12
  num.trees <- 500 # 500; 30min x 5
  classification <- NULL # instead of default regression, factor(response), df(predictors)
  # dependent.variable.name <- "CASE"
  importance <- "permutation"

  #### tuning separately, blockfor ####
  # run tuning in parallel using batchtools
  # blockfor does not work because we need the formula and/or it exceeds walltime
  source(file.path(PRJ_DIR, "task13_realdata", "src", "bf_tuning.R"))

  saveRDS(bf_tuning_res, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_%s_R%s.rds", do.this, "tuning", blocktype, ade, RV)))

  if (!(nrow(bf_tuning_res) == num.trees.pre)) {
    bf_tuning_res <- readRDS(file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_%s_R%s.rds", do.this, "tuning", blocktype, ade, RV)))
  }

  cvalues <- as.matrix(bf_tuning_res[which.min(err), c(-1, -2, -ncol(bf_tuning_res)), with = FALSE])[1, ]

  forest <- blockForest::blockForest(
    formula = data$model_formula,
    data = dat,
    num.trees = num.trees,
    mtry = mtry,
    blocks = data$blocks,
    block.weights = cvalues,
    keep.inbag = TRUE,
    block.method = block.method,
    splitrule = splitrule,
    classification = classification,
    save.memory = save.memory,
    probability = probability,
    importance = importance,
    scale.permutation.importance = TRUE,
    num.threads = num.threads)

  saveRDS(forest, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_%s_R%s.rds", do.this, "forest", blocktype, ade, RV)))

  # Roughly calculate block probabilities
  forest_prediction <- predict(object = forest, data = data$x_test) # matrix of predicted probabilies too
  prediction <- forest_prediction$predictions[, '1'] # vector
  prediction <- round(prediction)
  performance <- assess_results(prediction, data$y_test)
  saveRDS(performance, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_%s_R%s.rds", do.this, "performance", blocktype, ade, RV)))

  # A: variable importance
  # B: block importance
  names(cvalues) <- names(data$blocks)

  nzvars <- list(
    blocks = data.frame(cvalues),
    variables = data.frame(forest$variable.importance)
  )

  saveRDS(nzvars, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_%s_R%s.rds", do.this, "nzvars", blocktype, ade, RV)))
}
