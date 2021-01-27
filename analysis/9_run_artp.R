# to run:
# $ qsub -I -l mem=250gb,nodes=1:ppn=12 -N artpicrft
# $ R
ade <- "gastro" # or "icr"
blocktype <- "ft" # or "atcicd"
# > source('run_artp.R')

PRJ_DIR <- 

# test future apply on predict
DATA_DIR <- 
RES_DIR <- "results"
do.this <- "artp"
RV <- R.version$major
ADE <- c("icr", "gastro")
blocktypes <- c("ft", "atcicd")

pkgs <- c('data.table', 'precrec', 'future')
lapply(pkgs, require, character.only = TRUE)
devtools::load_all(file.path(PRJ_DIR, "ARTPredict"))

# set seed
set.seed(1)

# make a log file and write messages in it
logfname <- sprintf("%s_%s_%s", do.this, blocktype, ade)
logger <- file(get_log_file_name(file_name = logfname), open = "at")
sink(logger, type = "message")

# get ttd pathway vector data
all_targets_l <- readRDS(file.path(PRJ_DIR, PKG_DIR, "data", "ttd_clean_data", "ttd_pathways_doac.Rds"))

#### read and prep data ####
logit("Read and prep data")

subcohort <- readRDS(file.path(PRJ_DIR, DATA_DIR, paste0(ade, "_cc_cov.Rds")))

data <- prep_for_analytics(cohort = subcohort, keep_sparse_matrix_intercept = FALSE,
  nchar_gkz5 = 3, dummy_gkz5 = TRUE, sparse_matrix = FALSE, matrix = TRUE)

demographic_vars <- strsplit(as.character(data$model_formula), split = ' + ', fixed = TRUE)[[3]][-1]

# FTs
ft <- create_covars_groups(dt_colnames = colnames(data$x_train), demographic_vars = demographic_vars,
  ftargetsl = all_targets_l, singletons_as_one_group = FALSE, create_groups_for_single_covariates = FALSE)

# ATC/ICD
atcicd <- create_covars_groups(dt_colnames = colnames(data$x_train), demographic_vars = demographic_vars,
  ftargetsl = NULL, nchar_atc = 4, nchar_icd = 3, create_groups_for_single_covariates = FALSE, singletons_as_one_group = FALSE)

# data$y_train <- as.factor(data$y_train)
adjust_vars <- as.numeric(na.omit(match(demographic_vars, colnames(data$x_train))))

data$blocks <- get(blocktype)

options(future.globals.maxSize = 4000*1024^2)
logit("Start artp fit")
fit <- artp.fit(X = data$x_train, y = data$y_train, groups = data$blocks,
  trunc.point = 3,
  n.permutations = 50,
  adjust_vars = adjust_vars,
  single_covariates = TRUE,
  parallel = TRUE, nc = 11,
  verbose = FALSE)
logit("Finish artp fit")
saveRDS(fit, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_%s_R%s.rds", do.this, "model", blocktype, ade, RV)))

logit("Start artp predict")
prediction <- artp.predict(fit = fit, X.new = data$x_test, alpha = 0.05, predict.all = TRUE)
logit("Finish artp predict")
saveRDS(prediction, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_%s_R%s.rds", do.this, "prediction", blocktype, ade, RV)))

# fdr high, make alpha low, step-by-step down with alpha and look at the values (mariam diss), 0.05
performance <- assess_results(prediction$y.hat, data$y_test)
saveRDS(performance, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_%s_R%s2.rds", do.this, "performance", blocktype, ade, RV)))
sink()
