# to run:
# $ qsub -I -l mem=250gb,nodes=1:ppn=12 -N rfgastro
# $ R
# > ade <- "gastro" or "icr"
# source()

PRJ_DIR <- 

DATA_DIR <- 
RES_DIR <- "results"
do.this <- 'randomforest'
RV <- R.version$major

pkgs <- c('data.table', 'ranger', 'precrec')
lapply(pkgs, require, character.only = TRUE)


# set seed
set.seed(1)

#### read and prep data ####
  # read dcasted cov data
subcohort <- readRDS(file.path(PRJ_DIR, DATA_DIR, paste0(ade, "_cc_cov.Rds")))

data <- prep_for_analytics(cohort = subcohort, keep_sparse_matrix_intercept = TRUE,
  nchar_gkz5 = 3, dummy_gkz5 = FALSE, sparse_matrix = FALSE, matrix = FALSE)

# rf parameters
mtry <- (ncol(data$x_train) / 2)
splitrule <- "gini" # for binary is best
save.memory <- FALSE # saves memory, no sorting
probability <- TRUE # estimate probability instead of hard classification
num.threads <- 12  # get entire memory of node qsub -I -l mem=90gb,nodes=1:ppn=12
num.trees <- 500 # 500; 30min x 5
classification <- NULL # instead of default regression, factor(response), df(predictors)
# dependent.variable.name <- "CASE"
importance <- "permutation"

# recode for tuning
dat <- cbind(data$x_train, "CASE" = data$y_train)
dat[, CASE := as.factor(CASE)]
dat[, GKZ5 := as.factor(GKZ5)]

forest <- ranger::ranger(
  formula = data$model_formula,
  data = dat,
  num.trees = num.trees,
  mtry = mtry,
  keep.inbag = TRUE,
  splitrule = splitrule,
  importance = importance,
  scale.permutation.importance = TRUE,
  classification = classification,
  save.memory = save.memory,
  probability = probability,
  num.threads = num.threads)

saveRDS(forest, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_R%s.rds", do.this, "forest", ade, RV)))

forest_prediction <- predict(object = forest, data = data$x_test) # matrix of predicted probabilies too

# Roughly calculate block probabilities
prediction <- forest_prediction$predictions[, '1'] # vector
prediction <- round(prediction)
performance <- assess_results(prediction, data$y_test)
saveRDS(performance, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_R%s.rds", do.this, "performance", ade, RV)))

# variable importance
nzvars <- list(
  variables = data.frame(forest$variable.importance)
)

saveRDS(nzvars, file.path(PRJ_DIR, RES_DIR, sprintf("%s_%s_%s_R%s.rds", do.this, "nzvars", ade, RV)))
