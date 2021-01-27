# ----------------- explain group/pathways argument
# list of indices artp, bf, ogl (artp makes single groups for single covars, bf and ogp)
# ARTP: List of groups. Each item is a vector with the indices of the covariates that belong to that group
# grepregpverlap: group Different from that in grpreg, group here must be a list of vectors, each containing integer indices or character names of variables in the group. variables that not belong to any groups will be disgarded.
# bf: A list of length equal to the number of groups. Each element of the list is a vector of the variable index in that group

# ----------------- functions --------------------------------------------------

#' @export
assess_results <- function(response, labels) {
  P <- sum(labels)
  N <- sum(!labels)
  TP <- sum(response & labels)
  TN <- sum(!response & !labels)
  FP <- sum(response & !labels)
  FN <- sum(!response & labels)

  # recall or tpr or sensitivity
  recall <- TP / (TP + FN)
  precision <- TP / (TP + FP)
  specificity <- TN / (TN + FP)
  fpr <- FP / (TN + FP)
  fdr <- FP / (TP + FP)
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  f1score <- 2 * ((precision * recall) / (precision + recall))
  # f1scoreb <- 2*TP / (2*TP + FP + FN)
  prc_baseline <- P / (P + N)
  aucs <- precrec::auc(curves =
    precrec::evalmod(scores = response, labels = labels)
  )[[4]]
  roc_auc <- aucs[1]
  prc_auc <- aucs[2]

  ## Generate an sscurve object that contains ROC and Precision-Recall curves
  ## already calculated AUC for the whole graph
  # or
  # auc_prc <- subset(auc(curves = sscurves), curvetypes == "PRC")[[4]]

  data.frame(
    list(
      recall = recall,
      precision = precision,
      specificity = specificity,
      fpr = fpr,
      fdr = fdr,
      accuracy = accuracy,
      f1score = f1score,
      roc_auc = roc_auc,
      prc_baseline = prc_baseline,
      prc_auc = prc_auc
    )
  )
}

variable2block <- function(variable_names, data_colnames, blocks) {
  blocks_names <- lapply(seq_along(blocks),
    function(x) {
      data_colnames[blocks[[x]]]
    })

  names(blocks_names) <- names(blocks)

  grps <- vector(mode = "list")

  for (v in seq_along(variable_names)) {
    grps[[v]] <- unlist(sapply(seq_along(blocks_names), function(k) {
      if (variable_names[v] %in% blocks_names[[k]]) {
        names(blocks_names[k])
      }
    }))
  }

  # length(unlist(blocks_names)) == length(unlist(grps))

  names(grps) <- variable_names
  return(grps)
}

assess_selected_vars <- function(variable_importance) {
  # get top vars

  # option 1: add the groups each variable is in, can be multiple groups, write how many groups
  variable_imp_blocks <- variable2block(variable_names = rownames(variable_imp),
    data_colnames = colnames(data$x_train),
    blocks = data$blocks)
  variable_imp$blocks <- sapply(seq_along(variable_imp_blocks), function(x) {
    paste0(variable_imp_blocks[[x]], collapse = ', ')
  })

  # option 2: rank
  num_top <- 20
  variable_imp <- sort(forest$variable.importance, decreasing = TRUE)
  top_var_imp <- variable_imp[1:num_top]
  top_var_imp <- data.frame(Rank = 1:num_top, explain_table(names(top_var_imp))[, c("code", "short_desc")])
  colnames(top_var_imp) <- c("Rank", "ICD code", "Description")
  # could be complemented with here: ~/bips_devel/task13_realdata/lib/ftsim/inst/extdata/icd10_2017/icd10gm2017syst_kodes.txt


  # Ranks for DOACs
  anticoag <- c("B01AF01", "B01AF02", "B01AF03", "B01AE07", "B01AA03", "B01AA04")
  names(anticoag) <- c("Riva", "Api", "Edox", "Dabi", "Warf", "Phenpro")
  rank(variable_imp)

  # option 3: add the % of variables selected (with hight importance)
}

assess_selected_blocks <- function(block_importance) {
  # rank/sort the blocks
}

prep_for_analytics <- function(cohort, keep_sparse_matrix_intercept = TRUE,
  nchar_gkz5 = 3, dummy_gkz5 = FALSE, matrix = TRUE, sparse_matrix = TRUE) {
  # rename (just to shorten it)
  try(setnames(cohort, "TIME_TO_EVENT", "TIME", skip_absent = TRUE))

  # sex as binary
  cohort[, SEX := as.numeric(SEX == "F")]

  # gkz5 as numeric (logical error), integer (trims first zero), factor
  cohort[, GKZ5 := as.factor(GKZ5)]

  # remove those died before index date
  # subcohort <- subcohort[!p[DEATH_DAT < INDEX_DAT, ], on = "IDNUM"]

  if (isTRUE(dummy_gkz5)) {
    # create GKZ5 dummy variables
    gkz5 <- dcast_dt(dt = cohort[, .(IDNUM, GKZ5 = paste0('GK', substr(GKZ5, 1, nchar_gkz5)))],
      dcols = 'IDNUM', dvar = 'GKZ5')

    # drop a high arbitrary observation
    cohort <- merge(cohort[, !c('GKZ5'), with = FALSE], gkz5[, !c('GK100'), with = FALSE])
  } else {
    # just trim GKZ5 code
    cohort[, GKZ5 := substr(GKZ5, 1, nchar_gkz5)]
    # produces logical error in case gkz5 is factor
    # recode_cohort <- cohort[, .SD, .(GKZ5 = substr(GKZ5, 1, nchar_gkz5))][, colnames(cohort), with = FALSE]
  }

  # recode train and test
  xytrain <- cohort[TRAIN == 1, !c("IDNUM", "TRAIN"), with = FALSE]

  # remove 0 variance from train and test
  no_variance_vars <- names(xytrain)[(xytrain[, lapply(.SD, var)] == 0)]

  xytrain <- xytrain[, !no_variance_vars, with = FALSE]
  xytest <- cohort[TRAIN == 0, !c("IDNUM", "TRAIN"), with = FALSE][, !no_variance_vars, with = FALSE]
  y_train <- cohort[TRAIN == 1, CASE]
  y_test <- cohort[TRAIN == 0, CASE]

  # dim(cohort)[2] - (length(no_variance_vars) + 2) == dim(xytrain)[2]

  # factorcols <- c("CASE", "SEX", names(xytrain)[grepl(names(xytrain), pattern = '^GK')])
  # xytrain[, (factorcols) := lapply(.SD, as.factor), .SDcols = factorcols]
  # xytest[, (factorcols) := lapply(.SD, as.factor), .SDcols = factorcols]
  # Be aware that GK... has number 1 appened to them because of the factorization
  # No difference in model performance between converting gk to factor and leaving them as integer

  # formula
  gkz5f <- paste0(names(xytrain)[grepl(names(xytrain), pattern = '^GK')], collapse = ' + ')
  # trainf <- as.formula(paste0(c("CASE ~ . -1", "SEX", "AGE", "TIME"), collapse = ' + ')) # remove intercept

   # keep or remove intercept, keep for lasso and remove for group-based?
  keep_sparse_matrix_intercept <- ifelse(isTRUE(keep_sparse_matrix_intercept), "CASE ~ .", "CASE ~ . -1")

  model_formula <- as.formula(paste0(c(keep_sparse_matrix_intercept, "SEX", "AGE", "TIME", gkz5f), collapse = ' + '))

  if (isTRUE(matrix)) {
    if (isTRUE(sparse_matrix)) {
      # create sparse matrix mode, the only way for glmnet to see both the confoundings, and to use sparse
      x_train <- sparse.model.matrix(model_formula, xytrain)
      x_test <- sparse.model.matrix(model_formula, xytest)

      # TODO: in case of factor variables (eg gkz5, default is to create dummies/contrasts)
      # x_train <- sparse.model.matrix(model_formula, xytrain[1:1000, ], contrasts.arg = list(GKZ5 = contrasts(xytrain$GKZ5, contrasts = FALSE)))
    } else {
      x_train <- model.matrix(model_formula, xytrain)
      x_test <- model.matrix(model_formula, xytest)
    }
  } else {
    x_train <- xytrain[, !c("CASE"), with = FALSE]
    x_test <- xytest[, !c("CASE"), with = FALSE]
    # model_formula <- as.formula(paste0(c("CASE ~ .", "SEX", "AGE", "TIME", gkz5f), collapse = ' + ')) # for bf and rf
  }

  data <- list(
    x_train = x_train,
    y_train = y_train,
    x_test = x_test,
    y_test = y_test,
    model_formula = model_formula
  )

  return(data)
}

# lasso
#' @export
lasso_wrapper <- function(data, best_lambda = "lambda.min", type_measure = "mse") {
  #### run cv fit ####
  if(!is.na(parallel::detectCores())) {
    parallel <- TRUE
    cl <- parallel::makePSOCKcluster(10)
    parallel::clusterExport(cl, c("glmnet"))
    doParallel::registerDoParallel(cl)
    # registerDoMC(cores = 8)
  } else {
    parallel <- FALSE
  }

  cvfit <- cv.glmnet(x = data$x_train, y = as.factor(data$y_train), alpha = 1,
                     family = "binomial", type.measure = type_measure,
                     intercept = FALSE, parallel = parallel)

  # fit <- glmnet(x = X_train, y = y_train, lambda = cvfit$lambda.min) # LD

  if(dim(showConnections())[[1]] > 0) {
    parallel::stopCluster(cl)
  }

  prediction <- predict(cvfit,
                        type = "response",
                        newx = data$x_test,
                        s = cvfit[[best_lambda]])

  # predict returns a matrix of predicted probabilities
  prediction <- round(prediction)

  myCoefs <- coef(cvfit, s = cvfit[[best_lambda]])

  ## Asseble into a data.frame
  nzvars <- data.frame(
    features = myCoefs@Dimnames[[1]][ which(myCoefs != 0 ) ],
    coefs    = myCoefs              [ which(myCoefs != 0 ) ]
  )

  performance <- assess_results(prediction, data$y_test)

  return(
    list(performance = performance,
    nzvars = nzvars,
    model = cvfit)
  )
}

# ------------------------ Naive group-regression ------------------------------
#' @export
naive_grpl_wrapper <- function(data,
				                       method = "sum") {

  xytrain <- cbind(get_group_value(X = data$x_train,
    groups = data$blocks, method = method), "CASE" = data$y_train)

  xytest <- cbind(get_group_value(X = data$x_test,
    groups = data$blocks, method = method), "CASE" = data$y_test)

  data$x_train <- sparse.model.matrix(data$model_formula, xytrain)
  data$x_test <- sparse.model.matrix(data$model_formula, xytest)

  rm(xytrain, xytest); gc()

	rs <- lasso_wrapper(data)

  rs
}

get_group_value <- function(X, groups, method = c("sum", "any"), parallel = TRUE, nc = 3, create_groups_for_single_covariates = TRUE) {
  method <- match.arg(method)
  nr.groups <- length(groups)
  n.cov <- ncol(X)
  # weights
  weights <- sqrt(sapply(groups, length))
  
  # calculate group values
  if (isTRUE(parallel)) {
    group.data <- mcmapply(function(k) {
      g <- groups[[k]]
      sr <- rowSums(subset(X, select = g))
      switch(method,
             "sum" =	sr,
             "any" = sign(sr))
    }, mc.cores = 3, k = 1:nr.groups, SIMPLIFY = TRUE)
    
  } else {
    group.data <- matrix(nrow = nrow(X), ncol = nr.groups)
    
    for (k in 1:nr.groups) {
      g <- groups[[k]]
      sr <- rowSums(subset(X, select = g))
      group.data[, k] <- switch(method,
                                "sum" =	sr,
                                "any" = sign(sr))
    }
  }
  
  # transform group value to 1/sqrt(Kj)
  group.data.weights <- matrix(nrow = nrow(X), ncol = nr.groups)
  for (k in 1:nr.groups) {
    group.data.weights[, k] <- group.data[, k] * (1/weights[k])
  }
  
  if (isTRUE(create_groups_for_single_covariates)) {
    gnames <- names(groups)
    variables_in_groups <- do.call(c, groups)
    variables_not_in_groups <- setdiff(1:n.cov, variables_in_groups)
    
    # add covariates that are not in groups to group list
    groups <- c(groups, as.list(variables_not_in_groups))
    names(groups) <- c(gnames, colnames(X)[variables_not_in_groups])
    # names(groups) <- c(gnames, paste0("grp", nr.groups+1:length(variables_not_in_groups)))
    
    # add covariates that are not in groups to group.data
    group.data.weights <- cbind(group.data.weights, subset(X, select = variables_not_in_groups))
  }
  
  colnames(group.data.weights) <- names(groups)
  
  group.data.weights
}
