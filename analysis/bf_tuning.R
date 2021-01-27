set.seed <- 42

# Registry ----------------------------------------------------------------
reg_name <- paste0("blockforest_tuning_", ade)
reg_dir <- file.path(PRJ_DIR, "registries", reg_name)
unlink(reg_dir, recursive = TRUE)
makeRegistry(file.dir = reg_dir, packages = c("blockForest"))

# Function to apply --------------------------------------------------------
myfun <- function(set, data, formula, num.trees.pre, blocks, block.method, splitrule, mtry, ...) {
  M <- length(blocks)
  cvalues <- sample(c(sort(runif(M-1)), 1))

  forest <- blockForest::blockForest(
                        # dependent.variable.name = dependent.variable.name,
                        data = data,
                        formula = formula,
                        num.trees = num.trees.pre,
                        blocks = blocks,
                        block.weights = cvalues,
                        mtry = mtry, keep.inbag = TRUE,
                        block.method = block.method,
                        splitrule = splitrule, write.forest = FALSE,
                        ...)

  c(cvalues = cvalues,
    err = forest$prediction.error)
}

# Create jobs ----------------------------------------------------------------
batchMap(myfun,
         set = seq(1, nsets),
         more.args = list(data = dat,
                          formula = data$model_formula,
                          # dependent.variable.name = dependent.variable.name,
                          num.trees.pre = num.trees.pre,
                          blocks = data$blocks,
                          block.method = block.method,
                          splitrule = splitrule,
                          mtry = mtry,
                          classification = classification,
													save.memory = save.memory,
                          probability = probability,
													num.threads = num.threads))

# Test run -----------------------------------------------------------
# testJob(id = 1)

# Submit -----------------------------------------------------------
ids <- findNotDone()
ids[, chunk := chunk(job.id, chunk.size = 200)]
submitJobs(ids = ids, # walltime in seconds, 10 days max, memory in MB
           resources = list(name = reg_name, chunks.as.arrayjobs = TRUE,
                            ncpus = 12, memory = 90000, walltime = 10*24*3600,
                            max.concurrent.jobs = 200))


waitForJobs()
# showLog(1)

# Get results -------------------------------------------------------------
bf_tuning_res <- ijoin(flatten(getJobPars()), flatten(reduceResultsDataTable()))

