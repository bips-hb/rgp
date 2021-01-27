################################################################################
#
# This is a run file for getting the CASE information. It reads data from RDS
# and does not use gepardr package.
#
# Author: Mariam R. Rizkallah
# Date: 25.07.2019
# Run on a computing node: qsub -I -l mem=250gb,nodes=1:ppn=12
#
################################################################################

# Define project directory
PRJ_DIR <- 
PKG_DIR <-
  
  # Source config file (study variables file)
source(file.path(PRJ_DIR, PKG_DIR, "analysis", "1_define_study_variables.R"))
source(file.path(PRJ_DIR, PKG_DIR, "R", "cohort_helpers.R"))

# ADE file
ADE_FILE <- file.path(PRJ_DIR, PKG_DIR, "inst", "extdata", "Paper13_events.csv")

# outcome_vars
ade_icd_list <- get_ade_icd_list(ADE_FILE)
outcome_icds_pattern <- lapply(X = ade_icd_list,
                               FUN = function(x) paste(x, collapse = '|'))
outcome_hosp_diag_type <- c('HH', 'NH') # there is "NN", "NE" more
outcome_amb_diag_type <- 'G'

#### matching ####
cohort <- readRDS(file.path(PRJ_DIR, DATA_DIR, "cohort.Rds"))
setindex(cohort, BIRTHY)
cohort[, `:=` (BLEED_ICR.PAIR = 0, BLEED_GASTRO.PAIR = 0)]

icr_cohort <- cohort[BLEED_ICR.INK == 1, ]
names(icr_cohort) <- gsub(names(icr_cohort), pattern = 'BLEED_ICR.', replacement = '')
gastro_cohort <- cohort[BLEED_GASTRO.INK == 1, ]
names(gastro_cohort) <- gsub(names(gastro_cohort), pattern = 'BLEED_GASTRO.', replacement = '')
rm(cohort); gc()

ncores <- 8
set.seed(1)
cl <- parallel::makePSOCKcluster(ncores)
parallel::clusterExport(cl, c("icr_cohort", "create_matched_cc_cohort"))
doParallel::registerDoParallel(cl)
# registerDoMC(cores = 8)
# cl <- makeCluster(ncores)
# registerDoParallel(cl)
icr_cc <- create_matched_cc_cohort(icr_cohort, first_split_var = "BIRTHY",
  nmatches = 4, ttflag = NULL)
stopCluster(cl)
icr_cc[, AGE := (year(INDEX_DAT) - BIRTHY)]
icr_cc[, TRAIN := ifelse(as.numeric(RID) %% 2 == 0, 1, 0)]
setkey(icr_cc, IDNUM, RID)

# unique(table(icr_cc$PAIR)) == 5 # check number of matched pairs per case
# icr_cc[is.na(IDNUM), ]
# icr_cc[, .N, .(CASE, TRAIN)]
# icr_cc[duplicated(IDNUM), ]

saveRDS(icr_cc, file.path(PRJ_DIR, DATA_DIR, "icr_cc.Rds"))

cl <- parallel::makePSOCKcluster(ncores)
parallel::clusterExport(cl, c("gastro_cohort", "create_matched_cc_cohort"))
doParallel::registerDoParallel(cl)
gastro_cc <- create_matched_cc_cohort(gastro_cohort, first_split_var = "BIRTHY",
  nmatches = 4, ttflag = NULL)
stopCluster(cl)
gastro_cc[, AGE := (year(INDEX_DAT) - BIRTHY)]
gastro_cc[, TRAIN := ifelse(as.numeric(RID) %% 2 == 0, 1, 0)]
setkey(gastro_cc, IDNUM, RID)

unique(table(gastro_cc$PAIR)) == 5 # check number of matched pairs per case
gastro_cc[is.na(IDNUM), ]
gastro_cc[, .N, .(CASE, TRAIN)]
gastro_cc[duplicated(IDNUM), ]

saveRDS(gastro_cc, file.path(PRJ_DIR, DATA_DIR, "gastro_cc.Rds"))

#### disease and drug predictors ####
amb <- readRDS(file.path(PRJ_DIR, DATA_DIR, "amb_ade.Rds"))
diag <- readRDS(file.path(PRJ_DIR, DATA_DIR, "diag_ade.Rds"))
disp <- readRDS(file.path(PRJ_DIR, DATA_DIR, "disp.Rds"))

gastro_cc_ds <- get_ds_predictors(subcohort = gastro_cc, diag = diag, amb = amb, begin_yr = begin_yr)
gastro_cc_dr <- get_dr_predictors(subcohort = gastro_cc, disp = disp)
saveRDS(gastro_cc_ds, file.path(PRJ_DIR, DATA_DIR, "gastro_cc_ds.Rds"))
saveRDS(gastro_cc_dr, file.path(PRJ_DIR, DATA_DIR, "gastro_cc_dr.Rds"))

icr_cc_ds <- get_ds_predictors(subcohort = icr_cc, diag = diag, amb = amb, begin_yr = begin_yr)
icr_cc_dr <- get_dr_predictors(subcohort = icr_cc, disp = disp)
saveRDS(icr_cc_ds, file.path(PRJ_DIR, DATA_DIR, "icr_cc_ds.Rds"))
saveRDS(icr_cc_dr, file.path(PRJ_DIR, DATA_DIR, "icr_cc_dr.Rds"))

#### casting into a binary matrix ####
dcols <- c("IDNUM", "SEX", "AGE", "GKZ5", "TIME_TO_EVENT", "CASE", "TRAIN")
dvar <- 'COV'

icr_cc_cov <- create_cov_mat(subcohort = NULL,
  disease_covar = icr_cc_ds,
  drug_covar = icr_cc_dr,
  dcols, dvar,
  nchar_icd = 4)

saveRDS(icr_cc_cov, file.path(PRJ_DIR, DATA_DIR, "icr_cc_cov.Rds"))

gastro_cc_cov <- create_cov_mat(subcohort = NULL,
  disease_covar = gastro_cc_ds,
  drug_covar = gastro_cc_dr,
  dcols, dvar,
  nchar_icd = 4)

saveRDS(gastro_cc_cov, file.path(PRJ_DIR, DATA_DIR, "gastro_cc_cov.Rds"))

# x <- gastro_cc_cov
# has_more = apply(x[, !c("IDNUM", "SEX", "AGE", "CASE", "GKZ5", "TIME_TO_EVENT"),
# with = FALSE] > 1, 1, any); any(isTRUE(has_more))
# dim(x); x[1:10, 1:10]; sum(colSums(is.na(x)))
# x[rowSums(x[, !c("IDNUM", "SEX", "AGE", "CASE", "GKZ5", "TIME_TO_EVENT", "TRAIN"), with = FALSE]) == 0]
