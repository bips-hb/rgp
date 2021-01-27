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

# ------------------------- add outcome columns --------------------------------
# read amb and diag
diag <- readRDS(file.path(PRJ_DIR, DATA_DIR, "diag.Rds"))
setkey(diag, IDNUM, ICD)
setindex(diag, DIAG_TYPE, DIAG_DAT)
# make a column that contains type of ADE if any
diag <- add_ade_col(diag, outcome_hosp_diag_type, outcome_icds_pattern)
# check ADE in case NN or NE
diag <- add_sec_ade_col(diag, outcome_hosp_diag_type, outcome_icds_pattern)

amb <- readRDS(file.path(PRJ_DIR, DATA_DIR, "amb.Rds"))
setkey(amb, IDNUM, ICD)
setindex(amb, DIAG_TYPE, DIAG_DAT)
# make a column that contains type of ADE if any
amb <- add_ade_col(amb, outcome_amb_diag_type, outcome_icds_pattern)
# check ADE in case not G
amb <- add_sec_ade_col(amb, outcome_amb_diag_type, outcome_icds_pattern)

# head(diag); dim(diag); key(diag); indices(diag); tail(diag); diag[, .N, by = ADE]; diag[, .N, by = ADE2]
# head(amb); dim(amb); key(amb); indices(amb); tail(amb); amb[, .N, by = ADE]; amb[, .N, by = ADE2]

# save as backup
saveRDS(amb, file.path(PRJ_DIR, DATA_DIR, "amb_ade.Rds"))
saveRDS(diag, file.path(PRJ_DIR, DATA_DIR, "diag_ade.Rds"))

rm(amb); gc()

#### workflow in a nutshell ####
# select rows of the outcome (ADE)
# get the index date
# get time-to-event since cohort_beg

#### patient and outcome information: cases + eligibility ####
pat <- readRDS(file.path(PRJ_DIR, DATA_DIR, "pat_iperiods_death.Rds"))

# get first occurence
disease_outcome <- diag[ADE != '']
disease_outcome[, INDEX_DAT := min(DIAG_DAT), by = .(IDNUM, ADE)]
disease_outcome <- unique(disease_outcome[, list(IDNUM, INDEX_DAT, ADE)])
setkey(disease_outcome, IDNUM, ADE)

# inclusion criteria (outcome >= study begin + 90)
disease_outcome[, INK := 0]
disease_outcome[INDEX_DAT >= (as.IDate(paste0(begin_yr, '-01-01')) + 90),
  INK := 1]

# calc time-to-event as number of days
disease_outcome[, TIME_TO_EVENT := numeric()]
disease_outcome[INK == 1, TIME_TO_EVENT :=
  INDEX_DAT - (as.IDate(paste0(begin_yr, '-01-01')))]

# n in the baseline and first quarter
# disease_outcome[, .N, by = .(INK, ADE)]

cases <- disease_outcome[INK == 1, ]

# diagnostics
# duplicated: (> 1  ADE)
# uniqueN(cases, by = c('IDNUM', 'ADE')) -
#   uniqueN(cases, by = c('IDNUM'))
# cases[IDNUM %in% cases[duplicated(IDNUM), IDNUM],][duplicated(IDNUM), IDNUM]

# 2 ADE on the same day?
# cases[IDNUM %in% cases[duplicated(cases[, c("IDNUM", "INDEX_DAT")]), IDNUM], ]

for (i in seq_along(outcome_icds_pattern)) {
  ade <- names(outcome_icds_pattern)[i]
  cases_ade <- disease_outcome[ADE == ade]
  pat <- merge(x = pat, y = cases_ade, by = "IDNUM", all.x = TRUE)
  pat[is.na(INK), `:=`(INK = 1, CASE = 0)]
  pat[is.na(CASE), CASE := ifelse(INK == 0, NA, INK)]
  case_define_cols <- c(names(cases_ade)[-1], "CASE")
  setnames(pat, case_define_cols, paste0(ade, ".", case_define_cols))
}

saveRDS(pat, file.path(PRJ_DIR, DATA_DIR, "cohort.Rds"))

# diagnostics
# pat[, .N, by = .(BLEED_GASTRO.INK, BLEED_GASTRO.CASE, BLEED_ICR.INK, BLEED_ICR.CASE)]
# pat[(!is.na(DEATH_DAT) & BLEED_GASTRO.INK == 1 & DEATH_DAT == BLEED_GASTRO.INDEX_DAT), ]
# pat[(!is.na(DEATH_DAT) & BLEED_ICR.INK == 1 & DEATH_DAT == BLEED_ICR.INDEX_DAT), ]
# pat[(!is.na(DEATH_DAT) & (BLEED_ICR.INK == 1 | BLEED_GASTRO.INK == 1) & (DEATH_DAT < BLEED_ICR.INDEX_DAT | DEATH_DAT < BLEED_GASTRO.INDEX_DAT)), ]
