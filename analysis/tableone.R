library(data.table)
library(arsenal)
library(dplyr)

PRJ_DIR <- 

DATA_DIR <- 

#### read and prepare data ####
gastro_cc_dr <- readRDS(file.path(PRJ_DIR, DATA_DIR, "gastro_cc_dr.Rds"))
gastro_cc_ds <- readRDS(file.path(PRJ_DIR, DATA_DIR, "gastro_cc_ds.Rds"))
icr_cc_dr <- readRDS(file.path(PRJ_DIR, DATA_DIR, "icr_cc_dr.Rds"))
icr_cc_ds <- readRDS(file.path(PRJ_DIR, DATA_DIR, "icr_cc_ds.Rds"))

# add ADE
gastro <- rbind(gastro_cc_dr, gastro_cc_ds)
gastro[, SUBCOHORT := 'gastro']
icr <- rbind(icr_cc_dr, icr_cc_ds)
icr[, SUBCOHORT := 'icr']

# bind
p <- rbind(gastro, icr)

intcols <- c("AGE", "CASE", "TRAIN")
p[, (intcols) := lapply(.SD, as.integer), .SDcols = intcols]

# set keys
setkey(p, IDNUM)

p <- unique(p[, COV := NULL])

#### Descritpive statistics ####
ca <- factor(p$SUBCOHORT)

t1 <- tableby(data = p,
  formula = CASE ~ SEX + TIME_TO_EVENT + AGE,
  strata = ca,
  test = FALSE)
t1s <- summary(t1, text = TRUE,
  digits = 0, digits.p = 1, digits.pct = 3)
t1s

# for totals
p[, .N, by = c('CASE', 'SUBCOHORT')]
# write.csv(t1s, '/my/path/here/my_table.csv')

# add death data
d <- readRDS(file.path(PRJ_DIR, DATA_DIR, "pat_iperiods_death.Rds"))
p <- d[!is.na(DEATH_DAT), list(IDNUM, DEATH_DAT)][p, , on = 'IDNUM']
begin_yr <- 2015
p[, D_EVENT := ifelse(is.na(DEATH_DAT), 'ALIVE', 'DEAD')]
p[, D_TIME := ifelse(!is.na(DEATH_DAT), DEATH_DAT - (as.IDate(paste0(begin_yr, '-01-01'))), NA)]

# roughly convery TIME_TO_EVENT back to 'INDEX_DAT'
p[, INDEX_DAT := (as.IDate(paste0(begin_yr, '-01-01')) + TIME_TO_EVENT)]

# span between index date and death date
p[, D_I_SPAN := ifelse(!is.na(DEATH_DAT), (DEATH_DAT - INDEX_DAT), NA)]

# death date < index date? D_I_SPAN should be -
p[!is.na(DEATH_DAT), D_EVENT_BEFORE := ifelse(D_I_SPAN < 0, 'YES', 'NO')]
p[!is.na(DEATH_DAT), D_EVENT_ON := ifelse(DEATH_DAT == INDEX_DAT, 'YES', 'NO')]

t3 <- tableby(
  data = p,
  formula = CASE ~ D_EVENT + D_TIME + D_I_SPAN + D_EVENT_BEFORE + D_EVENT_ON,
  strata = ca,
  test = FALSE
)
summary(t3, text = TRUE,
  digits = 0, digits.p = 1, digits.pct = 3)


#### Covariables ####
rm(icr_cc_dr, icr_cc_ds, gastro_cc_dr, gastro_cc_ds, gastro, icr, p, d); gc()

#### read ####
ades <- c('gastro', 'icr')
pld <- lapply(ades, function(x) {readRDS(file.path(PRJ_DIR, DATA_DIR, paste0(x, "_cc_cov.Rds")))})
names(pld) <- toupper(ades)
dcols <- c("IDNUM", "SEX", "AGE", "GKZ5", "TIME_TO_EVENT", "CASE", "TRAIN")

#### describe cov matrix ####
ncol <- (length(ades) * 3)
t4 <- data.frame(matrix(ncol = ncol, nrow = 21))
colnames(t4) <- c("gib_ca", "gib_ctrl", "gib_tot",
                  "icb_ca", "icb_ctrl", "icb_tot"
)
rownames(t4) <- c("cov_tot", "cov_m", "cov_sd", "cov_min", "cov_max", "cov_zn", "cov_zp",
                  "dr_tot", "dr_m", "dr_sd", "dr_min", "dr_max", "dr_zn", "dr_zp",
                  "ds_tot", "ds_m", "ds_sd", "ds_min", "ds_max", "ds_zn", "ds_zp"
)

seg <- c(1:3)

for (i in 1:length(pld)) {
  s <- pld[[i]]

  # put prefix for col names
  try({
    setnames(s, colnames(s)[(length(dcols) + 1):which(colnames(s) == 'V90N')],
             paste0('dr', colnames(s)[(length(dcols) + 1):which(colnames(s) == 'V90N')]))

    setnames(s, colnames(s)[(which(colnames(s) == 'drV90N')+1):ncol(s)],
             paste0('ds', colnames(s)[(which(colnames(s) == 'drV90N')+1):ncol(s)]))
  })

  # n. covars
  start <- length(dcols) + 1
  end <- ncol(s)
  t4[1:7, seg] <- describe_cov_matrix(s, start = start, end = end)

  # n. drugs
  start <- length(dcols) + 1
  end <- which(colnames(s) == 'drV90N')
  t4[8:14, seg] <- describe_cov_matrix(s, start = start, end = end)

  # n. diseases
  start <- (which(colnames(s) == 'drV90N') + 1)
  end <- ncol(s)
  t4[15:21, seg] <- describe_cov_matrix(s, start = start, end = end)

  # shift by 3
  seg <- seg + 3
}

# to save
t4n <- data.frame(t(t4)) %>%
          summarise("N. covariables" = '',
                    "Total" = cov_tot,
                    "mean (SD)" = paste0(cov_m, ' (', cov_sd, ')'),
                    "range" = paste0('(', cov_min, ' - ', cov_max, ')'),
                    "N. non-informative (%)" = paste0(cov_zn, ' (', cov_zp, ')'),

                    "N. drugs" = '',
                    "DR Total" = dr_tot,
                    "DR mean (SD)" = paste0(dr_m, ' (', dr_sd, ')'),
                    "DR range" = paste0('(', dr_min, ' - ', dr_max, ')'),
                    "DR N. non-informative (%)" = paste0(dr_zn, ' (', dr_zp, ')'),

                    "N. comorbidities" = '',
                    "DS Total" = ds_tot,
                    "DS mean (SD)" = paste0(ds_m, ' (', ds_sd, ')'),
                    "DS range" = paste0('(', ds_min, ' - ', ds_max, ')'),
                    "DS N. non-informative (%)" = paste0(ds_zn, ' (', ds_zp, ')')
          ) %>%
          t() %>%
          data.frame()
names(t4n) <- colnames(t4)
write.table(t4n, file.path("reports", "t4.tsv"),
            sep = '\t', quote = FALSE)

#### describe grouping ####
ftargetsl <- readRDS(file.path(PRJ_DIR, PKG_DIR, "data", "ttd_clean_data", "ttd_pathways_doac.Rds"))
pld <- lapply(ades, function(x) {readRDS(file.path(PRJ_DIR, DATA_DIR, paste0(x, "_cc_cov.Rds")))})

t5 <- data.frame(matrix(ncol = ncol, nrow = 9))
colnames(t5) <- c("gib_ca", "gib_ctrl", "gib_tot",
                  "icb_ca", "icb_ctrl", "icb_tot"
)
rownames(t5) <- c("g_tot", "g_totp", "g_n", "g_m", "g_sd", "g_min", "g_max", "s_tot", "s_totp")

seg <- c(1:3)

for (i in 1:length(pld)) {
  s <- pld[[i]]

  # for the purposes of descriptive analysis, we only work with predictors names
  # i.e., we exclude confounding and 'split' variabes (CASE and TRAIN)

  start <- length(dcols) + 1
  end <- ncol(s)
  covnames <- colnames(s)[start:end]
  groups <- create_covars_groups(dt_colnames = covnames, ftargetsl, create_groups_for_single_covariates = TRUE,
    singletons_as_one_group = FALSE)
  t5[, seg] <- describe_group_matrix(s, groups, start, end)
  # shift by 3
  seg <- seg + 3
}

# formate
t5n <- data.frame(t(t5)) %>%
          summarise("Total has FTs (% of predictors)" = paste0(g_tot, ' (', g_totp, ')'),
                    "No. groups" = g_n,
                    "mean (SD)" = paste0(g_m, ' (', g_sd, ')'),
                    "range" = paste0('(', g_min, ' - ', g_max, ')'),
                    "No. singletons (% of predictors)" = paste0(s_tot, ' (', s_totp, ')'),
          ) %>%
          t() %>%
          data.frame()
names(t5n) <- colnames(t5)
write.table(t5n, file.path("reports", "t5.tsv"),
            sep = '\t', quote = FALSE)

# t7; atcicd
ncol <- (length(ades) * 3)
t7 <- data.frame(matrix(ncol = ncol, nrow = 9))
colnames(t7) <- c("gib_ca", "gib_ctrl", "gib_tot",
                  "icb_ca", "icb_ctrl", "icb_tot"
)
rownames(t7) <- c("g_tot", "g_totp", "g_n", "g_m", "g_sd", "g_min", "g_max", "s_tot", "s_totp")
seg <- c(1:3)

for (i in 1:length(pld)) {
  s <- pld[[i]]

  # for the purposes of descriptive analysis, we only work with predictors names
  # i.e., we exclude confounding and 'split' variabes (CASE and TRAIN)

  start <- length(dcols) + 1
  end <- ncol(s)
  covnames <- colnames(s)[start:end]
  groups <- create_covars_groups(dt_colnames = covnames, demographic_vars = NULL, ftargetsl = NULL,
    create_groups_for_single_covariates = TRUE,
    singletons_as_one_group = FALSE)
  t7[, seg] <- describe_group_matrix(s, groups, start, end)

  # get number of drug groups
  # dr_n <- length(grep(x = names(groups), pattern = '^dr'))
  # t7[10:11, seg] <- c(dr_n,
  #                     round((dr_n / length(covnames)) * 100, 1)
  #                   )
  # shift by 3
  seg <- seg + 3
}

# formate
t7n <- data.frame(t(t7)) %>%
          summarise("Total grouped (% of predictors)" = paste0(g_tot, ' (', g_totp, ')'),
                    "No. groups" = g_n,
                    "mean (SD)" = paste0(g_m, ' (', g_sd, ')'),
                    "range" = paste0(g_min, ' - ', g_max)
          ) %>%
          t() %>%
          data.frame()
names(t7n) <- colnames(t7)
write.table(t7n, file.path("reports", "t7.tsv"),
            sep = '\t', quote = FALSE)

# doac
t6 <- data.frame(matrix(nrow = ncol, ncol= 2))
rownames(t6) <- c("gib_ca", "gib_ctrl", "gib_tot",
                  "icb_ca", "icb_ctrl", "icb_tot"
)
colnames(t6) <- c("doacs_n", "doacs_p")

seg <- c(1:3)
doacs <- c("B01AF01", "B01AF02", "B01AF03", "B01AE07")
for (i in 1:length(pld)) {
  doac_tot <- cube(pld[[i]][, .SD, .SDcols = c('CASE', doacs)],
                   j = c(list(COUNT = .N), lapply(.SD, sum)),
                   .SDcols = doacs,
                   by = 'CASE')[c(2,1,3), ]
  doac_tot[, `:=` (doacs_n = apply(doac_tot[, .SD, .SDcols = doacs], 1, sum))]
  doac_tot[, `:=` (doacs_p = round((doacs_n / COUNT) * 100, 1))]

  t6[seg, ] <- doac_tot[, list(doacs_n, doacs_p)]
  seg <- seg + 3
}

t6n <- t6 %>%
          summarise("Yes (%)" = paste0(doacs_n, ' (', doacs_p, ')')) %>%
          t() %>%
          data.frame()
names(t6n) <- rownames(t6)
write.table(t6n, file.path("reports", "t6.tsv"),
            sep = '\t', quote = FALSE)

