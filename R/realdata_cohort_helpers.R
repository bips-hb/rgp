# provide list of ICDs as outcome of interest for cases (ADE)
get_ade_icd_list <- function(ADE_FILE) {
  # read csv
  ade <- read.csv(ADE_FILE, sep = ',')

  # split VAR_NAME by VAR_GRP
  ade_icd_list <- split(ade$VAR_NAME, factor(ade$VAR_GRP))

  return(ade_icd_list)
}

# functions to filter cases based on diagnosis
add_ade_col <- function(diag, outcome_diag_type, outcome_icds_pattern) {
  diag[, ADE := ""][]

  sapply(X = names(outcome_icds_pattern),
         FUN = function(x) {
           diag[((DIAG_TYPE %in% outcome_diag_type) &
                   grepl(x = ICD, pattern = outcome_icds_pattern[[x]])),
                ADE := x][]
         })

  return(diag)
}

add_sec_ade_col <- function(diag, outcome_diag_type, outcome_icds_pattern) {
  diag[, ADE2 := ""][]

  sapply(X = names(outcome_icds_pattern),
         FUN = function(x) {
           diag[((!(DIAG_TYPE %in% outcome_diag_type)) &
                   grepl(x = ICD, pattern = outcome_icds_pattern[[x]])),
                ADE2 := x][]
         })

  return(diag)
}

# function to split data table by birth year into chunks
# and to sample paired controls from those chunks
# https://stackoverflow.com/questions/8358098/how-to-set-seed-for-random-simulations-with-foreach-and-domc-packages
create_matched_cc_cohort <- function(pat,
                                     ttflag = NULL,
                                     first_split_var = "BIRTHY",
                                     nmatches = 4) {
  if (!is.null(ttflag)) {
    pat <- pat[TRAIN == ttflag]
  }

  expected <- split(pat, by = first_split_var)
  expected_full <- vector(mode = "list")

  # lapply(expected, function(x) {x[CASE == 1, .N]})
  for (j in 1:length(expected)) {
    if (dim(expected[[j]][CASE == 1])[1] != 0) {
    expected_full[[j]] <- expected[[j]]
    }
  }

  expected_full <- Filter(Negate(is.null), expected_full)

  rs <-  foreach(i = 1:length(expected_full),
                .packages = c('data.table'),
                .combine = rbind) %dopar% {
    set.seed(i)
    subs <- expected_full[[i]]
    subs_cases <- subs[CASE == 1]
    subs_ctrl <- subs[CASE == 0 & PAIR == 0]

    l <- vector(mode = "list", length = nrow(subs_cases))

    for (ca in 1:nrow(subs_cases)) {
      r <- subs_cases[ca, RID]
      s <- subs_cases[ca, SEX]
      b <- subs_cases[ca, BIRTHY]
      ade <- subs_cases[ca, ADE]
      # in case byear +/- 1 year: b <- seq(b-1, b+1)

      # find match
      gp <- subs_ctrl[SEX == s & BIRTHY == b & PAIR == 0 & (is.na(DEATH_DAT) | DEATH_DAT <= INDEX_DAT), ]

      # sample from them
      # in case nrow < nmatches required then all matches are stored
      # in case 0 matches found, an empty (NA) row will return
      # TODO: ifelse(nrow(gp) > 0, because it returns
      # one empty record if nrow(gp) == 0,
      # or nrow(gp) which can be < 4
      smr <- tryCatch(
        sample(nrow(gp), nmatches),
        error = function(e) {
          1:nrow(gp)
        })

      # TODO: read covars table and check if ncovar = 0
      # sm <- gp[smr, ][, INDEX_DAT := subs_cases[ca, INDEX_DAT]][] #PROD
      sm <- gp[smr, ][, c("INDEX_DAT", "TIME_TO_EVENT") := subs_cases[ca, list(INDEX_DAT, TIME_TO_EVENT)]][]
      rm(gp); gc()

      # flag and never remove sampled items
      # because cases from same year need ctrls too (without replacement)
      # delete row slows down the function by reallocating memory for the whole
      # data.table. Slow down from 1.5h to 24h
      subs_ctrl[sm, on = "IDNUM", PAIR := r] #PROD
      # subs_ctrl[sm, on = "IDNUM"]

      l[[ca]] <- rbind(subs_cases[ca, ], sm)
      l[[ca]][, c("PAIR", "SUBCOHORT") := list(r, ade)]
    }

    rbindlist(l)
  }

  return(rs)
}

dcast_dt <- function(dt, dcols, dvar) {
  dformula <- as.formula(
    paste0(paste0(dcols, collapse = ' + '), ' ~ ', dvar)
  )
  dcasted <- data.table::dcast.data.table(
    dt,
    dformula,
    value.var = dvar,
    fun.aggregate = function(x) length(x),
    drop = c(TRUE, FALSE),
    fill = 0
  )
  # remove the NA column (in case of straightforward dcast)
  suppressWarnings(dcasted[, "NA" := NULL][])

  data.table::setkeyv(dcasted, dcols)

  return(dcasted)
}

chunk_dcast_dt <- function(num_chunks = 50, dt, dcols, dvar) {
  indx <- split(seq(nrow(dt)), dt[[dvar]])
  n <- length(indx)
  chunks <- BBmisc::chunk(1:n, n.chunks = num_chunks)
  dcast_list <- vector(mode = 'list', length = num_chunks)

  for (i in 1:num_chunks) {
    ind <- unlist(indx[chunks[[i]]], use.names = FALSE)
    dcast_list[[i]] <- dcast_dt(dt[ind, ], dcols, dvar)
  }

  # merge chuncks
  # https://daranzolin.github.io/2016-12-10-join-list-dataframes/
  dcm <- dcast_list %>%
    purrr::reduce(dplyr::full_join) %>%
    dplyr::mutate_if(is.numeric, as.integer) %>%
    data.table(key = 'IDNUM')
    # mutate_if(is.numeric, replace_na, replace = 0)

  rm(dcast_list); gc()

  # replace NA rows with 0, they come from the Reduce merge
  # dcm[is.na(dcm)] <- 0

  # workaround to return patients with 0 covars
  if (nrow(dcm) != uniqueN(dt, by = 'IDNUM')) {
    dcasted <- rbind(dt[!dcm, on = 'IDNUM'][, ..dcols],
                     dcm,
                     fill = TRUE)
  } else {
    dcasted <- dcm
  }

  # replace NA rows with 0, they come from the Reduce merge AND rbind
  # dcasted[is.na(dcasted)] <- 0 # produces erros
  # workaround: https://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table
  for (j in seq_len(ncol(dcasted))) {
    set(dcasted, which(is.na(dcasted[[j]])), j, 0)
  }

  setkeyv(dcasted, dcols)

  return(dcasted)
}

call_dcast_dt <- function(...) {
  tryCatch(
    dcast_dt(...),
    error = function(e) { chunk_dcast_dt(...) }
  )
}

# create a matrix for each ADE sub-cohort
create_cov_mat <- function(subcohort = NULL,
                            disease_covar, drug_covar,
                            dcols, dvar, nchar_icd = 4) {

  if (is.null(subcohort)) {
    dr_dt <- drug_covar
    ds_dt <- disease_covar[, .(COV = unique(substr(COV, 1, nchar_icd))), by = dcols]
  }
  else {
    dr_dt <- drug_covar[subcohort][DAT <= INDEX_DAT, .(IDNUM, COV)][subcohort[, ..dcols], on = 'IDNUM']
    ds_dt <- disease_covar[subcohort][DAT <= INDEX_DAT, .(IDNUM, COV)
      ][, .(COV = unique(substr(COV, 1, nchar_icd))),
      by = 'IDNUM'][subcohort[, ..dcols], on = 'IDNUM']
  }

 # keep drugs before index date, and cast them
  pdc <- call_dcast_dt(
    dt = dr_dt,
    dcols = dcols,
    dvar = dvar
  )

 # same with diseases
 pdsc <- call_dcast_dt(
  dt = ds_dt,
  dcols = dcols,
  dvar = dvar
 )

  # TODO: convert to matrix?
  ddcov <- merge(pdc, pdsc, all = TRUE)
  # ddcov[, IDNUM := NULL]

  return(ddcov)
}

#### summary ####
summarize_cov_matrix <- function(sums) {
  list(mean = round(mean(sums), 1),
       sd = round(sd(sums), 1),
       range1 = min(sums),
       range2 = max(sums))
}

describe_cov_matrix <- function(s, start, end) {
  t4 <- data.frame(matrix(ncol = 3, nrow = 7))
  colnames(t4) <- c("ca", "ctrl", "tot")
  rownames(t4) <- c("tot", "m", "sd", "min", "max", "zn", "zp")

  tot <- cube(s, j = lapply(.SD, sum), .SDcols = start:end, by = 'CASE')[c(2,1,3), ]
  ncov <- length(start:end)
  zn <- apply(tot[, .SD == 0, .SDcols = 2:ncol(tot)], 1, function(x) length(which(x)))

  # total (dim exclue zero)
  t4[1, 1:3] <- ncov - zn

  t4[2:5, 1] <- unlist(summarize_cov_matrix(s[CASE == 1, rowSums(.SD), .SDcols = start:end]))
  t4[2:5, 2] <- unlist(summarize_cov_matrix(s[CASE == 0, rowSums(.SD), .SDcols = start:end]))
  t4[2:5, 3] <- unlist(summarize_cov_matrix(s[, rowSums(.SD), .SDcols = start:end]))

  t4[6, 1:3] <- zn
  t4[7, 1:3] <- round(zn / ncov * 100, 1)

  return(t4)
}

describe_groups_in_cohort <- function(groups) {
  tg <- vector(length = 9)
  names(tg) <- c("g_tot", "g_totp", "g_n", "g_m", "g_sd", "g_min", "g_max",
                 "s_tot", "s_totp")

  ncov <- length(unique(unlist(groups)))

  # the groups are both of FTs and dummy groups for singletons
  # is it possible that FTs can contain only one predictor?
  ft <- names(groups)[!(grepl(x = names(groups), pattern = '^grp'))]

  g_tot <- length(unique(unlist(groups[ft])))

  s_tot <- ncov - g_tot

  tg[1] <- g_tot
  tg[2] <- round((g_tot / ncov) * 100, 1)
  tg[3] <- length(ft)
  tg[4:7] <- unlist(summarize_cov_matrix(
                        sapply(groups[ft], length)
  ))
  tg[8] <- s_tot
  tg[9] <-  round((s_tot / ncov) * 100, 1)

  # the code below assumes that we want to describe groups > 1 even if they are
  # FTs of 1
  # g_tot <- length(unique(unlist(
  #   Filter(groups, f = function(x) {length(x) > 1})
  # )))
  #
  # s_tot <- length(unique(unlist(
  #   Filter(groups, f = function(x) {length(x) == 1})
  # )))
  #
  # tg[1] <- g_tot
  # tg[2] <- round((g_tot / ncov) * 100, 1)
  # tg[3] <- length(Filter(groups, f = function(x) {length(x) > 1}))
  # tg[4:7] <- unlist(summarize_cov_matrix(
  #                       sapply(Filter(groups, f = function(x) {length(x) > 1}),
  #                              length)
  # ))
  # tg[8] <- s_tot
  # tg[9] <-  round((s_tot / ncov) * 100, 1)

  return(tg)
}

describe_group_matrix <- function(s, groups, start, end) {
  tg <- data.frame(matrix(ncol = 3, nrow = 9))
  colnames(tg) <- c("ca", "ctrl", "tot")
  rownames(tg) <- c("g_tot", "g_totp", "g_n", "g_m", "g_sd", "g_min", "g_max",
                    "s_tot", "s_totp")

  tg[, 3] <- describe_groups_in_cohort(groups)

  # remove predictors that are not observed in CA/CTRLs
  tot <- cube(s, j = lapply(.SD, sum), .SDcols = start:end, by = 'CASE')[c(2,1), ]
  zero_list <- apply(tot[, .SD == 0, .SDcols = 2:ncol(tot)], 1, function(x) which(x))

  # cases
  groups_ca <- Filter(length, lapply(groups,
        function(x) {
          setdiff(x, zero_list[[1]])
        }
  ))

  tg[, 1] <- describe_groups_in_cohort(groups_ca)

  # ctrl
  groups_ctrl <- Filter(length, lapply(groups,
        function(x) {
          setdiff(x, zero_list[[2]])
        }
  ))

  tg[, 2] <- describe_groups_in_cohort(groups_ctrl)

  return(tg)
}

split_dt_by_group <- function(s, covnames, groups) {
  predictors <- s[, ..covnames]
  lapply(groups, function(x) {
    subset(predictors, select = x)
  })
}

# n.cov is the number of covariates, and groups are the groups as usual.
create_groups_for_single_covariates <- function(ncov, groups, dt_colnames, rename = FALSE) {
  # Any variable not in a group gets assigned their own group
  gnames <- names(groups)
  glen <- length(groups)
  variables_in_groups <- do.call(c, groups)
  variables_not_in_groups <- as.list(setdiff(1:ncov, variables_in_groups))

  if (isTRUE(rename)) {
    names(variables_not_in_groups) <- paste0("grp", glen+1:length(variables_not_in_groups))
  } else {
    names(variables_not_in_groups) <- dt_colnames[unlist(variables_not_in_groups)]
  }

  return(variables_not_in_groups)
}

# n.cov is the number of covariates, and groups are the groups as usual.
append_groups_for_single_covariates <- function(ncov, groups) {
  # Any variable not in a group gets assigned their own group
  variables_not_in_groups <- create_groups_for_single_covariates(ncov, groups)
  return(c(groups, variables_not_in_groups))
}

assign_to_atc_icd <- function(dt_colnames, demographic_vars, nchar_atc = 4,
  nchar_icd = 3) {
    # avoid using length function
    drug_disease_vars <- dt_colnames[!dt_colnames %in% demographic_vars]

    drugs_vars <- drug_disease_vars[1:which(drug_disease_vars == 'V90N')]
    diseases_vars <- drug_disease_vars[(which(drug_disease_vars == 'V90N') + 1):length(drug_disease_vars)]

    drug_grp_name <- paste0('dr', substr(drugs_vars, 1, nchar_atc))
    disease_grp_name <- paste0('ds',substr(diseases_vars, 1, nchar_icd))

    pathways <- data.table(var_name = drug_disease_vars, grp_name = c(drug_grp_name, disease_grp_name), index = as.numeric(na.omit(match(drug_disease_vars, dt_colnames))))
    pathways <- split(pathways[, list(index, grp_name)], by = 'grp_name')
    pathways <- lapply(pathways, function(x) {x[, index]})

    return(pathways)
}

# assign all columns of a data table to ft groups, return index
create_covars_groups <- function(dt_colnames, demographic_vars = NULL, ftargetsl = NULL,
  nchar_atc = 4,
  nchar_icd = 3,
  create_groups_for_single_covariates = TRUE,
  singletons_as_one_group = FALSE,
  rename = FALSE) {
  if (is.null(ftargetsl)) {
    pathways <- assign_to_atc_icd(dt_colnames = dt_colnames,
      demographic_vars = demographic_vars,
      nchar_atc = nchar_atc, nchar_icd = nchar_icd)

  } else {
    # assign a vector of covariate names to groups by string matching, return index
    pathways <- assign_to_ft(dt_colnames, ftargetsl)

    # get leftover covariates
    singletons <- setdiff(1:length(dt_colnames), unique(unlist(pathways)))

    # get leftover dr and ds
    extra_colnames <- dt_colnames[unlist(singletons, use.names = FALSE)]

    # remove drugs, here we use full drug code, no need to increase search time
    extra_colnames <- extra_colnames[(which(extra_colnames == 'V90N') + 1):length(extra_colnames)]

    # assign by grep
    extra_pathways <- assign_to_ft_by_grep_disease(dt_colnames,
      extra_colnames, ftargetsl)

    # merge the two lists
    pathways <- Map(function(val1, val2) union(val1, val2), pathways, extra_pathways)

    # remove empty pathways
    pathways <- Filter(length, pathways)
  }

  if (isTRUE(create_groups_for_single_covariates)) {
    # get singletons
    singletons <- create_groups_for_single_covariates(length(dt_colnames), pathways,
      dt_colnames, rename = rename)

    if (isTRUE(singletons_as_one_group)) {
      # gnames <- names(singletons)[1:2]
      demographics_grp <- as.numeric(na.omit(match(demographic_vars, dt_colnames)))
      predictors_grp <- setdiff(unlist(singletons, use.names = FALSE), demographics_grp)

      singletons <- list(
        demographics_grp,
        predictors_grp
      )

      # not assume that there were singleton predictors
      singletons <- Filter(length, singletons)

      names(singletons) <- paste0("grp", length(pathways)+1:length(singletons))
    }

    # length(unique(unlist(pathways))) + length(unlist(singletons)) == length(dt_colnames)
    pathways <- c(pathways, singletons)
  }

  return(pathways)
}

# add the diseases that start with the ICD
assign_to_ft_by_grep_disease <- function(dt_colnames, extra_colnames, ftargetsl) {
  extra_pathways <- lapply(ftargetsl, function(x) {
    matches <- names(Filter(length, sapply(extra_colnames, function(j) {
      x[grep(x = x, pattern = paste0("^", j))]
    })));
    as.numeric(na.omit(match(matches, dt_colnames)))
  })

  return(extra_pathways)
}

# assign a vector of covariate names to groups, return index
assign_to_ft <- function(dt_colnames, ftargetsl) {
  # try to clean colnames
  try(
    dt_colnames <- gsub(x = dt_colnames,  pattern = "^d(r|s)", replacement = "")
  )

  # create a list, fill it in with covars indices, and name it
  pathways <- lapply(ftargetsl, function(x) {
    as.numeric(na.omit(match(x, dt_colnames)))
    # exact equivalent to which(dt_colnames %in% x)
  })

  return(pathways)
}

# get the predictors for each person in an ADE subcohort
get_ds_predictors <- function(subcohort, diag, amb, begin_yr) {
  ade <- subcohort[1, SUBCOHORT]

  # diseases
  # ds_response <- diag[(IDNUM %in% subcohort[, IDNUM] & ADE == ade), ]

  ds_pred <- rbind(diag[(IDNUM %in% subcohort[, IDNUM] & ADE != ade & !is.na(ICD) & DIAG_DAT >= as.IDate(paste0(begin_yr, '-01-01'))), ],
amb[(IDNUM %in% subcohort[, IDNUM] & DIAG_DAT >= as.IDate(paste0(begin_yr, '-01-01'))), ])
  setnames(ds_pred, c("ICD", "DIAG_DAT"), c("COV", "DAT"))

  setkey(ds_pred, IDNUM)
  setindex(ds_pred, COV)
  ds_pred[, FIRST_DAT := min(DAT), by = .(IDNUM, COV)]
  ds_pred <- unique(ds_pred[, list(IDNUM, COV, FIRST_DAT)])
  setnames(ds_pred, "FIRST_DAT", "DAT")
  # uniqueN(ds_pred, by = 'IDNUM')

  # merge disp with cohort
  sc_ds <- merge(subcohort, ds_pred, all.x = TRUE)

  # keep cov before INDEX_DAT
  sc_ds <- sc_ds[DAT <= INDEX_DAT, list(IDNUM, SEX, GKZ5, AGE, TIME_TO_EVENT, COV, CASE, TRAIN)]

  if (uniqueN(sc_ds, by = 'IDNUM') != uniqueN(subcohort, by = 'IDNUM')) {
      sc_ds <- rbind(subcohort[!sc_ds, list(IDNUM, SEX, GKZ5, AGE, TIME_TO_EVENT, CASE, TRAIN), on = 'IDNUM'], sc_ds, fill = TRUE)
  }
  # uniqueN(sc_ds, by = 'IDNUM')

  return(sc_ds)
}

get_dr_predictors <- function(subcohort, disp) {
  ade <- subcohort[1, SUBCOHORT]
  setkey(disp, IDNUM)

  sc_dr <- disp[subcohort, ]
  setnames(sc_dr, c("ATC", "DEL_DAT"), c("COV", "DAT"))
  setkey(sc_dr, IDNUM)
  setindex(sc_dr, COV)
  # uniqueN(sc_dr, by = 'IDNUM')

  # keep cov before INDEX_DAT
  sc_dr <- sc_dr[DAT <= INDEX_DAT, list(IDNUM, SEX, GKZ5, AGE, TIME_TO_EVENT, COV, CASE, TRAIN)]

  if (uniqueN(sc_dr, by = 'IDNUM') != uniqueN(subcohort, by = 'IDNUM')) {
      sc_dr <- rbind(subcohort[!sc_dr, list(IDNUM, SEX, GKZ5, AGE, TIME_TO_EVENT, CASE, TRAIN), on = 'IDNUM'], sc_dr, fill = TRUE)
  }
  # uniqueN(sc_dr, by = 'IDNUM')

  return(sc_dr)
}

matrix_to_list <- function(boolean_matrix) {
	llist <- vector(mode = "list", length = dim(boolean_matrix)[1])
	names(llist) <- rownames(boolean_matrix)

	for (row_index in 1:dim(boolean_matrix)[1]) {
		m = dimnames(boolean_matrix[row_index, boolean_matrix[row_index,] == 1, drop = F])
		llist[[row_index]] = m[[2]]
	}
	return(llist)
}
