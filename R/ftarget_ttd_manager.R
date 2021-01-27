# ----------------------- data mangement ------------------------------
#' @import RSQLite
#' @export
connect_ttd <- function(sqlitedb = "TTD.db") {
  ttd_con <- dbConnect(SQLite(), sqlitedb)
  return(ttd_con)
}

#' @import RSQLite
#' @export
disconnect_ttd <- function(ttd_con) {
  dbDisconnect(ttd_con)
}

#' @import RSQLite
#' @export
process_query <- function(query, param = NULL) {
  rs <- dbSendQuery(ttd_con, query)
  if(!missing(param)) {
    dbBind(rs, param = param)
  }
  df <- dbFetch(rs)
  dbClearResult(rs)
  return(df)
}

#' @import RSQLite
#' @export
read_ttd <- function(ttd_con, ttd_table) {
  t <- data.table(dbReadTable(ttd_con, ttd_table))
  return(t)
}

ttd_target2kegg <- function() {
  # create ttd_target_id kegg_pathway_id map (only for successful targets)
  # drop table if exists
  if (dbExistsTable(ttd_con, 'target_kegg_map')) {
    dbRemoveTable(ttd_con, 'target_kegg_map')
  }
  create_kegg_map_q <- 'CREATE TABLE target_kegg_map AS SELECT ttd_target_id, GROUP_CONCAT(kegg_pathway_id, "; ") AS kegg_pathway_id
    FROM (SELECT DISTINCT ttd_target_id, kegg_pathway_id FROM ttd_target_kegg) GROUP BY ttd_target_id'
  nraffected <- dbExecute(ttd_con, create_kegg_map_q)
  # assign data type not possible
}

# ----------------------- drug2atc2path ------------------------------
ttd_drug2atc2path <- function() {
  # select all the drugs and join their successful targets (target can be NA, maybe experimental/trial)
  drug_target_q <- 'SELECT a.ttd_drug_id, a.drug_value, b.ttd_target_id FROM ttd_crossmatching a
  LEFT JOIN (SELECT ttd_target_id, target_value FROM ttd_target WHERE target_key = :x AND ttd_target_id IN
  (SELECT DISTINCT ttd_target_id FROM ttd_target WHERE target_key = :y AND target_value = :z)) b ON
  a.drug_value = b.target_value WHERE a.drug_key = :zz'
  drug_target_p <- list(x = "Drug(s)", y = "Type of target", z = "Successful target", zz = "DrugName")
  drug_target <- process_query(drug_target_q, drug_target_p)
  # full outer join example: http://www.sqlitetutorial.net/sqlite-full-outer-join/

  # drop table if exists
  if (dbExistsTable(ttd_con, 'drug_target')) {
    dbRemoveTable(ttd_con, 'drug_target')
  }
  # write it the sqlitedb
  dbWriteTable(ttd_con, "drug_target", drug_target)

  # and add atc column
  add_atc_column_q <- 'ALTER TABLE drug_target ADD atc TEXT;'
  nraffected <- dbExecute(ttd_con, add_atc_column_q) #statement, getRowsAffected, clearResults

  # should be only one ATC line per ttd_drug_id
  insert_atc_q <- 'UPDATE drug_target SET atc = (SELECT drug_value FROM ttd_crossmatching
  WHERE drug_target.ttd_drug_id = ttd_crossmatching.ttd_drug_id AND ttd_crossmatching.drug_key = "SuperDrug ATC")'
  nraffected <- dbExecute(ttd_con, insert_atc_q)

  # add the kegg pathway id
  # TTD targets can encompase >1 KEGG pathway
  add_kegg_column_q <- 'ALTER TABLE drug_target ADD kegg_pathway_id TEXT;'
  nraffected <- dbExecute(ttd_con, add_kegg_column_q)

  # subquery did not work in this condition
  # solutions: https://stackoverflow.com/questions/18285713/how-to-avoid-duplication-in-group-concat

  insert_kegg_q <- 'UPDATE drug_target SET kegg_pathway_id = (SELECT kegg_pathway_id FROM target_kegg_map
  WHERE drug_target.ttd_target_id = target_kegg_map.ttd_target_id)'
  nraffected <- dbExecute(ttd_con, insert_kegg_q)
}

# ----------------------- ds2icd2path ------------------------------
ttd_ds2icd2path <- function() {
  # # drop table if exists
  if (dbExistsTable(ttd_con, 'disease_target')) {
    dbRemoveTable(ttd_con, 'disease_target')
  }

  # select all the drugs and join their successful targets (target can be NA, maybe experimental/trial)
  create_ds_icd_q <- 'CREATE TABLE disease_target AS SELECT ttd_target_id, indication, icd10 FROM ttd_target_disease
  WHERE ttd_target_id IN
  (SELECT DISTINCT ttd_target_id FROM ttd_target WHERE target_key = "Type of target" AND target_value = "Successful target")'

  # only successful targets (this makes all diseases have targets, can be complemented by drug-disease pairs)
  # create_ds_icd_q <- 'CREATE TABLE disease_target AS SELECT ttd_target_id, indication, icd10 FROM target_disease WHERE ttd_target_id IN (SELECT ttd_target_id FROM target_kegg_map)'
  nraffected <- dbExecute(ttd_con, create_ds_icd_q)

  # add the kegg pathway id
  # TTD targets can encompase >1 KEGG pathway
  add_kegg_column_q <- 'ALTER TABLE disease_target ADD kegg_pathway_id TEXT;'
  nraffected <- dbExecute(ttd_con, add_kegg_column_q)

  insert_kegg_q <- 'UPDATE disease_target SET kegg_pathway_id = (SELECT kegg_pathway_id FROM target_kegg_map
  WHERE disease_target.ttd_target_id = target_kegg_map.ttd_target_id)'
  nraffected <- dbExecute(ttd_con, insert_kegg_q)

  # or add to target_disease then if not successful target, drop the target row
}
# ----------------------- dr2ds ------------------------------
ttd_dr2ds <- function() {
  # TODO: create a new table where combinations are split
  # drop table if exists
  if (dbExistsTable(ttd_con, 'drug_disease')) {
    dbRemoveTable(ttd_con, 'drug_disease')
  }
  # sqlite option: https://gist.github.com/dannguyen/9a2b1505bbe097b765a9e7c2e1f7a23c
  drug_disease <- read_ttd(ttd_con, ttd_table = 'ttd_drug_disease')

  # keep important columns
  drug_disease <- drug_disease[, c(1,3,5)]

  # split combinations
  drug_disease <- cSplit(drug_disease, splitCols = 'ttd_drug_id', sep = "-", direction = "long", fixed = TRUE)

  # write it the sqlitedb
  dbWriteTable(ttd_con, "drug_disease", drug_disease)

  # create index on ttd_drug_id
  dr_ds_index_q <- 'CREATE INDEX drug_disease_id_ix ON drug_disease (ttd_drug_id)'
  nraffected <- dbExecute(ttd_con, dr_ds_index_q)
  dr_ds_index_q <- "CREATE INDEX drug_icd10_ix ON drug_disease (ttd_drug_id, indication, icd10)"
  nraffected <- dbExecute(ttd_con, dr_ds_index_q)
  # and add atc column
  add_atc_column_q <- 'ALTER TABLE drug_disease ADD atc TEXT;'
  nraffected <- dbExecute(ttd_con, add_atc_column_q)

  # should be only one ATC line per ttd_drug_id
  insert_atc_q <- 'UPDATE drug_disease SET atc = (SELECT drug_value FROM ttd_crossmatching
  WHERE drug_disease.ttd_drug_id = ttd_crossmatching.ttd_drug_id AND ttd_crossmatching.drug_key = "SuperDrug ATC")'
  nraffected <- dbExecute(ttd_con, insert_atc_q)
}

# ------------------- work on the tables as data.tables ------------------------
# exclude ttd_ids that do not have a disease
get_drugs_to_keep <- function(ttd_con) {
  drug_disease <- read_ttd(ttd_con, ttd_table = 'drug_disease')
  # remove drugs that have no atc
  drug_disease <- drug_disease[!is.na(atc)]
  # exclude atcs that do not have icd10 e.g., ttd_drug_id D0P0HT
  drug_disease <- drug_disease[icd10 != '']

  drug_disease <- cSplit(drug_disease, splitCols = 'atc', sep = "; ", direction = "long", fixed = TRUE)

  atc_with_disease <- unique(drug_disease$atc)
  # drug combinations are excluded
  # SQLite have no restrictions on length of VARCHAR:
  # https://stackoverflow.com/questions/6109532/what-is-the-maximum-size-limit-of-varchar-data-type-in-sqlite

  # length(unique(drug_target$ttd_drug_id)) #1856
  # length(unique(drug_disease$ttd_drug_id)) #1709
  # drug_target[!(ttd_drug_id %in% drug_disease$ttd_drug_id) & is.na(ttd_target_id), ]
  # ttd targets are double those of kegg
  # length(unique(drug_target$ttd_target_id))
  # length(unique(na.omit(unlist(tstrsplit(drug_target$kegg_pathway_id, "; ", fixed=TRUE)))))

  # there are drugs with no ttd target id
  # remove drugs that have no target and no disease
  # remove drugs that have no disease
  # drugs disease-, target+
  # excl_drugs <- drug_target[!(ttd_drug_id %in% unique(drug_disease$ttd_drug_id)) & !is.na(kegg_pathway_id), ttd_drug_id]
  # # drugs disease-, target-
  # excl_drugs <- c(excl_drugs, drug_target[!(ttd_drug_id %in% unique(drug_disease$ttd_drug_id)) & is.na(kegg_pathway_id), ttd_drug_id])

  ## drugs disease+, target-
  # length(unique(drug_disease[(ttd_drug_id %in% drug_target[is.na(kegg_pathway_id), ttd_drug_id]), ttd_drug_id]))
  ## drugs disease+, target+
  # length(unique(drug_target[(ttd_drug_id %in% drug_disease$ttd_drug_id) & !is.na(kegg_pathway_id), ttd_drug_id]))
  return(atc_with_disease)
}

# ---------------------- make the matrices ---------------------
make_drug_path_matrix <- function(ttd_con) {
  atc_with_disease <- get_drugs_to_keep(ttd_con)
  drug_target <- read_ttd(ttd_con, ttd_table = 'drug_target')
  # remove drugs that have no atc
  drug_target <- drug_target[!is.na(atc)]

  # 1709 unique ttd drug ids
  # drugs are duplicated because a drug can have >1 ttd target

  # split by atc
  drug_target <- cSplit(drug_target, splitCols = 'atc', sep = "; ", direction = "long", fixed = TRUE)

  # keep only drugs with atc that have a disease
  drug_target <- drug_target[(atc %in% atc_with_disease), ]

  # split by kegg_pathway_id
  drug_target <- cSplit(drug_target, splitCols = 'kegg_pathway_id', sep = "; ", direction = "long", fixed = TRUE)

  drug_target_m <- data.frame(t(create_boolean_matrix(na.omit(drug_target[, c(4,5)]))))
  return(drug_target_m)

  # two questions:
  # dummy pathways [HOW TO MAKE DUMMY PATHWAYS, same name, x1..]
  # pathways that are not in humans (do not start with 'hsa') [CHECKED: KEEP]
  # double check
  # names(drug_target_m)[!(names(drug_target_m) %in% atc_with_disease)]
}

make_drug_disease_matrix <- function(ttd_con) {
  # excl_drugs <- get_drugs_to_keep(ttd_con)
  atc_with_disease <- get_drugs_to_keep(ttd_con)
  drug_disease <- read_ttd(ttd_con, ttd_table = 'drug_disease')
  # remove drugs that have no atc
  drug_disease <- drug_disease[!is.na(atc)]
  # exclude atcs that do not have icd10 e.g., ttd_drug_id D0P0HT
  drug_disease <- drug_disease[icd10 != '']

  # split by atc
  drug_disease <- cSplit(drug_disease, splitCols = 'atc', sep = "; ", direction = "long", fixed = TRUE)
  # keep only drugs with atc that have a disease
  drug_disease <- drug_disease[(atc %in% atc_with_disease), ]
  # split by icd10
  drug_disease <- cSplit(drug_disease, splitCols = 'icd10', sep = ", ", direction = "long", fixed = TRUE)

  # make table smaller
  # drug_disease <- drug_disease[, c(1,3,5,6)]

  # expand
  drug_disease_expanded <- expand_icd_2017(icd_dt = drug_disease, gm = T)

  drug_disease_m <- data.frame(t(create_boolean_matrix(na.omit(drug_disease_expanded[, c(4,5)]))))
  return(drug_disease_m)
}

make_disease_target_matrix <- function(ttd_con) {
  # what about disease_target
  disease_target <- read_ttd(ttd_con, ttd_table = 'disease_target')
  # all diseases that have icd
  # disease_target[is.na(icd10)]
  # drug_disease[is.na(icd10)]

  # keep diseases without target
  # disease drug-, target-
  # disease_target[(indication %in% unique(drug_disease$indication)) & is.na(kegg_pathway_id)]
  disease_target <- cSplit(disease_target, splitCols = 'icd10', sep = ", ", direction = "long", fixed = TRUE)
  # I triple checked them and they are correct
  disease_target[icd10 == 'K86.0K86.1', 'icd10'] <- 'K86.0, K86.1'
  # disease_target[icd10 == 'S00.0S09', 'icd10'] <- 'S00-S09' has no successful target
  disease_target <- cSplit(disease_target, splitCols = 'icd10', sep = ", ", direction = "long", fixed = TRUE)

  # split kegg pathway
  disease_target <- cSplit(disease_target, splitCols = 'kegg_pathway_id', sep = "; ", direction = "long", fixed = TRUE)

  # all diseases have either successful kegg targets or no targets at all

  # expand
  disease_target_expanded <- expand_icd_2017(icd_dt = disease_target, gm = T)

  disease_target_m <- data.frame(create_boolean_matrix(na.omit(disease_target_expanded[, c(4,5)])))

  return(disease_target_m)
}

# --------------------- expand icd ---------------------------------------

#' @export
read_icd10cm_2017 <- function(icd10cm_2017_file = paste(system.file("extdata", package = "RGP"),
                                                        "icd10_2017", "icd10cm_codes_2017.txt", sep = "/"),
                              level = 5) {
  # read CM ICDs file
  icd10cm_2017 <- readLines(icd10cm_2017_file)
  icd10cm_2017 <- data.table(icd10cm_2017)
  # split into two columns taking in the first 8 char
  icd10cm_2017$icd <- substr(icd10cm_2017$icd10cm_2017, 1, 8)
  icd10cm_2017$description <- substr(icd10cm_2017$icd10cm_2017, 9, nchar(icd10cm_2017$icd10cm_2017))
  # trim all trailing white spaces
  icd10cm_2017$icd <- trimws(icd10cm_2017$icd)
  # drop the first column
  icd10cm_2017 <- icd10cm_2017[, c(2,3)]
  # transform into . format
  icd10cm_2017$code_dot <- ifelse(nchar(icd10cm_2017$icd) > 3,
                                  paste(substr(icd10cm_2017$icd, 1, 3),
                                        substr(icd10cm_2017$icd, 4, level),
                                        # make it shorter than nchar(icd10cm_2017$icd)
                                        sep = '.'),
                                  icd10cm_2017$icd)
  # remove the too refined level
  icd10cm_2017_cmprsd <- icd10cm_2017[!duplicated(icd10cm_2017$code_dot), ]
  # icd-cm until 7
  # sort(unique(nchar(icd10cm_2017$icd)))
  # now compressed to 5 a little improved than GM, because GM has 5 levels for selected diseases

  # make a highest order column
  icd10cm_2017_cmprsd$chapter <- substr(icd10cm_2017_cmprsd$icd, 1, 3)

  return(icd10cm_2017_cmprsd)
}

read_icd10gm_2017 <- function(icd10gm_2017_file = paste(system.file("extdata", package = "RGP"), "icd10_2017", "icd10gm2017syst_kodes.txt", sep = "/")) {
  # read GM ICDs file
  icd10gm_2017 <- read.csv(icd10gm_2017_file, sep = ';', header = F)
  icd10gm_2017 <- data.table(icd10gm_2017)
  # make it smaller
  icd10gm_2017 <- icd10gm_2017[, c(1,6,7,8,11)]

  # no dot
  # code_dot
  # make a highest order column (chapter)
  # remove A00.- the highest level and keep the column
  icd10gm_2017 <- icd10gm_2017[grep(pattern = '-', fixed = T,
                                    icd10gm_2017$V6, invert = T)]
  icd10gm_2017$V6 <- substr(icd10gm_2017$V6, 1, 3)

  names(icd10gm_2017) <- c('l', 'chapter', 'code_dot', 'icd10', 'description')
  return(icd10gm_2017)
}


expand_icd_2017_start_stop <- function(icd10cm_2017_cmprsd, icd_pair) {
  icd_s <- unlist(strsplit(as.character(icd_pair), split = '-', fixed = TRUE))
  istart <- ifelse(grepl(x = icd_s[1], pattern = '.', fixed = T),
                   which(icd10cm_2017_cmprsd$code_dot == icd_s[1])[1],
                   which(icd10cm_2017_cmprsd$chapter == icd_s[1])[1])
  iend <- ifelse(grepl(x = icd_s[2], pattern = '.', fixed = T),
                 which(icd10cm_2017_cmprsd$code_dot == icd_s[2])[length(which(icd10cm_2017_cmprsd$code_dot == icd_s[2]))],
                 which(icd10cm_2017_cmprsd$chapter == icd_s[2])[length(which(icd10cm_2017_cmprsd$chapter == icd_s[2]))])
  expanded_icds <- icd10cm_2017_cmprsd$code_dot[istart:iend]
  return(expanded_icds)
}

#' @export
expand_icd_2017 <- function(icd_dt, gm = F) {
  if (isTRUE(gm)) {
    icd10cm_2017_cmprsd <- read_icd10gm_2017() # GM
  } else {
    icd10cm_2017_cmprsd <- read_icd10cm_2017() # CM
  }

  # make sure all icds are split
  icd_dt[icd10 == 'D37C33-C34', 'icd10'] <- 'D37, C33-C34'
  # "J12-18"
  icd_dt[icd10 == 'J12-18', 'icd10'] <- 'J12-J18'
  # "H00-H99"
  icd_dt[icd10 == 'H00-H99', 'icd10'] <- 'H00-H95'

  icd_dt <- cSplit(icd_dt, splitCols = 'icd10', sep = ', ', direction = 'long')

  # initialize new column
  icd_dt$expanded_icd <- icd_dt$icd10

  # handle (individuals not full classes)
  icd_wo_boarders <- unique(icd_dt$icd10[grep(icd_dt$icd10, pattern = '-', invert = T)])
  for (i in 1:length(icd_wo_boarders)) {
    # icd by icd
    x <- as.character(icd_wo_boarders[i])
    # get children
    expanded_icds <- icd10cm_2017_cmprsd[(grepl(x = icd10cm_2017_cmprsd$code_dot, pattern = x, fixed = T)), code_dot]

    icd_dt[icd10 == x, "expanded_icd"] <- ifelse(length(expanded_icds) > 0,
                                                 paste(expanded_icds, collapse=', '),
                                                 x)
  }

  if (isTRUE(gm)) {
    # handle these in GM
    # "I10-I16" # I10-I15
    # "E08-E13" # E10-E14
    # "T14.0-T14.1" # T14.00-T14

    # handle these in GM
    # K00-K95 #K00-K93
    # N39.3-N39.4 # N39.3-N39.48
    # C00-D49 # C00-D48
    icd_dt[icd10 == 'I10-I16', 'icd10'] <- 'I10-I15'
    icd_dt[icd10 == 'E08-E13', 'icd10'] <- 'E10-E14'
    icd_dt[icd10 == 'T14.0-T14.1', 'icd10'] <- 'T14.00-T14.1'

    icd_dt[icd10 == 'K00-K95', 'icd10'] <- 'K00-K93'
    icd_dt[icd10 == 'N39.3-N39.4', 'icd10'] <- 'N39.3-N39.48'
    icd_dt[icd10 == 'C00-D49', 'icd10'] <- 'C00-D48'
    # "O00-O9A"
    icd_dt[icd10 == 'O00-O9A', 'icd10'] <- 'O00-O99'
  } else {
    # handle these in CM
    "N39.3-N39.4" #"N39.3-N39.49 in CM
    icd_dt[icd10 == 'N39.3-N39.4', 'icd10'] <- 'N39.3-N39.49'
    # [1] "T14.0-T14.1" #only in GM
    icd_dt[icd10 == 'T14.0-T14.1', 'icd10'] <- 'T14.8'
    # [1] "R52.1-R52.2" #only in GM, CM R52 only
    icd_dt[icd10 == 'R52.1-R52.2', 'icd10'] <- 'R52'
    # [1] "E00-E90" E90 in GM, E89 in CM
    icd_dt[icd10 == 'E00-E90', 'icd10'] <- 'E00-E89'
    # [1] "B20-B24" # no B24 in CM, B24 in GM
    icd_dt[icd10 == 'B20-B24', 'icd10'] <- 'B20' # HIV
    # [1] "E70-E90" # E8989 in CM, E90 in GM
    icd_dt[icd10 == 'E70-E90', 'icd10'] <- 'E70-E89'
    # [1] "K00-K93" K93 is last in GM K92 in CM
    icd_dt[icd10 == 'K00-K93', 'icd10'] <- 'K00-K92'
    # [1] "T79.0-T79.1" T79.0X-T79.1X ##
    icd_dt[icd10 == 'T79.0-T79.1', 'icd10'] <- 'T79.0X-T79.1X'
    # [1] "F00-F99" F00 in GM, F01 in CM
    icd_dt[icd10 == 'F00-F99', 'icd10'] <- 'F01-F99'
    # [1] "C00-C97" C96 in CM, C97 in GM
    icd_dt[icd10 == 'C00-C97', 'icd10'] <- 'C00-C96'
    # [1] "K29.0-K29.7" K29.00-K29.71 ##
    icd_dt[icd10 == 'K29.0-K29.7', 'icd10'] <- 'K29.00-K29.71'
    # [1] "E80.0-E80.2" E80.0-E80.29 ##
    icd_dt[icd10 == 'E80.0-E80.2', 'icd10'] <- 'E80.0-E80.29'
    # [1] "F50.0-F50.1" in GM, in CM F50.00-F50.02 ##
    icd_dt[icd10 == 'F50.0-F50.1', 'icd10'] <- 'F50.00-F50.02'
    # [1] "S00-T98" in CM T88, in GM it is T98
    icd_dt[icd10 == 'S00-T98', 'icd10'] <- 'S00-T88'
    # [1] "E27.1-E27.4" E27.1-E2749 ##
    icd_dt[icd10 == 'E27.1-E27.4', 'icd10'] <- 'E27.1-E27.49'
  }

  icd_w_boarders <- unique(icd_dt$icd10[grep(icd_dt$icd10, pattern = '-')])

  # if classes
  if (length(icd_w_boarders) > 0) {

    for (i in 1:length(icd_w_boarders)) {
      # icd by icd
      x <- as.character(icd_w_boarders[i])
      rs <- try(expanded_icds <- expand_icd_2017_start_stop(icd10cm_2017_cmprsd, icd_pair = x))
      if ("try-error" %in% class(rs)) {
        expanded_icds <- x
        print(x) }

      icd_dt[icd10 == x, "expanded_icd"] <- paste(expanded_icds, collapse=', ')
    }

  }

  # expand the whole table
  icd_dt <- cSplit(icd_dt, splitCols = 'expanded_icd', sep = ', ', direction = 'long')
  return(icd_dt)
}

# ----- main function
etl_ttd <- function(create_tables = FALSE) {
  # required packages
  devtools::load_all('bips_devel/rgp')
  require(RSQLite)
  require(splitstackshape)

  sqlitedb <- paste(system.file("extdata", package = "RGP"), "ttd/TTD.db", sep = "/")
  if (!file.exists(sqlitedb)) {
    # create the TTD database if not created
    setwd(paste(system.file("extdata", package = "RGP"), "ttd", sep = "/"))
    system(paste(system.file("extdata", package = "RGP"), "ttd/load_data.sh", sep = "/"))
    setwd('~')
  }

  ttd_con <- connect_ttd(sqlitedb)
  dbListTables(ttd_con)

  if (isTRUE(create_tables)) {
    # make and modify the tables required for the matrices
    ttd_target2kegg()
    ttd_drug2atc2path()
    ttd_ds2icd2path()
    ttd_dr2ds()
  }

  # make the drug_path matrix
  ttd_drug_target_m <- make_drug_path_matrix(ttd_con)
  save(ttd_drug_target_m, file = paste(system.file("extdata", package = "RGP"), "ttd/ttd_drug_target_m.rda", sep = "/"))

  # make drug_disease matrix
  ttd_drug_disease_m <- make_drug_disease_matrix(ttd_con)
  save(ttd_drug_disease_m, file = paste(system.file("extdata", package = "RGP"), "ttd/ttd_drug_disease_m.rda", sep = "/"))

  # make target_disease matrix
  ttd_disease_target_m <- make_disease_target_matrix(ttd_con)
  save(ttd_disease_target_m, file = paste(system.file("extdata", package = "RGP"), "ttd/ttd_disease_target_m.rda", sep = "/"))

  # double check
  # names(ttd_drug_disease_m)[!(names(ttd_drug_disease_m) %in% atc_with_disease)]
  # names(ttd_drug_target_m)[!(names(ttd_drug_target_m) %in% names(ttd_drug_disease_m))]

  dbListTables(ttd_con)
  disconnect_ttd(ttd_con)
}
