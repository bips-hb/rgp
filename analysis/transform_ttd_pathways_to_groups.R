# Define project directory
PRJ_DIR 
PKG_DIR <- 
DATA_DIR <- 

# read ttd data
load(file.path(PRJ_DIR, PKG_DIR, DATA_DIR, "ttd_clean_data", "ttd_drug_target_m.rda"))
load(file.path(PRJ_DIR, PKG_DIR, DATA_DIR, "ttd_clean_data", "ttd_disease_target_m.rda"))

# add edoxaban and dabigatran: missing anticoag
# check
# anticoag <- c("B01AF01", "B01AF02", "B01AF03", "B01AE07", "B01AA03", "B01AA04")
# names(anticoag) <- c("Riva", "Api", "Edox", "Dabi", "Warf", "Phenpro")
# anticoag %in% colnames(ttd_drug_target_m)
# anticoagpath <- c("hsa04610", "hsa04080")
# anticoagpath %in% rownames(ttd_drug_target_m)

# From KEGG
# Edoxaban	(2011/2015)	ATC_WIDO = B01AF03 [FT: hsa04610]
# Dabigatran 	(2010)	ATC_WIDO = B01AE07 [FT: hsa04610, hsa04080]

n <- subset(ttd_drug_target_m, select = c("B01AF01", "B01AF02"))
names(n) <- c("B01AF03", "B01AE07")
n["hsa04080", "B01AE07"] <- 1
n["hsa04810", "B01AE07"] <- 1
ttd_drug_target_m <- cbind(ttd_drug_target_m, n)
ttd_drug_target_m <- ttd_drug_target_m[, order(colnames(ttd_drug_target_m))]

# remove . from icds (clean like gepard)
names(ttd_disease_target_m) <- gsub(names(ttd_disease_target_m), pattern = '.',
  fixed = TRUE, replacement = '')

# merge disease and drug targets
ttd_drug_target_m <- data.table(ttd_drug_target_m, keep.rownames = TRUE)
ttd_disease_target_m <- data.table(ttd_disease_target_m, keep.rownames = TRUE)
setkey(ttd_disease_target_m, rn)
setkey(ttd_drug_target_m, rn)

all_targets <- merge(ttd_disease_target_m, ttd_drug_target_m, all = TRUE)
all_targets[is.na(all_targets)] <- 0
# all_targets[, "NA"]

# transform the data table into matrix into a list of vectors
all_targets_m <- data.matrix(all_targets[, -1]) # remove the rn column
rownames(all_targets_m) <- all_targets$rn
all_targets_l <- matrix_to_list(all_targets_m)

saveRDS(all_targets_l, file.path(PRJ_DIR, PKG_DIR, DATA_DIR, "ttd_clean_data", "ttd_pathways_doac.Rds"))
