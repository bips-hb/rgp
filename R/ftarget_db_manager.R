#' @title Fetch KEGG Entries using KEGGREST
# KEGGREST: https://bioconductor.org/packages/release/bioc/vignettes/KEGGREST/inst/doc/KEGGREST-vignette.html
# link kegg databases: can work for drug-atc, disease-pathway and drug-disease, and drug-gene, and disease-gene, and pathway-gene.
#' @export
source2target <- function(kegg_source, kegg_target, file_name = NULL) {
	kegg_linked = keggLink(target = kegg_target, source = kegg_source)
	kegg_linked = data.table(kegg_source = names(kegg_linked), kegg_target = kegg_linked, keep.rownames = F)
	names(kegg_linked) = c(eval(kegg_source), eval(kegg_target))
	if (!missing(file_name)) {
	  write.table(kegg_linked, file = file_name, quote = F, sep = "\t", row.names = F)
	}
	return(kegg_linked)
}

#' @title Merge two data tables with the same key, usually kegg id
#' @importFrom utils write.table
#' @export
kegg_merge <- function(kegg_table_1, kegg_table_2, key_column, merged_file_name) {
	kegg_merged = merge(x = kegg_table_1, y = kegg_table_2, all = T)
	setkeyv(kegg_merged, key_column)
	if(!missing(merged_file_name)) {
		write.table(kegg_merged, file = merged_file_name, quote = F, sep = "\t", row.names = F)
	}
	return(kegg_merged)
}

#' @importFrom utils read.table
#' @export
read_kegg_merged <- function(merged_file_name, key_column) {
  kegg_merged = data.table(read.table(merged_file_name, header = T, sep = "\t"), key = key_column)
  return(kegg_merged)
}

#' @importFrom utils write.table
#' @export
drug2atc <- function(drug_atc_file, atc_level_4 = F) {
  drug_atc <- source2target(kegg_source = "drug", kegg_target = "atc", file_name = drug_atc_file)
  setkey(drug_atc, "drug")
  drug_atc[, atc := substr(atc, 5, nchar(atc))][] #remove 'atc:'
  if (isTRUE(atc_level_4)) {
    drug_atc[, atc := substr(atc, 1, nchar(atc)-2)][]
  }
  if (!missing(drug_atc_file)) {
    write.table(drug_atc, file = drug_atc_file, quote = F, sep = "\t", row.names = F)
  }
  return(drug_atc)
}

# ICD-10 International Classification of Diseases
## http://www.genome.jp/kegg-bin/get_htext?br08410.keg
# Human Diseases in ICD-10 Classification
## http://www.genome.jp/kegg-bin/get_htext?br08403
## using br08403+ is useless as I rely on ICD-10 classification because it is more detailed than the human disease one
#' @export
get_icd_by_ds <- function(ds_id) {
  # it is now 2019 and it is ICD-11!!!
	ds_entry_dblinks <- try(keggGet(dbentries = paste0("br08403_10+", ds_id))[[1]]$DBLINK) # br08403_10+ or br08410+
	ds_icd_entry = as.character(ds_entry_dblinks[grepl(x = ds_entry_dblinks, pattern = "ICD-10:")])
	if (length(ds_icd_entry) > 0) {
		if (length(ds_icd_entry) > 1) {
			ds_icd_entry = paste(unlist(ds_icd_entry), collapse = ' ')
		}
		ds_icd_entry = gsub(ds_icd_entry, pattern = 'ICD-10: ', replacement = '')
	}
	else {
		ds_icd_entry = NA
	}
	return(ds_icd_entry)
}

#' @title link disease to icd
#' @importFrom stats na.omit
#' @importFrom utils write.table
#' @export
ds2icd <- function(ds_list = NULL, file_name = NULL, icd_level_C = F) {
  if (missing(ds_list)){
    # compile diseases from pathway, drug and gene
    ds_pathway = source2target(kegg_source = "ds", kegg_target = "pathway")
    ds_pathway = ds_pathway[grepl(pattern = "hsa", x = ds_pathway$pathway, fixed = T), ]

    ds_gene = source2target(kegg_source = "ds", kegg_target = "hsa")
    ds_gene = ds_gene[grepl(pattern = "hsa", x = ds_gene$hsa, fixed = T), ]

    ds_drug = source2target(kegg_source = "ds", kegg_target = "drug")

    ds_list = c(ds_pathway$ds, ds_gene$ds, ds_drug$ds)
  }

	ds_icd = data.table(ds = unique(ds_list), key = "ds")
	# get ICDs
	ds_icd[, icd := sapply(ds, get_icd_by_ds, simplify = TRUE)][]
	ds_icd[, icd := gsub(x = icd, pattern = " ", replacement = ",")][]
	ds_icd = ds_icd[, .(ds, icd=unlist(strsplit(icd, split = ","))), by = seq_len(nrow(ds_icd))]
	ds_icd = na.omit(ds_icd[, 2:3])
	if (isTRUE(icd_level_C)) {
	  ds_icd[, icd := gsub(x = icd, pattern = "\\.[0-9]*", replacement = "")][]
	}

	if (!missing(file_name)) {
	  write.table(ds_icd, file = file_name, quote = F, sep = "\t", row.names = F)
	}
	return(ds_icd)
}

#' @title Find disease target: link disease to target
#' @export
ds2icd2target <- function(ds_icd_target_file, kegg_target) {
	ds_target = source2target(kegg_source = "ds", kegg_target = kegg_target)
	# remove map
	ds_target = ds_target[grepl(pattern = "hsa", x = ds_target[[kegg_target]], fixed = T), ]
	# ds_icd = ds2icd(ds_target$ds)
	ds_icd = read_kegg_merged(system.file("extdata", "ds2icd.txt", package = "RGP"), "ds")
	ds_icd_target = kegg_merge(kegg_table_1 = ds_icd, kegg_table_2 = ds_target,
	                           key_column = "icd", merged_file_name = ds_icd_target_file)
	return(ds_icd_target)
}

#' @export
get_path_by_dr <- function(dr_id) {
  dr_path_entry = try(keggGet(dbentries = dr_id)[[1]]$TARGET)
  if ("PATHWAY" %in% names(dr_path_entry)) {
    dr_path_entry = as.character(dr_path_entry$PATHWAY)
    dr_path_entry = sub("\\(.*", "", dr_path_entry)
    dr_path_entry = paste(paste0('path:', unlist(dr_path_entry)), collapse = ' ')
  }
  else {
    dr_path_entry = NA
  }
  return(dr_path_entry)
}

#' @title link drug to pathway
#' @importFrom stats na.omit
#' @importFrom utils write.table
#' @export
dr2path <- function(dr_list, file_name = NULL) {
  dr_path = data.table(drug = unique(dr_list), key = "drug")
  dr_path[, pathway := sapply(drug, get_path_by_dr, simplify = TRUE)][]

  dr_path[, pathway := gsub(x = pathway, pattern = " ", replacement = ",")][]
  dr_path = dr_path[, .(drug, pathway=unlist(strsplit(pathway, split = ","))), by = seq_len(nrow(dr_path))]
  dr_path = na.omit(dr_path[, 2:3])

  if (!missing(file_name)) {
    write.table(dr_path, file = file_name, quote = F, sep = "\t", row.names = F)
  }
  return(dr_path)
}

#' @title find drug target: link drug to target ("pathway" or gene: "hsa") + trim first 4 characters (atc:)
# can be br08303_target.keg
#' @export
drug2atc2target <- function(drug_atc_target_file, kegg_target) {
  drug_atc = read_kegg_merged(system.file("extdata", "drug2atc.txt", package = "RGP"), "drug")

  if (kegg_target == "pathway") {
    # drug_target = dr2path(atc_drug$drug)
     if (file.exists(system.file("extdata", "drug2path.txt", package = "RGP"))) {
       drug_target = read_kegg_merged(system.file("extdata", "drug2path.txt", package = "RGP"), "drug")
     } else {
       drug_target <- dr2path(dr_list = drug_atc$drug, file_name = system.file("extdata", "drug2path.txt", package = "RGP"))
     }
    drug_target = drug_target[grepl(x = pathway, pattern = ":hsa", fixed = T), ]
  }
  if (kegg_target == "hsa") {
    drug_target = source2target(kegg_source = "drug", kegg_target = kegg_target)
  }

	drug_atc_target = kegg_merge(kegg_table_1 = drug_atc,
	  kegg_table_2 = drug_target,
	  key_column = "atc",
	  merged_file_name = drug_atc_target_file
	)
	return(drug_atc_target)
}

#' @title get icd descriptions (all names)
#' @export
parse_br08410 <- function(br08410_file, icd_level) {
	br08410_file_con = file(br08410_file, open = "rt", blocking = F)
	icd = readLines(br08410_file_con, encoding = "UTF-8")
	close(br08410_file_con)
	icd = iconv(icd, "UTF-8", "UTF-8", sub = '') #grep -axv '.*' br08410.keg
	icd = gsub(x = icd, pattern = "^!", replacement = "#")
	icd_at_level = icd[grepl(x = icd, pattern = paste0("^", icd_level), fixed = F, perl = T)]
	icd_at_level = sub(icd_at_level, pattern = paste0("^", icd_level, "\\s+"), fixed = F, replacement = '')
	icd_at_level_dt = data.table(icd_at_level)
	setDT(icd_at_level_dt)[, c("icd_id", "icd_name") := tstrsplit(icd_at_level, "	", type.convert = TRUE, fixed = TRUE)]
	icd_at_level_dt[,icd_at_level:=NULL][]
	setkey(icd_at_level_dt, icd_id)
	return(icd_at_level_dt)
}

#' @importFrom stats na.omit
#' @export
create_kegg_drug_ds_dt <- function(icd_level_C, atc_level_4, target_type) {
	# https://stackoverflow.com/questions/15347282/split-delimited-strings-in-a-column-and-insert-as-new-rows
	if (target_type == "gene") {
		drug_file_name = "drug2atc2gene.txt"
	}
	if (target_type == "pathway") {
		drug_file_name = "drug2atc2path.txt"
	}
	drug_atc_target = read_kegg_merged(system.file("extdata", drug_file_name, package = "PGE"), key_column = "atc")
	drug_atc_target = na.omit(drug_atc_target[, 2:3])
	if (dr_path_entry == "path:ko00550  Peptidoglycan biosynthesis") {
	  dr_path_entry = NA
	}
	if (isTRUE(atc_level_4)) {
	  drug_atc_target[, atc:=substr(atc, 1, 5)][]
	}

	ds_icd_target = prep_ds_icd_target(icd_level_C, target_type)
	ds_icd_target = na.omit(ds_icd_target[, 2:3])

	drug_disease_target = rbindlist(list(drug_atc_target, ds_icd_target), use.names = F)
	drug_disease_target = drug_disease_target[, c(2,1)]
	names(drug_disease_target) = c("element", target_type)
	return(drug_disease_target)
}

#' @export
create_kegg_drug_ds_matrix <- function(icd_level_C, target_type) {
	drug_disease_target = create_kegg_drug_ds_dt(icd_level_C = icd_level_C, target_type = target_type)
	drug_disease_target_boolean_matrix = create_boolean_matrix(drug_disease_target)
	return(drug_disease_target_boolean_matrix)
}

#' @export
kegg_subset <- function(kegg_merged, id_list, id_col_name) {
	kegg_subsetted = kegg_merged[grepl(pattern = paste(id_list, collapse='|'), ignore.case = T, x = kegg_merged[[id_col_name]])]
	return(kegg_subsetted)
}

#' @export
prep_ds_icd_target <- function(icd_level_C, target_type) {
	if (target_type == "gene") {
		ds_file_name = "ds2icd2gene.txt"
	}
	if (target_type == "pathway") {
		ds_file_name = "ds2icd2path.txt"
	}
	ds_icd_target = read_kegg_merged(system.file("extdata", ds_file_name, package = "PGE"), "icd")
	ds_icd_target[, icd := gsub(x = icd, pattern = " ", replacement = ",")][]
	if (target_type == "gene") {
		ds_icd_target = ds_icd_target[, .(ds, icd=unlist(strsplit(icd, split = ",")), hsa), by = seq_len(nrow(ds_icd_target))]
	}
	if (target_type == "pathway") {
		ds_icd_target = ds_icd_target[, .(ds, icd=unlist(strsplit(icd, split = ",")), pathway), by = seq_len(nrow(ds_icd_target))]
	}
	ds_icd_target = ds_icd_target[, c(2:4)]
	if (isTRUE(icd_level_C)) {
		ds_icd_target[, icd := gsub(x = icd, pattern = "\\.[0-9]*", replacement = "")][]
	}
	return(ds_icd_target)
}

# ------------------------------- KEGG start with pathway ---------------------
get_all_hsa_path <- function() {
  # list all KEGG human pathways
  pathway_list <- keggList(database = "pathway", organism = "hsa")
  pathway_list <- names(pathway_list)
  return(pathway_list)
}

# TODO: functions required to handle pathways diseases and drugs and map back to ATC and ICD coding

# ------------------------------- KEGG ds2dr -----------------------------------
#
ds2dr_kegg <- function(ds_dr_file) {
  # get disease-drug pairs from KEGG (based on genetics probably)
  if (file.exists(ds_dr_file)) {
    ds_dr <- read_kegg_merged(paste(system.file("extdata", package = "RGP"), "ds2dr.txt", sep = "/"), "ds")
  } else {
    ds_dr = source2target(kegg_source = "ds", kegg_target = "drug", file_name = ds_dr_file)
  }
  # get the atc code of each drug
  drug_atc_map <- read_kegg_merged(paste(system.file("extdata", package = "RGP"), "drug2atc.txt", sep = "/"), "drug")
  ds_dr[drug_atc_map, atc := atc, on = "drug"]

  # get the icd code of each disease
  ds_icd_map <- read_kegg_merged(paste(system.file("extdata", package = "RGP"), "ds2icd.txt", sep = "/"), "ds")
  ds_dr[ds_icd_map, icd := icd, on = "ds"]

  # discard entries with icd/atc NA
  ds_dr <- na.omit(ds_dr)

  # split multiple icds into multiple rows
  # ds_dr <- cSplit(ds_dr, "icd", sep = " ", direction = "long")

  # discard kegg ds and dr ids
  write.table(ds_dr[, c(3, 4)],
              file = paste(system.file("extdata", package = "RGP"), "icd2atc.txt", sep = "/"),
              quote = F, sep = "\t", row.names = F)
  return(ds_dr[, c(3, 4)])
}


# ------------------------------- STITCH ---------------------------------------
# Fetch STITCH scores

#' @export
get_stitch_chemical_scores_from_tsv <- function(stitch_dir) {
	if(missing(stitch_dir)) {
		stitch_dir = '/bips/gruppen/pv-monitor/Auswertung/Task13/stitch/stitch_tsv/'
	}
	stitch_chemicals = "chemical.sources.v5.0.tsv"
	stitch_scores = "chemical_chemical.links.detailed.v5.0.tsv"

	stitch_chemical_ids = fread(sprintf('grep "KEGG\\|ATC" %s | grep -v "^#"', paste0(stitch_dir, stitch_chemicals)),
	  header = F, col.names = c("chemical", "alias", "source", "cross_id"))
	stitch_chemical_ids_merged = merge(stitch_chemical_ids[source=="ATC"],
	  stitch_chemical_ids[(source=="KEGG" & startsWith(cross_id, "D"))],
	  by = c("chemical", "alias"), all = T, suffixes = c(".atc", ".kegg"))
	stitch_chemical_ids_merged[,c("source.atc", "source.kegg"):=NULL][]

	stitch_chemical_scores = fread(paste0(stitch_dir, stitch_scores), header = T)

	stitch_chemical_scores1 = merge(stitch_chemical_ids_merged,
	  stitch_chemical_scores, by.x = "chemical", by.y = "chemical1", all.x = T)
	stitch_chemical_scores = merge(stitch_chemical_ids_merged,
	  stitch_chemical_scores1,
	  by.x = "chemical",
	  by.y = "chemical2",
	  all.x = T, suffixes = c(".1", ".2"))

	return(stitch_chemical_scores)
}

#' @export
create_stitch_score_matrix <- function(chemicals_data_table,
                                       score_type = "combined_score",
                                       score_cutoff = 400, atc = T) {
  if(missing(chemicals_data_table)) {
    chemicals_data_table = get_stitch_chemical_scores_from_tsv()
  }
	chemicals_data_table = chemicals_data_table[(!(is.na(cross_id.atc.1) | is.na(cross_id.atc.2))
	  & chemicals_data_table[[score_type]] >= score_cutoff),
	  c("cross_id.atc.1", "cross_id.atc.2", eval(score_type)), with=F]
	if (!isTRUE(atc)){
		chemicals_data_table = chemicals_data_table[(!(is.na(cross_id.kegg.1) | is.na(cross_id.kegg.2))
		  & chemicals_data_table[[score_type]] >= score_cutoff),
		  c("cross_id.kegg.1", "cross_id.kegg.2", eval(score_type)), with=F]
	}
	stitch_square_matrix = create_square_matrix(chemicals_data_table)
	return(stitch_square_matrix)
}

# ------------- sources of ICD -------------------------------------------------
# DIMDI: ICD-10-WHO-2016: https://www.dimdi.de/dynamic/.downloads/klassifikationen/icd-10-who/version2016/x1wmt2016.zip
# icd.data: ICD-10-CM-2016: icd.data::icd10cm2016 or icd::icd10_sources$`2016`$dx_flat
# KEGG: br08403_10
# WHO: ICD-10-WHO-2016: http://apps.who.int/classifications/apps/icd/ClassificationDownload/DLArea/Download.aspx
# Gitbug: ICD-10-CM-2018: https://github.com/kamillamagna/ICD-10-CSV
# ICD-10-CM-2017: https://www.cms.gov/Medicare/Coding/ICD10/2017-ICD-10-CM-and-GEMs.html
# ICD-10-CM-2017: https://github.com/ellessenne/comorbidity rda
# Parser: https://github.com/ultrah/ClaMLParser
# check if there are children in either the WHO-2016 or CM-2016
# using icd.data::
# children('I70', short_code = TRUE, defined = TRUE)
# mapping: http://www.nber.org/icd-10-cm-mappings/2017/
# start with WHO
#' @export
expand_icd_2016 <- function(icd10who2016_file = paste(system.file("extdata", package = "RGP"),
                                                 "icd10who2016syst_kodes.txt", sep = "/"),
                       icd_dt) {
  # initialize new column
  icd_dt$expanded_icd <- icd_dt$icd

  # read WHO ICDs file
  icd10who2016 <- data.table(read.csv(icd10who2016_file, sep = ";", header = F))

  # check if first 10 icds have a child
  for (i in 1:length(unique(icd_dt$icd))) {
    # icd by icd
    x <- icd_dt$icd[i]
    # get children
    expanded_icds <- icd10who2016[(V1 == 4 & grepl(x = icd10who2016$V6, pattern = x, fixed = T)), V6]
    # inflate
    icd_dt[icd == x, "expanded_icd"] <- ifelse(length(expanded_icds) > 0,
                                                  paste(expanded_icds, collapse=', '),
                                                  x)
  }

  # expand the whole table
  icd_dt <- cSplit(icd_dt, splitCols = 'expanded_icd', sep = ', ', direction = 'long')
  return(icd_dt)
}
