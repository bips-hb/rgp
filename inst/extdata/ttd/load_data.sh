#!/bin/bash

readonly DB=TTD.db
readonly SQLITE="sqlite3 $DB"

# define functions
download_data() {
  # At the time: latest version 6.1.01 (2017.10.04) from the latest paper (2018)
  wget https://db.idrblab.org/ttd/sites/default/files/ttd_database/P1-01-TTD_download.txt
  wget https://db.idrblab.org/ttd/sites/default/files/ttd_database/P1-02-TTD_crossmatching.txt
  wget https://db.idrblab.org/ttd/sites/default/files/ttd_database/P1-04-Drug_disease.txt
  wget https://db.idrblab.org/ttd/sites/default/files/ttd_database/P1-05-Target_disease.txt
  wget https://db.idrblab.org/ttd/sites/default/files/ttd_database/P4-02-Target-KEGGpathway_succ.txt
}

cut_first_n_lines() {
  # If the table already exists, the sqlite3 tool uses all the rows,
  # including the first row, in the CSV file as the actual data to import.
  # Therefore, you should delete the first row of the CSV file.
  local f="$1"
  local n="$2"
  tail -n +$n $f.txt > $f.tsv
}

prepare_files() {
  cut_first_n_lines P1-01-TTD_download 13
  cut_first_n_lines P1-02-TTD_crossmatching 13
  cut_first_n_lines P1-04-Drug_disease 14
  cut_first_n_lines P1-05-Target_disease 14
  cut_first_n_lines P4-02-Target-KEGGpathway_succ 14
}

drop_schema() {
  for t in ttd_target ttd_crossmatching ttd_drug_disease ttd_target_disease ttd_target_kegg; do
    $SQLITE "DROP TABLE IF EXISTS $t"
  done
}
prepare_schema() {
  $SQLITE "CREATE TABLE ttd_target (ttd_target_id varchar(6), target_key text, target_value text);"
  $SQLITE "CREATE INDEX ttd_target_ix ON ttd_target (ttd_target_id)"
  $SQLITE "CREATE INDEX target_value_ix ON ttd_target (ttd_target_id, target_key, target_value)"

  $SQLITE "CREATE TABLE ttd_crossmatching (ttd_drug_id varchar(6), col_id int, drug_key text, drug_value text);"
  $SQLITE "CREATE INDEX ttd_crossmatching_ix ON ttd_crossmatching (ttd_drug_id)"
  $SQLITE "CREATE INDEX drug_value_ix ON ttd_crossmatching (drug_key, drug_value)"

  $SQLITE "CREATE TABLE ttd_drug_disease (ttd_drug_id varchar(6), lnm text, indication text, icd9 text, icd10 text);"
  $SQLITE "CREATE INDEX ttd_drug_disease_ix ON ttd_drug_disease (ttd_drug_id)"
  $SQLITE "CREATE INDEX ttd_drug_icd10_ix ON ttd_drug_disease (ttd_drug_id, indication, icd10)"

  $SQLITE "CREATE TABLE ttd_target_disease (ttd_target_id varchar(6), target_name text, indication text, icd9 text, icd10 text);"
  $SQLITE "CREATE INDEX ttd_target_disease_ix ON ttd_target_disease (ttd_target_id, indication, icd10)"

  $SQLITE "CREATE TABLE ttd_target_kegg (ttd_target_id varchar(6), kegg_pathway_id varchar(8), kegg_pathway_name text);"
  $SQLITE "CREATE INDEX ttd_target_kegg_ix ON ttd_target_kegg (ttd_target_id, kegg_pathway_id)"
}

load_file() {
  # can be using RSQLite: https://stackoverflow.com/questions/25194568/how-to-import-tab-delimited-data-to-sqlite-using-rsqlite
  local f="$1"
  local table="$2"
  $SQLITE <<EOF
.mode tabs
.import $f.tsv $table
EOF
}

load_files() {
  load_file P1-01-TTD_download ttd_target
  load_file P1-02-TTD_crossmatching ttd_crossmatching
  load_file P1-04-Drug_disease ttd_drug_disease
  load_file P1-05-Target_disease ttd_target_disease
  load_file P4-02-Target-KEGGpathway_succ ttd_target_kegg
}

clean_up_files() {
  rm P1-01-TTD_download.tsv P1-02-TTD_crossmatching.tsv P1-04-Drug_disease.tsv P1-05-Target_disease.tsv P4-02-Target-KEGGpathway_succ.tsv
}

# execute
prepare_files
drop_schema
prepare_schema
load_files
clean_up_files
