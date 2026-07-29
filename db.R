library(DBI)
library(RSQLite)
library(rhandsontable)

EPITO_DB <- "epitoscope_results.sqlite"   # lives next to the app; swap for a server DSN later

.db_con <- function(path = EPITO_DB) DBI::dbConnect(RSQLite::SQLite(), path)

# ---------save analysis--------------
save_analysis_to_db <- function(lst, meta_table, quantity_cols, col_map = NULL,
                                software = "", submitted_by = "", description = "",
                                path = EPITO_DB) {
  con <- .db_con(path); on.exit(DBI::dbDisconnect(con))
  aid <- paste0("A_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  
  ## 1) analyses
  analyses <- data.frame(analysis_id = aid, timestamp = as.character(Sys.time()),
                         submitted_by = submitted_by, description = description,
                         n_samples = length(lst), stringsAsFactors = FALSE)   # no `software` column here
  
  ## map annotation measurement -> data column, so peptides & metadata share one key
  meas2col <- if (!is.null(col_map))
    vapply(col_map, `[[`, character(1), "col") else NULL   # names = measurement, values = data col
  
  ## 2) sample_metadata (long) — keyed on the data column so it joins peptides
  mt <- meta_table
  mt$data_col <- if (!is.null(meas2col)) meas2col[mt$measurement] else mt$measurement
  idc <- intersect(c("Sample","measurement","data_col"), names(mt))
  sample_metadata <- tidyr::pivot_longer(mt, cols = setdiff(names(mt), idc),
                                         names_to = "field_name", values_to = "field_value")
  sample_metadata <- data.frame(analysis_id = aid,
                                Sample      = as.character(sample_metadata$Sample),
                                measurement = as.character(sample_metadata$data_col),          # join key = data column
                                field_name  = sample_metadata$field_name,
                                field_value = as.character(sample_metadata$field_value), stringsAsFactors = FALSE)
  
  ## 3) peptides — pivot measurement columns to long; measurement = data column name
  pep <- c("STRIPPED","PEPTIDE","PROTEIN","LENGTH","MASS","MZ","RT","K0","PPM","SCORE","CHARGE","PTM")
  peptides <- dplyr::bind_rows(lapply(names(lst), function(s) {
    d  <- lst[[s]]; mc <- intersect(quantity_cols, names(d)); if (!length(mc)) return(NULL)
    tidyr::pivot_longer(d[, c(intersect(pep, names(d)), mc), drop = FALSE],
                        cols = dplyr::all_of(mc), names_to = "measurement", values_to = "quantity") |>
      dplyr::mutate(analysis_id = aid, Sample = s)
  }))
  
  DBI::dbWriteTable(con, "analyses",        analyses,        append = TRUE)
  DBI::dbWriteTable(con, "sample_metadata", sample_metadata, append = TRUE)
  DBI::dbWriteTable(con, "peptides",        peptides,        append = TRUE)
  aid
}

build_meta_table <- function(lst, quantity_cols, annotation_df, col_map, software_map) {
  has_ann <- !is.null(annotation_df) && isTRUE(attr(annotation_df, "has_measurement")) && !is.null(col_map)
  if (has_ann) {
    ann <- annotation_df
    cond_cols <- attr(annotation_df, "condition_cols"); if (is.null(cond_cols)) cond_cols <- character(0)
    rep_cols  <- intersect(c("biological_replicate","technical_replicate"), names(ann))
    tbl <- ann[, unique(c("name","measurement", cond_cols, rep_cols)), drop = FALSE]
    names(tbl)[names(tbl) == "name"] <- "Sample"
  } else {
    tbl <- dplyr::bind_rows(lapply(names(lst), function(s) {
      mc <- intersect(quantity_cols, names(lst[[s]]))
      if (!length(mc)) return(NULL)
      data.frame(Sample = s, measurement = mc, stringsAsFactors = FALSE)
    }))
  }
  # always-present editable columns, pre-filled
  tbl$software   <- unname(software_map[tbl$Sample])   # per-sample software
  tbl$instrument <- ""
  tbl[] <- lapply(tbl, as.character)
  tbl
}

list_analyses <- function(path = EPITO_DB) {
  if (!file.exists(path)) return(data.frame())
  con <- .db_con(path); on.exit(DBI::dbDisconnect(con))
  if (!DBI::dbExistsTable(con, "analyses")) return(data.frame())
  DBI::dbGetQuery(con, "SELECT analysis_id, timestamp, submitted_by, description, n_samples
                        FROM analyses ORDER BY timestamp DESC")
}

load_run_data <- function(run_id, path = EPITO_DB) {
  con <- .db_con(path); on.exit(DBI::dbDisconnect(con))
  r <- DBI::dbGetQuery(con, "SELECT data_blob FROM run_data WHERE run_id = ?", params = list(run_id))
  if (!nrow(r)) return(NULL)
  unserialize(memDecompress(r$data_blob[[1]], type = "gzip"))   # full data list, all columns
}

# ---------browse database--------------
db_tables <- function(path = EPITO_DB) {
  if (!file.exists(path)) return(character(0))
  con <- .db_con(path); on.exit(DBI::dbDisconnect(con))
  DBI::dbListTables(con)
}

db_schema <- function(path = EPITO_DB) {
  con <- .db_con(path); on.exit(DBI::dbDisconnect(con))
  tabs <- DBI::dbListTables(con)
  setNames(lapply(tabs, function(t) DBI::dbListFields(con, t)), tabs)   # table -> columns
}