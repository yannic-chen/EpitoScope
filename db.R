library(DBI)
library(RSQLite)
library(rhandsontable)

# metadata fields that should ALWAYS be offered (edit to taste)
BASE_META_FIELDS <- c("instrument", "biological source", "cell line","condition",
                      "biological replicate", "technical replicate", "species", "reference database","notes")

if (!exists("EPITO_DB")) EPITO_DB <- "epitoscope_results.sqlite"


.db_con <- function(path = EPITO_DB) DBI::dbConnect(RSQLite::SQLite(), path)

# ---------save analysis--------------
save_analysis_to_db <- function(lst, meta_table, quantity_cols, col_map = NULL,
                                spectra_cols = NULL, data_info = NULL, mod_map = NULL,
                                software = "", submitted_by = "", description = "",
                                allow_duplicate = FALSE, path = EPITO_DB) {
  
  con <- .db_con(path); on.exit(DBI::dbDisconnect(con))
  aid <- paste0("A_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  
  ## ---- split binding predictions out of the analysis data ----
  allele_cols <- unique(unlist(lapply(lst, mhc_allele_cols)))
  
  if (length(allele_cols)) {
    pred_long <- dplyr::bind_rows(lapply(lst, function(df) {
      ac <- intersect(allele_cols, colnames(df))
      df[, c("STRIPPED", ac), drop = FALSE]
    })) %>%
      tidyr::pivot_longer(cols = dplyr::any_of(allele_cols),
                          names_to = "allele", values_to = "rank") %>%
      dplyr::filter(!is.na(rank)) %>%
      dplyr::transmute(peptide = STRIPPED,
                       allele  = sub(paste0("^", MHC_PREFIX), "", allele),
                       rank    = as.numeric(rank)) %>%
      dplyr::distinct(peptide, allele, .keep_all = TRUE)
    
    # remove allele columns from EVERY sample frame
    lst <- lapply(lst, function(df) df[, setdiff(names(df), allele_cols), drop = FALSE])
  } else {
    pred_long <- NULL
  }
  
  ## ---- fingerprint the analysis WITHOUT predictions ---
  fp <- analysis_fingerprint(lst, quantity_cols, spectra_cols) #compute hash

  ##---- analysis row ----
  analyses <- data.frame(analysis_id = aid, timestamp = as.character(Sys.time()),
                         submitted_by = submitted_by, description = description,
                         n_samples = length(lst), content_hash = fp, stringsAsFactors = FALSE)
  
  meas2col <- if (!is.null(col_map)) vapply(col_map, `[[`, character(1), "col") else NULL
  
  ##---- metadata (long) ----
  mt <- meta_table
  mt$data_col <- if (!is.null(meas2col) && "measurement" %in% names(mt))
    meas2col[mt$measurement] else mt$measurement
  idc <- intersect(c("Sample", "measurement", "data_col"), names(mt))
  sample_metadata <- tidyr::pivot_longer(mt, cols = setdiff(names(mt), idc),
                                         names_to = "field_name", values_to = "field_value")
  sample_metadata <- data.frame(
    analysis_id = aid,
    Sample      = as.character(sample_metadata$Sample),
    measurement = as.character(sample_metadata$data_col),
    field_name  = as.character(sample_metadata$field_name),
    field_value = as.character(sample_metadata$field_value),
    stringsAsFactors = FALSE)
  keep <- !is.na(sample_metadata$field_value) & nzchar(trimws(sample_metadata$field_value))
  sample_metadata <- sample_metadata[keep, , drop = FALSE]
  
  ##---- check duplicate ----
  if (!allow_duplicate &&
      "analyses" %in% DBI::dbListTables(con) &&
      "content_hash" %in% DBI::dbListFields(con, "analyses")) {
    dup <- DBI::dbGetQuery(con,
                           "SELECT analysis_id, timestamp FROM analyses WHERE content_hash = ? LIMIT 1",
                           params = list(fp))
    if (nrow(dup))
      stop(structure(
        class = c("epito_duplicate", "error", "condition"),
        list(message = paste0("Identical to ", dup$analysis_id[1], " (saved ", dup$timestamp[1], ")"),
             call = NULL, duplicate_of = dup$analysis_id[1], timestamp = dup$timestamp[1])))
  }
  
  ##---- normalized raw data (identity + quantities, normalized schema) ----
  pep_keys <- setdiff(unique(unlist(lapply(column_schema, names))), c("QUANTITY", "SPECTRA"))
  
  DBI::dbBegin(con)
  tryCatch({
    ensure_condition_terms(con)
    ensure_peptide_tables(con)
    
    # canonicalize vocab metadata values in place
    vrow <- vapply(sample_metadata$field_name, is_vocab_field, logical(1))
    if (any(vrow)) {
      sample_metadata$field_value[vrow] <- mapply(
        function(f, v) canonicalize_condition(con, f, v),
        sample_metadata$field_name[vrow], sample_metadata$field_value[vrow],
        USE.NAMES = FALSE)
    }
    
    # integer peptide_id, unique across analyses
    start_id <- if ("peptides" %in% DBI::dbListTables(con))
      DBI::dbGetQuery(con, "SELECT COALESCE(MAX(peptide_id),0) AS m FROM peptides")$m else 0
    pid <- as.integer(start_id)
    
    pep_rows <- list(); qty_rows <- list()
    for (s in names(lst)) {
      d  <- lst[[s]]
      mc <- intersect(quantity_cols, names(d)); if (!length(mc)) next
      sc <- intersect(spectra_cols,  names(d))          # <-- defines sc
      measurement_cols <- unique(c(mc, sc))             # <-- defines measurement_cols
      n  <- nrow(d)
      ids <- pid + seq_len(n); pid <- pid + n
      
      id_df <- d[, setdiff(names(d), measurement_cols), drop = FALSE]   # uses it here
      id_df$peptide_id  <- ids
      id_df$analysis_id <- aid
      id_df$Sample      <- s
      pep_rows[[s]] <- id_df
      
      q <- d[, mc, drop = FALSE]; q$peptide_id <- ids
      qty_rows[[s]] <- tidyr::pivot_longer(q, dplyr::all_of(mc),
                                           names_to = "measurement", values_to = "quantity")
    }
    peptides   <- dplyr::bind_rows(pep_rows)
    quantities <- dplyr::bind_rows(qty_rows)
    
    add_missing_columns(con, "peptides", peptides)     # auto-widen for new columns
    
    if ("analyses" %in% DBI::dbListTables(con)) add_missing_columns(con, "analyses", analyses)
    DBI::dbWriteTable(con, "analyses",        analyses,        append = TRUE)
    DBI::dbWriteTable(con, "sample_metadata", sample_metadata, append = TRUE)
    DBI::dbAppendTable(con, "peptides",   peptides)
    DBI::dbAppendTable(con, "quantities", quantities)
    cmap <- build_column_map_long(aid, data_info)
    if (!is.null(cmap)) DBI::dbWriteTable(con, "column_map", cmap, append = TRUE)
    mm <- build_mod_map_long(aid, mod_map)
    if (!is.null(mm)) DBI::dbWriteTable(con, "mod_map", mm, append = TRUE) 
    
    # ---- shared binding-prediction cache ----
    if (!is.null(pred_long) && nrow(pred_long)) {
      DBI::dbExecute(con, "CREATE TABLE IF NOT EXISTS predictions
                           (peptide TEXT, allele TEXT, rank REAL, PRIMARY KEY (peptide, allele))")
      DBI::dbWriteTable(con, "predictions_stg", pred_long, temporary = TRUE, overwrite = TRUE)
      DBI::dbExecute(con, "INSERT OR IGNORE INTO predictions
                           SELECT peptide, allele, rank FROM predictions_stg")
      DBI::dbExecute(con, "DROP TABLE predictions_stg")
    }
    
    DBI::dbCommit(con)
  }, error = function(e) { DBI::dbRollback(con); stop(e) })

  aid
}

sql_affinity <- function(x) {
  if (is.integer(x)) "INTEGER" else if (is.numeric(x)) "REAL" else "TEXT"
}

# add any df column the table doesn't have yet (old rows get NULL)
add_missing_columns <- function(con, table, df) {
  existing <- DBI::dbListFields(con, table)
  for (nm in setdiff(names(df), existing)) {
    DBI::dbExecute(con, sprintf('ALTER TABLE "%s" ADD COLUMN "%s" %s',
                                table, nm, sql_affinity(df[[nm]])))
  }
}

build_mod_map_long <- function(analysis_id, mod_map) {
  if (is.null(mod_map) || length(mod_map) == 0 ||
      (length(mod_map) == 1 && is.na(mod_map))) return(NULL)
  data.frame(analysis_id = analysis_id, token = names(mod_map),
             symbol = unname(as.character(mod_map)), stringsAsFactors = FALSE)
}
reconstruct_mod_map <- function(mm) {
  if (is.null(mm) || !nrow(mm)) return(NULL)
  setNames(mm$symbol, mm$token)
}

# base tables — the identity column set is DERIVED from column_schema, so
# extending column_schema automatically extends fresh databases.
ensure_peptide_tables <- function(con) {
  pep_keys  <- setdiff(unique(unlist(lapply(column_schema, names))),
                       c("QUANTITY", "SPECTRA"))          # per-peptide identity/spectral cols
  int_keys  <- c("LENGTH", "CHARGE")
  real_keys <- c("MASS", "MZ", "RT", "K0", "PPM", "SCORE")
  col_ddl <- vapply(pep_keys, function(k) {
    ty <- if (k %in% int_keys) "INTEGER" else if (k %in% real_keys) "REAL" else "TEXT"
    sprintf('"%s" %s', k, ty)
  }, character(1))
  
  DBI::dbExecute(con, sprintf('
    CREATE TABLE IF NOT EXISTS peptides (
      peptide_id INTEGER PRIMARY KEY,
      analysis_id TEXT, Sample TEXT, %s )', paste(col_ddl, collapse = ", ")))
  DBI::dbExecute(con, "CREATE TABLE IF NOT EXISTS quantities (
      peptide_id INTEGER, measurement TEXT, quantity REAL )")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_qty_pid ON quantities(peptide_id)")
}

# flatten data_info_r() (merged_info) -> long rows for storage
build_column_map_long <- function(analysis_id, data_info) {
  if (is.null(data_info) || !nrow(data_info)) return(NULL)
  sample_cols <- setdiff(names(data_info), "final_name")
  rows <- list()
  for (i in seq_len(nrow(data_info))) {
    fn <- as.character(data_info$final_name[i])
    for (sc in sample_cols) {
      vals <- data_info[[sc]][[i]]          # list-cell -> character vector of original names
      vals <- vals[!is.na(vals)]
      for (v in vals)
        rows[[length(rows) + 1]] <- data.frame(
          analysis_id = analysis_id, Sample = sc,
          final_name = fn, original_name = as.character(v),
          stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# Generate hash to "remember" saved analysis to avoid duplicates
# THe hash is made of column names, row count, per column NA, distinct count, and sum-of-squares
analysis_fingerprint <- function(lst, quantity_cols, spectra_cols = NULL) {
  parts <- vapply(sort(names(lst)), function(s) {
    d  <- lst[[s]]
    sc <- intersect(spectra_cols, names(d))
    d  <- d[, sort(setdiff(names(d), sc)), drop = FALSE]
    
    col_sig <- lapply(d, function(col) {
      miss  <- sum(is.na(col) | (is.character(col) & col %in% ""))
      nuniq <- length(unique(col))
      if (is.numeric(col)) {
        v <- col[!is.na(col)]
        list(miss = miss, nuniq = nuniq,
             sumsq = sum(v * v))
      } else {
        list(miss = miss, nuniq = nuniq)
      }
    })
    digest::digest(list(cols = names(d), nrow = nrow(d), col_sig = col_sig),
                   algo = "xxhash64")
  }, character(1))
  digest::digest(parts, algo = "xxhash64")
}

# ---- Harmonize vocabulary ----

# Which meta columns get controlled dropdowns. condition_* (from annotation)
# always qualifies; the rest are the fixed categorical BASE fields.
VOCAB_FIELDS <- c("instrument", "biological source", "cell line", "condition", "species", "reference database")

is_vocab_field <- function(name) {
  name %in% VOCAB_FIELDS || grepl("^condition", name, ignore.case = TRUE)
}

# normalize a value to a match key: trim, collapse whitespace, lowercase
.norm_term <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", " ", x)
  tolower(x)
}

# `%||%` guard in case it isn't already defined app-wide
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

ensure_condition_terms <- function(con) {
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS condition_terms (
      term_id     INTEGER PRIMARY KEY AUTOINCREMENT,
      field_name  TEXT NOT NULL,
      value       TEXT NOT NULL,           -- canonical display form
      value_norm  TEXT NOT NULL,           -- match key
      n_uses      INTEGER NOT NULL DEFAULT 0,
      first_seen  TEXT DEFAULT CURRENT_TIMESTAMP,
      UNIQUE(field_name, value_norm)
    )")
  invisible(TRUE)
}

# distinct known values for a field, most-used first — feeds the dropdown
condition_term_source <- function(field_name, con = NULL) {
  own <- is.null(con); if (own) { con <- .db_con(); on.exit(DBI::dbDisconnect(con)) }
  ensure_condition_terms(con)
  DBI::dbGetQuery(con,
                  "SELECT value FROM condition_terms WHERE field_name = ? ORDER BY n_uses DESC, value",
                  params = list(field_name))$value
}

# canonicalize one value; register it if new. Returns canonical spelling.
canonicalize_condition <- function(con, field_name, raw) {
  vn <- .norm_term(raw)
  if (is.na(vn) || vn == "") return(NA_character_)
  hit <- DBI::dbGetQuery(con,
                         "SELECT term_id, value FROM condition_terms WHERE field_name=? AND value_norm=?",
                         params = list(field_name, vn))
  if (nrow(hit)) {
    DBI::dbExecute(con, "UPDATE condition_terms SET n_uses = n_uses + 1 WHERE term_id = ?",
                   params = list(hit$term_id[1]))
    return(hit$value[1])                       # reuse the existing spelling
  }
  DBI::dbExecute(con,
                 "INSERT INTO condition_terms(field_name, value, value_norm, n_uses) VALUES(?,?,?,1)",
                 params = list(field_name, trimws(as.character(raw)), vn))
  trimws(as.character(raw))
}

# load a seed file once (idempotent); seeded terms start at n_uses = 0
seed_condition_terms <- function(con, seed_path = "condition_seed.csv") {
  ensure_condition_terms(con)
  if (!file.exists(seed_path)) return(invisible(FALSE))
  seed <- utils::read.csv(seed_path, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(seed))) {
    vn <- .norm_term(seed$value[i])
    ex <- DBI::dbGetQuery(con,
                          "SELECT 1 FROM condition_terms WHERE field_name=? AND value_norm=?",
                          params = list(seed$field_name[i], vn))
    if (!nrow(ex))
      DBI::dbExecute(con,
                     "INSERT INTO condition_terms(field_name, value, value_norm, n_uses) VALUES(?,?,?,0)",
                     params = list(seed$field_name[i], trimws(seed$value[i]), vn))
  }
  invisible(TRUE)
}

## ---- low-information detector ------

# TRUE for values that carry no meaning on their own (yes/no, single chars,
# bare numbers, na/nd/etc.). Empty/NA is NOT flagged — that's "missing", not "junk".
.lowinfo <- function(v) {
  s <- trimws(tolower(as.character(v)))
  present <- !is.na(s) & nzchar(s)
  junk <- grepl("^(yes|no|y|n|true|false|t|f|na|n/?a|n\\.?d\\.?|none|null|x|\\?|\\+|-|[0-9]+(\\.[0-9]+)?)$", s)
  present & (junk | nchar(s) == 1)
}

# scan a meta table for flagged condition-ish cells
scan_lowinfo_conditions <- function(meta_table) {
  cols <- names(meta_table)[vapply(names(meta_table), is_vocab_field, logical(1))]
  cols <- setdiff(cols, "instrument")
  sample_col <- if ("name" %in% names(meta_table)) "name" else "Sample"
  out <- list()
  for (cc in cols) {
    for (i in which(.lowinfo(meta_table[[cc]]))) {
      out[[length(out) + 1]] <- data.frame(
        row = i, sample = as.character(meta_table[[sample_col]][i]),
        field = cc, value = as.character(meta_table[[cc]][i]),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out))
    return(data.frame(row = integer(), sample = character(),
                      field = character(), value = character()))
  do.call(rbind, out)
}

## ---- flattened preview (your comma-list, on demand, never stored) -------

condition_preview <- function(analysis_id, con = NULL) {
  own <- is.null(con); if (own) { con <- .db_con(); on.exit(DBI::dbDisconnect(con)) }
  DBI::dbGetQuery(con, "
    SELECT Sample,
           GROUP_CONCAT(field_name || '=' || field_value, ', ') AS conditions
    FROM sample_metadata
    WHERE analysis_id = ? AND field_name LIKE 'condition%'
    GROUP BY Sample", params = list(analysis_id))
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
  
  # software from detection — always available
  tbl$software <- unname(software_map[tbl$Sample])
  
  # guarantee the baseline columns exist (empty if not already provided)
  base <- BASE_META_FIELDS
  if (any(grepl("^condition", names(tbl), ignore.case = TRUE)))
    base <- setdiff(base, "condition")            # annotation already supplies condition_* cols
  for (f in base) if (!f %in% names(tbl)) tbl[[f]] <- ""
  
  # keys first, then software, then the rest
  key <- c("Sample", "measurement", "software")
  tbl <- tbl[, c(key, setdiff(names(tbl), key)), drop = FALSE]
  tbl[] <- lapply(tbl, as.character)
  tbl
}


# ------- check data ----
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


# ---------load from database--------------
load_analysis_from_db <- function(aid, path = EPITO_DB) {
  con <- .db_con(path); on.exit(DBI::dbDisconnect(con))
  
  pep  <- DBI::dbGetQuery(con, "SELECT * FROM peptides WHERE analysis_id = ?", params = list(aid))
  if (!nrow(pep)) stop("No peptide data found for ", aid)
  qty  <- DBI::dbGetQuery(con,
                          "SELECT q.peptide_id, q.measurement, q.quantity
       FROM quantities q JOIN peptides p USING(peptide_id)
      WHERE p.analysis_id = ?", params = list(aid))
  meta <- DBI::dbGetQuery(con, "SELECT * FROM sample_metadata WHERE analysis_id = ?",
                          params = list(aid))
  
  cmap <- if ("column_map" %in% DBI::dbListTables(con))
    DBI::dbGetQuery(con, "SELECT * FROM column_map WHERE analysis_id = ?", params = list(aid)) else NULL
  mm <- if ("mod_map" %in% DBI::dbListTables(con))
    DBI::dbGetQuery(con, "SELECT token, symbol FROM mod_map WHERE analysis_id = ?", params = list(aid)) else NULL
  
  # peptides is a shared wide table; keep only the columns this analysis actually used
  analysis_finals <- if (!is.null(cmap) && nrow(cmap)) unique(cmap$final_name) else character(0)
  book      <- c("peptide_id", "analysis_id", "Sample")
  keep_cols <- names(pep)[
    names(pep) %in% book            |
      names(pep) %in% analysis_finals |
      vapply(pep, function(col) any(!is.na(col)), logical(1))
  ]
  pep <- pep[, keep_cols, drop = FALSE]
  
  drop_book <- c("peptide_id", "analysis_id", "Sample")
  samples   <- unique(pep$Sample)
  
  data_list <- lapply(samples, function(s) {
    id_s <- pep[pep$Sample == s, , drop = FALSE]
    q_s  <- qty[qty$peptide_id %in% id_s$peptide_id, , drop = FALSE]
    d <- if (nrow(q_s)) {
      qw <- tidyr::pivot_wider(q_s, id_cols = peptide_id,
                               names_from = "measurement", values_from = "quantity")
      dplyr::left_join(id_s, qw, by = "peptide_id")
    } else id_s
    as.data.frame(d[, setdiff(names(d), drop_book), drop = FALSE], stringsAsFactors = FALSE)
  })
  names(data_list) <- samples
  
  meas    <- unique(qty$measurement)
  col_map <- setNames(lapply(meas, function(m) list(name = m, col = m)), meas)
  sw <- meta[meta$field_name == "software", c("Sample", "field_value")]
  software_map <- setNames(sw$field_value, sw$Sample)
  
  list(data = data_list, col_map = col_map, software = software_map,
       meta = meta, annotation = reconstruct_annotation(meta),
       data_info = reconstruct_data_info(cmap),
       data_mod_map = reconstruct_mod_map(mm))
}

reconstruct_annotation <- function(meta) {
  if (!nrow(meta)) return(NULL)
  wide <- tidyr::pivot_wider(
    meta[, c("Sample", "measurement", "field_name", "field_value")],
    names_from = "field_name", values_from = "field_value")
  wide <- as.data.frame(wide, stringsAsFactors = FALSE)
  names(wide)[names(wide) == "Sample"] <- "name"

  if (!"source" %in% names(wide))
    wide$source <- if ("measurement" %in% names(wide)) wide$measurement else wide$name
  wide
}

# rebuild merged_info (final_name + one list-column per sample) from the long table
reconstruct_data_info <- function(cmap) {
  if (is.null(cmap) || !nrow(cmap)) return(NULL)
  samples <- unique(cmap$Sample)
  finals  <- unique(cmap$final_name)
  di <- data.frame(final_name = finals, stringsAsFactors = FALSE)
  for (s in samples) {
    di[[s]] <- lapply(finals, function(fn)
      cmap$original_name[cmap$Sample == s & cmap$final_name == fn])
  }
  di
}

#----------source order problem-------------

local({
  con <- .db_con(); on.exit(DBI::dbDisconnect(con))
  ensure_condition_terms(con)
  seed_condition_terms(con, "condition_seed.csv")
})