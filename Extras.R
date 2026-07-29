# Split each sample df into per-measurement frames; the quantity column name IS the measurement name
split_into_measurements <- function(lst, quantity_cols) {
  pep <- c("STRIPPED","PEPTIDE","PROTEIN","LENGTH","CHARGE")
  out <- list()
  for (s in names(lst)) {
    d <- lst[[s]]
    for (col in intersect(quantity_cols, names(d))) {
      md <- d[, intersect(pep, names(d)), drop = FALSE]
      md$quantity   <- d[[col]]
      md$Sample     <- s
      md$measurement <- col                       # column name → measurement name
      out[[col]] <- if (is.null(out[[col]])) md else rbind(out[[col]], md)
    }
  }
  out   # named list keyed by measurement (= the quantity column name)
}