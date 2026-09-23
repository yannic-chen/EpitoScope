#' @keywords internal
"_PACKAGE"

#' @rawNamespace import(shiny, except = c(column, actionButton, tabsetPanel, insertTab, navbarMenu))
#' @import bs4Dash
#' @import dplyr
#' @import tidyr
#' @rawNamespace import(ggplot2, except = c(last_plot))
#' @import plotly
#' @importFrom rlang .data .env sym
#' @importFrom data.table := as.data.table
#' @importFrom purrr imap imap_dfr map
#' @importFrom stringr str_replace_all str_trim
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggseqlogo ggseqlogo
#' @importFrom patchwork wrap_plots
#' @importFrom bslib layout_column_wrap
#' @importFrom ggVennDiagram ggVennDiagram
#' @importFrom viridisLite viridis
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom grid gpar textGrob convertWidth stringWidth
#' @importFrom grDevices adjustcolor col2rgb colorRampPalette
#' @importFrom graphics hist par
#' @importFrom stats approx cor density dist hclust model.matrix p.adjust prcomp setNames t.test
#' @importFrom utils capture.output combn head read.csv read.delim read.table tail write.csv
NULL

utils::globalVariables(c(
  ".", "A", "Allele", "Binder", "Class", "Count", "EG.ModifiedPeptide",
  "EG.Qvalue", "FG.Quantity", "Group", "Intensity", "Item", "Length_bin",
  "Mapped.Proteins", "Max", "Mean", "Mean_G1", "Mean_G2", "Measurement",
  "Mid", "Min", "Missing", "Modified.Peptide", "Modified.Sequence", "Non",
  "Non_NA_Count", "Peptide", "Percent_Non_NA", "Percent_Rank",
  "Precursor.Quantity", "Protein", "Protein_Mapped.Proteins", "Q.Value",
  "Quantity", "R.FileName", "Rank", "RawFile", "Run", "Sample", "Set",
  "Significance", "Spectrum", "Strong", "Total", "Weak", "X",
  "adj_pval_BH", "adj_pval_Bonf", "allele", "coalesced", "final_name",
  "hl", "log2FC", "name", "negLog10AdjP_BH", "netMHCpan",
  "netmhcpan_use_wsl", "original_name", "peptide", "peptide_id", "perc",
  "pval", "score", "selected", "value", "y",
  "CHARGE", "K0", "LENGTH", "MAX_QUANTITY", "MZ", "PEPTIDE", "PTM", "STRIPPED"
))
