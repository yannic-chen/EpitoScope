# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Script Name:        global.R
# Purpose:            This script creates a shiny app for easy immunopeptidomics analysis
# Author:             Yannic Chen
# Date Created:       2025-11-19
# Last Modified:      2026-02-02
# Version:            0.1
# R Version:          4.5.2
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Details:
# This script contains most of the functions to generate plots for the Shiny.R app.
# These are also called when generating the report in the report.Rmd script.
#
# To Do:
# For package management, switch to pkg::fun() style (i.e. dplyr::mutate()).
# create a required_pkgs list and install these when not already installed:
#        missing <- setdiff(required_pkgs, rownames(installed.packages()))
#       if (length(missing)) {
#          install.packages(missing)
#       }
# Perhaps check our renv::init()?

#These are for Shiny UI
library(shiny)
library(shinyBS)
library(bs4Dash)
library(bslib)
library(DT)
library(shinyjs) #Only used to grey out buttons.
#These are for general data handling
library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(stringr)
library(reshape2)
library(curl)
library(anticlust) #this is used for same size clustering to speed up heatmap generation for large datasets by clustering datapoints together.
library(limma) #only for grouped comparison when only 1 measurement exists in the group. Limma is used to infer p-value.
#These are for plotting
library(viridis)
library(ggplot2)
library(ggrepel)
library(ggseqlogo)
library(ggVennDiagram)
library(plotly)
library(ComplexUpset)
library(ComplexHeatmap)
library(TSP) #required for heatmaply
library(registry) #required for heatmaply
library(ca) #required for heatmaply
library(heatmaply)
library(circlize)
library(grid)
library(patchwork) #This is only used for the 1/k0 vs m/z plot. Could remove this by changing the code.
#These are exclusive for GO-term.
library(clusterProfiler) #this one masks a lot of dplyr and other package functions
library(org.Hs.eg.db)
#These are exclusive for STRING-DB
library(httr)
library(jsonlite)
library(igraph)
library(ggraph)
library(tidyverse)

if (requireNamespace("fastcluster", quietly = TRUE)) {
  hclust_fn <- fastcluster::hclust
} else {
  hclust_fn <- stats::hclust
}

#Suppress warnings from everywhere, including plot_ly
globalCallingHandlers(warning = function(w) {
  if (grepl("group_hm", conditionMessage(w), fixed = TRUE))
    invokeRestart("muffleWarning")
})

#-----------Column extraction------------
# Here we initiate all the possible column names important for us from all different input formats
# QUANTITY and SPECTRA values are treated as column name string to be searchs, since one column exists for each measurement in the sample.
# It is likely that they are a combination of a word and the measurement name, mostly a suffix or prefix + measurement name.

column_schema <- list(
  #If multiple column for a stat is detected, an error will occur.
  PEAKS = list( #Only checked for peptide.tsv
    PEPTIDE        = c("Peptide"),                       
    STRIPPED       = c(),                       # not present in peptide.tsv
    LENGTH         = c("Length"),
    MASS           = c("Mass"),
    MZ             = c("m.z"),                  #Could also switch to Raw.M.z.
    SCORE          = c("X.10LgP"),
    CHARGE         = c("z"),
    RT             = c("RT"),                   #Could also switch to Raw.RT.
    K0             = c("X1.k0.Range"),          #PEAKS 12 and 13 uses "X1.k0.Start" and "X1.k0.End", for which the K0 needs to be calculated from the middle value. PEAKS Online returns "X1.k0.Range"
    PPM            = c("ppm"),
    PROTEIN        = c("Accession"),
    QUANTITY       = c("area."),                 #Prefix. Some versions also has intensity, but we prefer area.
    SPECTRA        = c("x.feature.", "X.Spec."),  #Prefix. X.Spec for PEAKS 11 Online. X.Feature for PEAKS 12 studio. PEAKS 13 returns both. Prefer x.spec over x.feature. PEAKS is unique in that peptides can be found but still have 0 or NA quantity, hence using spectra here is crucial
    PTM            = c("PTM")
  ),
  
  Fragpipe = list( #Checked psm.tsv and combined.peptide.tsv
    PEPTIDE        = c("Modified.Peptide", "Modified.Sequence"),
    STRIPPED       = c("Peptide", "Peptide.Sequence"),           #not present in combined_peptide.tsv
    LENGTH         = c("Peptide.Length"),                        #not present in combined_modified_peptide.tsv, combined_peptide.tsv
    MASS           = c("Observed.Mass"),                         #can also switch to Calculated.Peptide.Mass, #not present in combined_modified_peptide.tsv, combined_peptide.tsv
    MZ             = c("Observed.M.Z"),       #can also switch to Calibrated.Observed.M.Z or Calculated.M.Z , #not present in combined_modified_peptide.tsv, combined_peptide.tsv
    SCORE          = c("Probability", "PeptideProphet.Probability"), #not present in combined_modified_peptide.tsv, combined_peptide.tsv
    CHARGE         = c("Charge", "Charges"),
    RT             = c("Retention"),                             #not present in combined_modified_peptide.tsv, combined_peptide.tsv
    K0             = c("ion.mobility"),                          #not present in combined_modified_peptide.tsv, combined_peptide.tsv
    PPM            = c("Delta.Mass"),                            #not present in combined_modified_peptide.tsv, combined_peptide.tsv
    PROTEIN        = c("Protein_Mapped.Proteins"),               #This column is created later from Protein and Mapped.Protein column. #combined_modified_peptide.tsv, combined_peptide.tsv only has Protein column.
    QUANTITY       = c("maxlfq.intensity", "Intensity"),         #suffix. prefer maxlfq.intensity (exist in peptide.tsv, but not psm.tsv).
    SPECTRA        = c("Spectral.Count"),                        #suffix. Cannot use spectrum in psm.csv file since these are only the names and are lost
    PTM            = c("Assigned.Modifications")                 #not present in combined_modified_peptide.tsv, combined_peptide.tsv
  ),
  
  DIANN = list( #This is for report.pr_matrix.tsv from the DIANN of the fragpipe pipeline (may be the same as original DIANN)
    PEPTIDE        = c("modified.sequence"),
    STRIPPED       = c("stripped.sequence"),
    LENGTH         = c(),                  # not present in report.pr_matrix.tsv
    MASS           = c(),                  # not present in report.pr_matrix.tsv
    MZ             = c(),                  # not present in report.pr_matrix.tsv
    SCORE          = c(),                  # not present in report.pr_matrix.tsv
    CHARGE         = c(),                  # not present in report.pr_matrix.tsv
    RT             = c(),                  # not present in report.pr_matrix.tsv
    K0             = c(),                  # not present in report.pr_matrix.tsv
    PPM            = c(),                  # not present in report.pr_matrix.tsv
    PROTEIN        = c("Protein.names"),   #Can also switch with Protein.IDs, Protein.Group or Genes (although some proteins lack gene, like the CONTAs)
    QUANTITY       = c("D..data"),         #This is actually the path to the file. So its a prefix.
    SPECTRA        = c("D..data"),         #Can use the same as quantity.
    PTM            = c("")                 # not present in report.pr_matrix.tsv
  ),
  
  DIANN_parquet = list( #This is for DIANN parquet file, which is in long format.
    PEPTIDE        = c("modified.sequence"),
    STRIPPED       = c("stripped.sequence"),
    LENGTH         = c(),                  # not present in parquet
    MASS           = c(),                  # not present in parquet
    MZ             = c("Precursor.Mz"),
    SCORE          = c("Q.Value"),         # this is used for FDR. But Global.Q.Value can also be used. PEP does not replace Q.value, but can be used as additional filter.
    CHARGE         = c("Precursor.Charge"),   
    RT             = c("RT"),              # could also use iRT, predicted.RT or predicted.iRT
    K0             = c("IM"),              # There is also indexed and predictet and predicted.iIM
    PPM            = c(),                  # not present in parquet
    PROTEIN        = c("Protein.names"),     #Can also switch with Protein.IDs, Protein.Group or Genes (although some proteins lack gene, like the CONTAs)
    QUANTITY       = c("Precursor.quantity"), #Could Use Precursor.Normalized (but values look the same to me). Transformed from long format. Why not use MS1.Area/Ms1.normalised
    SPECTRA        = c(),                   #Transformed from long format.
    PTM            = c()                   # not present in parquet
  )
)

#these columns are used to identify software
signature <- list(
  PEAKS    = c("X.10LgP"),
  Fragpipe = c("prev.aa"),
  DIANN    = c("First.Protein.Description"), #this is for report.pr_matrix.tsv.
  DIANN_parquet    = c("Run.Index")
)

#read netMHCpan allelenames.
hla_alleles <- tryCatch({
  lines <- readLines(file.path(getwd(), "allelenames.netMHCpan"))
  lines <- trimws(lines)
  lines <- lines[nchar(lines) > 0 & !startsWith(lines, "#")]
  # First column is the allele name passed to netMHCpan
  alleles <- sapply(strsplit(lines, "\\s+"), `[`, 1)
  # Group by species prefix (HLA-A, HLA-B, BoLA, etc.) for organised dropdown
  prefixes <- sub("([^-]+-[^:0-9]*).*", "\\1", alleles)
  split(alleles, prefixes)
}, error = function(e) {
  warning("allelenames.netMHCpan not found, using default allele list.")
  list(
    "HLA-A" = c(
      "HLA-A01:01",
      "HLA-A02:01", "HLA-A02:03", "HLA-A02:06",
      "HLA-A03:01",
      "HLA-A11:01",
      "HLA-A23:01",
      "HLA-A24:02",
      "HLA-A26:01",
      "HLA-A29:02",
      "HLA-A30:01", "HLA-A30:02",
      "HLA-A31:01",
      "HLA-A32:01",
      "HLA-A33:01",
      "HLA-A68:01", "HLA-A68:02"
    ),
    "HLA-B" = c(
      "HLA-B07:02",
      "HLA-B08:01",
      "HLA-B13:01",
      "HLA-B15:01",
      "HLA-B18:01",
      "HLA-B27:05",
      "HLA-B35:01",
      "HLA-B39:01",
      "HLA-B40:01",
      "HLA-B44:02", "HLA-B44:03",
      "HLA-B51:01",
      "HLA-B57:01",
      "HLA-B58:01"
    ),
    "HLA-C" = c(
      "HLA-C03:03", "HLA-C03:04",
      "HLA-C04:01",
      "HLA-C05:01",
      "HLA-C06:02",
      "HLA-C07:01", "HLA-C07:02",
      "HLA-C08:02",
      "HLA-C12:03"
    )
  )
  
})

#-----------Helper functions------------------
safe_reactive <- function(x) {
  tryCatch(
    x(),
    error = function(e) NULL
  )
}

safe_validate <- function(cond, message) {
  if (!cond) {
    if (isTRUE(getOption("knitr.in.progress"))) {
      # In R Markdown → DO NOT STOP
      return(list(.skip = TRUE, message = message))
    } else {
      # In Shiny → show validation message
      shiny::validate(shiny::need(FALSE, message))
    }
  }
  return(NULL)
}

#This is for avoiding grid viewport stack corruption
safe_draw <- function(ht, ...) {
  tryCatch(
    ComplexHeatmap::draw(ht, ...),
    error = function(e) {
      try(grid::upViewport(0), silent = TRUE)
      shiny::validate(shiny::need(FALSE, conditionMessage(e)))
    }
  )
}

remove_ptms <- function(x) {
  gsub("\\(.*?\\)|\\[.*?\\]|\\{.*?\\}", "", x)
}

convert_with_mod_map <- function(sequences, mod_map) {
  
  convert_single <- function(seq) {
    out <- character(0)
    i <- 1
    
    while (i <= nchar(seq)) {
      rest <- substr(seq, i, nchar(seq))
      
      matched <- FALSE
      for (token in names(mod_map)) {
        if (startsWith(rest, token)) {
          out <- c(out, mod_map[[token]])
          i <- i + nchar(token)
          matched <- TRUE
          break
        }
      }
      
      if (!matched) {
        out <- c(out, substr(seq, i, i))
        i <- i + 1
      }
    }
    
    paste0(out, collapse = "")
  }
  
  unname(vapply(sequences, convert_single, character(1)))
}

extract_mod_tokens <- function(sequences) {
  tokens <- character(0)
  
  for (seq in sequences) {
    rest <- seq
    while (nchar(rest) > 0) {
      
      if (grepl("^n\\[[0-9.]+\\][A-Z]", rest)) {
        token <- regmatches(rest, regexpr("^n\\[[0-9.]+\\][A-Z]", rest))
        tokens <- c(tokens, token)
        rest <- substr(rest, nchar(token) + 1, nchar(rest))
        next
      }
      
      if (grepl("^[A-Z]\\[[0-9.]+\\]", rest)) {
        token <- regmatches(rest, regexpr("^[A-Z]\\[[0-9.]+\\]", rest))
        tokens <- c(tokens, token)
        rest <- substr(rest, nchar(token) + 1, nchar(rest))
        next
      }
      
      if (grepl("^[A-Z]\\(\\+[0-9.]+\\)", rest)) {
        token <- regmatches(rest, regexpr("^[A-Z]\\(\\+[0-9.]+\\)", rest))
        tokens <- c(tokens, token)
        rest <- substr(rest, nchar(token) + 1, nchar(rest))
        next
      }
      
      rest <- substr(rest, 2, nchar(rest))
    }
  }
  
  unique(tokens)
}

make_symbol <- function(i) {
  pool <- c(as.character(1:9), letters)
  if (i > length(pool))
    stop("Too many PTMs to encode as single characters")
  pool[i]
}

build_mod_map <- function(tokens) {
  setNames(
    vapply(seq_along(tokens), make_symbol, character(1)),
    tokens
  )
}

extract_protein_prefixes <- function(accessions) {
  all_acc <- unlist(strsplit(accessions, ";"))
  proteins <- sub(" .*", "", trimws(all_acc))  # truncate at first space
  unique(proteins[nchar(proteins) > 0])
}

aa_comp_from_peptides <- function(peptides) {
  aa <- unlist(strsplit(peptides, ""))
  tab <- table(aa)
  
  freq <- as.numeric(tab)
  names(freq) <- names(tab)
  
  100 * freq / sum(freq)
}

check_data_error <- function(data, required_cols = NULL, na_policy = c("any", "all", "ignore")) {
  #'na_policy determines how strict we allow NAs:
  #'         ignore: the columns can have any number of NA
  #'         allL: only if all values in column are NA, break operation
  #'         any: if even just one NA exist, break it.
  #'
  #'
  #'
  na_policy <- match.arg(na_policy)
  
  stop_with_msg <- function(msg) {
    shiny::validate(shiny::need(FALSE, msg))
  }
  
  if (is.null(data)) stop_with_msg("Data not available")
  
  # Single data.frame case
  if (is.data.frame(data)) {
    missing_cols <- dplyr::setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) stop_with_msg(
      paste("Missing column(s):", paste(missing_cols, collapse = ", "))
    )
    
    if (na_policy != "ignore") {
      for (col in required_cols) {
        if (na_policy == "any" && any(is.na(data[[col]]))) stop_with_msg(paste("NA in column:", col))
        if (na_policy == "all" && all(is.na(data[[col]]))) stop_with_msg(paste("All NA in column:", col))
      }
    }
    
    return(TRUE)
  }
  
  #named list case
  if (is.list(data)) {
    # Keep only data frames
    valid_dfs <- data[vapply(data, is.data.frame, logical(1))]
    if (length(valid_dfs) == 0) stop_with_msg("No valid data.frames in the list")
    
    # Check which dfs have all required columns
    dfs_with_cols <- valid_dfs[vapply(valid_dfs, function(df) all(required_cols %in% colnames(df)), logical(1))]
    
    if (length(dfs_with_cols) == 0) stop_with_msg(
      paste("None of the data.frames contain required column(s):", paste(required_cols, collapse = ", "))
    )
    
    # NA check: only on dfs that have the columns
    if (na_policy != "ignore") {
      for (df in dfs_with_cols) {
        for (col in required_cols) {
          if (na_policy == "any" && any(is.na(df[[col]]))) {
            stop_with_msg(paste("NA values found in column:", col))
          }
          if (na_policy == "all" && all(is.na(df[[col]]))) {
            stop_with_msg(paste("All values are NA in column:", col))
          }
        }
      }
    }
    
    return(TRUE)
  }
  
  stop_with_msg("Unsupported data type")
}

compute_padding <- function(labels, fontsize = 10, rot = 90) {
  
  # Create text grob
  tg <- textGrob(labels,
                 gp = gpar(fontsize = fontsize),
                 rot = rot)
  
  # Measure width
  max_width <- max(convertWidth(stringWidth(labels), "mm", valueOnly = TRUE))
  
  # Add a little buffer (important)
  padding_mm <- max_width + 1
  
  padding_mm
}

#for report.Rmd
plot_per_sample_grid <- function(lst, plot_fn, ncol = 2, ...) {
  plots <- lapply(names(lst), function(nm) {
    tryCatch(
      plot_fn(df = lst[[nm]], title_name = nm, ...),
      error = function(e) {
        ggplot() +
          annotate("text", x = .5, y = .5, label = paste("Error:", e$message),
                   size = 4, color = "red") +
          theme_void() +
          ggtitle(nm)
      }
    )
  })
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) return(invisible(NULL))
  
  pw <- patchwork::wrap_plots(plots, ncol = ncol)
  tryCatch(
    print(pw),
    error = function(e) {
      print(
        ggplot() +
          annotate("text", x = .5, y = .5,
                   label = "Plot window too small — resize the Plots pane and regenerate.",
                   size = 4, color = "grey40") +
          theme_void()
      )
    }
  )
  invisible(NULL)
}

#for report.Rmd
plot_motif_grid <- function(lst, lengths, ncol = 3, 
                            ptm_only = FALSE, 
                            namespace = NULL) {
  length_panels <- lapply(lengths, function(L) {
    
    sample_plots <- lapply(names(lst), function(nm) {
      df <- lst[[nm]]
      
      # PTM filter — only modified peptides
      if (ptm_only) {
        df <- df[!is.na(df$PTM) & df$PTM != "", ]
        peptide_col <- "PTM_Pseudo"
      } else {
        peptide_col <- "STRIPPED"
      }
      
      df_L     <- df[df$LENGTH == L, ]
      peptides <- unique(df_L[[peptide_col]])
      peptides <- peptides[!is.na(peptides) & nchar(peptides) == L]
      
      if (length(peptides) < 5) {
        return(
          ggplot() +
            annotate("text", x = .5, y = .5,
                     label = "Not enough peptides",
                     size = 3, color = "grey50") +
            theme_void() +
            ggtitle(paste(nm, "• Length", L))
        )
      }
      
      par(mar = c(1.5, 1.5, 2, 0.5))
      plot_seqlogo(
        peptides,
        title     = paste(nm, "• Length", L),
        namespace = namespace
      )
    })
    
    patchwork::wrap_plots(
      c(sample_plots),
      ncol = length(lst)
    )
  })
  
  patchwork::wrap_plots(length_panels, ncol = 1)
}

get_quantity_cols <- function(data_info) {
  data_info %>%
    dplyr::filter(final_name == "QUANTITY") %>%
    dplyr::select(-final_name) %>%
    unlist(recursive = TRUE, use.names = FALSE)
}

#-----------Data handling/transformation functions------------------
build_generic_schema <- function(schema) {
  
  # Get all field names (peptide, peptidoform, etc.)
  fields <- unique(unlist(lapply(schema, names)))
  
  generic <- list()
  
  for (field in fields) {
    # Collect all candidates across software types
    combined <- unique(unlist(lapply(schema, function(x) x[[field]])))
    generic[[field]] <- combined
  }
  
  return(generic)
}

register_custom_schema <- function(custom_schema    = NULL, custom_signature = NULL,replace_schema   = FALSE) {
  # Capture all_keys before any nulling — needed for filling missing keys
  all_keys <- names(column_schema[[1]])
  
  if (!is.null(custom_schema)) {
    
    if (replace_schema) {
      column_schema <<- list()  # empty list rather than NULL — safer for [[<- assignment
      message("Replaced default schema.")
    }
    
    for (nm in names(custom_schema)) {
      
      has_peptide  <- !is.null(custom_schema[[nm]]$PEPTIDE)  && length(custom_schema[[nm]]$PEPTIDE)  > 0
      has_stripped <- !is.null(custom_schema[[nm]]$STRIPPED) && length(custom_schema[[nm]]$STRIPPED) > 0
      
      if (!has_peptide && !has_stripped) {
        stop(paste0(
          "Custom schema '", nm, "' must specify at least one of PEPTIDE or STRIPPED."
        ))
      }
      
      # Fill missing keys silently with empty vectors
      complete_entry <- setNames(
        lapply(all_keys, function(k) c()),  # all keys start empty
        all_keys
      )
      
      # Overwrite with user-provided values
      for (k in names(custom_schema[[nm]])) {
        complete_entry[[k]] <- custom_schema[[nm]][[k]]
      }
      
      column_schema[[nm]] <<- complete_entry
      message(paste0("Added custom schema: '", nm, "'"))
    }
  }
  
  if (!is.null(custom_signature)) {
    
    if (replace_schema) {
      signature <<- list()  # same — empty list rather than NULL
      message("Replaced default signature.")
    }
    
    for (nm in names(custom_signature)) {
      signature[[nm]] <<- custom_signature[[nm]]
      message(paste0("Registered custom signature: '", nm, "'"))
    }
  }
}

check_annotation_table <- function(df) {
  #check the file headers
  message("Annotation file given.")
  df <- distinct(df)
  colnames(df) <- tolower(colnames(df))
  header <- colnames(df)
  
  # ── Required columns ──────────────────────────────────────────────────────
  if(!("name" %in% header && "source" %in% header)) {
    stop("Missing either Name or Source in Annotation Table")
  }
  
  # The same source file cannot apply to multiple names, but multiple source file can apply to the same name (e.g. fragpipe psm.tsv)
  distinctmap <- all(tapply(df$name, df$source, function(x) length(unique(x)) == 1))

  if(!all(distinctmap)) {
    stop("Some 'source' files map to multiple 'name' values — each source must only map to one name.")
  }
    
  if(ncol(df) == 2) {
    # if ncol is 2, we know it doesnt contain any optional info
    message("No condition, replicate or measurement info given. Proceed as default")
    return(df)
  }

  #now we need to check. It is possible to allow conditions based only on the sample name. However, in that case, one can just use the default grouping mechanis.
  #Only when one has multiple conditions associated to the same sample name does it become important to have the measurement.
  
  has_measurement <- "measurement" %in% header
  condition_cols  <- grep("^condition", header, value = TRUE)
  replicate_cols  <- dplyr::intersect(c("biological_replicate", "technical_replicate"), header)
  
  
  # ── Type coercion ─────────────────────────────────────────────────────────
  if (length(condition_cols) > 0) {
    df[condition_cols] <- lapply(df[condition_cols], as.character)
  }
  if (length(replicate_cols) > 0) {
    df[replicate_cols] <- lapply(df[replicate_cols], as.character)
  }
  
  # ── Measurement column ────────────────────────────────────────────────────
  
  if(has_measurement) {
    #measurement must be distinct. 
    if(any(duplicated(df$measurement))){
      #if we find duplicates, then we can try to make them unique, by combining them with names.
      #This is under the assumption that for analysis of different samples, one can have the same measurement name. In that case, the sample name will make them distinct.
      warning("There are duplicates in the measurement names. Attempting to make unique by prepending sample name.")
      df$measurement <- paste(df$name, df$measurement, sep = "_")
      if(any(duplicated(df$measurement))){
        stop("Duplicate measurement names remain after prepending sample name.")
        }
    }
    
  } else {
    message("no measurement column found.")

    df %>%
      dplyr::select(all_of(c("name", condition_cols, replicate_cols))) %>%
      dplyr::group_by(name) %>%
      dplyr::summarise(n = n(), .groups = "drop") %>%
      { 
        if (any(.$n > 1)) {
          stop("Multiple condition + replicate rows identified for the same sample name. 
               Without 'measurement' column, each name can only be assignet to one condition + replicate combination.
               Add a 'measurement' column to distinguish individual sources.")
        }
      }
  }
  
  # ── Attach metadata as attributes for downstream use ──────────────────────
  attr(df, "replicate_cols") <- replicate_cols
  attr(df, "condition_cols") <- condition_cols
  attr(df, "has_measurement") <- has_measurement
  
  df
  
}

load_from_annotation <- function(annotation_df) {
  
  annotation_df <- check_annotation_table(annotation_df)
  
  split_df <- split(annotation_df, annotation_df$name)
  
  data_list <- lapply(names(split_df), function(nm) {
    
    sub_df <- split_df[[nm]]
    
    sub_df <- unique(sub_df[, c("name", "source")])
    
    dfs <- lapply(sub_df$source, function(path) {
      
      # Normalise path
      path <- trimws(path)
      path <- gsub("\\\\", "/", path)
      path <- gsub("/+", "/", path)
      path <- normalizePath(path, winslash = "/", mustWork = FALSE)
      
      if (!file.exists(path)) {
        stop(paste0("File not found for '", nm, "': ", path))
      }
      
      ext <- tolower(tools::file_ext(path))
      
      message(paste0("Loading '", nm, "' from: ", path, " (", ext, ")"))
      
      switch(ext,
             csv     = read.csv(path,  stringsAsFactors = FALSE),
             tsv     = read.delim(path, stringsAsFactors = FALSE),
             txt     = read.delim(path, stringsAsFactors = FALSE),
             parquet = arrow::read_parquet(path),
             stop(paste0("Unsupported file type '.", ext, "' for '", nm, "'."))
      )
    })
    
    # Combine all sources for this name
    dfs <- dplyr::bind_rows(dfs)
  })
  
  names(data_list) <- names(split_df)
  
  attr(data_list, "annotation") <- annotation_df
  
  data_list
}

build_measurement_col_map <- function(data_list, annotation_df) {
  if (!isTRUE(attr(annotation_df, "has_measurement"))) return(NULL)
  
  col_map <- list()
  
  for (nm in unique(annotation_df$name)) {
    sub_ann <- annotation_df[annotation_df$name == nm, ]
    df_cols <- colnames(data_list[[nm]])
    
    for (meas in sub_ann$measurement) {
      hits <- grep(tolower(meas), tolower(df_cols), value = TRUE, fixed = TRUE)
      
      if (length(hits) == 0) {
        stop(paste0("Measurement '", meas, "' (sample '", nm, "') ",
                    "not found in any column of the data. ",
                    "Ensure the measurement name appears in the quantity column names."))
      }
      if (length(hits) > 1) {
        stop(paste0("Measurement '", meas, "' (sample '", nm, "') ",
                    "matches multiple columns: ",
                    paste(hits, collapse = ", "), ". ",
                    "Measurement names must be unique enough to identify exactly one column."))
      }
      
      col_map[[meas]] <- list(name = nm, col = hits)
    }
  }
  
  col_map
}


# Evaluate a condition expression against an annotation data.frame.
# expr_terms: list of list(col, val, op)
#   col — condition column name
#   val — character vector of selected values
#   op  — NULL for first term; "AND", "OR", or "NOT" for subsequent terms
# Returns a named logical vector: sample_name → TRUE/FALSE
eval_condition_expr <- function(expr_terms, ann_df) {
  if (length(expr_terms) == 0) return(setNames(logical(0), character(0)))
  
  # Keep all rows, no collapsing to unique sample
  sample_names <- unique(ann_df$name)
  
  # Base mask for the first term
  t1 <- expr_terms[[1]]
  
  # Evaluate value mask
  val_mask <- ann_df[[t1$col]] %in% t1$val
  
  # Evaluate measurement mask if provided (case-insensitive)
  if (isTRUE(attr(ann_df, "has_measurement")) && !is.null(t1$measurement)) {
    meas_mask <- tolower(ann_df$measurement) %in% tolower(t1$measurement)
    val_mask <- val_mask & meas_mask
  }
  
  # Create a result column per row instead of collapsing yet
  ann_df$result <- val_mask
  
  # Process remaining terms
  if (length(expr_terms) > 1) {
    for (i in 2:length(expr_terms)) {
      t <- expr_terms[[i]]
      tmask <- ann_df[[t$col]] %in% t$val
      
      if (isTRUE(attr(ann_df, "has_measurement")) && !is.null(t$measurement)) {
        meas_mask <- tolower(ann_df$measurement) %in% tolower(t$measurement)
        tmask <- tmask & meas_mask
      }
      
      # Combine row-wise
      ann_df$result <- switch(t$op,
                              "AND" = ann_df$result & tmask,
                              "OR"  = ann_df$result | tmask,
                              "NOT" = ann_df$result & !tmask,
                              stop(paste("Unknown operator:", t$op))
      )
    }
  }
  
  out_cols <- c("name",
                if (isTRUE(attr(ann_df, "has_measurement"))) "measurement",
                "result")
  ann_df[, out_cols, drop = FALSE]
}

detect_software <- function(df, signature, fallback = "Generic") {
  
  df_cols <- tolower(colnames(df))
  
  matches <- sapply(names(signature), function(soft) {
    any(tolower(signature[[soft]]) %in% df_cols)
  })
  
  n_matches <- sum(matches)
  
  if (n_matches == 1) {
    detected <- names(matches)[matches]
    message("Detected software: ", detected)
    return(detected)
  }
  
  if (n_matches > 1) {
    stop("Multiple software signatures detected: ",
         paste(names(matches)[matches], collapse = ", "))
  }
  
  message("No signature detected. Using fallback: ", fallback)
  return(fallback)
}

find_transform_column <- function(df, candidates, name, exclude = character(0)) {
  for (cand in candidates) {
    # Exact case-sensitive match, excluding already-renamed columns
    exact_match <- which(colnames(df) == cand & !colnames(df) %in% exclude)
    if (length(exact_match) == 1) {
      original_name <- colnames(df)[exact_match]
      colnames(df)[exact_match] <- name
      return(list(df = df, matched_column = original_name))
    }
    
    # Case-insensitive fallback, excluding already-renamed columns
    match <- which(tolower(colnames(df)) == tolower(cand) & !colnames(df) %in% exclude)
    if (length(match) == 1) {
      original_name <- colnames(df)[match]
      colnames(df)[match] <- name
      return(list(df = df, matched_column = original_name))
    } else if (length(match) > 1) {
      stop("Ambiguous: multiple columns match '", cand, "' for ", name, ": ",
           paste(colnames(df)[match], collapse = ", "))
    }
  }
  
  message(paste0("No matching column found for: ", name))
  return(NULL)
}


#!!!important, the order for stripped must be before peptide, since peptide is a rather ambiguous name, which can result in false assignment.
transform_columns <- function(df, schema, software, targets = c("STRIPPED", "PEPTIDE", "LENGTH", "MASS", "MZ", "SCORE", "CHARGE", "RT", "PPM", "PROTEIN", "PTM")) {
  mapping_log <- data.frame(
    final_name = character(),
    original_name = character(),
    stringsAsFactors = FALSE
  )
  
  # Get the software-specific column map
  col_map <- schema[[software]]
  
  # Limit to requested targets if provided
  if (!is.null(targets)) {
    col_map <- col_map[names(col_map) %in% targets]
  }
  
  already_renamed <- character(0)
  
  for (target_name in names(col_map)) {
    candidates <- col_map[[target_name]]
    
    # Skip if no candidates defined
    if (length(candidates) == 0) next
    
    res <- find_transform_column(df, candidates, target_name, exclude = already_renamed)
    if (is.null(res)) next
    df <- res$df
    
    already_renamed <- c(already_renamed, target_name)
    
    # Log the mapping
    mapping_log <- rbind(mapping_log,
                         data.frame(final_name = target_name,
                                    original_name = res$matched_column,
                                    stringsAsFactors = FALSE))
  }
  
  return(list(df = df, log = mapping_log))
}

log_rename <- function(log_df, res, new_name) {
  rbind(log_df,
        data.frame(final_name = new_name,
                   original_name = res$matched_column,
                   stringsAsFactors = FALSE))
}

get_midpoint <- function(x) { #this function retrieves the middle value given a character like this: "15-35"
  sapply(x, function(val) {
    # If already numeric, keep as-is
    if (is.numeric(val)) return(val)
    
    # Convert to character
    val <- as.character(val)
    
    # Check if it contains a "-"
    if (grepl("-", val)) {
      parts <- as.numeric(strsplit(val, "-")[[1]])
      return(mean(parts, na.rm = TRUE))
    } else {
      return(as.numeric(val))
    }
  })
}

aggregate_fragpipe_psm <- function(df) {
  #' We groupby the fragpipe psm.tsv according to file and modified.peptide
  #' We collapse by taking the row with maximum PeptideProphet.Probability or Probability.
  
  df <- df %>%
    dplyr::mutate(
      RawFile = sub("\\.\\d+\\..*$", "", Spectrum)
    )
  
  num_cols <- df %>%
    dplyr::select(where(is.numeric)) %>%
    names() %>%
    dplyr::setdiff(c("RawFile", "Modified.Peptide"))
  
  char_cols <- df %>%
    dplyr::select(where(is.character)) %>%
    names() %>%
    dplyr::setdiff(c("RawFile", "Modified.Peptide", "Spectrum", "Spectrum.File"))
  
  df %>%
    dplyr::group_by(RawFile, Modified.Peptide) %>%
    dplyr::slice_max(
      .data[[dplyr::intersect(c("PeptideProphet.Probability", "Probability"), names(df))[1]]], #Try "PeptideProphet.Probability", then "Probability" column
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup()
}

column_schema$Generic <- build_generic_schema(column_schema)

normalize_df <- function(df) {
  #Here we transform the data to only contain the minimum columns required for analysis and rename column names to work with code.
  #To unify columns, I make them all capitalized
  
  #Initiate the dataframe used for mapping the original columns to the ones we need
  original <- data.frame(
    final_name = character(),
    original_name = character(),
    stringsAsFactors = FALSE
  )
  
  software <- detect_software(df, signature)
  
  #Since Fragpipe and PEAKS both contain "Peptide" column, but they mean different things, we need to differentiate them.
  if(software == "Fragpipe") {
    #Since we know it is a Fragpipe file, we can select only the column we are interested in improve speed of subsequent steps.
    columns_to_keep <- c("Spectrum", "Peptide", "Peptide.Sequence", "Modified.Peptide", "Modified.Sequence", "Charge", "Charges", "Retention", "Observed.Mass", "Observed.M.Z", 
                         "PeptideProphet.Probability", "Probability", "Intensity", "Ion.Mobility", "Protein", "Mapped.Proteins","Delta.Mass",
                         "Length", "Peptide.Length", "Assigned.Modifications") #Here include all Columns in peptide.tsv and psm.tsv that could be interesting
    
    if(any(tolower(colnames(df)) == "peptide.sequence")) { #for combined_peptide.tsv, it follows a wide format similar to PEAKS. For Intensity, we use MaxLFQ.Intensity columns.
      columns_to_keep <- c(columns_to_keep, colnames(df)[grepl("maxlfq.intensity", tolower(colnames(df)))])
    }
    
    df <- df %>% dplyr::select(any_of(columns_to_keep)) #Get the columns that actually exist.
    
    #Here we detect the long format used in Frapipe psm.tsv file. So that we convert to wide format
    if(any(tolower(colnames(df)) == "modified.peptide")) { # Only psm.tsv has this column. in combined_modified_peptide.tsv this is called modified.sequence.
      #For the psm.tsv file, we need to convert from long format to wide format for the intensities and filter unique.
      df <- df %>% dplyr::mutate(Modified.Peptide = dplyr::coalesce(na_if(Modified.Peptide, ""), na_if(Peptide, "")))
      
      #This function group_by on rawfile and modified.peptide.
      df <- aggregate_fragpipe_psm(df)
      
      #Now we can get the intensity per sample in wide format
      df_wide <- df %>%
        dplyr::select(Modified.Peptide, RawFile, Intensity) %>%
        pivot_wider(
          names_from  = RawFile,
          values_from = Intensity,
          names_prefix = "Intensity." #Add prefix, so that these are easy to select for Area/Quantity calculations.
        )
      
      df_wide$X.Spec <- rowSums(!is.na(df_wide[,-1])) #Add the info how many samples the peptide was found.
      
      #Finally unique for the representative peptide, using the row with the best PeptideProphet.probability.
      df<- df %>%
        dplyr::group_by(Modified.Peptide) %>%
        dplyr::slice_max(
          .data[[dplyr::intersect(c("PeptideProphet.Probability", "Probability"), names(df))[1]]], #Try "PeptideProphet.Probability", then "Probability" column
          n = 1,
          with_ties = FALSE
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Intensity)
      
      if (nrow(df) != nrow(df_wide)) {
        stop("Error: The number of rows in df and df_wide do not match!")
      }
      
      df <- left_join(df, df_wide, by = "Modified.Peptide")
      
    }
    
    #Here we combine the Proteins and Mapped.Proteins together, to get all Accessions
    df <- df %>%
      dplyr::mutate(
        # Replace commas in Mapped.Protein with ;
        Mapped.Proteins = str_replace_all(Mapped.Proteins, ", ", ";"),
        
        # Combine columns with ";" (remove empty strings to avoid ;;)
        Protein_Mapped.Proteins = str_trim(paste(Protein, Mapped.Proteins, sep = ";")),
        
        # Remove duplicate ;;
        Protein_Mapped.Proteins = str_replace_all(Protein_Mapped.Proteins, ";+", ";"),
        
        # Remove leading/trailing ; if any
        Protein_Mapped.Proteins = str_replace_all(Protein_Mapped.Proteins, "^;|;$", ""),
        
        # Remove sp| from protein names
        Protein_Mapped.Proteins = str_replace_all(Protein_Mapped.Proteins, "sp\\|", "")
      )
    
    if(any(tolower(colnames(df)) == "assigned.modifications")) { #detect fragpipe
      res <- find_transform_column(df, column_schema[[software]][["PTM"]], "PTM")
      df <- res$df
      #For Fragpipe, we need to extract the numbers
      df$PTM <- sapply(df$PTM, function(x) {
        ptms <- stringr::str_extract_all(x, "\\(\\d+\\.\\d+\\)")[[1]]
        ptms <- gsub("[()]", "", ptms)
        if (length(ptms) == 0) return("")
        paste(ptms, collapse = ";")
      })
      original <- log_rename(original, res, "PTM")
    } else { #for peptide.tsv and combined_peptide.tsv we dont have modification, thus we need to assign as NA.
      df$PTM <- NA
    }
    
  }
  
  if(software == "DIANN_parquet") {
    columns_to_keep <- c("Run", "Modified.Sequence", "Stripped.Sequence", "Precursor.Charge", "Precursor.Mz", "Protein.Names", "RT", "IM", "Precursor.Quantity", "Q.Value")
    
    df <- df %>% dplyr::select(any_of(columns_to_keep))
    
    df_wide <- df %>% dplyr::select(Modified.Sequence, Run, Precursor.Quantity) %>%
      pivot_wider(
        names_from  = Run,
        values_from = Precursor.Quantity,
        names_prefix = "Precursor.Quantity.",
        values_fn = ~ if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE),
        values_fill = NA_real_
        )
    
    df_top <- df %>%
      dplyr::group_by(Modified.Sequence) %>%
      dplyr::arrange(Q.Value, .by_group = TRUE) %>%
      dplyr::slice_head(n = 1) %>%
      ungroup()
    
    df <- left_join(df_top %>% dplyr::select(-Run, -any_of("Precursor.Quantity")), df_wide, by = "Modified.Sequence")
    
  }
  
  res <- transform_columns(df, column_schema, software)
  
  df <- res$df
  original <- res$log
  
  #-------------Derive missing columns if possible----------------
  #CHARGE
  if (!"CHARGE" %in% colnames(df)) {
    df$CHARGE <- 0
    message("Charge column missing → set to 0")
    original <- rbind(original, data.frame(final_name = "CHARGE", original_name = "[no CHARGE column]", stringsAsFactors = FALSE))
  } else {
    df$CHARGE <- sapply(as.character(df$CHARGE), function(x) {
      nums <- suppressWarnings(as.integer(unlist(regmatches(x, gregexpr("[0-9]+", x)))))
      nums <- nums[!is.na(nums)]
      if (length(nums) == 0) NA_integer_ else min(nums) #takes the minimum digit, if multiple are given.
    }, USE.NAMES = FALSE)
  }
  #STRIPPED
  if (!"STRIPPED" %in% colnames(df)) {
    df$STRIPPED <- remove_ptms(df$PEPTIDE)
    message("Stripped column missing: manual remove PTMs")
    original <- rbind(original, data.frame(final_name = "STRIPPED", original_name = "[generated from PEPTIDE]", stringsAsFactors = FALSE))
    }
  #PEPTIDE
  if (!"PEPTIDE" %in% colnames(df)) {
    df$PEPTIDE <- df$STRIPPED
    message("Peptidoform column missing: STRIPPED → PEPTIDE")
    original <- rbind(original, data.frame(final_name = "PEPTIDE", original_name = "[= STRIPPED]", stringsAsFactors = FALSE))
    }
  #LENGTH
  if (!"LENGTH" %in% colnames(df)) {
    df$LENGTH <- nchar(df$STRIPPED)
    message("Length column missing → calculated from peptide")
    original <- rbind(original, data.frame(final_name = "LENGTH", original_name = "[Calculated]", stringsAsFactors = FALSE))
  }
  #M/Z
  if (!"MASS" %in% colnames(df) & "CHARGE" %in% colnames(df) & "MZ" %in% colnames(df)) {
    df$MASS <- df$MZ * df$CHARGE
    message("Mass column missing → calculated from m/z and charge")
    original <- rbind(original, data.frame(final_name = "MASS", original_name = "[MZ * CHARGE]", stringsAsFactors = FALSE))
  }

  
  if (all(c("X1.k0.Start", "X1.k0.End") %in% colnames(df))) { #Generate the X1.k0 from two columns in PEAKS 12 and 13 Studio.
    df$K0 <- rowMeans(df[, c("X1.k0.Start", "X1.k0.End")], na.rm = TRUE)
    res <- find_transform_column(df, "K0", "K0") #PEAKS is fine, but if it is FragPipe, 
    df <- res$df
    original <- bind_rows(original, data.frame(final_name="K0", original_name="X1.k0.Start, X1.k0.End"))
  } else {
    res <- find_transform_column(df, column_schema[[software]][["K0"]], "K0")
    if(!is.null(res)) {
      df <- res$df
      original <- log_rename(original, res, "K0")
      if (is.character(df$K0)){ #This is only for PEAKS, since it returns a range.
        df <- df %>% dplyr::mutate(K0 = get_midpoint(K0))
      }
    } else {
      df$K0 <- 0
      message("Ion mobility column missing → set to 0")
      original <- rbind(original, data.frame(final_name = "K0", original_name = "[no Ion Mobility column]", stringsAsFactors = FALSE))
    }
  }
  
  #These two columns are the minimum required. If there is no PTM, PEPTIDE will simply be the same as STRIPPED.
  if (!all(c("PEPTIDE", "STRIPPED") %in% colnames(df))) { 
    stop(sprintf("Either PEPTIDE or STRIPPED column missing in one or more dataframes."))
  }
  
  #-------------Quantity -----------------
  # 1. Check schema-defined candidates first (prefix match to capture multi-sample columns)
  schema_qty <- column_schema[[software]][["QUANTITY"]]
  sample <- integer(0)
  
  if (length(schema_qty) > 0) {
    for (cand in schema_qty) {
      hits <- which(grepl(tolower(cand), tolower(colnames(df)), fixed = TRUE))
      if (length(hits) > 0) { sample <- hits; break }
    }
  }
  
  # 2. Hard-coded fallbacks for DIANN pr/pr_matrix
  if (software %in% c("DIANN")) {
    sample <- which(grepl("^[a-z]\\.\\.", tolower(colnames(df)))) #this is the path to the file.
  }
  
  if (length(sample) == 0) {
    warning("No quantity/intensity columns found. Define QUANTITY in your custom_schema.")
  }
  
  original <- bind_rows(original, data.frame(final_name = "QUANTITY",
                                             original_name = list(colnames(df[, sample, drop = FALSE]))))
  
  df$MAX_QUANTITY <- apply(df[, sample, drop = FALSE], 1, function(x) {
    if (all(is.na(x))) NA else max(x, na.rm = TRUE)
  })
  
  #-------------Spectra -----------------
  schema_spec <- column_schema[[software]][["SPECTRA"]]
  spec <- integer(0)
  
  if (length(schema_spec) > 0 && !identical(sort(schema_spec), sort(schema_qty))) {
    for (cand in schema_spec) {
      hits <- which(grepl(tolower(cand), tolower(colnames(df)), fixed = TRUE) &
                      !seq_along(colnames(df)) %in% sample)
      if (length(hits) > 0) { spec <- hits; break }
    }
  }
  
  if (software == "PEAKS") {
    if (any(startsWith(tolower(colnames(df)),"x.spec."))) {
      spec <- startsWith(tolower(colnames(df)),"x.spec.")
    } else if(any(startsWith(tolower(colnames(df)),"x.feature."))) {
      spec <- startsWith(tolower(colnames(df)),"x.feature.")
    } else {
      warning("PEAKS format detected, but no spectra column found.")
      spec <- integer(0)
    }
  }

  #In case we can use Quantity columns (sample) as indicator of Spectral match
  if (!any(spec)) {
    original <- bind_rows(original, data.frame(final_name="SPECTRA", original_name=list(colnames(df[,sample, drop = FALSE])))) #We use the same columns as Intensity for Spectra found per sample in Fragpipe
  } else {
    df[, spec][df[, spec] == 0] <- NA #For spectra column, the 0 actually means not identified. We need these to be converted to 0 to work with FragPipe format for data completeness plot.
    original <- bind_rows(original, data.frame(final_name="SPECTRA", original_name=colnames(df[, spec, drop = FALSE]), stringsAsFactors=FALSE))
  }
  
  #-------------Finalize dataframe and generate mapping table-----------------
  
  #Only keep cols that were found
  keep_cols <- original %>%
    dplyr::filter(!is.na(original_name)) %>%
    dplyr::filter(purrr::map_lgl(final_name, ~ is.character(.x) && length(.x) == 1)) %>% #This removes entries where the original column names are a list of strings. These are added separately.
    dplyr::pull(final_name)
  keep_cols <- c(keep_cols, colnames(df)[sample], colnames(df)[spec], "STRIPPED", "MAX_QUANTITY")
  
  #Now we can prepare the summary table since we have the columns to keep.
  original <- original %>% dplyr::mutate(coalesced = do.call(coalesce, across(-1))) %>% dplyr::select(1, coalesced) %>% dplyr::group_by(final_name) %>% dplyr::summarise(coalesced_list = list(coalesced), .groups = "drop")
  
  return(list(df = df %>% dplyr::select(any_of(unique(keep_cols))), table = original))
}

# @param lst      Named list of data.frames (active_data_list())
# @param filters  Named list of slider values from input$*
# @return         Filtered named list of data.frames
apply_filters <- function(lst, filters) {
  lapply(lst, function(df) {
    
    # Length filter
    if (!is.null(filters$length_range) && "LENGTH" %in% names(df)) {
      df <- df[df$LENGTH >= filters$length_range[1] & df$LENGTH <= filters$length_range[2], ]
    }
    
    # Quantity filter
    if (!is.null(filters$quantity_range) && "MAX_QUANTITY" %in% names(df) &&
        length(filters$quantity_range) == 2) {
      df <- df[df$MAX_QUANTITY >= filters$quantity_range[1] &
                 df$MAX_QUANTITY <= filters$quantity_range[2], ]
    }
    
    # Score filter
    if (!is.null(filters$score_range) && "SCORE" %in% names(df)) {
      df <- df[df$SCORE >= filters$score_range[1] & df$SCORE <= filters$score_range[2], ]
    }
    
    # Charge filter
    if (!is.null(filters$charge_range) && "CHARGE" %in% names(df)) {
      df <- df[df$CHARGE >= filters$charge_range[1] & df$CHARGE <= filters$charge_range[2], ]
    }
    
    # Mass filter
    if (!is.null(filters$mass_range) && "MASS" %in% names(df)) {
      df <- df[df$MASS >= filters$mass_range[1] & df$MASS <= filters$mass_range[2], ]
    }
    
    # RT filter
    if (!is.null(filters$RT_range) && "RT" %in% names(df)) {
      df <- df[df$RT >= filters$RT_range[1] & df$RT <= filters$RT_range[2], ]
    }
    
    # Remove all-NA rows
    df[rowSums(!is.na(df)) > 0, , drop = FALSE]
  })
}

# @param lst        Named list of data.frames
# @param col        Column name to compute range from
# @param input_id   Shiny inputId for the sliderInput
# @param label      Display label (can be tagList with icon)
# @param step       Slider step (default NULL = continuous)
# @return           sliderInput widget or a styled error div
make_range_slider_ui <- function(lst, col, input_id, label, step = NULL, digits = NULL) {
  if (is.null(lst) || !any(sapply(lst, function(df) col %in% colnames(df)))) {
    return(tags$div(
      style = "color: #b30000; font-style: italic;",
      paste0("No ", col, " column in data.")
    ))
  }
  
  min_val <- min(sapply(lst, function(df) if (col %in% names(df)) min(df[[col]], na.rm = TRUE) else Inf),  na.rm = TRUE)
  max_val <- max(sapply(lst, function(df) if (col %in% names(df)) max(df[[col]], na.rm = TRUE) else -Inf), na.rm = TRUE)
  
  # Guard against degenerate range (e.g. single unique value)
  if (!is.finite(min_val) || !is.finite(max_val) || min_val == max_val) {
    max_val <- min_val + 1
  }
  
  if (!is.null(digits)) {
    min_val <- floor(min_val * 10^digits) / 10^digits
    max_val <- ceiling(max_val * 10^digits) / 10^digits
    if (is.null(step)) step <- 10^(-digits)
  }
  
  sliderInput(input_id, label,
              min   = min_val,
              max   = max_val,
              value = c(min_val, max_val),
              step  = step)
}

#This is for creating a matrix used for the group-based peptide heatmap
prepare_peptide_matrix <- function(lst, groups, quantity_cols, group_peptide_sets, col_map = NULL) {
  
  peptides_all <- unique(unlist(group_peptide_sets))
  pep_mat <- matrix(NA_real_,
                    nrow = length(peptides_all),
                    ncol = length(groups),
                    dimnames = list(peptides_all, names(groups)))
  
  for (g in names(groups)) {
    group_items      <- groups[[g]]
    allowed_peptides <- group_peptide_sets[[g]]
    
    if (!is.null(col_map) && is.list(group_items) && !is.null(group_items$measurement)) {
      
      
      # Measurement mode: grp_items = list(name=c(...), measurement=c(...))
      measurements <- group_items$measurement
      names_vec    <- group_items$name
      df_group <- dplyr::bind_rows(lapply(seq_along(measurements), function(i) {
        nm   <- names_vec[i]
        meas <- measurements[i]
        entry <- col_map[[meas]]
        if (is.null(entry)) return(NULL)
        df <- lst[[nm]]
        if (is.null(df) || !"STRIPPED" %in% colnames(df)) return(NULL)
        actual_col <- colnames(df)[tolower(colnames(df)) == tolower(entry$col)][1]
        if (is.na(actual_col)) return(NULL)
        df_sub <- df[df$STRIPPED %in% allowed_peptides, c("STRIPPED", actual_col), drop = FALSE]
        names(df_sub)[names(df_sub) == actual_col] <- "Quantity"
        df_sub$Quantity[df_sub$Quantity == 0] <- NA #set 0 to NA
        df_sub
      }))
      
      if (is.null(df_group) || nrow(df_group) == 0) next
      
      agg <- df_group %>%
        dplyr::group_by(STRIPPED) %>%
        dplyr::summarise(value = max(Quantity, na.rm = TRUE), .groups = "drop")
    } else {
      # Sample mode: each item is a sample name
      df_group <- dplyr::bind_rows(lapply(group_items, function(s) {
        df <- lst[[s]]
        cols <- dplyr::intersect(c("STRIPPED", quantity_cols), colnames(df))
        if (length(cols) < 2) return(NULL)
        df[, cols, drop = FALSE]
      }))
      
      if (is.null(df_group) || nrow(df_group) == 0) next
      
      df_group <- df_group[df_group$STRIPPED %in% allowed_peptides, ]
      
      # Vectorized row-wise max — avoids c_across row-by-row overhead
      qty_present  <- intersect(quantity_cols, colnames(df_group))
      qty_mat      <- as.matrix(df_group[, qty_present, drop = FALSE])
      df_group$value <- do.call(pmax, c(as.data.frame(qty_mat), list(na.rm = TRUE)))
      
      agg <- df_group[, c("STRIPPED", "value")] |>
        dplyr::group_by(STRIPPED) |>
        dplyr::summarise(value = max(value, na.rm = TRUE), .groups = "drop")
    }
    
    pep_mat[agg$STRIPPED, g] <- agg$value
  }
  
  # Drop peptides not identified in any group (i.e. NA in both groups)
  pep_mat <- pep_mat[rowSums(!is.na(pep_mat)) > 0, , drop = FALSE]
  
  pep_mat
}

## This is for heatmap generation of the grouped comparison where each column is a peptide. 
## In cases of very large number of peptides, this is not possible to show, thus we "bin" them to improve visualization.
## The bins sizes are evenly distributed.
bin_heatmap_columns <- function(mat, max_cols = 1000) {
  if (ncol(mat) <= max_cols) return(mat)
  n <- ncol(mat)
  k <- min(max_cols, n)
  mat_imp <- mat
  mat_imp[is.na(mat_imp)] <- 0
  bin_map <- list() #this is for the interactive heatmap
  
  presence_key   <- apply(mat_imp > 0, 2, function(x) paste(as.integer(x), collapse = ""))
  patterns       <- unique(presence_key)
  pattern_counts <- setNames(sapply(patterns, function(p) sum(presence_key == p)), patterns)
  
  # Proportional allocation, immediately capped at actual peptide count
  pattern_bins <- pmax(1L, floor(k * pattern_counts / n))
  pattern_bins <- pmin(pattern_bins, pattern_counts)
  names(pattern_bins) <- patterns
  
  # Redistribute remaining bins only to patterns that still have room
  remaining <- k - sum(pattern_bins)
  while (remaining > 0) {
    can_expand <- patterns[pattern_counts[patterns] > pattern_bins[patterns]]
    if (length(can_expand) == 0) break
    n_add  <- min(remaining, length(can_expand))
    add_to <- can_expand[order(pattern_counts[can_expand] - pattern_bins[can_expand],
                               decreasing = TRUE)][seq_len(n_add)]
    pattern_bins[add_to] <- pattern_bins[add_to] + 1L
    remaining <- remaining - n_add
  }
  
  all_bins <- unlist(lapply(patterns, function(pat) {
    idx   <- which(presence_key == pat)
    n_pat <- length(idx)
    k_pat <- min(as.integer(pattern_bins[[pat]]), n_pat)
    if (is.na(k_pat) || k_pat < 1L) k_pat <- 1L
    
    pep_names <- colnames(mat)[idx]
    mat_pat <- mat[, idx, drop = FALSE]
    
    if (n_pat == 1 || k_pat == 1) {
      result <- rowMeans(mat_pat, na.rm = TRUE)
      result[is.nan(result)] <- NA
      bin_name <- paste0("bin_", pat, "_1")
      bin_map[[bin_name]] <<- pep_names
      return(list(matrix(result, ncol = 1,
                         dimnames = list(rownames(mat), paste0("bin_", pat, "_1")))))
    }
    
    pc1 <- tryCatch(
      prcomp(t(mat_imp[, idx, drop = FALSE]), center = TRUE, scale. = FALSE)$x[, 1],
      error = function(e) seq_len(n_pat)
    )
    
    ord       <- order(pc1)
    mat_pat <- mat_pat[, order(pc1), drop = FALSE]
    pep_names <- pep_names[ord]
    
    bin_ids <- ceiling(seq_len(n_pat) * k_pat / n_pat)
    lapply(seq_len(k_pat), function(i) {
      cols <- which(bin_ids == i)
      if (length(cols) == 1) return(mat_pat[, cols, drop = FALSE])
      result <- rowMeans(mat_pat[, cols, drop = FALSE], na.rm = TRUE)
      result[is.nan(result)] <- NA
      matrix(result, ncol = 1,
             dimnames = list(rownames(mat), paste0("bin_", pat, "_", i)))
    })
  }), recursive = FALSE)
  
  result <- do.call(cbind, all_bins)
  attr(result, "bin_map") <- bin_map
  result
}

prepare_measurement_matrix <- function(lst, quantity_cols) {
  
  peptides_all <- unique(unlist(lapply(lst, function(df) df$PEPTIDE)))
  
  col_names <- unlist(lapply(names(lst), function(nm) {
    cols <- dplyr::intersect(quantity_cols, colnames(lst[[nm]]))
    paste0(nm, " | ", cols)
  }))
  
  pep_mat <- matrix(NA_real_,
                    nrow = length(peptides_all),
                    ncol = length(col_names),
                    dimnames = list(peptides_all, col_names))
  
  safe_max <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else max(x)
  }
  
  for (nm in names(lst)) {
    df   <- lst[[nm]]
    cols <- dplyr::intersect(quantity_cols, colnames(df))
    if (length(cols) == 0 || !"PEPTIDE" %in% colnames(df)) next
    
    df_sub <- df[, c("PEPTIDE", cols), drop = FALSE]
    
    for (col in cols) {
      agg <- df_sub %>%
        dplyr::group_by(PEPTIDE) %>%
        dplyr::summarise(value = safe_max(.data[[col]]), .groups = "drop")
      
      pep_mat[agg$PEPTIDE, paste0(nm, " | ", col)] <- agg$value
    }
  }
  
  # Remove peptides with no finite measurements across any column
  pep_mat <- pep_mat[rowSums(is.finite(pep_mat)) > 0, , drop = FALSE]
  
  pep_mat
}

#-----------Plotting functions------------------
plot_unique_counts <- function(lst, column, y_label, transform_fn = identity, color = "default", quantity_cols = NULL) {
  
  stats <- dplyr::bind_rows(lapply(names(lst), function(s) {
    df <- lst[[s]]
    if (!column %in% colnames(df)) return(NULL)
    
    if (!is.null(quantity_cols)) {
      meas_cols <- intersect(quantity_cols, colnames(df))
      
      if (length(meas_cols) == 0) {
        # No measurement columns in this sample — fall back to total unique count
        return(data.frame(Sample = s, Mean = length(unique(transform_fn(df[[column]]))),
                          Min = NA_real_, Max = NA_real_, stringsAsFactors = FALSE))
      }
      
      counts <- sapply(meas_cols, function(m) {
        col <- df[[m]]
        detected <- if (any(is.na(col))) {
          !is.na(col)       # NA = not found; 0 and >0 both mean found
        } else {
          col > 0           # no NAs → 0 means not found
        }
        length(unique(transform_fn(df[[column]][detected])))
      })
      
      data.frame(Sample = s, Mean = mean(counts), Min = min(counts), Max = max(counts),
                 stringsAsFactors = FALSE)
    } else {
      data.frame(Sample = s, Mean = length(unique(transform_fn(df[[column]]))),
                 Min = NA_real_, Max = NA_real_, stringsAsFactors = FALSE)
    }
  }))
  
  if (is.null(stats) || nrow(stats) == 0) return(NULL)
  
  fill_colors <- if (color == "default") NULL else viridis(nrow(stats), option = color)
  
  has_range <- any(!is.na(stats$Min) & stats$Min != stats$Max)
  
  p <- ggplot(stats, aes(x = Sample, y = Mean, fill = Sample)) +
    geom_col() +
    labs(x = "Sample", y = y_label) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  if (has_range) {
    p <- p +
      geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.25, linewidth = 0.8) +
      geom_text(aes(label = sprintf("%.0f\n[%d\u2013%d]", Mean, Min, Max)),
                vjust = -0.3, size = 3.5)
  } else {
    p <- p + geom_text(aes(label = round(Mean)), vjust = -0.3, size = 5)
  }
  
  if (!is.null(fill_colors)) p <- p + scale_fill_manual(values = fill_colors)
  p
}


plot_histogram <- function(df, column, x_label, title_name = "RT plot", color = "default") {
  
  # Create 1-minute bins (change to 60 if your data is in seconds)
  vals   <- df[[column]][is.finite(df[[column]])]
  breaks <- pretty(vals, n = 50)

  
  if (color == "default") {
    color <- "grey35"  # ggplot will use default fill colors
  } else {
    color <- viridis(1, option = color)
  }
  
  p <- ggplot(df, aes_string(x = column)) +   # <- use column, not data_col
    geom_histogram(
      breaks = breaks,
      fill = color,
      color = "black"
    ) +
    labs(
      title = title_name,
      x = x_label,
      y = "Count"
    ) +
    theme_minimal()
  
  return(p)
}

plot_stacked_bar <- function(lst, column, fill_label = NULL, rev_levels = TRUE, percentage = FALSE, color = "default") {
  combined <- bind_rows(lapply(names(lst), function(name) {
    df <- lst[[name]]
    if (!column %in% colnames(df)) return(NULL)  # skip samples missing the column
    df <- df[, column, drop = FALSE]    # only the column needed
    if (nrow(df) == 0) {
      return(data.frame(Sample = character(0), df[0, , drop = FALSE]))
    }
    df$Sample <- name
    df
  }))
  
  if (is.null(combined) || nrow(combined) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = paste("No data for", column),
                 size = 5, color = "grey50") +
        theme_void()
    )
  }
  
  # Convert to factor and optionally reverse levels
  if (rev_levels) {
    combined[[column]] <- factor(combined[[column]], levels = rev(sort(unique(combined[[column]]))))
  } else {
    combined[[column]] <- factor(combined[[column]], levels = sort(unique(combined[[column]])))
  }
  
  if (color == "default") {
    color <- NULL  # ggplot will use default fill colors
  } else {
    color <- viridis(length(unique(combined[[column]])), option = color)
  }

  # Compute counts or percentages
  if (percentage) {
    combined <- combined %>%
      dplyr::group_by(Sample) %>%
      dplyr::mutate(perc = 100 * (n() / n())) %>%  # initially gives 100%, will correct below
      ungroup()
    
    # Actually, we need % per factor level per sample
    perc_df <- combined %>%
      dplyr::group_by(Sample, !!sym(column)) %>%
      dplyr::summarise(count = n(), .groups = "drop") %>%
      dplyr::group_by(Sample) %>%
      dplyr::mutate(perc = 100 * count / sum(count)) %>%
      ungroup()
    
    p <- ggplot(perc_df, aes(x = Sample, y = perc, fill = !!sym(column))) +
      geom_bar(stat = "identity", position = "stack") +
      geom_text_repel(aes(label = round(perc, 1)), position=position_stack(vjust = 0.5), direction="y") +
      labs(x = "Sample", y = "Percentage", fill = fill_label) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    
  } else {
    # Default: count
    p <- ggplot(combined, aes_string(x = "Sample", fill = column)) +
      geom_bar(position = "stack") +
      geom_text_repel(stat = "count", aes(label = after_stat(count)), position=position_stack(vjust = 0.5), direction="y") +
      labs(x = "Sample", y = "Count", fill = fill_label) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  if (!is.null(color)) {
    p <- p + scale_fill_manual(values = color)
  }
  
  return(p)
}

plot_density <- function(lst, column, transform = NULL, x_label = NULL, alpha = 0.3, color = "default") {
  if (is.null(x_label)) x_label <- column
  n          <- length(lst)
  plotly_pal <- c('#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A',
                  '#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52')
  cols <- if (color == "default") {
    setNames(plotly_pal[((seq_len(n) - 1L) %% length(plotly_pal)) + 1L], names(lst))
  } else {
    setNames(viridis(n, option = color), names(lst))
  }
  
  densities <- lapply(names(lst), function(sname) {
    vals <- lst[[sname]][[column]]
    if (!is.null(transform)) vals <- transform(vals)
    vals <- vals[!is.na(vals) & is.finite(vals)]
    if (length(vals) < 2L) return(NULL)
    list(vals = vals, dens = density(vals))
  })
  names(densities) <- names(lst)
  
  all_y_max <- max(sapply(densities, function(d) if (is.null(d)) 0 else max(d$dens$y)), na.rm = TRUE)
  row_h     <- all_y_max * 0.04
  rug_ys    <- setNames(-row_h * seq_len(n), names(lst))
  
  p <- plotly::plot_ly()
  for (sname in names(lst)) {
    d <- densities[[sname]]
    if (is.null(d)) next
    p <- p %>%
      plotly::add_lines(
        x = d$dens$x, y = d$dens$y, name = sname,
        line = list(color = cols[[sname]], width = 2),
        legendgroup = sname, showlegend = TRUE,
        hovertemplate = paste0("<b>", sname, "</b>: %{y:.4g}<extra></extra>")
      ) %>%
      plotly::add_trace(
        x = d$vals, y = rep(rug_ys[[sname]], length(d$vals)),
        type = "scatter", mode = "markers", name = sname,
        marker = list(symbol = "line-ns-open", size = 8,
                      color = cols[[sname]], opacity = 0.5),
        legendgroup = sname, showlegend = FALSE,
        hoverinfo = "none"
      )
  }
  
  p %>% plotly::layout(
    xaxis = list(title = x_label),
    yaxis = list(
      title    = "Density",
      range    = c(-(n + 0.5) * row_h, all_y_max * 1.05),
      zeroline = TRUE
    ),
    hovermode = "x unified",
    legend = list(title = list(text = "Sample"))
  )
}

plot_violin <- function(lst, column, x_label = "Sample", y_label = NULL, title = NULL, add_boxplot = TRUE, color = "default") {
  if (is.null(y_label)) y_label <- column
  n          <- length(lst)
  plotly_pal <- c('#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A',
                  '#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52')
  cols <- if (color == "default") {
    setNames(plotly_pal[((seq_len(n) - 1L) %% length(plotly_pal)) + 1L], names(lst))
  } else {
    setNames(viridis(n, option = color), names(lst))
  }
  
  all_vals <- unlist(lapply(lst, function(df) {
    vals <- df[[column]]
    vals[!is.na(vals) & is.finite(vals)]
  }))
  y_range <- range(all_vals, na.rm = TRUE)
  
  p <- plotly::plot_ly()
  for (i in seq_along(lst)) {
    sname <- names(lst)[i]
    vals  <- lst[[sname]][[column]]
    vals  <- vals[!is.na(vals) & is.finite(vals)]
    if (length(vals) == 0L) next
    p <- p %>% plotly::add_trace(
      x = rep(sname, length(vals)),
      y = vals, type = "violin", name = sname,
      fillcolor = adjustcolor(cols[[i]], alpha.f = 0.5),
      line = list(color = cols[[i]]),
      box = list(visible = add_boxplot, fillcolor = "white",
                 line = list(color = cols[[i]])),
      meanline = list(visible = TRUE, color = cols[[i]]),
      points = FALSE
    )
  }
  
  p %>% plotly::layout(
    xaxis = list(title = x_label, tickangle = -45),
    yaxis = list(title = y_label, range = y_range),
    showlegend = FALSE
  )
}

plot_length_distribution <- function(lst, color = "default") {
  req(lst)
  
  if (color == "default") {
    color <- NULL  # ggplot will use default fill colors
  } else {
    cols <- viridis(length(lst), option = color)
    color <- setNames(cols, names(lst))
  }
  
  # Combine datasets into one long dataframe
  combined <- bind_rows(lapply(names(lst), function(name) {
    df <- lst[[name]]
    if (!"LENGTH" %in% colnames(df)) return(NULL)
    df <- df[, "LENGTH", drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df$Sample <- name
    df
  }))
  req(nrow(combined) > 0)
  
  combined$LENGTH <- as.numeric(combined$LENGTH)
  combined <- combined %>% dplyr::mutate(Length_bin = factor(LENGTH))
  
  # Compute counts and percentages
  count_df <- combined %>% dplyr::count(Length_bin, Sample, name = "Count")
  percent_df <- count_df %>% dplyr::group_by(Sample) %>%
    dplyr::mutate(Percent = Count / sum(Count) * 100)
  
  samples <- unique(count_df$Sample)
  n <- length(samples)
  
  fig <- plot_ly()
  
  # --- Add COUNT traces first ---
  for (i in seq_along(samples)) {
    s <- samples[i]
    sub_c <- count_df[count_df$Sample == s, ]
    sub_p <- percent_df[percent_df$Sample == s, ]
    if (nrow(sub_c) == 0 || nrow(sub_p) == 0) next #checkpoint for when one sample has 0 data.
    
    fig <- fig %>%
      add_bars(
        x = sub_c$Length_bin,
        y = sub_c$Count,
        name = s,
        legendgroup = s,
        showlegend = TRUE,
        visible = TRUE,
        customdata = round(sub_p$Percent, 2),
        meta = s,
        marker = if (!is.null(color)) list(color = color[s]) else NULL,
        hovertemplate =
          "Sample: %{meta}<br>Length: %{x}<br>Count: %{y}<br>Percent: %{customdata}%<extra></extra>"
      )
  }
  
  # --- Add PERCENT traces next ---
  for (i in seq_along(samples)) {
    s <- samples[i]
    sub_c <- count_df[count_df$Sample == s, ]
    sub_p <- percent_df[percent_df$Sample == s, ]
    if (nrow(sub_c) == 0 || nrow(sub_p) == 0) next #checkpoint for when one sample has 0 data.
    
    fig <- fig %>%
      add_bars(
        x = sub_p$Length_bin,
        y = sub_p$Percent,
        name = s,
        legendgroup = s,
        showlegend = TRUE,
        visible = FALSE,
        customdata = sub_c$Count,
        meta = s,
        marker = if (!is.null(color)) list(color = color[s]) else NULL,
        hovertemplate =
          "Sample: %{meta}<br>Length: %{x}<br>Percent: %{y:.2f}%<br>Count: %{customdata}<extra></extra>"
      )
  }
  
  # --- Layout and buttons ---
  fig <- fig %>%
    layout(
      barmode = "group",
      xaxis = list(title = "Peptide Length"),
      yaxis = list(title = "Count"),
      updatemenus = list(
        list(
          type = "buttons",
          direction = "right",
          x = 1.05, y = 1.15,
          buttons = list(
            list(
              label = "Count",
              method = "update",
              args = list(
                list(visible = c(rep(TRUE, n), rep(FALSE, n))),
                list(yaxis = list(title = "Count"))
              )
            ),
            list(
              label = "Percent",
              method = "update",
              args = list(
                list(visible = c(rep(FALSE, n), rep(TRUE, n))),
                list(yaxis = list(title = "Percent (%)"))
              )
            )
          )
        )
      )
    )
  fig
}

plot_length_range_per_measurement <- function(lst, quantity_cols, lo = 8, hi = 13, color = "default") {
  
  stats <- dplyr::bind_rows(lapply(names(lst), function(s) {
    df <- lst[[s]]
    if (!"LENGTH" %in% colnames(df)) return(NULL)
    
    meas_cols <- intersect(quantity_cols, colnames(df))
    
    if (length(meas_cols) == 0) {
      pct <- mean(df$LENGTH >= lo & df$LENGTH <= hi, na.rm = TRUE) * 100
      return(data.frame(Sample = s, Mean = pct, Min = NA_real_, Max = NA_real_,
                        stringsAsFactors = FALSE))
    }
    
    pcts <- sapply(meas_cols, function(m) {
      col      <- df[[m]]
      detected <- if (any(is.na(col))) !is.na(col) else col > 0
      rows     <- df[detected, ]
      if (nrow(rows) == 0) return(NA_real_)
      mean(rows$LENGTH >= lo & rows$LENGTH <= hi, na.rm = TRUE) * 100
    })
    
    pcts <- pcts[!is.na(pcts)]
    if (length(pcts) == 0) return(NULL)
    
    data.frame(Sample = s, Mean = mean(pcts), Min = min(pcts), Max = max(pcts),
               stringsAsFactors = FALSE)
  }))
  
  if (is.null(stats) || nrow(stats) == 0) return(NULL)
  
  has_range <- any(!is.na(stats$Min) & stats$Min != stats$Max)
  fill_colors <- if (color == "default") NULL else viridis(nrow(stats), option = color)
  
  p <- ggplot(stats, aes(x = Sample, y = Mean, fill = Sample)) +
    geom_col() +
    labs(x = "Sample", y = paste0("% peptides ", lo, "\u2013", hi, "mer")) +
    ylim(0, 115) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  
  if (has_range) {
    p <- p +
      geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.25, linewidth = 0.8) +
      geom_text(aes(y = Max, label = sprintf("%.1f%%\n[%.1f\u2013%.1f]", Mean, Min, Max)),
                vjust = -0.3, size = 3.5)
  } else {
    p <- p + geom_text(aes(label = sprintf("%.1f%%", Mean)), vjust = -0.3, size = 5)
  }
  
  if (!is.null(fill_colors)) p <- p + scale_fill_manual(values = fill_colors)
  p
}


plot_seqlogo <- function(peptides, title = "Motif", namespace = NULL) {
  peptides <- peptides[!is.na(peptides)]
  if (length(peptides) == 0) return(NULL)
  
  L <- unique(nchar(peptides))
  if (length(L) != 1) return(NULL)  # must be same length
  
  df <- data.frame(seq = peptides)
  
  if (is.null(namespace)) {
    suppressWarnings(ggseqlogo(peptides, method = 'bits') +
      ggtitle(title) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
      ))
  } else {
    suppressWarnings(ggseqlogo(peptides, method = 'bits', namespace = namespace, seq_type = "other") +
      ggtitle(title) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
      ))
  }
}

extract_legend <- function(p) {
  tmp <- ggplot2::ggplotGrob(p)
  leg <- gtable::gtable_filter(tmp, "guide-box")
  grid::grid.newpage()
  leg
}

plot_charge_per_measurement <- function(lst, quantity_cols, column = "CHARGE",
                                        fill_label = "Charge", percentage = FALSE,
                                        color = "default") {
  rows <- list()
  for (sname in names(lst)) {
    df        <- lst[[sname]]
    if (!column %in% colnames(df)) next
    meas_cols <- intersect(quantity_cols, colnames(df))
    if (length(meas_cols) == 0) next
    for (m in meas_cols) {
      col_vals <- df[[m]]
      detected <- if (any(is.na(col_vals))) !is.na(col_vals) else col_vals > 0
      sub_df   <- df[detected, column, drop = FALSE]
      if (nrow(sub_df) == 0) next
      counts <- as.data.frame(table(sub_df[[column]]), stringsAsFactors = FALSE)
      colnames(counts) <- c("Charge", "Count")
      counts$Sample      <- sname
      counts$Measurement <- m
      rows[[length(rows) + 1]] <- counts
    }
  }
  if (length(rows) == 0) return(NULL)
  
  df_long <- do.call(rbind, rows)
  if (percentage) {
    df_long <- df_long %>%
      dplyr::group_by(Sample, Measurement) %>%
      dplyr::mutate(Count = 100 * Count / sum(Count)) %>%
      dplyr::ungroup()
  }
  
  # Preserve sample order from lst, measurements in the order they appear
  meas_order <- unlist(lapply(names(lst), function(s) {
    intersect(quantity_cols, colnames(lst[[s]]))
  }))
  meas_order <- meas_order[meas_order %in% df_long$Measurement]
  
  charges    <- sort(unique(df_long$Charge), decreasing = FALSE)
  n_ch       <- length(charges)
  plotly_pal <- c('#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A',
                  '#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52')
  charge_cols <- if (color == "default") {
    setNames(plotly_pal[((seq_len(n_ch) - 1L) %% length(plotly_pal)) + 1L], charges)
  } else {
    setNames(viridis(n_ch, option = color), charges)
  }
  
  meas_df <- data.frame(Measurement = meas_order, stringsAsFactors = FALSE)
  p <- plotly::plot_ly()
  for (ch in charges) {
    d <- merge(meas_df,
               df_long[df_long$Charge == ch, c("Measurement", "Count")],
               by = "Measurement", all.x = TRUE)
    d$Count[is.na(d$Count)] <- 0
    d <- d[match(meas_order, d$Measurement), ]
    p <- p %>% plotly::add_trace(
      type = "bar",
      x    = meas_order,
      y    = d$Count,
      name = paste(fill_label, ch),
      marker = list(color = charge_cols[[ch]]),
      hovertemplate = paste0(fill_label, " ", ch, ": %{y}<extra>%{x}</extra>")
    )
  }
  
  # Sample group labels + dividers
  annotations <- list()
  shapes      <- list()
  for (i in seq_along(names(lst))) {
    sname     <- names(lst)[i]
    meas_cols <- intersect(quantity_cols, colnames(lst[[sname]]))
    meas_cols <- meas_cols[meas_cols %in% meas_order]
    if (length(meas_cols) == 0) next
    positions <- which(meas_order %in% meas_cols) - 1L  # 0-based for plotly
    annotations[[length(annotations) + 1]] <- list(
      x = mean(positions), y = -0.18,
      xref = "x", yref = "paper",
      text = paste0("<b>", sname, "</b>"),
      showarrow = FALSE, xanchor = "center",
      font = list(size = 11)
    )
    if (i > 1) {
      shapes[[length(shapes) + 1]] <- list(
        type = "line",
        x0 = min(positions) - 0.5, x1 = min(positions) - 0.5,
        y0 = 0, y1 = 1,
        xref = "x", yref = "paper",
        line = list(color = "grey60", width = 1, dash = "dot")
      )
    }
  }
  
  p %>% plotly::layout(
    barmode = "stack",
    xaxis   = list(title = "", tickangle = -45, tickfont = list(size = 8)),
    yaxis   = list(title = if (percentage) "Percentage (%)" else "Count"),
    legend  = list(title = list(text = fill_label)),
    margin  = list(b = 100),
    annotations = annotations,
    shapes      = shapes
  )
}

plot_density_envelope <- function(lst, quantity_cols, column, x_label, color = "default") {
  result <- dplyr::bind_rows(lapply(names(lst), function(s) {
    df        <- lst[[s]]
    meas_cols <- intersect(quantity_cols, colnames(df))
    if (length(meas_cols) == 0 || !column %in% colnames(df)) return(NULL)
    all_vals <- df[[column]][is.finite(df[[column]])]
    if (length(all_vals) < 2) return(NULL)
    x_grid <- seq(min(all_vals), max(all_vals), length.out = 512)
    dens_mat <- do.call(rbind, Filter(Negate(is.null), lapply(meas_cols, function(m) {
      col_vals <- df[[m]]
      detected <- if (any(is.na(col_vals))) !is.na(col_vals) else col_vals > 0
      vals     <- df[[column]][detected & is.finite(df[[column]])]
      if (length(vals) < 2) return(NULL)
      d <- density(vals, from = min(x_grid), to = max(x_grid), n = 512)
      approx(d$x, d$y, xout = x_grid)$y
    })))
    if (is.null(dens_mat) || nrow(dens_mat) == 0) return(NULL)
    data.frame(
      Sample = s, x = x_grid,
      Mean   = colMeans(dens_mat, na.rm = TRUE),
      ymin   = apply(dens_mat, 2, min, na.rm = TRUE),
      ymax   = apply(dens_mat, 2, max, na.rm = TRUE)
    )
  }))
  if (is.null(result) || nrow(result) == 0) return(NULL)
  
  samples    <- unique(result$Sample)
  n          <- length(samples)
  plotly_pal <- c('#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A',
                  '#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52')
  cols <- if (color == "default") {
    setNames(plotly_pal[((seq_len(n) - 1L) %% length(plotly_pal)) + 1L], samples)
  } else {
    setNames(viridis(n, option = color), samples)
  }
  
  fmt <- function(x) formatC(x, digits = 3, format = "g")
  
  p <- plotly::plot_ly()
  for (s in samples) {
    d     <- result[result$Sample == s, ]
    col   <- cols[[s]]
    rgb_v <- col2rgb(col)[, 1]
    fill  <- sprintf("rgba(%d,%d,%d,0.15)", rgb_v[1], rgb_v[2], rgb_v[3])
    d$range_text <- paste0(fmt(d$ymin), " \u2013 ", fmt(d$ymax))
    
    p <- p %>%
      plotly::add_ribbons(
        data = d, x = ~x, ymin = ~ymin, ymax = ~ymax,
        customdata = ~range_text,
        name = s, legendgroup = s, showlegend = TRUE,
        fillcolor = fill,
        line = list(color = "rgba(0,0,0,0)"),
        hovertemplate = paste0("<b>", s, "</b> range: %{customdata}<extra></extra>")
      ) %>%
      plotly::add_lines(
        data = d, x = ~x, y = ~Mean,
        name = s, legendgroup = s, showlegend = FALSE,
        line = list(color = col, width = 1.5),
        hovertemplate = paste0("<b>", s, "</b> mean: %{y:.4g}<extra></extra>")
      )
  }
  
  p %>% plotly::layout(
    xaxis = list(title = x_label),
    yaxis = list(title = "Density", rangemode = "tozero"),
    hovermode = "x unified",
    legend = list(title = list(text = "Sample"))
  )
}

plot_violin_envelope <- function(lst, quantity_cols, column, x_label = "Sample", y_label = NULL, color = "default") {
  if (is.null(y_label)) y_label <- column
  n          <- length(lst)
  plotly_pal <- c('#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A',
                  '#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52')
  cols <- if (color == "default") {
    setNames(plotly_pal[((seq_len(n) - 1L) %% length(plotly_pal)) + 1L], names(lst))
  } else {
    setNames(viridis(n, option = color), names(lst))
  }
  
  p <- plotly::plot_ly()
  
  for (i in seq_along(lst)) {
    sname     <- names(lst)[i]
    df        <- lst[[sname]]
    meas_cols <- intersect(quantity_cols, colnames(df))
    if (length(meas_cols) == 0 || !column %in% colnames(df)) next
    
    all_vals <- df[[column]][is.finite(df[[column]])]
    if (length(all_vals) < 2) next
    y_grid <- seq(min(all_vals), max(all_vals), length.out = 256)
    
    dens_list <- Filter(Negate(is.null), lapply(meas_cols, function(m) {
      col_vals <- df[[m]]
      detected <- if (any(is.na(col_vals))) !is.na(col_vals) else col_vals > 0
      vals     <- df[[column]][detected & is.finite(df[[column]])]
      if (length(vals) < 2) return(NULL)
      d <- density(vals, from = min(y_grid), to = max(y_grid), n = 256)
      approx(d$x, d$y, xout = y_grid)$y
    }))
    if (length(dens_list) == 0) next
    
    dens_mat <- do.call(rbind, dens_list)
    mean_d   <- colMeans(dens_mat, na.rm = TRUE)
    min_d    <- apply(dens_mat, 2, min, na.rm = TRUE)
    max_d    <- apply(dens_mat, 2, max, na.rm = TRUE)
    
    scale_f  <- 0.4 / max(max_d, na.rm = TRUE)
    mean_d   <- mean_d * scale_f
    min_d    <- min_d  * scale_f
    max_d    <- max_d  * scale_f
    
    col   <- cols[[sname]]
    rgb_v <- col2rgb(col)[, 1]
    fill  <- sprintf("rgba(%d,%d,%d,0.2)", rgb_v[1], rgb_v[2], rgb_v[3])
    
    # Mirrored polygon helpers
    violin_x <- function(d) c(i + d, i - rev(d))
    violin_y <- function()  c(y_grid, rev(y_grid))
    
    p <- p %>%
      # Shaded envelope (max extent across measurements)
      plotly::add_trace(
        type = "scatter", mode = "lines",
        x = violin_x(max_d), y = violin_y(),
        fill = "toself", fillcolor = fill,
        line = list(color = "rgba(0,0,0,0)"),
        name = sname, legendgroup = sname, showlegend = TRUE,
        hovertemplate = "<extra></extra>"
      ) %>%
      # Mean violin outline
      plotly::add_trace(
        type = "scatter", mode = "lines",
        x = violin_x(mean_d), y = violin_y(),
        fill = "toself", fillcolor = "rgba(0,0,0,0)",
        line = list(color = col, width = 2),
        name = sname, legendgroup = sname, showlegend = FALSE,
        hovertemplate = "%{y:.4g}<extra></extra>"
      ) %>%
      # Min violin outline (dashed)
      plotly::add_trace(
        type = "scatter", mode = "lines",
        x = violin_x(min_d), y = violin_y(),
        fill = "toself", fillcolor = "rgba(0,0,0,0)",
        line = list(color = col, width = 1, dash = "dot"),
        name = sname, legendgroup = sname, showlegend = FALSE,
        hovertemplate = "<extra></extra>"
      )
  }
  
  p %>% plotly::layout(
    xaxis = list(
      title    = x_label,
      tickvals = seq_len(n),
      ticktext = names(lst),
      tickangle = -45,
      range    = c(0.5, n + 0.5)
    ),
    yaxis     = list(title = y_label),
    hovermode = "closest",
    showlegend= FALSE
  )
}

plot_rt_histogram_range <- function(df, sample_name, quantity_cols, color = "steelblue") {
  
  meas_cols <- intersect(quantity_cols, colnames(df))
  
  all_rt <- df$RT[is.finite(df$RT)]
  if (length(all_rt) == 0) return(NULL)
  
  vals   <- df[["RT"]][is.finite(df[["RT"]])]
  breaks <- pretty(vals, n = 50)
  
  if (length(meas_cols) == 0) {
    counts <- list(data.frame(
      Mid   = hist(all_rt, breaks = breaks, plot = FALSE)$mids,
      Count = hist(all_rt, breaks = breaks, plot = FALSE)$counts
    ))
  } else {
    counts <- Filter(Negate(is.null), lapply(meas_cols, function(m) {
      col      <- df[[m]]
      detected <- if (any(is.na(col))) !is.na(col) else col > 0
      vals     <- df$RT[detected & is.finite(df$RT)]
      if (length(vals) == 0) return(NULL)
      h <- hist(vals, breaks = breaks, plot = FALSE)
      data.frame(Mid = h$mids, Count = h$counts)
    }))
  }
  
  if (length(counts) == 0) return(NULL)
  
  stats <- dplyr::bind_rows(counts) %>%
    dplyr::group_by(Mid) %>%
    dplyr::summarise(Mean = mean(Count), Min = min(Count), Max = max(Count), .groups = "drop")
  
  ggplot(stats, aes(x = Mid)) +
    geom_col(aes(y = Mean), fill = color, alpha = 0.6, width = diff(breaks)[1] * 0.9) +
    geom_ribbon(aes(ymin = Min, ymax = Max), fill = color, alpha = 0.3) +
    geom_line(aes(y = Max), color = color, linewidth = 0.4) +
    geom_line(aes(y = Min), color = color, linewidth = 0.4) +
    labs(x = "RT", y = "Count", title = sample_name) +
    theme_minimal()
}

#This is for subplotting dynamic range to a single plot so it only counts as a single context.
dynamic_range_subplot <- function(lst, data_col = "MAX_QUANTITY", ncol = 3) {
  samples <- names(lst)[vapply(lst, function(df)
    nrow(df) >= 10 && data_col %in% colnames(df), logical(1))]
  
  if (length(samples) == 0)
    return(plotly::plot_ly() %>% plotly::layout(
      title = "Not enough peptides to plot dynamic range",
      xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
  
  figs <- lapply(samples, function(sname) {
    df <- lst[[sname]]
    df[[data_col]][df[[data_col]] == 0] <- NA
    df <- df[!is.na(df[[data_col]]), , drop = FALSE]
    
    if (nrow(df) == 0) {
      rank_v <- numeric(0); y_v <- numeric(0); hover <- character(0)
    } else {
      prot_vec <- if ("PROTEIN" %in% names(df)) df$PROTEIN else rep("", nrow(df))
      pep_vec  <- if ("PEPTIDE" %in% names(df)) df$PEPTIDE else rep("", nrow(df))
      rank_v  <- rank(-df[[data_col]], ties.method = "first")
      y_v     <- log2(df[[data_col]])
      gene_nm <- sub(".*?GN=([0-9A-Z/\\-]+).*", "\\1", prot_vec, perl = TRUE)
      disp_nm <- ifelse(gene_nm == prot_vec, prot_vec, paste0(gene_nm, "<br>", prot_vec))
      hover   <- paste0("<b>", disp_nm, "</b><br>", pep_vec, "<br>",
                        "log2(", data_col, "): ", round(y_v, 2), "<br>Rank: ", rank_v)
    }
    
    plotly::plot_ly() %>%
      plotly::add_trace(type = "scattergl", mode = "markers",           # base
                        x = rank_v, y = y_v, text = hover, hoverinfo = "text",
                        marker = list(color = "steelblue", size = 4, opacity = 0.5), showlegend = FALSE) %>%
      plotly::add_trace(type = "scattergl", mode = "markers",           # protein match
                        x = numeric(0), y = numeric(0), text = character(0), hoverinfo = "text",
                        marker = list(color = "red", size = 8, symbol = "square"), showlegend = FALSE) %>%
      plotly::add_trace(type = "scattergl", mode = "markers",           # peptide match
                        x = numeric(0), y = numeric(0), text = character(0), hoverinfo = "text",
                        marker = list(color = "#2CA02C", size = 8, symbol = "square"), showlegend = FALSE) %>%
      plotly::layout(
        xaxis = list(title = ""), yaxis = list(title = ""),
        annotations = list(list(text = sname, x = 0.5, y = 0.98,
                                xref = "paper", yref = "paper", xanchor = "center", yanchor = "top",
                                showarrow = FALSE, font = list(size = 11))))
  })
  
  nrows <- ceiling(length(figs) / ncol)
  plotly::subplot(figs, nrows = nrows, shareX = FALSE, shareY = FALSE,
                  titleX = FALSE, titleY = FALSE,
                  heights = rep(1 / nrows, nrows)) %>%
    plotly::layout(hovermode = "closest", hoverdistance = 30)
}

dynamic_range_plot_combined <- function(df_list, data_col = "MAX_QUANTITY", color = "default", rank_mode = "absolute") {
  # Same inclusion rule the observer relies on for trace-index alignment
  df_list <- df_list[sapply(df_list, function(df)
    nrow(df) >= 10 && data_col %in% colnames(df))]
  
  empty_msg <- function(msg)
    plotly::plot_ly() %>% plotly::layout(
      title = msg,
      xaxis = list(visible = FALSE), yaxis = list(visible = FALSE))
  
  if (length(df_list) == 0)
    return(empty_msg("Not enough peptides to plot combined dynamic range"))
  
  n_samples  <- length(df_list)
  plotly_pal <- c('#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A',
                  '#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52')
  cols <- if (color == "default") {
    setNames(plotly_pal[((seq_len(n_samples) - 1L) %% length(plotly_pal)) + 1L], names(df_list))
  } else setNames(viridis(n_samples, option = color), names(df_list))
  
  p <- plotly::plot_ly()
  any_points <- FALSE
  for (sname in names(df_list)) {
    df <- df_list[[sname]]
    df[[data_col]][df[[data_col]] == 0] <- NA
    df <- df[!is.na(df[[data_col]]), , drop = FALSE]
    
    if (nrow(df) == 0) {
      # keep an (empty) trace so indices stay aligned with the observer
      p <- p %>% plotly::add_trace(
        type = "scattergl", mode = "markers",
        x = numeric(0), y = numeric(0), text = character(0), hoverinfo = "text",
        marker = list(color = cols[[sname]], size = 4, opacity = 0.5),
        name = sname, legendgroup = sname)
      next
    }
    any_points <- TRUE
    
    prot_vec <- if ("PROTEIN" %in% names(df)) df$PROTEIN else rep("", nrow(df))
    pep_vec  <- if ("PEPTIDE" %in% names(df)) df$PEPTIDE else rep("", nrow(df))
    rank_v   <- rank(-df[[data_col]], ties.method = "first")
    n_pts    <- length(rank_v)
    rank_x   <- if (rank_mode == "relative") 100 * rank_v / n_pts else rank_v
    y_v      <- log2(df[[data_col]])
    gene_nm  <- sub(".*?GN=([0-9A-Z/\\-]+).*", "\\1", prot_vec, perl = TRUE)
    disp_nm  <- ifelse(gene_nm == prot_vec, prot_vec, paste0(gene_nm, "<br>", prot_vec))
    hover    <- paste0("<b>", sname, "</b><br>", disp_nm, "<br>", pep_vec, "<br>",
                       "log2(", data_col, "): ", round(y_v, 2), "<br>Rank: ", rank_v)
    
    p <- p %>% plotly::add_trace(
      type = "scattergl", mode = "markers",
      x = rank_v, y = y_v, text = hover, hoverinfo = "text",
      marker = list(color = cols[[sname]], size = 4, opacity = 0.5),
      name = sname, legendgroup = sname)
  }
  if (!any_points) return(empty_msg("No data points available to display"))
  
  # Highlight traces — appended AFTER all sample traces (indices n_samples, n_samples+1)
  p <- p %>%
    plotly::add_trace(type = "scattergl", mode = "markers",
                      x = numeric(0), y = numeric(0), text = character(0), hoverinfo = "text",
                      marker = list(color = "red", size = 8, symbol = "square"), name = "Protein match") %>%
    plotly::add_trace(type = "scattergl", mode = "markers",
                      x = numeric(0), y = numeric(0), text = character(0), hoverinfo = "text",
                      marker = list(color = "#2CA02C", size = 8, symbol = "square"), name = "Peptide match")
  
  x_title <- if (rank_mode == "relative") "Rank (percentile %)" else "Rank"
  p %>% plotly::layout(
    title = list(text = "Combined Dynamic Range for All Samples"),
    xaxis = list(title = x_title),
    yaxis = list(title = paste0("log2(", data_col, ")")),
    hovermode = "closest", hoverdistance = 40,
    legend = list(title = list(text = "Sample")))
}

plot_scatter_mz_k0_plotly <- function(lst, color = "default", ncol = 3) {
  ok <- vapply(lst, function(df)
    !is.null(df) && nrow(df) > 0 && all(c("MZ", "K0") %in% colnames(df)), logical(1))
  samples <- names(lst)[ok]
  if (length(samples) == 0)
    return(plotly::plotly_empty() %>% plotly::layout(title = "No m/z / 1/K0 data"))
  
  charges <- sort(unique(unlist(lapply(lst[samples], function(df) df$CHARGE))))
  charges <- charges[!is.na(charges)]
  if (length(charges) == 0) charges <- NA
  
  plotly_pal <- c('#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A',
                  '#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52')
  ccols <- if (color == "default")
    setNames(plotly_pal[((seq_along(charges) - 1) %% length(plotly_pal)) + 1], as.character(charges))
  else setNames(viridisLite::viridis(length(charges), option = color), as.character(charges))
  
  figs <- lapply(seq_along(samples), function(si) {
    df <- lst[[samples[si]]]
    df <- df[is.finite(df$MZ) & is.finite(df$K0), , drop = FALSE]
    p <- plotly::plot_ly()
    for (ch in charges) {
      d <- df[!is.na(df$CHARGE) & df$CHARGE == ch, , drop = FALSE]
      if (nrow(d) == 0) next
      p <- p %>% plotly::add_trace(
        x = d$MZ, y = d$K0, type = "scattergl", mode = "markers",
        name = paste0("z = ", ch), legendgroup = paste0("z", ch),  # shared toggle
        showlegend = (si == 1),                                    # one legend only
        marker = list(color = ccols[[as.character(ch)]], size = 3, opacity = 0.5),
        hoverinfo = "skip")                                        # no per-dot info
    }
    p %>% plotly::layout(annotations = list(list(
      text = samples[si], x = 0.5, y = 1.03, xref = "paper", yref = "paper",
      xanchor = "center", yanchor = "bottom", showarrow = FALSE, font = list(size = 11))))
  })
  
  fig <- plotly::subplot(figs, nrows = ceiling(length(figs) / ncol),
                         shareX = FALSE, shareY = FALSE,       # aligned axes + fewer ticks to draw
                         titleX = FALSE, titleY = FALSE, margin = 0.03)
  # shared axis labels
  fig$x$layout$annotations <- c(fig$x$layout$annotations, list(
    list(text = "m/z",  x = 0.5,  y = -0.04, xref = "paper", yref = "paper",
         xanchor = "center", yanchor = "top",    showarrow = FALSE),
    list(text = "1/K0", x = -0.04, y = 0.5,  xref = "paper", yref = "paper",
         xanchor = "right",  yanchor = "middle", textangle = -90, showarrow = FALSE)))
  fig %>% plotly::layout(legend = list(title = list(text = "Charge")),
                         margin = list(l = 50, b = 45))
}

plot_completeness <- function(lst, spectra_cols, percent = FALSE, title = "Data Completeness",color = "default") {
  
  make_completeness <- function(df, spectra_cols) {
    
    spectra_cols <- dplyr::intersect(spectra_cols, colnames(df))
    
    df %>%
      transmute(
        Non_NA_Count = rowSums(across(all_of(spectra_cols), ~ {
          if (any(is.na(.x))) !is.na(.x) else .x > 0
        }))
      ) %>%
      arrange(desc(Non_NA_Count)) %>%
      dplyr::mutate(
        Rank           = row_number(),
        Percent_Non_NA = 100 * Non_NA_Count / length(spectra_cols),
        Percent_Rank   = 100 * Rank / max(Rank)
      )
  }
  
  ## color handling
  if (color == "default") {
    color <- NULL
  } else {
    color <- viridis(length(lst), option = color)
  }
  
  df <- map(lst, make_completeness, spectra_cols = spectra_cols) %>%
    bind_rows(.id = "Sample")
  
  if (percent) {
    p <- ggplot(df, aes(Percent_Rank, Percent_Non_NA, color = Sample)) +
      geom_line(linewidth = 1) +
      labs(
        title = title,
        x = "Percent rank (by completeness)",
        y = "Percentage of samples",
        color = "Sample"
      ) +
      theme_minimal() +
      ylim(0,100)
      
  } else {
    p <- ggplot(df, aes(Rank, Non_NA_Count, color = Sample)) +
      geom_line(linewidth = 1) +
      labs(
        title = title,
        x = "Peptide rank (by completeness)",
        y = "Number of samples",
        color = "Sample"
      ) +
      theme_minimal()
  }
  
  if (!is.null(color)) {
    p <- p + scale_color_manual(values = color)
  }
  
  return(p)
}

plot_upset <- function(data_list, min_size = 2, min_size_is_percent = TRUE, title = "Upset Plot", stripped = TRUE, min_degree = 1, n_intersections = 40) {
  
  # Extract peptides per dataset
  l <- purrr::imap(data_list, function(df, name) {
    peptides <- df[["PEPTIDE"]]
    
    # Strip PTMs
    if (stripped) {
      peptides <- df[["STRIPPED"]]
    } else {
      peptides <- df[["PEPTIDE"]]
    }
    
    unique(peptides)
  })
  
  # Combined set of all peptides
  all_items <- unique(unlist(l))
  
  # Binary membership matrix
  df <- lapply(l, function(x) all_items %in% x) |>
    as.data.frame()
  df$Item <- all_items
  df_upset <- dplyr::select(df, -Item)
  
  # Convert min_size percent to actual number
  if (min_size_is_percent) {
    max_size <- nrow(df)
    min_size <- ceiling((min_size / 100) * max_size)
  }
  
  # Upset plot
  p <- suppressWarnings(ComplexUpset::upset(
    df_upset,
    intersect             = names(df_upset),
    name                  = "Sample",
    min_size              = min_size,
    min_degree            = min_degree,
    n_intersections       = n_intersections,
    width_ratio           = 0.1,
    sort_intersections_by = "cardinality"
  ) + ggtitle(title))
  
  return(p)
}

compute_upset_intersections <- function(data_list, stripped = TRUE,
                                        min_size = 2, min_degree = 1,
                                        n_intersections = 40) {
  col  <- if (stripped) "STRIPPED" else "PEPTIDE"
  sets <- lapply(data_list, function(df) {
    if (is.null(df) || !col %in% names(df)) return(NULL)
    unique(df[[col]])
  })
  sets <- Filter(function(s) !is.null(s) && length(s) > 0, sets)
  if (length(sets) < 2) return(list())
  
  all_items  <- unique(unlist(sets))
  if (length(all_items) == 0) return(list())
  sample_nms <- names(sets)
  
  membership <- matrix(FALSE, nrow = length(all_items), ncol = length(sample_nms),
                       dimnames = list(all_items, sample_nms))
  for (s in sample_nms) membership[all_items %in% sets[[s]], s] <- TRUE
  
  profiles       <- apply(membership, 1, function(row) paste(sample_nms[row], collapse = " & "))
  profile_groups <- split(all_items, profiles)
  
  # Assign names BEFORE Filter so they survive filtering
  candidates        <- lapply(names(profile_groups), function(profile) {
    peps   <- profile_groups[[profile]]
    degree <- length(strsplit(profile, " & ")[[1]])
    if (length(peps) >= min_size && degree >= min_degree) peps else NULL
  })
  names(candidates) <- names(profile_groups)
  results           <- Filter(Negate(is.null), candidates)
  
  results <- results[order(sapply(results, length), decreasing = TRUE)]
  if (length(results) > n_intersections) results <- results[seq_len(n_intersections)]
  results
}

plot_shared_peptide_interactive <- function(lst, color = "default",
                                            mode = c("count", "percent"),
                                            percent_type = c("min", "union"),
                                            cluster_mode = "both") {
  mode <- match.arg(mode); percent_type <- match.arg(percent_type)
  peptides <- lapply(lst, function(df) unique(as.character(df$STRIPPED)))
  samples <- names(peptides); n <- length(samples)
  m <- matrix(0, n, n, dimnames = list(samples, samples))
  for (i in samples) for (j in samples) {
    shared <- length(dplyr::intersect(peptides[[i]], peptides[[j]]))
    if (mode == "count") {
      m[i, j] <- shared
    } else {
      denom <- switch(percent_type,
                      min   = min(length(peptides[[i]]), length(peptides[[j]])),
                      union = length(union(peptides[[i]], peptides[[j]])))
      m[i, j] <- if (denom > 0) 100 * shared / denom else NA
    }
  }
  lab  <- if (mode == "count") "# shared peptides" else "% shared peptides"
  lims <- if (mode == "count") c(0, max(m, na.rm = TRUE)) else c(0, 100)
  if (diff(lims) == 0) lims <- c(0, 1)                      # guard degenerate scale
  plot_correlation_heatmap_interactive(m, color = color, cluster_mode = cluster_mode,
                                       label = lab, limits = lims)
}

plot_pairwise_quant_correlation_interactive <- function(lst, method = "pearson",
                                                        min_shared = 2, color = "default",
                                                        cluster_mode = "both") {
  peptides <- lapply(lst, function(df)
    df %>% dplyr::select(STRIPPED, MAX_QUANTITY) %>%
      dplyr::group_by(STRIPPED) %>%
      dplyr::summarise(MAX_QUANTITY = max(MAX_QUANTITY, na.rm = TRUE), .groups = "drop"))
  samples <- names(peptides); n <- length(samples)
  m <- matrix(NA_real_, n, n, dimnames = list(samples, samples))
  for (i in samples) for (j in samples) {
    merged <- merge(peptides[[i]], peptides[[j]], by = "STRIPPED", suffixes = c("_i", "_j"))
    if (nrow(merged) >= min_shared)
      m[i, j] <- suppressWarnings(cor(merged$MAX_QUANTITY_i, merged$MAX_QUANTITY_j,
                                      method = method, use = "complete.obs"))
  }
  plot_correlation_heatmap_interactive(m, color = color, cluster_mode = cluster_mode,
                                       label = paste(method, "correlation"),
                                       limits = NULL)        # correlation → auto range (diag = 1)
}
bin_heatmap_columns <- function(mat, max_cols = 1000) {
  if (ncol(mat) <= max_cols) return(mat)
  n <- ncol(mat); k <- min(max_cols, n)
  mat_imp <- mat; mat_imp[is.na(mat_imp)] <- 0
  bin_map <- list()
  presence_key <- apply(mat_imp > 0, 2, function(x) paste(as.integer(x), collapse = ""))
  patterns <- unique(presence_key)
  pattern_counts <- setNames(sapply(patterns, function(p) sum(presence_key == p)), patterns)
  pattern_bins   <- pmax(1L, floor(k * pattern_counts / n))
  pattern_bins   <- pmin(pattern_bins, pattern_counts)
  names(pattern_bins) <- patterns
  remaining <- k - sum(pattern_bins)
  while (remaining > 0) {
    can_expand <- patterns[pattern_counts[patterns] > pattern_bins[patterns]]
    if (length(can_expand) == 0) break
    n_add  <- min(remaining, length(can_expand))
    add_to <- can_expand[order(pattern_counts[can_expand] - pattern_bins[can_expand],
                               decreasing = TRUE)][seq_len(n_add)]
    pattern_bins[add_to] <- pattern_bins[add_to] + 1L
    remaining <- remaining - n_add
  }
  all_bins <- unlist(lapply(patterns, function(pat) {
    idx    <- which(presence_key == pat)
    n_pat  <- length(idx)
    k_pat  <- min(as.integer(pattern_bins[[pat]]), n_pat)
    if (is.na(k_pat) || k_pat < 1L) k_pat <- 1L
    pep_names <- colnames(mat)[idx]
    mat_pat   <- mat[, idx, drop = FALSE]
    if (n_pat == 1 || k_pat == 1) {
      result   <- rowMeans(mat_pat, na.rm = TRUE)
      result[is.nan(result)] <- NA
      bin_name <- paste0("bin_", pat, "_1")
      bin_map[[bin_name]] <<- pep_names
      return(list(matrix(result, ncol = 1, dimnames = list(rownames(mat), bin_name))))
    }
    pc1 <- tryCatch(
      prcomp(t(mat_imp[, idx, drop = FALSE]), center = TRUE, scale. = FALSE)$x[, 1],
      error = function(e) seq_len(n_pat)
    )
    ord       <- order(pc1)
    mat_pat   <- mat_pat[, ord, drop = FALSE]
    pep_names <- pep_names[ord]
    bin_ids   <- ceiling(seq_len(n_pat) * k_pat / n_pat)
    lapply(seq_len(k_pat), function(i) {
      cols     <- which(bin_ids == i)
      bin_name <- paste0("bin_", pat, "_", i)
      bin_map[[bin_name]] <<- pep_names[cols]
      if (length(cols) == 1) return(mat_pat[, cols, drop = FALSE])
      result <- rowMeans(mat_pat[, cols, drop = FALSE], na.rm = TRUE)
      result[is.nan(result)] <- NA
      matrix(result, ncol = 1, dimnames = list(rownames(mat), bin_name))
    })
  }), recursive = FALSE)
  result <- do.call(cbind, all_bins)
  attr(result, "bin_map") <- bin_map
  result
}

plot_heatmap_interactive <- function(pep_mat, color = "default") {
  mat <- log10(pep_mat + 1)
  mat <- t(mat)
  
  bin_map  <- NULL
  max_cols <- 1000L
  if (ncol(mat) > max_cols) {
    mat     <- bin_heatmap_columns(mat, max_cols = max_cols)
    bin_map <- attr(mat, "bin_map")
  }
  
  # Row clustering
  clust_na0 <- function(x) { x2 <- x; x2[!is.finite(x2)] <- 0; dist(x2) }
  row_ord <- tryCatch({
    d  <- clust_na0(mat)
    hc <- hclust(d, method = "complete")
    hc$order
  }, error = function(e) seq_len(nrow(mat)))
  mat <- mat[row_ord, , drop = FALSE]
  
  # Build hover text matrix
  hover_mat <- matrix("", nrow = nrow(mat), ncol = ncol(mat),
                      dimnames = dimnames(mat))
  for (j in seq_len(ncol(mat))) {
    col_name <- colnames(mat)[j]
    if (!is.null(bin_map) && col_name %in% names(bin_map)) {
      peps  <- bin_map[[col_name]]
      n     <- length(peps)
      shown <- paste(head(peps, 30), collapse = "<br>")
      extra <- if (n > 30) paste0("<br>... +", n - 30, " more") else ""
      label <- paste0("<b>", n, " peptides in bin</b><br>", shown, extra)
    } else {
      label <- paste0("<b>", col_name, "</b>")
    }
    hover_mat[, j] <- label
  }
  
  # Color scale
  if (color == "default") {
    colorscale <- list(list(0, "lightyellow"), list(1, "red"))
  } else {
    cols <- viridis(10, option = color)
    colorscale <- lapply(seq_along(cols) - 1,
                         function(i) list(i / (length(cols) - 1), cols[i + 1]))
  }
  
  plotly::plot_ly(
    x         = colnames(mat),
    y         = rownames(mat),
    z         = mat,
    text      = hover_mat,
    type      = "heatmap",
    colorscale = colorscale,
    hovertemplate = "%{text}<extra></extra>",
    showscale = TRUE
  ) %>%
    plotly::layout(
      xaxis = list(title = "", showticklabels = ncol(mat) <= 50,
                   tickfont = list(size = 9)),
      yaxis = list(title = "", tickfont = list(size = 8),
                   autorange = "reversed"),
      margin = list(l = 120, b = if (ncol(mat) <= 50) 120 else 40)
    )
}

plot_correlation_heatmap_interactive <- function(cor_mat, color = "default",
                                                 cluster_mode = "both", groups = NULL,
                                                 label = "Pearson", tick_limit = 60, limits = NULL) {
  n <- ncol(cor_mat)
  if (is.null(rownames(cor_mat))) rownames(cor_mat) <- colnames(cor_mat)
  
  anchors <- if (color == "default") c("lightyellow", "orange", "red")
  else viridisLite::viridis(3, option = color)
  pal_vec   <- colorRampPalette(anchors)(256)
  show_tick <- n <= tick_limit
  side <- if (!is.null(groups)) data.frame(Sample = as.factor(groups)) else NULL
  
  if (cluster_mode == "sample" && !is.null(groups)) {
    # order by sample; side-bar shows the blocks (Colv/Rowv FALSE → order preserved)
    ord     <- order(groups)
    cor_mat <- cor_mat[ord, ord, drop = FALSE]
    if (!is.null(side)) side <- side[ord, , drop = FALSE]
    dend <- "none"; Rowv <- FALSE; Colv <- FALSE
    
  } else {
    dend <- switch(cluster_mode, both = "both", rows = "row",
                   columns = "column", none = "none", "both")
    Rowv <- dend %in% c("both", "row"); Colv <- dend %in% c("both", "column")
    # side stays in original order — heatmaply reorders it with the dendrogram
  }
  
  heatmaply::heatmaply(
    cor_mat,
    dendrogram      = dend, Rowv = Rowv, Colv = Colv,
    seriate         = "none",
    hclustfun       = fastcluster::hclust,
    colors          = pal_vec, 
    limits          = if (is.null(limits)) range(cor_mat, na.rm = TRUE) else limits,
    col_side_colors = side,
    showticklabels  = c(show_tick, show_tick),
    key.title       = label)
}

pca_scatter_plotly <- function(pca, labels, groups = NULL, color = "default",
                               show_labels = TRUE, dim = "2d") {
  imp <- summary(pca)$importance[2, ] * 100
  npc <- ncol(pca$x)
  use3d <- (dim == "3d") && npc >= 3
  
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                  PC3 = if (npc >= 3) pca$x[, 3] else 0,
                  Label = labels,
                  Group = if (is.null(groups)) "All" else as.character(groups),
                  stringsAsFactors = FALSE)
  glev <- unique(d$Group); n_grp <- length(glev)
  plotly_pal <- c('#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A',
                  '#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52')
  cols <- if (color == "default")
    setNames(plotly_pal[((seq_len(n_grp) - 1) %% length(plotly_pal)) + 1], glev)
  else setNames(viridisLite::viridis(n_grp, option = color), glev)
  
  mode <- if (show_labels) "markers+text" else "markers"
  p <- plotly::plot_ly()
  for (g in glev) {
    dg <- d[d$Group == g, , drop = FALSE]
    common <- list(data = dg, name = g, mode = mode, text = ~Label,
                   textposition = "top center", textfont = list(size = 9),
                   customdata = ~Group,
                   hovertemplate = "<b>%{text}</b><br>Group: %{customdata}<extra></extra>")
    if (use3d) {
      p <- do.call(plotly::add_trace, c(common, list(type = "scatter3d",
                                                     x = ~PC1, y = ~PC2, z = ~PC3, marker = list(size = 4, color = cols[[g]])), list(p = p)))
    } else {
      p <- do.call(plotly::add_trace, c(common, list(type = "scatter",
                                                     x = ~PC1, y = ~PC2, marker = list(size = 10, color = cols[[g]])), list(p = p)))
    }
  }
  if (use3d) {
    p %>% plotly::layout(scene = list(
      xaxis = list(title = sprintf("PC1 (%.1f%%)", imp[1])),
      yaxis = list(title = sprintf("PC2 (%.1f%%)", imp[2])),
      zaxis = list(title = sprintf("PC3 (%.1f%%)", imp[3]))),
      legend = list(title = list(text = "Group")))
  } else {
    p %>% plotly::layout(
      xaxis = list(title = sprintf("PC1 (%.1f%%)", imp[1]), zeroline = FALSE),
      yaxis = list(title = sprintf("PC2 (%.1f%%)", imp[2]), zeroline = FALSE),
      legend = list(title = list(text = "Group")), hovermode = "closest")
  }
}

plot_pca_variance <- function(pca, max_pc = 10) {
  ve  <- summary(pca)$importance[2, ] * 100
  k   <- min(length(ve), max_pc)
  ve  <- ve[seq_len(k)]; cum <- cumsum(ve); pcs <- paste0("PC", seq_len(k))
  plotly::plot_ly() %>%
    plotly::add_bars(x = pcs, y = ve, name = "Individual",
                     marker = list(color = "#636EFA"),
                     hovertemplate = "%{x}: %{y:.1f}%<extra></extra>") %>%
    plotly::add_trace(x = pcs, y = cum, type = "scatter", mode = "lines+markers",
                      name = "Cumulative", line = list(color = "#EF553B"),
                      hovertemplate = "%{x}: %{y:.1f}% cumulative<extra></extra>") %>%
    plotly::layout(
      xaxis = list(title = "Principal component", categoryorder = "array", categoryarray = pcs),
      yaxis = list(title = "Explained variance (%)", range = c(0, 100)),
      legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.15),
      hovermode = "x unified")
}

# --- best across alleles: one bar per sample (min rank) ---
plot_binders_plotly <- function(df, color = "default", percent = TRUE, alleles = NULL) {
  allele_cols <- grep("^HLA", colnames(df), value = TRUE)
  if (!is.null(alleles) && length(alleles) > 0) allele_cols <- intersect(allele_cols, alleles)
  shiny::validate(shiny::need(length(allele_cols) > 0, "No alleles selected."))
  
  minrank <- suppressWarnings(do.call(pmin, c(df[, allele_cols, drop = FALSE], na.rm = TRUE)))
  d <- data.frame(Set = df$Set, min = minrank) %>%
    dplyr::mutate(class = dplyr::case_when(
      !is.finite(min) ~ "NA",          # all alleles unpredicted -> NA (see note)
      min <= 0.5      ~ "Strong",
      min <= 2        ~ "Weak",
      TRUE            ~ "Non-binder")) %>%
    dplyr::count(Set, class, name = "n") %>%
    dplyr::group_by(Set) %>% dplyr::mutate(percent = 100 * n / sum(n)) %>% dplyr::ungroup()
  
  build_binder_stack(d, catcol = "Set", percent = percent, orientation = "h",)
}

# --- per allele: stacked bars, one subplot panel per sample ---
plot_binders_per_allele_plotly <- function(df, color = "default", percent = TRUE,
                                           alleles = NULL,
                                           facet_by = c("sample", "allele")) {
  facet_by <- match.arg(facet_by)
  allele_cols <- grep("^HLA", colnames(df), value = TRUE)
  if (!is.null(alleles) && length(alleles) > 0) allele_cols <- intersect(allele_cols, alleles)
  shiny::validate(shiny::need(length(allele_cols) > 0, "No alleles selected."))
  
  long <- df %>%
    tidyr::pivot_longer(cols = dplyr::all_of(allele_cols),
                        names_to = "Allele", values_to = "Rank") %>%
    dplyr::mutate(class = dplyr::case_when(
      is.na(Rank) ~ "NA",
      Rank <= 0.5 ~ "Strong",
      Rank <= 2   ~ "Weak",
      TRUE        ~ "Non-binder")) %>%
    dplyr::count(Set, Allele, class, name = "n") %>%
    tidyr::complete(Set, Allele, class, fill = list(n = 0)) %>%
    dplyr::group_by(Set, Allele) %>%
    dplyr::mutate(percent = if (sum(n) > 0) 100 * n / sum(n) else 0) %>%
    dplyr::ungroup()
  
  # facet variable vs. y-axis category
  if (facet_by == "sample") { facet_col <- "Set";    cat_col <- "Allele" }
  else                      { facet_col <- "Allele"; cat_col <- "Set"    }
  
  facets <- sort(unique(long[[facet_col]]))
  figs <- lapply(seq_along(facets), function(i)
    build_binder_stack(long[long[[facet_col]] == facets[i], , drop = FALSE],
                       catcol = cat_col, percent = percent, orientation = "h",
                       show_legend = (i == 1), title = facets[i]))
  
  plotly::subplot(figs, nrows = length(figs), shareX = TRUE,
                  titleY = FALSE, margin = 0.03)
}

## --- shared: build one stacked-bar figure from a class-count data.frame ---
build_binder_stack <- function(d, catcol, percent = TRUE, show_legend = TRUE,
                               title = NULL, orientation = "v") {
  class_levels <- c("Strong", "Weak", "Non-binder", "NA")
  class_cols   <- c("Strong" = "#31a354", "Weak" = "#fee08b",
                    "Non-binder" = "#969696", "NA" = "#e0e0e0")
  vcol   <- if (percent) "percent" else "n"
  horiz  <- orientation == "h"
  valfmt <- if (percent) ".1f}%" else ".0f}"
  vtitle <- if (percent) "% of peptides" else "# peptides"
  
  p <- plotly::plot_ly()
  for (cl in class_levels) {
    dc <- d[d$class == cl, , drop = FALSE]
    if (nrow(dc) == 0) next
    if (horiz) {
      p <- p %>% plotly::add_bars(
        y = dc[[catcol]], x = dc[[vcol]], orientation = "h", name = cl,
        marker = list(color = unname(class_cols[cl])),
        legendgroup = cl, showlegend = show_legend,
        hovertemplate = paste0("%{y}<br>", cl, ": %{x:", valfmt, "<extra></extra>"))
    } else {
      p <- p %>% plotly::add_bars(
        x = dc[[catcol]], y = dc[[vcol]], name = cl,
        marker = list(color = unname(class_cols[cl])),
        legendgroup = cl, showlegend = show_legend,
        hovertemplate = paste0("%{x}<br>", cl, ": %{y:", valfmt, "<extra></extra>"))
    }
  }
  
  ax_cat <- list(title = "", tickangle = if (horiz) 0 else -45)
  ax_val <- list(title = vtitle)
  p %>% plotly::layout(
    barmode = "stack",
    xaxis = if (horiz) ax_val else ax_cat,
    yaxis = if (horiz) c(ax_cat, list(autorange = "reversed")) else ax_val,
    legend = list(title = list(text = "Binding class")),
    annotations = if (is.null(title)) NULL else list(list(
      text = title, x = 0.5, y = 1.03, xref = "paper", yref = "paper",
      xanchor = "center", yanchor = "bottom", showarrow = FALSE,
      font = list(size = 12))))
}

#-----Table function------

Colum_mapping_table <- function(df, max_chars = 20) {
  #' This function is to generate the table for column mapping. i.e. which columns are used for what info
  #'
  #' @param df data.frame to display
  #' @param max_chars integer, maximum number of characters to show before truncating. Avoid wide columns
  #' @return data.frame with HTML <span> for truncated display
  #' 
  
  # Convert all columns to character, handling lists
  df_display <- lapply(df, function(col) {
    if (is.list(col)) {
      sapply(col, function(x) {
        paste(unlist(x), collapse = ",")
      })
    } else {
      as.character(col)
    }
  })
  df_display <- as.data.frame(df_display, stringsAsFactors = FALSE)
  
  # Truncate and add hover function
  df_display[] <- lapply(df_display, function(col) {
    sapply(col, function(x) {
      short <- if (nchar(x, type = "width") > max_chars) {
        paste0(substr(x, 1, max_chars), "…")
      } else {
        x
      }
      full <- gsub('"', "&quot;", x)  # escape quotes for HTML title
      sprintf('<span title="%s">%s</span>', full, short)
    }, USE.NAMES = FALSE)
  })
  
  df_display
}

render_summary_pre <- function(lst, height = "800px", font_size = "10px") {
  summary_text <- lapply(lst, summary)
  summary_lines <- lapply(summary_text, function(x) paste(capture.output(x), collapse = "\n"))
  full_text <- paste(summary_lines, collapse = "\n\n")
  
  tags$pre(full_text, style = paste0(
    "height:", height, ";",
    "overflow-y:auto;",
    "overflow-x:auto;",
    "font-size:", font_size, ";"
  ))
}

summarize_binders <- function(df) {
  
  # Select only HLA columns
  hla_cols <- grep("^HLA", names(df), value = TRUE)
  
  df %>%
    dplyr::select(all_of(hla_cols)) %>%
    pivot_longer(cols = everything(), names_to = "Allele", values_to = "Affinity") %>%
    dplyr::mutate(
      Binder = case_when(
        is.na(Affinity) ~ "NA",
        Affinity < 0.5  ~ "Strong",
        Affinity >= 0.5 & Affinity <= 1 ~ "Weak",
        Affinity > 1 ~ "Non"
      )
    ) %>%
    dplyr::group_by(Allele, Binder) %>%
    dplyr::summarise(Count = n(), .groups = "drop") %>%
    pivot_wider(names_from = Binder, values_from = Count, values_fill = 0)
}

#------Statistical caluclations-------
#calcualte group comparison statistics
compute_group_comp_stats <- function(lst, groups, allowed_peptides_g1, allowed_peptides_g2,
                                     g1, g2, pep_col, quantity_cols, col_map = NULL, use_measurements = FALSE) {
  
  grp_g1 <- groups[[g1]]
  grp_g2 <- groups[[g2]]

 
  if (!use_measurements) {
    n_g1 <- sum(sapply(grp_g1, function(s) length(intersect(quantity_cols, colnames(lst[[s]])))))
    n_g2 <- sum(sapply(grp_g2, function(s) length(intersect(quantity_cols, colnames(lst[[s]])))))
  } else {
    n_g1 <- length(grp_g1)
    n_g2 <- length(grp_g2)
  }
  
  use_limma <- FALSE
  #use_limma <- n_g1 == 1 || n_g2 == 1
  
  # ── Build long-format data ──────────────────────────────────────────────────
  build_long <- function(grp_items, grp_name, allowed) {
    #this is the standard sample mode: if no annotation table measurement data is given.
    if (!use_measurements) {
      dplyr::bind_rows(lapply(grp_items, function(s) {
        df <- lst[[s]]
        cols <- intersect(c(pep_col, quantity_cols, "PROTEIN"), colnames(df))
        df[df[[pep_col]] %in% allowed, cols, drop = FALSE]
      })) %>%
        tidyr::pivot_longer(cols = dplyr::any_of(quantity_cols),
                            names_to = "Sample", values_to = "Quantity") %>%
        dplyr::mutate(Quantity = replace(Quantity, Quantity == 0, NA),
                      Group = grp_name)
    } else if (is.list(grp_items) && !is.null(grp_items$measurement)) {
      #this is if the annotation table is given and also has measurement column.
      # grp_items = list(name = c(...), measurement = c(...)), parallel vectors
      measurements <- grp_items$measurement
      names_vec    <- grp_items$name
      dplyr::bind_rows(lapply(seq_along(measurements), function(i) {
        nm   <- names_vec[i]
        meas <- measurements[i]
        entry <- col_map[[meas]]
        if (is.null(entry)) {
          warning(paste0("[build_long] No col_map entry for measurement '", meas, "'"))
          return(NULL)
        }
        df <- lst[[nm]]
        if (is.null(df)) {
          warning(paste0("[build_long] No data in lst for sample '", nm, "'"))
          return(NULL)
        }
        actual_col <- colnames(df)[tolower(colnames(df)) == tolower(entry$col)][1]
        if (is.na(actual_col)) {
          warning(paste0("[build_long] Column '", entry$col, "' not found (case-insensitive) in '",
                         nm, "'. Available cols: ", paste(head(colnames(df), 30), collapse = ", ")))
          return(NULL)
        }
        keep_cols <- intersect(c(pep_col, "PROTEIN", actual_col), colnames(df))
        df <- df[df[[pep_col]] %in% allowed, keep_cols, drop = FALSE]
        names(df)[names(df) == actual_col] <- "Quantity"
        df$Quantity[df$Quantity == 0] <- NA
        if (!"PROTEIN" %in% colnames(df)) df$PROTEIN <- NA_character_
        df$Sample <- meas
        df$Group  <- grp_name
        df
      }))
    } else {
      #fallback method, in unforseen circumstances.
      print("[compute_group_comp_stats - build_long] Fallback branch. investigate why.")
      dplyr::bind_rows(lapply(grp_items, function(s) {
        df <- lst[[s]]
        cols <- intersect(c(pep_col, quantity_cols, "PROTEIN"), colnames(df))
        df[df[[pep_col]] %in% allowed, cols, drop = FALSE]
      })) %>%
        tidyr::pivot_longer(cols = dplyr::any_of(quantity_cols),
                            names_to = "Sample", values_to = "Quantity") %>%
        dplyr::mutate(Quantity = replace(Quantity, Quantity == 0, NA),
                      Group = grp_name)
    }
    
  }
  
  df_long <- dplyr::bind_rows(
    build_long(grp_g1, g1, allowed_peptides_g1),
    build_long(grp_g2, g2, allowed_peptides_g2)
  )
  
  if (nrow(df_long) == 0) return(NULL)
  if (!"Quantity" %in% colnames(df_long)) {
    warning(paste0("[compute_group_comp_stats] 'Quantity' column missing from df_long. ",
                   "use_measurements=", use_measurements,
                   "; colnames=[", paste(colnames(df_long), collapse = ", "), "]"))
    return(NULL)
  }
  
  # Bail out if either group has no quantifiable data (all NA/zero after imputation)
  if (!any(is.finite(df_long$Quantity[df_long$Group == g1]))) {
    message(paste0("[compute_group_comp_stats] Group '", g1, "' has no finite quantities — skipping comparison."))
    return(NULL)
  }
  if (!any(is.finite(df_long$Quantity[df_long$Group == g2]))) {
    message(paste0("[compute_group_comp_stats] Group '", g2, "' has no finite quantities — skipping comparison."))
    return(NULL)
  }
  
  # ── Means and FC ─────────────────────────────────────────────────────────────
  summary_df <- df_long %>%
    dplyr::group_by(.data[[pep_col]]) %>%
    dplyr::summarise(
      PROTEIN = dplyr::first(.data[["PROTEIN"]]),
      Mean_G1 = mean(.data[["Quantity"]][.data[["Group"]] == g1], na.rm = TRUE),
      Mean_G2 = mean(.data[["Quantity"]][.data[["Group"]] == g2], na.rm = TRUE),
      log2FC  = log2(Mean_G2 + 1) - log2(Mean_G1 + 1),
      A       = 0.5 * (log2(Mean_G1 + 1) + log2(Mean_G2 + 1)),
      .groups = "drop"
    )
  
  
  # ── P-values ──────────────────────────────────────────────────────────────────
  safe_limma_pvals <- function(df_g1, df_g2, qty_g1, qty_g2, pep_col) {
    if (!requireNamespace("limma", quietly = TRUE))
      stop("Package 'limma' is required for n=1 testing. Install via BiocManager::install('limma')")
    
    mat_g1 <- as.matrix(df_g1 %>% dplyr::select(any_of(qty_g1)) %>% log2())
    mat_g2 <- as.matrix(df_g2 %>% dplyr::select(any_of(qty_g2)) %>% log2())
    
    # align peptide rows
    peps   <- intersect(df_g1[[pep_col]], df_g2[[pep_col]])
    mat_g1 <- mat_g1[df_g1[[pep_col]] %in% peps, , drop = FALSE]
    mat_g2 <- mat_g2[df_g2[[pep_col]] %in% peps, , drop = FALSE]
    
    combined <- cbind(mat_g1, mat_g2)
    group    <- factor(c(rep("g1", ncol(mat_g1)), rep("g2", ncol(mat_g2))))
    design   <- model.matrix(~group)
    
    fit  <- limma::lmFit(combined, design)
    fit  <- limma::eBayes(fit)
    limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")$P.Value
  }
  
  if (use_limma) {
    df_g1 <- df_long[df_long$Group == g1, ]
    df_g2 <- df_long[df_long$Group == g2, ]
    pval_df <- safe_limma_pvals(df_g1, df_g2, "Quantity", "Quantity", pep_col)
    test_method <- "limma (n=1 fallback)"
  } else {
    pval_df <- df_long %>%
      dplyr::group_by(.data[[pep_col]]) %>%
      dplyr::summarise(
        pval = {
          x <- Quantity[Group == g1 & is.finite(Quantity)]
          y <- Quantity[Group == g2 & is.finite(Quantity)]
          if (length(x) < 2 || length(y) < 2) NA_real_
          else tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_)
        },
        .groups = "drop"
      )
    test_method <- "t-test"
  }
  
  # ── Combine and adjust ────────────────────────────────────────────────────────
  dplyr::left_join(summary_df, pval_df, by = pep_col) %>%
    dplyr::mutate(
      adj_pval_BH       = p.adjust(pval, method = "BH"),
      adj_pval_Bonf     = p.adjust(pval, method = "bonferroni"),
      negLog10P         = -log10(pval),
      negLog10AdjP_BH   = -log10(adj_pval_BH),
      negLog10AdjP_Bonf = -log10(adj_pval_Bonf),
      test_method       = test_method
    )
}


#binding prediction summary
compute_binder_summary <- function(df, alleles = NULL) {
  allele_cols <- grep("^HLA", colnames(df), value = TRUE)
  if (!is.null(alleles) && length(alleles) > 0){
    allele_cols <- intersect(allele_cols, alleles)
  }
  
  df %>%
    tidyr::pivot_longer(cols = dplyr::all_of(allele_cols), names_to = "Allele", values_to = "Rank") %>%
    dplyr::mutate(
      Class = dplyr::case_when(
        is.na(Rank) ~ "Missing",
        Rank <= 0.5 ~ "Strong",
        Rank <= 2   ~ "Weak",
        TRUE        ~ "Non"
      ),
      Class = factor(Class, levels = c("Strong", "Weak", "Non", "Missing"))
    ) %>%
    dplyr::count(Set, Allele, Class) %>%
    tidyr::complete(Set, Allele, Class, fill = list(n = 0)) %>%
    tidyr::pivot_wider(names_from = Class, values_from = n) %>%
    dplyr::mutate(
      Total      = Strong + Weak + Non + Missing,
      Strong_pct = round(100 * Strong  / Total, 2),
      Weak_pct   = round(100 * Weak    / Total, 2),
      Non_pct    = round(100 * Non     / Total, 2),
      NA_pct     = round(100 * Missing / Total, 2)
    )
}

group_venn_plotting <- function(sets, color = "default") {
  n_groups <- length(sets)
  
  p <- ggVennDiagram::ggVennDiagram(
    sets,
    label_alpha = 0
  ) +
    ggplot2::theme_void()  +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = margin(10,10,10,10)
    )
  
  if (color != "default") {
    p <- p + ggplot2::scale_fill_viridis_c(
      option = color
    )  
  } else {
    p <- p + scale_fill_distiller(palette = "RdBu")
  }
  
  return(p)
}

group_volcano_plot <- function(df,sel,plot_name) {
  # Check if groups have data
  shiny::validate(shiny::need(any(!is.na(df$log2FC)) && any(!is.na(df$negLog10AdjP_BH)), "have data, but no log2FC and p-values."))
  
  if (is.null(sel)) sel <- NA_character_  # <- avoids length 0
  
  # Add 'selected' column
  volc_df <- df %>%
    dplyr::mutate(selected = !is.na(sel) & PEPTIDE == sel)
  
  # Split by significance for plotting
  ns_df        <- volc_df %>% dplyr::filter(Significance == "Not significant")
  sig_df       <- volc_df %>% dplyr::filter(Significance == "Significant", !selected)
  selected_df  <- volc_df %>% dplyr::filter(selected)
  
  #Axis setting
  safe_max_neglogp <- suppressWarnings(max(df$negLog10AdjP_BH, na.rm = TRUE))
  safe_min_fc      <- suppressWarnings(min(df$log2FC,   na.rm = TRUE))
  safe_max_fc      <- suppressWarnings(max(df$log2FC,   na.rm = TRUE))
  
  if (!is.finite(safe_max_neglogp)) safe_max_neglogp <- 1
  if (!is.finite(safe_min_fc))      safe_min_fc <- -1
  if (!is.finite(safe_max_fc))      safe_max_fc <- 1
  
  #Plotting
  plot_ly(source = "group_diff") %>%
    add_trace(
      data = ns_df, x = ~log2FC, y = ~negLog10AdjP_BH,
      type = "scatter", mode = "markers",
      name = "Not significant",
      marker = list(color = "grey40", size = 6),
      hoverinfo = "none",
      inherit = FALSE
    ) %>%
    add_trace(
      data = sig_df, x = ~log2FC, y = ~negLog10AdjP_BH,
      type = "scatter", mode = "markers",
      name = "Significant",
      marker = list(color = "red", size = 8),
      text = ~PEPTIDE,
      key = ~PEPTIDE,
      hoverinfo = "text+x+y"
    ) %>%
    add_trace(
      data = selected_df,
      x = ~log2FC, y = ~negLog10AdjP_BH,
      type = "scatter", mode = "markers",
      name = "Selected peptide",
      marker = list(color = "lightgreen", size = 12),
      text = ~PEPTIDE,
      key = ~PEPTIDE,
      hoverinfo = "text+x+y"
    ) %>%
    layout(
      title = paste0("volcano plot: ",plot_name),
      xaxis = list(title = "log2 Fold Change"),
      yaxis = list(title = "-log10 p-value"),
      shapes = list(
        # fold change vertical lines
        list(type = "line", x0 = -1, x1 = -1, y0 = 0, y1 = safe_max_neglogp,
             line = list(dash = "dash", color = "grey")),
        list(type = "line", x0 =  1, x1 =  1, y0 = 0, y1 = safe_max_neglogp,
             line = list(dash = "dash", color = "grey")),
        # p-value cutoff
        list(type = "line", x0 = safe_min_fc, x1 = safe_max_fc,
             y0 = 1.3, y1 = 1.3,
             line = list(dash = "dash", color = "grey"))
      )
    )
}

#group comparison: P-value Histogram
group_p_histogram <- function(df,plot_name) {
  shiny::validate(shiny::need(any(!is.na(df$negLog10AdjP_BH)), "have data, but no p-values."))
  
  pvals <- 10^(-df$negLog10AdjP_BH)
  pvals <- pvals[is.finite(pvals) & pvals >= 0 & pvals <= 1]
  
  plot_ly(
    x = pvals,
    type = "histogram",
    nbinsx = 50,
    marker = list(color = "grey40"),
    hoverinfo = "x+y"
  ) %>%
    layout(
      title = paste("P-value distribution:", plot_name),
      xaxis = list(title = "Adjusted p-value", range = c(0, 1)),
      yaxis = list(title = "Count", type = "log"),
      shapes = list(
        list(
          type = "line",
          x0 = 0.05, x1 = 0.05,
          y0 = 0, y1 = 1,
          xref = "x",
          yref = "paper",
          line = list(color = "red", dash = "dash")
        )
      ),
      annotations = list(
        list(
          x = 0.05, y = 1,
          xref = "x", yref = "paper",
          text = "0.05",
          showarrow = FALSE,
          xanchor = "left",
          font = list(color = "red")
        )
      )
    )
  
}

#group comparison: MA plot
group_MA_plot <- function(df,sel,plot_name) {
  shiny::validate(shiny::need(any(!is.na(df$negLog10AdjP_BH)), "have data, but no Average values."))

  if (is.null(sel)) sel <- NA_character_  # <- avoids length 0
  
  # Add 'selected' column
  ma_df <- df %>%
    dplyr::mutate(selected = !is.na(sel) & PEPTIDE == sel)
  
  # Split by significance for plotting
  ns_df        <- ma_df %>% dplyr::filter(Significance == "Not significant")
  sig_df       <- ma_df %>% dplyr::filter(Significance == "Significant", !selected)
  selected_df  <- ma_df %>% dplyr::filter(selected)
  
  plot_ly(source = "group_diff") %>%
    add_trace(
      data = ns_df,
      x = ~A, y = ~log2FC,
      type = "scatter", mode = "markers",
      name = "Not significant",
      marker = list(color = "grey40", size = 6),
      hoverinfo = "none",
      inherit = FALSE
    ) %>%
    add_trace(
      data = sig_df,
      x = ~A, y = ~log2FC,
      type = "scatter", mode = "markers",
      name = "Significant",
      marker = list(color = "red", size = 8),
      text = ~PEPTIDE,
      key = ~PEPTIDE,
      hoverinfo = "text+x+y"
    ) %>%
    add_trace(
      data = selected_df,
      x = ~A, y = ~log2FC,
      type = "scatter", mode = "markers",
      name = "Selected peptide",
      marker = list(color = "lightgreen", size = 12),
      text = ~PEPTIDE,
      key = ~PEPTIDE,
      hoverinfo = "text+x+y"
    ) %>%
    layout(
      title = paste("MA plot:", plot_name),
      xaxis = list(title = "Mean abundance (A)"),
      yaxis = list(title = "log2 Fold Change (M)"),
      shapes = list(
        list(type = "line", x0 = min(df$A, na.rm = TRUE),
             x1 = max(df$A, na.rm = TRUE),
             y0 = 0, y1 = 0,
             line = list(dash = "dash", color = "grey"))
      )
    )
}

#group comparison: Ranked Fold Change
group_rank_FC <- function(df,plot_name,sel) {
  shiny::validate(shiny::need(any(!is.na(df$negLog10AdjP_BH)), "have data, but no log2FC values."))
  
  if (is.null(sel)) sel <- NA_character_
  
  rank_df <- df %>%
    dplyr::filter(!is.na(log2FC)) %>%
    arrange(log2FC) %>%          # ascending; use desc(log2FC) if you prefer
    dplyr::mutate(
      rank = row_number(),
      selected = !is.na(sel) & PEPTIDE == sel
    )
  
  selected_df <- rank_df %>% dplyr::filter(selected)
  rest_df     <- rank_df %>% dplyr::filter(!selected)
  
  plot_ly(source = "group_diff") %>%
    
    # All peptides
    add_trace(
      data = rest_df,
      x = ~rank,
      y = ~log2FC,
      type = "scatter",
      mode = "markers",
      name = "Peptides",
      marker = list(color = "grey40", size = 6),
      hoverinfo = "none"
    ) %>%
    
    # Selected peptide
    add_trace(
      data = selected_df,
      x = ~rank,
      y = ~log2FC,
      type = "scatter",
      mode = "markers",
      name = "Selected peptide",
      marker = list(color = "lightgreen", size = 12),
      text = ~PEPTIDE,
      key = ~PEPTIDE,
      hoverinfo = "text+x+y"
    ) %>%
    
    layout(
      title = paste("Ranked fold change:", plot_name),
      xaxis = list(title = "Rank (based on Fold Change)"),
      yaxis = list(title = "log2 Fold Change"),
      shapes = list(
        list(
          type = "line",
          x0 = 0,
          x1 = max(rank_df$rank),
          y0 = 0,
          y1 = 0,
          line = list(dash = "dash", color = "grey")
        )
      )
    )
}

#group comparison: Fold Change table
group_FC_table <- function(df) {
  table_df <- df %>%
    dplyr::filter(!is.na(log2FC)) %>%
    dplyr::select(
      PEPTIDE,
      log2FC,
      adj_pval_BH,
      adj_pval_Bonf,
      Significance
    ) %>%
    arrange(adj_pval_BH)
  
  DT::datatable(
    table_df,
    rownames = FALSE,
    selection = "none",
    options = list(
      pageLength = 10,
      scrollX = TRUE,
      order = list(list(2, "asc"))
    )
  ) %>%
    DT::formatRound(c("log2FC", "adj_pval_BH","adj_pval_Bonf"), digits = 3)
}

##----static versions-----
group_volcano_plot_static <- function(df, plot_name) {
  ggplot(df, aes(x = log2FC, y = negLog10AdjP_BH, colour = Significance)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_colour_manual(values = c(
      "Significant"     = "red",
      "Not significant" = "grey40",
      "Missing"         = "grey80"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey60") +
    geom_hline(yintercept = 1.3,      linetype = "dashed", colour = "grey60") +
    labs(title = paste("Volcano:", plot_name),
         x = "log2 Fold Change", y = "-log10 adjusted p-value") +
    theme_minimal()
}

group_MA_plot_static <- function(df, plot_name) {
  ggplot(df, aes(x = A, y = log2FC, colour = Significance)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_colour_manual(values = c(
      "Significant"     = "red",
      "Not significant" = "grey40",
      "Missing"         = "grey80"
    )) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    labs(title = paste("MA plot:", plot_name),
         x = "Mean abundance (A)", y = "log2 Fold Change (M)") +
    theme_minimal()
}

group_p_histogram_static <- function(df, plot_name) {
  pvals <- 10^(-df$negLog10AdjP_BH)
  pvals <- pvals[is.finite(pvals) & pvals >= 0 & pvals <= 1]
  
  if (length(pvals) == 0) {
    return(ggplot() + annotate("text", x=.5, y=.5, label="No p-values available",
                               size=4, colour="grey50") + theme_void())
  }
  
  ggplot(data.frame(pval = pvals), aes(x = pval)) +
    geom_histogram(bins = 50, fill = "grey40", colour = "white") +
    geom_vline(xintercept = 0.05, linetype = "dashed", colour = "red") +
    scale_y_log10() +
    labs(title = paste("P-value distribution:", plot_name),
         x = "Adjusted p-value", y = "Count (log10)") +
    theme_minimal()
}

group_rank_FC_static <- function(df, plot_name) {
  rank_df <- df %>%
    dplyr::filter(!is.na(log2FC)) %>%
    arrange(log2FC) %>%
    dplyr::mutate(rank = row_number())
  
  ggplot(rank_df, aes(x = rank, y = log2FC, colour = Significance)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_colour_manual(values = c(
      "Significant"     = "red",
      "Not significant" = "grey40",
      "Missing"         = "grey80"
    )) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    labs(title = paste("Ranked fold change:", plot_name),
         x = "Rank", y = "log2 Fold Change") +
    theme_minimal()
}
#------Grids and precomputation for Markdown report---------
#The functions here are specifically for the Markdown report, where all plots in the grid tabs need to be precomputed. Opposite of Shiny reactive.
generate_motif_grid <- function(lst, lengths = 7:20) {
  library(ggplot2)
  library(ggseqlogo)
  library(patchwork)  # for arranging plots
  
  # global legend
  global_legend_plot <- ggseqlogo::ggseqlogo("ACDEFGHIKLMNPQRSTVWY") +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("Amino Acid Colors")
  
  length_plots <- list()
  
  for (L in lengths) {
    plots_per_length <- list()
    
    # Add the legend at the top
    plots_per_length[[1]] <- global_legend_plot
    
    for (sample_name in names(lst)) {
      df <- lst[[sample_name]]
      
      df_L <- df[df$LENGTH == L, ]
      peptides <- unique(df_L[["STRIPPED"]])
      peptides <- peptides[nchar(peptides) == L]
      
      if (length(peptides) < 5) {
        p <- ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = "Not enough peptides",
                   size = 6, color = "red") +
          theme_void()
      } else {
        p <- plot_seqlogo(peptides, title = paste(sample_name, "• Length", L))
      }
      
      plots_per_length[[length(plots_per_length) + 1]] <- p
    }
    
    # Combine horizontally with patchwork
    length_plots[[paste0("Length_", L)]] <- wrap_plots(plots_per_length, ncol = 3)
  }
  
  length_plots
}

#----Running external stuff------

run_netmhcpan <- function(peptides, allele, netmhcpan_path) {
  peptide_file <- tempfile(fileext = ".txt")
  output_file  <- tempfile(fileext = ".txt")
  writeLines(peptides, peptide_file)
  
  peptide_wsl <- trimws(system2("wsl", c("wslpath", "-a", shQuote(peptide_file)), stdout = TRUE))
  out_wsl     <- trimws(system2("wsl", c("wslpath", "-a", shQuote(output_file)),  stdout = TRUE))
  
  cmd <- paste(
    shQuote(netmhcpan_path),
    "-p", shQuote(peptide_wsl),
    "-a", shQuote(allele),
    "-l 8,9,10,11",
    "-xls",
    "-xlsfile", shQuote(out_wsl)
  )
  status <- system2("wsl", c("bash", "--login", "-c", shQuote(cmd)), stdout = NULL)
  
  if (status != 0)
    stop("netMHCpan exited with status ", status,
         ". Check that the path is correct and the executable exists: ", netmhcpan_path)
  if (!file.exists(output_file))
    stop("netMHCpan ran but produced no output file. Check the netMHCpan installation.")
  
  output_file
}

parse_netmhc_output <- function(output_file) {
  res     <- read.table(output_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  header1 <- colnames(res)
  header2 <- as.character(unlist(res[1, ]))
  
  colnames_new <- header1
  colnames_new[2] <- "Peptide"
  current_hla  <- NULL
  for (i in seq_along(colnames_new)) {
    if (grepl("^HLA", header1[i])) current_hla <- header1[i]
    if (!is.null(current_hla) && header2[i] != "") colnames_new[i] <- paste0(current_hla, "_", header2[i])
  }
  colnames(res) <- colnames_new
  res <- res[-1, ]
  
  keep_cols    <- c("Peptide", grep("_Rank$", colnames(res), value = TRUE, ignore.case = TRUE)) #This handles case insensitive _rank.
  res          <- res[, keep_cols, drop = FALSE]
  colnames(res) <- gsub("(_[A-Za-z]+)?_rank$", "", colnames(res), ignore.case = TRUE) #this now also handles _EL_rank as well as just _Rank.
  colnames(res)[-1] <- sub("^([^.]+\\.[^.]+)\\.", "\\1", colnames(res)[-1])
  res[-1] <- lapply(res[-1], as.numeric)
  res
}

run_go_enrichment <- function(df) {
  df_sig <- df %>%
    dplyr::filter(!is.na(log2FC), Significance == "Significant") %>%
    dplyr::select(PEPTIDE, Significance, PROTEIN)
  
  uni_ids <- df_sig$PROTEIN %>%
    strsplit(";") %>%
    unlist() %>%
    trimws() %>%
    sapply(function(prot) {
      parts <- strsplit(prot, "\\|")[[1]]
      if (length(parts) >= 3) parts[2] else parts[1]  # "sp|P04439|HLA_A" → "P04439"
    }) %>%
    unique()
  
  check <- safe_validate(!is.null(uni_ids) && nrow(as.data.frame(uni_ids)) > 0, "No significant IDs")
  if (!is.null(check)) return(check)
  
  gene_map <- tryCatch(
    suppressWarnings(clusterProfiler::bitr(uni_ids, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db)),
    error = function(e) data.frame(UNIPROT = character(), ENTREZID = character())
  )
  
  check <- safe_validate(nrow(gene_map) > 0, "No valid UniProt→Entrez mapping")
  if (!is.null(check)) return(check)
  
  ego <- clusterProfiler::enrichGO(
    gene          = gene_map$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    universe      = NULL,
    readable      = TRUE
  )
  
  check <- safe_validate(!is.null(ego) && nrow(as.data.frame(ego)) > 0, "No significant GO terms")
  if (!is.null(check)) return(check)
  
  suppressWarnings(barplot(ego, showCategory = 10))
  
}

run_string <- function(df) {
  df <- df %>%
    dplyr::filter(!is.na(log2FC)) %>%
    dplyr::filter(Significance == "Significant") %>%
    dplyr::select(PEPTIDE, Significance, PROTEIN)
  
  uni_ids <- df$PROTEIN %>%
    strsplit(";") %>%
    unlist() %>%
    trimws() %>%
    sapply(function(prot) {
      parts <- strsplit(prot, "\\|")[[1]]
      if (length(parts) >= 3) parts[2] else parts[1]
    }) %>%
    unique()
  
  #Check for empty uni_ids
  check <- safe_validate(!is.null(uni_ids) && nrow(as.data.frame(uni_ids)) > 0, "No significant IDs")
  
  url <- paste0(
    "https://string-db.org/api/json/network?",
    "identifiers=", paste(uni_ids, collapse = "%0d"),
    "&species=", 9606 # human
  )
  
  res <- GET(url)
  
  #Check HTTP response
  check <- safe_validate(http_status(res)$category == "Success", paste0("STRING request failed (HTTP ", res$status_code, ")"))
  if (!is.null(check)) return(check)
  
  # Try to parse JSON safely
  data <- tryCatch(
    {
      fromJSON(content(res, "text", encoding = "UTF-8"))
    },
    error = function(e) {
      safe_validate(FALSE, paste0("JSON parse error:\n", e$message))
      NULL
    }
  )
  
  # Check that parsing returned something
  check <- safe_validate(!is.null(data) && length(data) > 0,paste0("No STRING-DB result (HTTP ", res$status_code, ")"))
  if (!is.null(check)) return(check)
  
  required_cols <- c("preferredName_A", "preferredName_B", "score")
  if (!all(required_cols %in% colnames(data))) {
    return(list(.skip = TRUE, message = "STRING-DB result missing required columns"))
  }
  
  g <- graph_from_data_frame(
    data[, required_cols],
    directed = FALSE
  )
  
  ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = score), alpha = 0.8) +
    geom_node_point(size = 5, color = "steelblue") +
    geom_node_text(aes(label = name), repel = TRUE) +
    theme_void()
}


