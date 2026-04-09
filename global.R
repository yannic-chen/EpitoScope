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
    PEPTIDE        = c("Modified.Peptide", "modified.sequence"),
    STRIPPED       = c("Peptide", "peptide.sequence"),           #not present in combined_peptide.tsv
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
    SPECTRA        = c("Spectrum", "Spectral.Count"),            #suffix.
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
  # Split by semicolon
  all_acc <- unlist(strsplit(accessions, ";"))
  # Extract the part before first "|" in each entry
  proteins <- sapply(all_acc, function(x) {
    strsplit(x, "\\|")[[1]][2]  # 2nd element is usually the accession ID
  })
  proteins <- proteins[!is.na(proteins)]  # remove NAs if any
  unique(proteins)
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
  if (length(plots) == 0) return(NULL)
  patchwork::wrap_plots(plots, ncol = ncol)
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
      ncol = length(lst) + 1
    )
  })
  
  patchwork::wrap_plots(length_panels, ncol = 1)
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
      dplyr::distinct(name) %>%
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

find_transform_column <- function(df, candidates, name) {
  for (cand in candidates) {
    # Prefer exact case-sensitive match to avoid ambiguity when a prior rename
    # produced a column that differs only in case (e.g. "Peptide" vs "PEPTIDE")
    exact_match <- which(colnames(df) == cand)
    if (length(exact_match) == 1) {
      original_name       <- colnames(df)[exact_match]
      colnames(df)[exact_match] <- name
      return(list(df = df, matched_column = original_name))
    }

    # Fall back to case-insensitive when no exact match exists
    match <- which(tolower(colnames(df)) == tolower(cand))
    if (length(match) == 1) {
      original_name    <- colnames(df)[match]
      colnames(df)[match] <- name
      return(list(df = df, matched_column = original_name))
    } else if (length(match) > 1) {
      stop("Ambiguous: multiple columns match '", cand, "' for ", name, ": ",
           paste(colnames(df)[match], collapse = ", "))
    }
    # length == 0 → try next candidate
  }
  
  message(paste0("No matching column found for: ", name))
  return(NULL)
}

transform_columns <- function(df, schema, software, targets = c("PEPTIDE", "STRIPPED", "LENGTH", "MASS", "MZ", "SCORE", "CHARGE", "RT", "PPM", "PROTEIN", "PTM")) {
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
  
  for (target_name in names(col_map)) {
    candidates <- col_map[[target_name]]
    
    # Skip if no candidates defined
    if (length(candidates) == 0) next
    
    res <- find_transform_column(df, candidates, target_name)
    if (is.null(res)) next
    df <- res$df
    
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
        names_prefix = "Intensity.",
        values_fn = ~ if (all(is.null(.x))) NA_real_ else max(.x, na.rm = TRUE),
        values_fill = NA_real_
        )
    
    df_top <- df %>%
      dplyr::group_by(Modified.Sequence) %>%
      dplyr::arrange(Q.Value, .by_group = TRUE) %>%
      dplyr::slice_head(n = 1) %>%
      ungroup()
    
    df <- left_join(df_top %>% dplyr::select(-Run), df_wide, by = "Modified.Sequence")
    
  }
  res <- transform_columns(df, column_schema, software)
  
  df <- res$df
  original <- res$log
  
  #These two columns are the minimum required. If there is no PTM, PEPTIDE will simply be the same as STRIPPED.
  if (any(sapply(df, function(data) c("PEPTIDE", "STRIPPED") %in% colnames(data)))) { 
    stop(sprintf("Either PEPTIDE or STRIPPED column missing in one or more dataframes."))
    }
  
  #-------------Derive missing columns if possible----------------
  #CHARGE
  if (!"CHARGE" %in% colnames(df)) {
    df$CHARGE <- 0
    message("Charge column missing → set to 0")
    original <- rbind(original, data.frame(final_name = "CHARGE", original_name = "[no CHARGE column]", stringsAsFactors = FALSE))
  } else {
    df$CHARGE <- as.integer(df$CHARGE)
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
  
  # 2. Hard-coded fallbacks (existing behaviour for known software)
  if (software %in% c("DIANN", "DIANN_parquet")) {
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
  
  #-------------Quantity -----------------
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
  if (length(spec) == 0) {
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
make_range_slider_ui <- function(lst, col, input_id, label, step = NULL) {
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
    group_items    <- groups[[g]]
    # Normalise: strip $expr, keep list(name, measurement) or plain char vector
    if (!is.null(col_map) && is.list(group_items) && !is.null(group_items$measurement)) {
      group_items <- list(name = group_items$name, measurement = group_items$measurement)
    }
    allowed_peptides <- group_peptide_sets[[g]]

    if (!is.null(col_map)) {
      # Measurement mode: grp_items = list(name=c(...), measurement=c(...))
      measurements <- group_items$measurement
      names_vec    <- group_items$name
      df_group <- dplyr::bind_rows(lapply(seq_along(measurements), function(i) {
        nm   <- names_vec[i]
        meas <- measurements[i]
        entry <- col_map[[meas]]
        if (is.null(entry)) return(NULL)
        df <- lst[[nm]]
        if (is.null(df) || !"PEPTIDE" %in% colnames(df)) return(NULL)
        actual_col <- colnames(df)[tolower(colnames(df)) == tolower(entry$col)][1]
        if (is.na(actual_col)) return(NULL)
        df_sub <- df[df$PEPTIDE %in% allowed_peptides, c("PEPTIDE", actual_col), drop = FALSE]
        names(df_sub)[names(df_sub) == actual_col] <- "Quantity"
        df_sub
      }))

      if (is.null(df_group) || nrow(df_group) == 0) next

      agg <- df_group %>%
        dplyr::group_by(PEPTIDE) %>%
        dplyr::summarise(value = max(Quantity, na.rm = TRUE), .groups = "drop")
    } else {
      # Sample mode: each item is a sample name
      df_group <- dplyr::bind_rows(lapply(group_items, function(s) {
        df <- lst[[s]]
        cols <- dplyr::intersect(c("PEPTIDE", quantity_cols), colnames(df))
        if (length(cols) < 2) return(NULL)
        df[, cols, drop = FALSE]
      }))

      if (is.null(df_group) || nrow(df_group) == 0) next

      df_group <- df_group[df_group$PEPTIDE %in% allowed_peptides, ]

      agg <- df_group %>%
        dplyr::group_by(PEPTIDE) %>%
        dplyr::summarise(value = max(dplyr::c_across(dplyr::any_of(quantity_cols)), na.rm = TRUE), .groups = "drop")
    }

    pep_mat[agg$PEPTIDE, g] <- agg$value
  }

  pep_mat
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
plot_unique_counts <- function(lst, column, y_label, transform_fn = identity, color = "default") {
  stats <- lapply(names(lst), function(sample_name) {
    df <- lst[[sample_name]]
    if (!column %in% colnames(df)) return(NULL)
    values <- df[[column]]
    
    # Apply optional transformation (e.g., extract protein prefixes)
    values <- transform_fn(values)
    
    data.frame(
      Sample = sample_name,
      Count = length(unique(values)),
      stringsAsFactors = FALSE
    )
  }) %>% do.call(rbind, .)
  
  if (color == "default") {
    color <- NULL  # ggplot will use default fill colors
  } else {
    color <- viridis(length(lst), option = color)
  }
  
  p <- ggplot(stats, aes(x = Sample, y = Count, fill = Sample)) +
    geom_col() +
    geom_text(aes(label = Count), vjust = -0.3, size = 5) +
    labs(x = "Sample", y = y_label) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  if (!is.null(color)) {
    p <- p + scale_fill_manual(values = color)
  }
  p
}

plot_histogram <- function(df, column, x_label, title_name = "RT plot", color = "default") {
  
  # Create 1-minute bins (change to 60 if your data is in seconds)
  breaks <- seq(
    floor(min(df[[column]], na.rm = TRUE)),
    ceiling(max(df[[column]], na.rm = TRUE)),
    by = 1
  )
  
  if (color == "default") {
    color <- NULL  # ggplot will use default fill colors
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
  
  if (!is.null(color)) {
    p <- p + scale_fill_manual(values = color)
  }
  
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
      geom_text_repel(stat = "count", aes(label = ..count..), position=position_stack(vjust = 0.5), direction="y") +
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
  combined <- bind_rows(lapply(names(lst), function(name) {
    df <- lst[[name]]
    if (!column %in% colnames(df)) return(NULL)  # skip samples missing the column
    df <- df[, column, drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df$Sample <- name
    
    # Apply optional transformation
    if (!is.null(transform)) {
      df[[column]] <- transform(df[[column]])
    }
    
    df
  }))
  
  if (color == "default") {
    color <- NULL  # ggplot will use default fill colors
  } else {
    color <- viridis(length(lst), option = color)
  }
  
  # Automatic x-axis label if not provided
  if (is.null(x_label)) x_label <- column
  
  p <- ggplot(combined, aes_string(x = column, color = "Sample", fill = "Sample")) +
    geom_density(alpha = alpha) +
    labs(x = x_label, y = "Density") +
    theme_minimal()
  
  if (!is.null(color)) {
    p <- p + scale_fill_manual(values = color)
  }
  
  return(p)
}

plot_violin <- function(lst, column, x_label = "Sample", y_label = NULL, title = NULL, add_boxplot = TRUE, color = "default") {
  # Combine all samples
  df_all <- bind_rows(lapply(names(lst), function(s_name) {
    df <- lst[[s_name]]
    if (!(column %in% colnames(df)) || nrow(df) == 0) return(NULL)
    df$Sample <- s_name
    df[, c(column, "Sample"), drop = FALSE]
  }))
  
  req(nrow(df_all) > 0)
  
  if (is.null(y_label)) y_label <- column
  if (is.null(title)) title <- paste(column, "Distribution Across Samples")
  
  if (color == "default") {
    color <- NULL  # ggplot will use default fill colors
  } else {
    color <- viridis(length(lst), option = color)
  }
  
  p <- ggplot(df_all, aes(x = Sample, y = .data[[column]], fill = Sample)) +
    geom_violin(trim = TRUE, alpha = 0.6)
  
  if (add_boxplot) {
    p <- p + geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA)
  }
  
  p <- p + theme_minimal() +
    labs(title = title, x = x_label, y = y_label) +
    theme(legend.position = "none")
  
  if (!is.null(color)) {
    p <- p + scale_fill_manual(values = color)
  }
  
  return(p)
}

plot_length_distribution <- function(lst, color = "default") {
  req(lst)
  
  if (color == "default") {
    color <- NULL  # ggplot will use default fill colors
  } else {
    cols <- viridis(length(lst), option = color)
    color <- setNames(cols, length(lst))
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

dynamic_range_plot <- function(df, data_col = "MAX_QUANTITY", title_name = "Dynamic range plot", name_col = "PROTEIN", 
                               gene = "", rev_rank = TRUE, gene_regex = "(?<=GN=)[0-9A-Z//-]+") {
  
  # --- Safety ---
  if (!data_col %in% colnames(df)) {
    stop(paste("Column", data_col, "not found in dataframe"))
  }
  
  
  df[[data_col]][df[[data_col]] == 0] <- NA #convert 0 to NA, to avoid log2(0) = -inf error.
  if (all(is.na(df[[data_col]]))) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No data points available to display", size = 5) +
        theme_void()
    )
  }
    
  # --- Ranking ---
  df$Rank <- if (rev_rank) {
    rank(-df[[data_col]], ties.method = "first")
  } else {
    rank(df[[data_col]], ties.method = "first")
  }
  
  # --- Base plot (always present) ---
  p <- ggplot(df, aes(x = Rank, y = log2(.data[[data_col]]))) +
    geom_point() +
    labs(title = title_name)
  
  # --- Highlighting logic ---
  if (!is.null(name_col) && gene != "") {
    
    df$Gene <- stringr::str_extract(df[[name_col]], gene_regex)
    
    highlight_df <- df %>%
      dplyr::filter(!is.na(Gene),
                    stringr::str_detect(Gene, gene))
    
    if (nrow(highlight_df) > 0) {
      p <- p +
        geom_point(
          data = highlight_df,
          aes(x = Rank, y = log2(.data[[data_col]])),
          color = "red",
          size = 2
        ) +
        ggrepel::geom_text_repel(
          data = highlight_df,
          aes(x = Rank, y = log2(.data[[data_col]]), label = Gene),
          color = "red",
          box.padding = 0.5,
          max.overlaps = Inf
        )
    }
  } 
  
  p + theme_minimal()
}

dynamic_range_plot_combined <- function(df_list, data_col = "MAX_QUANTITY", color = "default") {
  
  # Filter out empty or insufficient samples
  df_list <- df_list[sapply(df_list, function(df) nrow(df) >= 10 && data_col %in% colnames(df))]
  
  if (length(df_list) == 0) {
    plot.new()
    text(0.5, 0.5, "Not enough peptides to plot combined dynamic range",
         cex = 1.2, col = "red")
    return(invisible(NULL))
  }
  
  
  if (color == "default") {
    color <- NULL  # ggplot will use default fill colors
  } else {
    color <- viridis(length(df_list), option = color)
  }
  
  df_all <- bind_rows(lapply(names(df_list), function(nm) {
    df <- df_list[[nm]]
    df <- df[, data_col, drop = FALSE]
    df <- df[!is.na(df[[data_col]]), , drop = FALSE]  # remove NAs
    if (nrow(df) == 0) return(NULL)
    df$Sample <- nm
    df
  }))
  
  # If bind_rows produced nothing
  if (is.null(df_all) || nrow(df_all) == 0 || all(is.na(df_all[[data_col]]))) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No data points available to display", size = 5) +
        theme_void()
    )
  }
  
  df_all <- df_all %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(Rank = rank(-.data[[data_col]], ties.method = "first")) %>%
    ungroup()
  
  p <- ggplot(df_all, aes(x = Rank, y = log2(.data[[data_col]]), color = Sample)) +
    geom_point(alpha = 0.5, size = 1.5) +
    theme_minimal() +
    labs(
      title = "Combined Dynamic Range for All Samples",
      x = "Rank",
      y = "log2(Intensity)"
    )
  
  if (!is.null(color)) {
    p <- p + scale_color_manual(values = color)
  }
  return(p)
}

generate_scatterplots <- function(lst, color = "default") {
  lapply(names(lst), function(s_name) {
    df <- lst[[s_name]]
    
    if (nrow(df) == 0) {
      p <- ggplot() + 
        annotate("text", x = 0.5, y = 0.5,
                 label = paste("No k0/mz data for sample:", s_name),
                 size = 6, color = "red") +
        theme_void()
      return(p)
    }
    
    charges <- sort(unique(df$CHARGE))
    
    if (color == "default") {
      cols <- NULL
    } else {
      cols <- viridis(length(unique(charges)), option = color)
      names(cols) <- charges
    }
    
    p <- ggplot(df, aes(x = MZ, y = K0, color = as.factor(CHARGE))) +
      geom_point(size = 1.5, alpha = 0.7) +
      theme_minimal() +
      labs(title = paste("1/k0 vs m/z:", s_name),
           x = "m/z",
           y = "1/k0",
           color = "Charge")
    if (!is.null(cols)) {
      p <- p + scale_color_manual(values = cols)
    }
    return(p)
  })
}

plot_completeness <- function(lst, spectra_cols, percent = FALSE, title = "Data Completeness",color = "default") {
  
  make_completeness <- function(df, spectra_cols) {
    
    spectra_cols <- dplyr::intersect(spectra_cols, colnames(df))
    
    df %>%
      transmute(
        Non_NA_Count = rowSums(!is.na(across(all_of(spectra_cols))))
      ) %>%
      arrange(desc(Non_NA_Count)) %>%
      dplyr::mutate(
        Rank = row_number(),
        Percent_Non_NA = 100 * Non_NA_Count / length(spectra_cols),
        Percent_Rank = 100 * Rank / max(Rank)
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

plot_upset <- function(data_list, min_size = 2, title = "Upset Plot",stripped = TRUE, min_size_is_percent = TRUE) {
  
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
    intersect = names(df_upset),
    name = "Sample",
    min_size = min_size,
    width_ratio=0.1
  ) +
    ggtitle(title))
  
  return(p)
}

plot_shared_peptide <- function(lst, color = "default", mode = c("count", "percent"), percent_type = c("min", "union")) {
  
  mode <- match.arg(mode)
  percent_type <- match.arg(percent_type)
  
  peptides <- lapply(lst, function(df) unique(as.character(df$STRIPPED)))
  samples <- names(peptides)
  n <- length(samples)
  
  shared_mat <- matrix(0, nrow = n, ncol = n,
                       dimnames = list(samples, samples))
  
  for (i in samples) {
    for (j in samples) {
      shared <- length(dplyr::intersect(peptides[[i]], peptides[[j]]))
      
      if (mode == "count") {
        shared_mat[i, j] <- shared
        
      } else {
        denom <- switch(
          percent_type,
          min   = min(length(peptides[[i]]), length(peptides[[j]])),
          union = length(union(peptides[[i]], peptides[[j]]))
        )
        
        shared_mat[i, j] <- ifelse(denom > 0, 100 * shared / denom, NA)
      }
    }
  }
  
  # Color mapping
  if (color == "default") {
    col_fun <- if (mode == "count") {
      colorRamp2(c(min(shared_mat), max(shared_mat, na.rm = TRUE)), c("white", "red"))
    } else {
      colorRamp2(c(0, 100), c("white", "red"))
    }
  } else {
    cols <- viridis(100, option = color)
    col_fun <- colorRamp2(
      range(shared_mat, na.rm = TRUE),
      c(cols[1], cols[100])
    )
  }
  
  Heatmap(
    shared_mat,
    name = if (mode == "count") "# shared peptides" else "% shared peptides",
    col = col_fun,
    na_col = "grey90"
  )
}

plot_pairwise_peptide_quant_correlation <- function(lst, method = "pearson", min_shared = 2, color = "default", cluster = c("none", "rows", "columns", "both")) {
 
  # Keep only peptide -> quantity mapping per sample
  peptides <- lapply(lst, function(df) {
    df %>%
      dplyr::select(STRIPPED, MAX_QUANTITY) %>%
      dplyr::group_by(STRIPPED) %>%
      dplyr::summarise(MAX_QUANTITY = max(MAX_QUANTITY, na.rm = TRUE), .groups = "drop") #this collapses for cases of PTM peptidoforms.
  })
  
  samples <- names(peptides)
  n <- length(samples)
  cluster <- match.arg(cluster)
  
  cor_mat <- matrix(NA_real_, nrow = n, ncol = n,
                    dimnames = list(samples, samples))
  
  for (i in samples) {
    for (j in samples) {
      
      merged <- merge(
        peptides[[i]],
        peptides[[j]],
        by = "STRIPPED",
        suffixes = c("_i", "_j")
      )
      
      if (nrow(merged) >= min_shared) {
        cor_mat[i, j] <- cor(
          merged$MAX_QUANTITY_i,
          merged$MAX_QUANTITY_j,
          method = method,
          use = "complete.obs"
        )
      }
    }
  }
  
  if (color == "default") {
    col_fun  <- colorRamp2(c(min(cor_mat), max(cor_mat)), c("white", "red"))
  } else {
    cols <- viridis(100, option = color)
    col_fun <- colorRamp2(range(cor_mat, na.rm = TRUE), c(cols[1], cols[100]))
  }
  
  Heatmap(
    cor_mat,
    name = paste(method, "correlation"),
    col = col_fun,
    na_col = "grey90",
    
    cluster_rows    = cluster %in% c("rows", "both"),
    cluster_columns = cluster %in% c("columns", "both"),
    
    clustering_distance_rows    = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows      = "complete",
    clustering_method_columns   = "complete"
  )
}

plot_heatmap <- function(mat, color = "default", log_transform = FALSE, cluster = "both", transpose = FALSE, row_groups = NULL, col_groups = NULL, label = "QUANTITY") {
  
  mat[is.na(mat)] <- 0
  
  # Log-transform if requested
  if (log_transform) {
    mat <- log10(mat)
    mat[is.infinite(mat)] <- 0
  }
  
  if (transpose) {
    mat <- t(mat)
  }
  
  # Define color function
  mat_range <- range(mat, na.rm = TRUE)
  if (!is.finite(mat_range[1]) || mat_range[1] == mat_range[2]) {
    mat_range <- c(mat_range[1] - 0.5, mat_range[1] + 0.5)
  }
  if (color == "default") {
    col_fun <- colorRamp2(mat_range, c("white", "red"))
  } else {
    cols <- viridis(100, option = color)
    col_fun <- colorRamp2(mat_range, c(cols[1], cols[100]))
  }
  
  # ---- Handle grouping vs clustering ----
  
  # Rows
  if (!is.null(row_groups)) {
    row_groups <- as.factor(row_groups)
    cluster_rows <- FALSE
  } else {
    cluster_rows <- cluster %in% c("rows", "both")
  }
  
  # Columns
  if (!is.null(col_groups)) {
    col_groups <- as.factor(col_groups)
    cluster_cols <- FALSE
  } else {
    cluster_cols <- cluster %in% c("columns", "both")
  }
  
  # ---- Build heatmap ----
  
  Heatmap(
    mat,
    name = if (log_transform) paste("log10(",label, ")") else label,
    col = col_fun,
    na_col = "grey90",
    
    # clustering
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    
    clustering_distance_rows    = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows      = "complete",
    clustering_method_columns   = "complete",
    
    # grouping (splitting)
    row_split = row_groups,
    column_split = col_groups,
    
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10)
  )
}

plot_PCA <- function(lst, color = "default") {
  peptides <- lapply(lst, function(df) {
    df %>%
      dplyr::select(STRIPPED, MAX_QUANTITY) %>%
      dplyr::group_by(STRIPPED) %>%
      dplyr::summarise(MAX_QUANTITY = max(MAX_QUANTITY, na.rm = TRUE), .groups = "drop") #this collapses for cases of PTM peptidoforms.
  })
  
  # Merge all samples into one data frame
  df_merged <- Reduce(function(x, y) full_join(x, y, by = "STRIPPED"), peptides)
  
  # Assign column names
  colnames(df_merged) <- c("STRIPPED", names(lst))
  
  # Convert to matrix: rows = peptides, columns = samples
  mat <- as.matrix(df_merged[,-1])
  rownames(mat) <- df_merged$STRIPPED
  
  #WIP: need to check how to solve missing values.
  mat[is.na(mat)] <- 0
  
  pca <- prcomp(t(mat), scale. = TRUE) 
  
  df <- as.data.frame(pca$x)
  df$Sample <- rownames(df)
  
  p <- ggplot(df, aes(x = PC1, y = PC2, label = Sample)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5) +
    xlab(paste0("PC1 (", round(summary(pca)$importance[2,1]*100,1), "%)")) +
    ylab(paste0("PC2 (", round(summary(pca)$importance[2,2]*100,1), "%)")) +
    theme_minimal()
  
  if (color == "default") {
    return(p)
  } else {
    p <- p + scale_color_manual(values = color)
    return(p)
  }
}

plot_binders <- function(df, color = "default", percent = TRUE, alleles = NULL) {
  allele_cols <- grep("^HLA", colnames(df), value = TRUE)
  
  if (!is.null(alleles) && length(alleles) > 0) {
    allele_cols <- intersect(allele_cols, alleles)
  }
  
  df$min <- do.call(pmin, c(df[,allele_cols], na.rm=TRUE))
  
  df <- df %>%
    dplyr::mutate(
      class = dplyr::case_when(
        is.na(min)       ~ "NA",
        min <= 0.5       ~ "Strong",
        min <= 2         ~ "Weak",
        TRUE             ~ "Non-binder"
      )
    ) %>%
    dplyr::count(Set, class, name = "n") %>%
    dplyr::group_by(Set) %>%
    dplyr::mutate(
      percent = 100 * n / sum(n),
      class = factor(
        class,
        levels = c("Strong", "Weak", "Non-binder", "NA")
      )
    ) %>%
    dplyr::ungroup()
  
  if(percent) {
    p <- ggplot(df, aes(x = Set, y = percent, fill = class)) +
      geom_col() +
      labs(
        x = "Set",
        y = "Percent of peptides",
        fill = "Binding class"
      ) + 
      geom_text(
        aes(label = ifelse(percent > 5, sprintf("%.1f%%", percent), "")),
        position = position_stack(vjust = 0.5),
        size = 3
      ) +
      scale_y_continuous(labels = scales::percent_format(scale = 1))
  } else {
    p <- ggplot(df, aes(x = Set, y = n, fill = class)) +
      geom_col() +
      labs(
        x = "Set",
        y = "Number of peptides",
        fill = "Binding class"
      )
  }
  
  p + theme_minimal()

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

plot_aa_composition <- function(lst, color = "default",show_numbers = TRUE) {
  
  aa_order <- c(
    "A","C","D","E","F","G","H","I","K","L",
    "M","N","P","Q","R","S","T","V","W","Y",
    "U","X"   # U and X are those rare amino acids that are often not found at all.
  )
  
  build_aa_matrix <- function(lst) {
    
    comp_list <- lapply(lst, function(df) {
      if (!"STRIPPED" %in% colnames(df))
        stop("STRIPPED column missing")
      
      aa_comp_from_peptides(df$STRIPPED)
    })
    
    mat <- matrix(0,
                  nrow = length(comp_list),
                  ncol = length(aa_order),
                  dimnames = list(names(comp_list), aa_order))
    
    for (i in seq_along(comp_list)) {
      mat[i, names(comp_list[[i]])] <- comp_list[[i]]  # missing AAs stay 0
      }
    mat
  }
  
  mat <- build_aa_matrix(lst)
  
  bg <- unlist(lapply(lst, function(df) df$STRIPPED))
  bg <- aa_comp_from_peptides(bg)
  bg_full <- numeric(length(aa_order))
  names(bg_full) <- aa_order
  bg_full[names(bg)] <- bg
  bg <- bg_full
  
  mat_diff <- sweep(mat, 2, bg, FUN = "-")
  
  #This was calculated separately from human_20586_2023-6-29.fasta
  human_ref <- c(
    A = 7.0114394770,
    C = 2.3020849101,
    D = 4.7350149446,
    E = 7.1070651632,
    F = 3.6501492529,
    G = 6.5774142723,
    H = 2.6236727595,
    I = 4.3322364677,
    K = 5.7304201058,
    L = 9.9669958082,
    M = 2.1321177950,
    N = 3.5837093961,
    P = 6.3150838295,
    Q = 4.7674806718,
    R = 5.6396265687,
    S = 8.3358604917,
    T = 5.3513901944,
    U = 0.0003157121,
    V = 5.9612758066,
    W = 1.2143339606,
    X = 0.0005875753,
    Y = 2.6617248369
  )
  
  if (color == "default") {
    col_fun <- colorRamp2(
      c(min(mat_diff), 0, max(mat_diff)),
      c("blue", "white", "red")
    )
  } else {
    cols <- viridis(100, option = color)
    col_fun <- colorRamp2(
      c(min(mat_diff), 0, max(mat_diff)),
      c(cols[1], cols[50], cols[100])
    )
  }
  
  # ---- numbers in cells ----
  cell_fun <- if (show_numbers) {
    function(j, i, x, y, w, h, fill) {
      grid::grid.text(
        sprintf("%.1f", mat[i, j]),
        x, y,
        gp = grid::gpar(fontsize = 8)
      )
    }
  } else NULL
  
  # ---- background bars (top annotation) ----
  top_anno <- NULL
  top_anno <- HeatmapAnnotation(
    "Uniprot Human" = anno_barplot(
      human_ref,
      gp = grid::gpar(fill = "grey70"),
      border = FALSE,
      height = unit(3, "cm"),
      cell_fun = function(j, x, y, width, height, fill) {
        grid::grid.rect(x, y, width, height, gp = grid::gpar(fill = fill))
        grid::grid.text(sprintf("%.1f", human_ref[j]),
                        x, y + height/2 + unit(1, "mm"),
                        gp = grid::gpar(fontsize = 8, col = "white"))
      }
    ),
    "Background" = anno_barplot(
      bg,
      gp = grid::gpar(fill = "grey70"),
      border = FALSE,
      height = unit(3, "cm"),
      cell_fun = function(j, x, y, width, height, fill) {
        grid::grid.rect(x, y, width, height, gp = grid::gpar(fill = fill))
        grid::grid.text(sprintf("%.1f", bg[j]),
                        x, y + height/2 + unit(1, "mm"),
                        gp = grid::gpar(fontsize = 8))
      }
    ),
    show_annotation_name = TRUE
  )
  
  Heatmap(
    mat_diff,
    name = "Abs. diff. to background",
    col = col_fun,
    top_annotation = top_anno,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    rect_gp = grid::gpar(col = "white"),
    cell_fun = cell_fun,
    row_names_side = "left",
    column_names_rot = 45
  )
}

#------Statistical caluclations-------
#calcualte group comparison statistics
compute_group_comp_stats <- function(lst, groups, allowed_peptides_g1, allowed_peptides_g2,
                                     g1, g2, pep_col, quantity_cols, col_map = NULL) {
  
  use_measurements <- !is.null(col_map)
  grp_g1 <- groups[[g1]]
  grp_g2 <- groups[[g2]]

  # Strip $expr so grp is consistently list(name, measurement) or a plain char vector
  if (use_measurements) {
    if (is.list(grp_g1) && !is.null(grp_g1$measurement))
      grp_g1 <- list(name = grp_g1$name, measurement = grp_g1$measurement)
    if (is.list(grp_g2) && !is.null(grp_g2$measurement))
      grp_g2 <- list(name = grp_g2$name, measurement = grp_g2$measurement)
  }

  if (!use_measurements) {
    n_g1 <- sum(sapply(grp_g1, function(s) length(intersect(quantity_cols, colnames(lst[[s]])))))
    n_g2 <- sum(sapply(grp_g2, function(s) length(intersect(quantity_cols, colnames(lst[[s]])))))
  } else {
    n_g1 <- if (is.list(grp_g1)) length(grp_g1$measurement) else length(grp_g1)
    n_g2 <- if (is.list(grp_g2)) length(grp_g2$measurement) else length(grp_g2)
  }
  use_limma <- n_g1 == 1 || n_g2 == 1

  # ── Build long-format data ──────────────────────────────────────────────────
  build_long <- function(grp_items, grp_name, allowed) {
    if (!use_measurements) {
      dplyr::bind_rows(lapply(grp_items, function(s) {
        df <- lst[[s]]
        cols <- intersect(c(pep_col, quantity_cols, "PROTEIN"), colnames(df))
        df[df[[pep_col]] %in% allowed, cols, drop = FALSE]
      })) %>%
        tidyr::pivot_longer(cols = dplyr::any_of(quantity_cols),
                            names_to = "Sample", values_to = "Quantity") %>%
        dplyr::mutate(Group = grp_name)
    } else {
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
        if (!"PROTEIN" %in% colnames(df)) df$PROTEIN <- NA_character_
        df$Sample <- meas
        df$Group  <- grp_name
        df
      }))
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
    pval_df     <- safe_limma_pvals(df_long, g1, g2, pep_col)
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

generate_dynrange_grid <- function(lst, data_col = "MAX_QUANTITY") {
  
  lapply(names(lst), function(sample_name) {
    df <- lst[[sample_name]]
    if (nrow(df) < 10) {
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = "Not enough peptides", size = 6) + 
        theme_void()
    } else {
      dynamic_range_plot(df, data_col, paste("Dynamic Range –", sample_name))
    }
  })
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
  system2("wsl", c("bash", "--login", "-c", shQuote(cmd)), stdout = NULL)
  
  output_file  # return path; caller reads it
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
  
  keep_cols    <- c("Peptide", grep("_Rank$", colnames(res), value = TRUE))
  res          <- res[, keep_cols, drop = FALSE]
  colnames(res) <- gsub("_Rank$", "", colnames(res))
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
    lapply(function(x) sapply(strsplit(x, "\\|"), `[`, 1)) %>%
    unlist() %>%
    unique()
  
  check <- safe_validate(!is.null(uni_ids) && nrow(as.data.frame(uni_ids)) > 0, "No significant IDs")
  if (!is.null(check)) return(check)
  
  gene_map <- suppressWarnings(clusterProfiler::bitr(uni_ids, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  
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
    strsplit(";") %>%                # split multiple proteins
    lapply(function(x) sapply(strsplit(x, "\\|"), `[`, 1)) %>%  # take first part of each
    unlist() %>%
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