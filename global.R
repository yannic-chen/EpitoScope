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
#Here we initiate all the possible column names important for us from all different input formats
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
    QUANTITY       = c("area"),                 # prefer area later
    SPECTRA        = c("x.feature", "X.Spec"),  #X.Spec for PEAKS 11 Online. X.Feature for PEAKS 12 studio. PEAKS 13 returns both. Prefer x.spec over x.feature.
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
    PROTEIN      = c("Protein_Mapped.Proteins"),               #This column is created later from Protein and Mapped.Protein column. #combined_modified_peptide.tsv, combined_peptide.tsv only has Protein column.
    QUANTITY       = c("maxlfq.intensity", "Intensity"),         #prefer maxlfq.intensity (exist in peptide.tsv, but not psm.tsv).
    SPECTRA        = c("Spectrum", "Spectral.Count"),            #However, this will be transformed anyway.
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
    PROTEIN      = c("Protein.names"),   #Can also switch with Protein.IDs, Protein.Group or Genes (although some proteins lack gene, like the CONTAs)
    QUANTITY       = c("D..data"), #Here NA means not found I guess
    SPECTRA        = c("D..data"), #However, this will be transformed anyway.
    PTM            = c("")                  # not present in report.pr_matrix.tsv
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

#-----------Helper functions------------------

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
  
  # Check if data exists
  if (is.null(data)) {
    stop_with_msg("Error: Data not available.")
  }
  
  #HANDLE LIST OF DATAFRAMES
  if (is.list(data)) {
    
    # Check if at least one dataframe exists
    if (length(data) == 0) {
      stop_with_msg("Error: Empty data list.")
    }
    
    # Check required columns across list
    if (!is.null(required_cols)) {
      has_col <- vapply(data, function(df) {
        is.data.frame(df) && all(required_cols %in% colnames(df))
      }, logical(1))
      
      if (!any(has_col)) {
        stop_with_msg(paste("Error: None of the datasets contain column(s):",
                            paste(required_cols, collapse = ", ")))
      }
      
      # NA checks (only on valid dfs)
      if (na_policy != "ignore") {
        valid_dfs <- data[vapply(data, function(df) {
          is.data.frame(df) && all(required_cols %in% colnames(df))
        }, logical(1))]
        
        if (length(valid_dfs) == 0) {
          stop_with_msg("Error: No valid dataframes for NA check.")
        }
        
        if (na_policy == "any") {
          if (any(vapply(valid_dfs, function(df) {
            any(is.na(df[[required_cols]]))
          }, logical(1)))) {
            stop_with_msg(paste("Error: NA values found in column:", required_cols))
          }
        }
        
        if (na_policy == "all") {
          if (any(vapply(valid_dfs, function(df) {
            all(is.na(df[[required_cols]]))
          }, logical(1)))) {
            stop_with_msg(paste("Error: Column all NA:", required_cols))
          }
        }
      }
    }
    return(TRUE)
  }
  
  #SINGLE DATAFRAME (fallback)
  if (is.data.frame(data)) {
    
    if (!is.null(required_cols)) {
      missing_cols <- setdiff(required_cols, names(data))
      if (length(missing_cols) > 0) {
        stop_with_msg(paste("Error: Missing column(s):", paste(missing_cols, collapse = ", ")))
      }
      
      if (na_policy != "ignore") {
        if (na_policy == "any" && any(is.na(data[[required_cols]]))) {
          stop_with_msg(paste("Error: NA values in column:", required_cols))
        }
        
        if (na_policy == "all" && all(is.na(data[[required_cols]]))) {
          stop_with_msg(paste("Error: Column all NA:", required_cols))
        }
      }
    }
    
    return(TRUE)
  }
  
  stop_with_msg("Error: Unsupported data type.")
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
column_schema$Generic <- build_generic_schema(column_schema)

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

find_transform_column <- function(df, list, name) {
  # Match independent of capitalization
  match <- which(tolower(colnames(df)) %in% tolower(list))
  
  if (length(match) == 1) {
    original_name <- colnames(df)[match]
    colnames(df)[match] <- name
  } else if (length(match) > 1) {
    stop("Multiple columns detected: ", name, ": ",
         paste(colnames(df)[match], collapse = ", "))
  } else {
    # No match found
    print(paste0("No matching column found for: ", name))
    return(NULL)
  }
  
  return(list(
    df = df,
    matched_column = original_name
  ))
}

transform_columns <- function(df, schema, software, targets = c("PEPTIDE", "STRIPPED", "LENGTH", "MASS", "MZ", "SCORE", "CHARGE", "RT", "PPM", "PROTEIN")) {
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
    mutate(
      RawFile = sub("\\.\\d+\\..*$", "", Spectrum)
    )
  
  num_cols <- df %>%
    dplyr::select(where(is.numeric)) %>%
    names() %>%
    setdiff(c("RawFile", "Modified.Peptide"))
  
  char_cols <- df %>%
    dplyr::select(where(is.character)) %>%
    names() %>%
    setdiff(c("RawFile", "Modified.Peptide", "Spectrum", "Spectrum.File"))
  
  df %>%
    group_by(RawFile, Modified.Peptide) %>%
    slice_max(
      .data[[intersect(c("PeptideProphet.Probability", "Probability"), names(df))[1]]], #Try "PeptideProphet.Probability", then "Probability" column
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup()
}

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
      df <- df %>% mutate(Modified.Peptide = coalesce(na_if(Modified.Peptide, ""), na_if(Peptide, "")))
      
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
        group_by(Modified.Peptide) %>%
        slice_max(
          .data[[intersect(c("PeptideProphet.Probability", "Probability"), names(df))[1]]], #Try "PeptideProphet.Probability", then "Probability" column
          n = 1,
          with_ties = FALSE
        ) %>%
        ungroup() %>%
        dplyr::select(-Intensity)
      
      if (nrow(df) != nrow(df_wide)) {
        stop("Error: The number of rows in df and df_wide do not match!")
      }
      
      df <- left_join(df, df_wide, by = "Modified.Peptide")
      
    }
    
    #Here we combine the Proteins and Mapped.Proteins together, to get all Accessions
    df <- df %>%
      mutate(
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
      group_by(Modified.Sequence) %>%
      arrange(Q.Value, .by_group = TRUE) %>%
      slice_head(n = 1) %>%
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
  
  #Fill in some missing columns.
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
        df <- df %>% mutate(K0 = get_midpoint(K0))
      }
    } else {
      df$K0 <- 0
      message("Ion mobility column missing → set to 0")
      original <- rbind(original, data.frame(final_name = "K0", original_name = "[no Ion Mobility column]", stringsAsFactors = FALSE))
    }
  }
  
  #Area/Intensity <- make this more elegant.
  sample  <- which(startsWith(tolower(colnames(df)), "area"))
  if (length(sample) == 0) {
    sample  <- which(startsWith(tolower(colnames(df)), "maxlfq.intensity"))
  } 
  if (length(sample) == 0) {
    sample  <- grep("intensity", tolower(colnames(df)))
  }
  if (length(sample) == 0) { #For DIANN report.pr_matrix.tsv this assumes that the file name is the absolute path of the file. Meaning it starts with the Drive letter e.g. D:\ -> d..
    sample  <-  which(grepl("^[a-z]\\.\\.", tolower(colnames(df)))) 
  }
  if (length(sample) == 0) {
    stop("No Area or Intensity columns")
  }
  
  original <- bind_rows(original, data.frame(final_name="QUANTITY", original_name=list(colnames(df[,sample, drop = FALSE]))))
  
  df$MAX_QUANTITY <- apply(df[, sample, drop = FALSE], 1, function(x) {
    if (all(is.na(x))) {
      NA
    } else {
      max(x, na.rm = TRUE)
    }
  })
  
  #We need another column to know how many samples the peptide was found in.
  #In X.Spec for PEAKS 11 Online. X.Feature for PEAKS 12 studio. PEAKS 13 returns both. Prefer x.spec over x.feature.
  if (any(tolower(colnames(df)) == "x.spec")) {
    spec <- grep("^x\\.spec", tolower(colnames(df)))
    spec <- spec[colnames(df)[spec] != "X.Spec"] #This column in PEAKS sums all Spectra found by all samples, however, one sample can find the peptide multiple times.
  } else if (any(tolower(colnames(df)) == "x.feature")) {
    spec <- grep("x.feature", tolower(colnames(df)))
    spec <- spec[colnames(df)[spec] != "X.Feature"]
  } else {
    spec <- integer(0)
  }

  #In case we can use Quantity columns (sample) as indicator of Spectral match
  if (length(spec) == 0) {
    original <- bind_rows(original, data.frame(final_name="SPECTRA", original_name=list(colnames(df[,sample, drop = FALSE])))) #We use the same columns as Intensity for Spectra found per sample in Fragpipe
  } else {
    df[, spec][df[, spec] == 0] <- NA #For spectra column, the 0 actually means not identified. We need these to be converted to 0 to work with FragPipe format for data completeness plot.
    original <- bind_rows(original, data.frame(final_name="SPECTRA", original_name=colnames(df[, spec, drop = FALSE]), stringsAsFactors=FALSE))
  }
  
  #Only keep cols that were found
  keep_cols <- original %>%
    filter(!is.na(original_name)) %>%
    filter(purrr::map_lgl(final_name, ~ is.character(.x) && length(.x) == 1)) %>% #This removes entries where the original column names are a list of strings. These are added separately.
    pull(final_name)
  keep_cols <- c(keep_cols, colnames(df)[sample], colnames(df)[spec], "STRIPPED", "MAX_QUANTITY")
  
  #Now we can prepare the summary table since we have the columns to keep.
  original <- original %>% mutate(coalesced = do.call(coalesce, across(-1))) %>% dplyr::select(1, coalesced) %>% group_by(final_name) %>% summarise(coalesced_list = list(coalesced), .groups = "drop")
  
  return(list(df = df %>% dplyr::select(any_of(unique(keep_cols))), table = original))
}

#This is for creating a matrix used for the group-based peptide heatmap
prepare_peptide_matrix <- function(lst, groups, quantity_cols, group_peptide_sets) {
  
  peptides_all <- unique(unlist(group_peptide_sets))
  pep_mat <- matrix(NA_real_, 
                    nrow = length(peptides_all), 
                    ncol = length(groups),
                    dimnames = list(peptides_all, names(groups)))
  
  for (g in names(groups)) {
    samples <- groups[[g]]
    allowed_peptides <- group_peptide_sets[[g]]
    
    # Combine all samples in group
    df_group <- bind_rows(lapply(samples, function(s) {
      df <- lst[[s]]
      cols <- intersect(c("PEPTIDE", quantity_cols), colnames(df))
      if (length(cols) < 2) return(NULL)
      df[, cols, drop = FALSE]
    }))
    
    if (is.null(df_group) || nrow(df_group) == 0) next
    
    # Filter to allowed peptides
    df_group <- df_group[df_group$PEPTIDE %in% allowed_peptides, ]
    
    # Aggregate: take max (or mean) per peptide across samples
    agg <- df_group %>%
      group_by(PEPTIDE) %>%
      summarise(value = max(c_across(any_of(quantity_cols)), na.rm = TRUE), .groups = "drop")
    
    pep_mat[agg$PEPTIDE, g] <- agg$value
  }
  
  pep_mat
}

#-----------Plotting functions------------------
plot_unique_counts <- function(lst, column, y_label, transform_fn = identity, color = "default") {
  stats <- lapply(names(lst), function(sample_name) {
    df <- lst[[sample_name]]
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
  
  if (color == "default") color <- "steelblue"
  
  ggplot(df, aes_string(x = column)) +   # <- use column, not data_col
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
}

plot_stacked_bar <- function(lst, column, fill_label = NULL, rev_levels = TRUE, percentage = FALSE, color = "default") {
  combined <- bind_rows(lapply(names(lst), function(name) {
    df <- lst[[name]]
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
      group_by(Sample) %>%
      mutate(perc = 100 * (n() / n())) %>%  # initially gives 100%, will correct below
      ungroup()
    
    # Actually, we need % per factor level per sample
    perc_df <- combined %>%
      group_by(Sample, !!sym(column)) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(Sample) %>%
      mutate(perc = 100 * count / sum(count)) %>%
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
  combined <- combined %>% mutate(Length_bin = factor(LENGTH))
  
  # Compute counts and percentages
  count_df <- combined %>% count(Length_bin, Sample, name = "Count")
  percent_df <- count_df %>% group_by(Sample) %>%
    mutate(Percent = Count / sum(Count) * 100)
  
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
    ggseqlogo(peptides, method = 'bits') +
      ggtitle(title) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
      )
  } else {
    ggseqlogo(peptides, method = 'bits', namespace = namespace, seq_type = "other") +
      ggtitle(title) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
      )
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
          box.padding = 0.5
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
    group_by(Sample) %>%
    mutate(Rank = rank(-.data[[data_col]], ties.method = "first")) %>%
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
    
    spectra_cols <- intersect(spectra_cols, colnames(df))
    
    df %>%
      transmute(
        Non_NA_Count = rowSums(!is.na(across(all_of(spectra_cols))))
      ) %>%
      arrange(desc(Non_NA_Count)) %>%
      mutate(
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
  p <- ComplexUpset::upset(
    df_upset,
    intersect = names(df_upset),
    name = "Sample",
    min_size = min_size,
    width_ratio=0.1
  ) +
    ggtitle(title)
  
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
      shared <- length(intersect(peptides[[i]], peptides[[j]]))
      
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

plot_heatmap <- function(mat, color = "default", log_transform = FALSE, cluster = "both", transpose = FALSE) {
  
  mat[is.na(mat)] <- 0
  
  # Log-transform if requested
  if (log_transform) {
    mat <- log10(mat)  # avoid log10(0)
    mat[is.infinite(mat)] <- 0   # replace -Inf / Inf with 0
  }
  
  if (transpose) {
    mat <- t(mat)
  }
  
  # Define color function
  if (color == "default") {
    col_fun <- colorRamp2(c(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE)), c("white", "red"))
  } else {
    cols <- viridis(100, option = color)
    col_fun <- colorRamp2(range(mat, na.rm = TRUE), c(cols[1], cols[100]))
  }
  
  # Cluster options
  cluster_rows <- cluster %in% c("rows", "both")
  cluster_cols <- cluster %in% c("columns", "both")
  
  Heatmap(
    mat,
    name = if(log_transform) "log10(QUANTITY)" else "QUANTITY",
    col = col_fun,
    na_col = "grey90",
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    clustering_distance_rows    = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows      = "complete",
    clustering_method_columns   = "complete",
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

plot_binders <- function(df, color = "default", percent = TRUE) {
  allele_cols <- grep("^HLA", colnames(df), value = TRUE)
  
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
    mutate(
      Binder = case_when(
        is.na(Affinity) ~ "NA",
        Affinity < 0.5  ~ "Strong",
        Affinity >= 0.5 & Affinity <= 1 ~ "Weak",
        Affinity > 1 ~ "Non"
      )
    ) %>%
    group_by(Allele, Binder) %>%
    summarise(Count = n(), .groups = "drop") %>%
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
