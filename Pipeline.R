# ============================================================================
# Script Name:        EDA_pipeline.R
# Purpose:            This script outputs various plots for immunopeptidomics data
# Author:             Yannic Chen
# Date Created:       2025-11-14
# Last Modified:      2025-10-14
# Version:            0.1
# R Version:          4.5.2
# ============================================================================
# Details:
# This script contains functions to generate various plots for exploratory data analysis (EDA) of immunopeptidomics MS data.
# It should take in files from PEAKS and Fragpipe.
# 
# The script is structured into the following sections:
# 1. Loading required libraries
# 2. Defining plotting functions
# 3. Helper functions
# 4. Loading data
# 5. Performing EDA on raw data
# 6. Data cleaning
# 7. Performing EDA on cleaned data
#
# Design Philosophy:
# The data should be kept as close to their original state as possible. 
# For ease, all data should be collected within a list. And each function should take in the data list as input and extract the relevat info themselves.
# Transformations to obtain the relevant data should only happen within the plotting function.
# This ensures that the original data remains intact and allows flexibility and modularity.
# Plotting functions should do all the heavy lifting, so that they only a 1-liner is required to plot from the data list.
# Plotting functions should contain options of different styles (e.g., percent vs absolute) if it makes sense.
#
# To Do:
# Create a function to standardize data sources (PEAKS, Fragpipe).
# Think of a way to group samples (e.g., by patient, by condition).
# Implement statistical tests for group comparisons.

#--------------Loading Libraries----------------


library(ggplot2)
library(dplyr)
library(ggpointdensity)
library(viridis)
library(purrr)
library(shiny)
library(plotly)
library(ComplexUpset)
library(GGally)
library(ggrepel)
library(data.table)


#--------------Plotly Main----------------
#--------------Plotting Functions----------------
theme <- theme(
  # add border
  panel.border = element_rect(colour = "gray", fill = NA, linetype = 1),
  # modify text, axis and colour
  axis.text = element_text(colour = "black", family = "Times New Roman", size = 10),
  axis.title = element_text(colour = "black", face = "bold", family = "Times New Roman", size = 11),
  # title
  plot.title = element_text(hjust = 0.5, face = "bold", family = "Times New Roman", size = 14) #put in the middle
)

plot_barchart <- function(data_list, x_col, y_col, title, x_label, y_label) {
  data <- purrr::imap_dfr(
    data_list,
    ~ {
      # Extract vector from column
      values <- .x[[y_col]]
      tibble(
        !!x_col := .y,
        Count = dplyr::n_distinct(values)
      )
    }
  )
  
  ggplot_object <- ggplot(data, aes_string(x = x_col, y = "Count", fill = x_col)) +
    geom_bar(stat = "identity") +
    theme +
    labs(title = title, x = x_label, y = y_label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(aes(label = round(Count, 1)),  # Add labels
              vjust = -0.5,                 # slightly above the bar
              size = 3)
  
  return(ggplot_object)
}

plot_stacked_barchart <- function(data_list, x_col, fill_col, title, x_label, y_label, percent = FALSE) {
  # Combine data from list into long format
  data <- purrr::imap_dfr(
    data_list,
    ~ tibble(
      !!x_col := .y,
      !!fill_col := unlist(.x[[fill_col]])   # <-- unlist ensures it's a vector
    )
  )
  
  # Count occurrences per (x, fill)
  data_summary <- data %>%
    count(!!sym(x_col), !!sym(fill_col), name = "Count")
  
  # Convert to percentage if requested
  if (percent) {
    data_summary <- data_summary %>%
      group_by(!!sym(x_col)) %>%
      mutate(Count = 100 * Count / sum(Count)) %>%
      ungroup()
  }
  
  # Plot
  ggplot_object <- ggplot(data_summary, aes_string(x = x_col, y = "Count", fill = fill_col)) +
    geom_bar(stat = "identity", position = "stack") +
    theme +
    labs(title = title, x = x_label, y = y_label, fill = fill_col) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(ggplot_object)
}

plot_groupedbarchart <- function(data_list, x_col, title, x_label, y_label) {
  
  # Combine all dataframes and assign a Sample column from the list names
  data <- purrr::imap_dfr(
    data_list,
    ~ tibble(
      Sample = .y,
      !!x_col := .x[[x_col]]
    )
  )
  
  # Count occurrences
  data_summary <- data %>%
    count(Sample, !!sym(x_col), name = "Count")
  
  # Plot: grouped bar chart
  ggplot_object <- ggplot(data_summary, aes(
    x = !!sym(x_col),
    y = Count,
    fill = Sample
  )) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    theme +
    labs(title = title, x = x_label, y = y_label, fill = "Sample")
  return(ggplot_object)
}

plot_pointdensity <- function(data, x_col, y_col, title, x_label, y_label, adjust_val = 0.001) {
  ggplot_object <- ggplot(data, aes_string(x = x_col, y = y_col)) +
        geom_pointdensity(adjust = adjust_val, aes(color = after_stat(log(n_neighbors)))) +
        scale_color_viridis() +
        theme +
        labs(title = title, x = x_label, y = y_label, color = "Log Density")
  return(ggplot_object)
}

plot_violinplot <- function(data_list, x_col, y_col, title, x_label, y_label) {
  data <- purrr::imap_dfr(
    data_list,
    ~ tibble(
      !!x_col := .y,
      !!y_col := .x[[y_col]]
    )
  )
  ggplot_object <- ggplot(data, aes_string(x = x_col, y = y_col, fill = x_col)) +
        geom_violin(trim = FALSE) +
        theme +
        labs(title = title, x = x_label, y = y_label) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(ggplot_object)
}

plot_completeness <- function(data_named_list, percent = FALSE, title) {
  make_completeness <- function(df) {
    base <- df %>%
      transmute(
        Non_NA_Count = rowSums(select(., matches("^X\\.Spec\\.")) > 0)
      ) %>%
      count(Non_NA_Count, name = "Frequency")
    
    # sort descending to compute cumulative "at least X samples"
    base %>%
      arrange(desc(Non_NA_Count)) %>%
      mutate(CumFreq = cumsum(Frequency)) %>%
      arrange(Non_NA_Count)
  }
  
    
    df <- map(data_named_list, make_completeness) %>%
        bind_rows(.id = "Sample") %>%
        group_by(Sample) %>% 
        mutate(
        Percent_Non_NA_Count = 100 * Non_NA_Count / max(Non_NA_Count),
        Percent_Frequency    = 100 * Frequency / sum(Frequency)
        )
        
  if (percent) {
    ggplot_object <- ggplot(df, aes(Percent_Non_NA_Count, Percent_Frequency, color = Sample)) +
      geom_line() +
      theme +
      labs(
        title = title,
        x = "Percentage of samples",
        y = "Percentage of peptides",
        color = "Sample"
      )
  } else {
    ggplot_object <- ggplot(df, aes(Non_NA_Count, Frequency, color = Sample)) +
      geom_line() +
      theme +
      labs(
        title = title,
        x = "Number of samples",
        y = "Number of peptides",
        color = "Sample"
      )
  }
    return(ggplot_object)
}

plot_upset <- function(data_list, col, min_size = 100,title) {
  # Convert list of vectors to binary membership matrix
  l <- purrr::imap(data_list, ~ unique(.x[[col]]))
  all_items <- unique(unlist(l))
  
  df <- sapply(l, function(x) all_items %in% x) %>% as.data.frame()
  df$Item <- all_items

  df_upset <- df %>% select(-Item)
  
  # ComplexUpset plot
  p <- upset(df_upset, 
             intersect = names(df_upset),
             name = "Sample",
             min_size=min_size) +
    ggtitle(title %||% paste("UpSet Plot:", col))
  
  return(p)
}

plot_ggpair <- function(data_list, score_col, peptide_col, samples = NULL, title) {
  
  # Use all samples if none specified
  if (is.null(samples)) samples <- names(data_list)
  
  # Preprocess each dataframe: ensure PTMless column exists
  df_list <- lapply(data_list[samples], function(df) {
    # Create PTMless if requested
    if (peptide_col == "PTMless" && !"PTMless" %in% names(df)) {
      df$PTMless <- gsub("\\(.*?\\)|\\[.*?\\]|\\{.*?\\}", "", df$Peptide)
    }
    
    # Keep only peptide and score columns
    df <- df %>% select(all_of(c(peptide_col, score_col)))
    
    # For peptidoforms (duplicates), keep max score per peptide
    df <- df %>%
      group_by(.data[[peptide_col]]) %>%
      summarise(!!score_col := max(.data[[score_col]], na.rm = TRUE), .groups = "drop")
    
    df
  })
  
  # Name the score column for each sample
  df_list <- purrr::imap(df_list, ~ rename(.x, !!.y := !!sym(score_col)))
  
  # Merge all dataframes by peptide_col (inner join = only shared peptides)
  combined_df <- reduce(df_list, full_join, by = peptide_col)
  
  # Remove peptides with missing scores in any sample
  combined_df <- combined_df %>% drop_na()
  
  # Remove peptide_col for ggpairs, keep as rownames for reference
  rownames(combined_df) <- combined_df[[peptide_col]]
  combined_df <- combined_df %>% select(-all_of(peptide_col))
  
  # Create ggpairs plot
  p <- ggpairs(combined_df) +
    ggtitle(title %||% "Score comparison of shared peptides across samples")
  
  return(p)
}

plot_dynamic_range <- function(data_list, quantity_col, peptide_col, accession_col = NULL, title_name, gene_regex = "(?<=GN=)[0-9A-Z//-]+", 
                               highlight_regex = "", log_scale = TRUE, rev_rank = TRUE) {
  
  # Combine all samples into one long dataframe
  df <- purrr::imap_dfr(data_list, ~ {
    d <- .x
    
    # Create PTMless column if requested
    if(peptide_col == "PTMless" && !"PTMless" %in% names(d)) {
      d$PTMless <- gsub("\\(.*?\\)|\\[.*?\\]|\\{.*?\\}", "", d$Peptide)
    }
    
    # Keep only relevant columns
    keep_cols <- c(peptide_col, quantity_col, accession_col)
    d <- d %>% select(any_of(keep_cols)) %>%
      filter(!is.na(.data[[quantity_col]]))
    
    # Rank peptides by descending quantity
    if(rev_rank) {
      d <- d %>% arrange(desc(.data[[quantity_col]]))
    } else {
      d <- d %>% arrange(.data[[quantity_col]])
    }
    
    d <- d %>% mutate(Rank = row_number())
    
    tibble(
      Sample = .y,
      Peptide = d[[peptide_col]],
      Quantity = d[[quantity_col]],
      Rank = d$Rank,
      Accession = d[[accession_col]]
    )
  })
  
  # Log-transform quantity if requested
  if(log_scale) df <- df %>% mutate(Quantity = log10(Quantity + 1e-6))
  
  # Identify highlighted peptides by regex on Accession
  if(highlight_regex != "") {
    df$Highlight <- str_detect(df[[accession_col]], highlight_regex)
  } else {
    df$Highlight <- FALSE
  }
  
  # Base scatter plot
  p <- ggplot(df, aes(x = Rank, y = Quantity, color = Sample)) +
    geom_point(alpha = 0.7, size = 2) +
    theme_minimal() +
    labs(
      title = title_name,
      x = "Peptide Rank (highest abundance = rank 1)",
      y = ifelse(log_scale, paste0("log10(", quantity_col, ")"), quantity_col)
    ) +
    theme(axis.text.x = element_blank())
  
  # Highlight peptides using geom_text_repel
  if(any(df$Highlight)) {
    p <- p + 
      ggrepel::geom_text_repel(
        data = df[df$Highlight, ],
        aes(label = Peptide),
        color = "red",
        box.padding = 0.5,
        max.overlaps = Inf
      )
  }
  
  # Return interactive Plotly
  ggplotly(p, tooltip = c("Peptide", "y", "x", "color", "Accession"))
}

#--------------Helper functions----------------
remove_ptms <- function(x) {
  gsub("\\(.*?\\)|\\[.*?\\]|\\{.*?\\}", "", x)
}

#--------------Loading and transforming data------------------
BT15_tumor <- read.csv("C:/PostDoc/Poschke/BT15/PEAKS/[COGNATE]BT15.[SAREKv0.99b]_BT15_tumor/db.peptides.csv")
BT15_plasma <- read.csv("C:/PostDoc/Poschke/BT15/PEAKS/[COGNATE]BT15.[SAREKv0.99b]_BT15_plasma/db.peptides.csv")

EMB95_tumor <- read.csv("C:/PostDoc/Poschke/EMB95_01-009/PEAKS/[COGNATE]EMB95_01-009.[SAREKv0.99b]EMB95_tumor/db.peptides.csv")
EMB95_plasma <- read.csv("C:/PostDoc/Poschke/EMB95_01-009/PEAKS/[COGNATE]EMB95_01-009.[SAREKv0.99b]EMB95_plasma/db.peptides.csv")
EMB95_CSF <- read.csv("C:/PostDoc/Poschke/EMB95_01-009/PEAKS/[COGNATE]EMB95_01-009.[SAREKv0.99b]EMB95_CSF/db.peptides.csv")

NOA01_01_tumor <- read.csv("C:/PostDoc/Poschke/NOA21_01-001/PEAKS/[COGNATE]NOA21-01-001.[SAREKv0.99b]NOA21_01_001_tumor/db.peptides.csv")
NOA01_01_plasma <- read.csv("C:/PostDoc/Poschke/NOA21_01-001/PEAKS/[COGNATE]NOA21-01-001.[SAREKv0.99b]NOA21_01_001_plasma/db.peptides.csv")
NOA01_01_serum <- read.csv("C:/PostDoc/Poschke/NOA21_01-001/PEAKS/[COGNATE]NOA21-01-001.[SAREKv0.99b]NOA21_01_001_serum/db.peptides.csv")           

NOA01_05_tumor <- read.csv("C:/PostDoc/Poschke/NOA21_01-005/PEAKS/[COGNATE]NOA21_01-005.[SAREKv0.99b]NOA21-01-005_tumor/db.peptides.csv")
NOA01_05_serum <- read.csv("C:/PostDoc/Poschke/NOA21_01-005/PEAKS/[COGNATE]NOA21_01-005.[SAREKv0.99b]NOA21-01-005_serum/db.peptides.csv")
NOA01_05_plasma <- read.csv("C:/PostDoc/Poschke/NOA21_01-005/PEAKS/[COGNATE]NOA21_01-005.[SAREKv0.99b]NOA21-01-005_plasma/db.peptides.csv")

x1088_tumor <- read.csv("C:/PostDoc/Poschke/1088/PEAKS/[COGNATE]1088.1088_tumor/db.peptides.csv")
x1088_CSF <- read.csv("C:/PostDoc/Poschke/1088/PEAKS/[COGNATE]1088.1088_CSF/db.peptides.csv")
x1088_serum <- read.csv("C:/PostDoc/Poschke/1088/PEAKS/[COGNATE]1088.1088_serum/db.peptides.csv")

D170_04_tumor <- read.csv("C:/PostDoc/Poschke/D170_04/PEAKS/[COGNATE]D170-04.D170_04_tumor/db.peptides.csv")
D170_04_serum <- read.csv("C:/PostDoc/Poschke/D170_04/PEAKS/[COGNATE]D170-04.D170_04_serum/db.peptides.csv")

D170_05_tumor <- read.csv("C:/PostDoc/Poschke/D170_05/PEAKS/[COGNATE]D170_05.D170_05_tumor/db.peptides.csv")
D170_05_serum <- read.csv("C:/PostDoc/Poschke/D170_05/PEAKS/[COGNATE]D170_05.D170_05_serum/db.peptides.csv")

D170_15_tumor <- read.csv("C:/PostDoc/Poschke/D170_15/PEAKS/[COGNATE]D170_15.D170_15_tumor/db.peptides.csv")
D170_15_serum <- read.csv("C:/PostDoc/Poschke/D170_15/PEAKS/[COGNATE]D170_15.D170_15_serum/db.peptides.csv")

D170_20_tumor <- read.csv("C:/PostDoc/Poschke/D170_20/PEAKS/[COGNATE]D170_20.D170_20_tumor/db.peptides.csv")
D170_20_serum <- read.csv("C:/PostDoc/Poschke/D170_20/PEAKS/[COGNATE]D170_20.D170_20_serum/db.peptides.csv")

D170_48_tumor <- read.csv("C:/PostDoc/Poschke/D170_48/PEAKS/[COGNATE]D170_48.D170_48_tumor/db.peptides.csv")
D170_48_serum <- read.csv("C:/PostDoc/Poschke/D170_48/PEAKS/[COGNATE]D170_48.D170_48_serum/db.peptides.csv")
D170_48_plasma <- read.csv("C:/PostDoc/Poschke/D170_48/PEAKS/[COGNATE]D170_48.D170_48_plasma/db.peptides.csv")
D170_48_CSF <- read.csv("C:/PostDoc/Poschke/D170_48/PEAKS/[COGNATE]D170_48.D170_48_CSF/db.peptides.csv")

peaks12 <- read.csv("C:/Users/Yannic/Downloads/PEAKS12.db.peptides.csv")
peaks13 <- read.csv("C:/Users/Yannic/Downloads/PEAKS13.peptideDb.peptides.csv")
  
common_cols <- intersect(colnames(NOA01_05_plasma), colnames(NOA01_01_tumor))

data_list <- list(
  NOA01_05_plasma = NOA01_05_plasma,
  NOA01_05_serum  = NOA01_05_serum,
  NOA01_05_tumor  = NOA01_05_tumor,
  NOA01_01_plasma = NOA01_01_plasma,
  NOA01_01_serum  = NOA01_01_serum,
  NOA01_01_tumor  = NOA01_01_tumor,
  BT15_tumor = BT15_tumor,
  BT15_plasma = BT15_plasma,
  EMB95_tumor = EMB95_tumor,
  EMB95_plasma = EMB95_plasma,
  EMB95_CSF = EMB95_CSF,
  x1088_tumor = x1088_tumor,
  x1088_CSF = x1088_CSF,
  x1088_serum = x1088_serum,
  D170_04_tumor = D170_04_tumor,
  D170_04_serum = D170_04_serum,
  D170_05_tumor = D170_05_tumor,
  D170_05_serum = D170_05_serum,
  D170_15_tumor = D170_15_tumor,
  D170_15_serum = D170_15_serum,
  D170_20_tumor = D170_20_tumor,
  D170_20_serum = D170_20_serum,
  D170_48_tumor = D170_48_tumor,
  D170_48_serum = D170_48_serum,
  D170_48_plasma = D170_48_plasma,
  D170_48_CSF = D170_48_CSF
)

data_list <- purrr::map(
  data_list,
  ~ dplyr::mutate(.x, PTMless = remove_ptms(Peptide)) #Get PTMless sequence
)

data_list <- lapply(data_list, function(df) {
  df %>%
    mutate(Max_Area = do.call(pmax, c(select(., starts_with("Area.")), na.rm = TRUE)))
})

netMHcpan <- as.data.table(bind_rows(
  read.csv("8mer_netMHCpan4.2_human_20586_2023-6-29.csv"),
  read.csv("9mer_netMHCpan4.2_human_20586_2023-6-29.csv"),
  read.csv("10mer_netMHCpan4.2_human_20586_2023-6-29.csv"),
  read.csv("11mer_netMHCpan4.2_human_20586_2023-6-29.csv")
))
netMHCpan <- netMHcpan[!duplicated(netMHcpan$Peptide), ]

#--------------Statistic on Raw------------
# The following statistics are performed on the raw data before any filtering or cleaning.
# This is mainly for quality checking purpose.
  
# Compare number of IDs.
plot_barchart(data_list, x_col = "Sample", y_col = "Peptide", title = "[RAW] Number of Identified Peptidoform per Sample",
   x_label = "Sample", y_label = "Number of Peptidoform")

plot_barchart(data_list, x_col = "Sample", y_col = "PTMless", title = "[RAW] Number of Identified Peptides per Sample",
              x_label = "Sample", y_label = "Number of Peptides")

# Compare score distribution
plot_violinplot(data_list, x_col = "Sample", y_col = "X.10LgP", title = "[RAW] Score Distribution per Sample",
   x_label = "Sample Type", y_label = "Score (10LgP)")

# Data completeness curve
plot_completeness(data_list, percent = FALSE, title = "[RAW] Data Completeness Curves (absolute)")
plot_completeness(data_list, percent = TRUE, title = "[RAW] Data Completeness Curves (relative)")

# charge distribution
plot_stacked_barchart(data_list, x_col = "Sample", fill_col = "z", title = "[RAW] Charge Distribution per Sample", 
                      x_label = "Sample Type", y_label = "Count")

plot_stacked_barchart(data_list, x_col = "Sample", fill_col = "z", title = "[RAW] Charge Distribution per Sample", 
                      x_label = "Sample Type", y_label = "Fraction", percent = TRUE)

# Length distribution
plot_groupedbarchart(data_list, x_col = "Tag.Length",title = "Peptide Length Distribution",
                     x_label = "Length",y_label = "Count")

plot_upset(data_list, col = "PTMless", title = "Unique Peptides", min_size = 500)

plot_ggpair(data_list, score_col = "X.10LgP", peptide_col = "Peptide", title = "Score comparison between shared peptidoform")
plot_ggpair(data_list, score_col = "X.10LgP", peptide_col = "PTMless", title = "Score comparison between shared peptides")

plot_dynamic_range(data_list, quantity_col = "Max_Area", peptide_col = "PTMless", accession_col = "Accession", title = "Peptide Dynamic Range", highlight_regex = "HLA")


#--------------Data Cleaning----------------



#--------------Collecting plots----------------
