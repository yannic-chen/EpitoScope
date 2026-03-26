# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Script Name:        Shiny.R
# Purpose:            This script creates a shiny app for easy immunopeptidomics analysis
# Author:             Yannic Chen
# Date Created:       2026-01-06
# Last Modified:      2026-01-06
# Version:            0.1
# R Version:          4.5.2
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Details:
# This script contains information on how to preload data for use in the shiny app.
#
# Important:
# The named list variable must be called "preloaded_data".
# The netMHCpan pre-generated variable must be called "netMHCpan".
#


#Load data
sample1_tumor <- read.csv("path to file")
sample2_plasma <- read.delim("path to file")
library(arrow)
DIANN_file <- read_parquet("path to file")

#generate a named list. You can name your data however you want.
preloaded_data <- list(
  NAME_tumor = sample1_tumor,
  NAME_plasma  = sample2_plasma
)

#Load netMHCpan pre-generated data
netMHCpan <- as.data.table(bind_rows(
  read.csv("8mer_netMHCpan4.2_human_20586_2023-6-29.csv"),
  read.csv("9mer_netMHCpan4.2_human_20586_2023-6-29.csv"),
  read.csv("10mer_netMHCpan4.2_human_20586_2023-6-29.csv"),
  read.csv("11mer_netMHCpan4.2_human_20586_2023-6-29.csv")
))

#Deduplicating as a savety measure
netMHCpan <- netMHCpan[!duplicated(netMHCpan$Peptide), ]

preloaded_data <- list(
  seph_well = combined_modified_peptide[,c(1:16,grep("seph_well",colnames(combined_modified_peptide)))],
  seph_filter = combined_modified_peptide[,c(1:16,grep("seph_filter",colnames(combined_modified_peptide)))], 
  mag_well = combined_modified_peptide[,c(1:16,grep("mag_well",colnames(combined_modified_peptide)))]
  )

preloaded_data <- lapply(preloaded_data, function(df) {
  maxlfq_cols <- grep("MaxLFQ", colnames(df))
  df[rowSums(df[, maxlfq_cols, drop = FALSE] != 0) > 0, ]
})

custom_schema <- list(
  Software1 = list(PEPTIDE = c("output_sequence"), 
                   QUANTITY = c("relative_intensity")), 
  Software2 = list(STRIPPED = c("string"), 
                   RT = c("retention")))

custom_signature <- list(
  Software1 = c("output_sequence"), 
  Software2 = c("retention"))

