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

D170_manual_combined <- read_parquet("C:/Users/Yannic/Desktop/report.parquet")
D170_tumor_only <- read_parquet("C:/Users/Yannic/Desktop/report2.parquet")
D170_fragpipe_combined <- read_parquet("C:/Users/Yannic/Desktop/report3.parquet")

preloaded_data <- list(
  D170_manual_combined = D170_manual_combined,
  D170_tumor_only = D170_tumor_only,
  D170_fragpipe_combined = D170_fragpipe_combined
)

preloaded_data <- list(
  DIANN_parquet = DIANN_parquet,
  NOA01_05_serum  = NOA01_05_serum,
  fragpipe_combined_psm  = read.delim("C:/Users/Yannic/Downloads/combined_psm.tsv"),
  DIANN_report.pr_matrix = read.delim("C:/Users/Yannic/Downloads/report.pr_matrix.tsv"),
  fragpipe_combined_modified_peptide = read.delim("C:/Users/Yannic/Downloads/combined_modified_peptide.tsv"),
  NOA01_01_tumor  = NOA01_01_tumor)

preloaded_data <- read.csv("C:/Users/Yannic/Downloads/test_annotation.csv", sep=";")

custom_schema <- list(
  Software1 = list(PEPTIDE = c("output_sequence"), 
                   QUANTITY = c("relative_intensity")), 
  Software2 = list(STRIPPED = c("string"), 
                   RT = c("retention")))

custom_signature <- list(
  Software1 = c("output_sequence"), 
  Software2 = c("retention"))


invisible(mapply(function(name, source) {
  assign(name, read_parquet(source), envir = .GlobalEnv)
}, annotation_COGNATE_MM$Name, annotation_COGNATE_MM$Source))

data_list <- list(
  NOA05_plasma = NOA01_05_plasma,
  NOA05_serum  = NOA01_05_serum,
  NOA05_tumor  = NOA01_05_tumor,
  NOA01_plasma = NOA01_01_plasma,
  NOA01_serum  = NOA01_01_serum,
  NOA01_tumor  = NOA01_01_tumor,
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

#DIANN new version actually now reports the decoy hits.
data_list <- lapply(data_list, function(df) {
  df %>%
    mutate(type = if_else(startsWith(Protein.Ids, "rev_"), "decoy", "target"))
})

#remove decoy
data_list <- lapply(data_list, function(df) {
  df[df$type == "target",]
})

# Find sequences present in ALL dataframes
shared_sequences <- Reduce(
  intersect,
  lapply(data_list , function(df) df$Stripped.Sequence)
)

# Remove them from every dataframe
data_list_without_universals <- lapply(data_list, function(df) {
  df[!df$Stripped.Sequence %in% shared_sequences, ]})

invisible(mapply(function(name, source) {
  assign(name, read_parquet(source), envir = .GlobalEnv)
}, annotation_COGNATE_MM_MHC2$Name, annotation_COGNATE_MM_MHC2$Source))

data_list2 <- list(
  NOA05_tumor_MHC2  = NOA01_05_tumor_MHC2,
  NOA01_tumor_MHC2  = NOA01_01_tumor_MHC2,
  BT15_tumor_MHC2 = BT15_tumor_MHC2,
  EMB95_tumor_MHC2 = EMB95_tumor_MHC2,
  x1088_tumor_MHC2 = x1088_tumor_MHC2,
  D170_04_tumor_MHC2 = D170_04_tumor_MHC2,
  D170_05_tumor_MHC2 = D170_05_tumor_MHC2,
  D170_15_tumor_MHC2 = D170_15_tumor_MHC2,
  D170_20_tumor_MHC2 = D170_20_tumor_MHC2,
  D170_48_tumor_MHC2 = D170_48_tumor_MHC2
)

data_list2 <- lapply(data_list2, function(df) {
  df %>%
    mutate(type = if_else(startsWith(Protein.Ids, "rev_"), "decoy", "target"))
})
data_list2 <- lapply(data_list2, function(df) {
  df[df$type == "target",]
})

annotation_COGNATE_MM_noCEDAR_5 <- read.delim("C:/PostDoc/Poschke/Data_analysis_MM_noCEDAR_5/annotation_COGNATE_MM_noCEDAR_5.tsv")
preloaded_data <- annotation_COGNATE_MM_noCEDAR_5

invisible(mapply(function(name, source) {
  assign(name, read_parquet(source), envir = .GlobalEnv)
}, annotation_COGNATE_MM_noCEDAR_5$Name, annotation_COGNATE_MM_noCEDAR_5$Source))

data_list <- list(
  NOA05_plasma = NOA01_05_plasma,
  NOA05_serum  = NOA01_05_serum,
  NOA05_tumor  = NOA01_05_tumor,
  NOA01_plasma = NOA01_01_plasma,
  NOA01_serum  = NOA01_01_serum,
  NOA01_tumor  = NOA01_01_tumor,
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
  D170_48_CSF = D170_48_CSF,
  D170_14_plasma = D170_14_plasma,
  D170_14_CSF = D170_14_CSF,
  D170_18_plasma = D170_18_plasma,
  D170_18_CSF = D170_18_CSF,
  JY = JY,
  plasmaQC = plasmaQC
)

#DIANN new version actually now reports the decoy hits.
data_list <- lapply(data_list, function(df) {
  df %>%
    mutate(type = if_else(startsWith(Protein.Ids, "rev_"), "decoy", "target"))
})

#remove decoy
data_list <- lapply(data_list, function(df) {
  df[df$type == "target",]
})

# Find sequences present in ALL dataframes
shared_sequences <- Reduce(
  intersect,
  lapply(data_list , function(df) df$Stripped.Sequence)
)

# Remove them from every dataframe
data_list_without_universals <- lapply(data_list, function(df) {
  df[!df$Stripped.Sequence %in% shared_sequences, ]})
