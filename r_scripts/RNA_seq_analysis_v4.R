install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("impute")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggpubr")
BiocManager::install("pathview")
BiocManager::install("GSEABase")
BiocManager::install("msigdbr")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

library("limma")
library("edgeR")
library("impute")
library(ggplot2)
library(ggrepel)  # For better text label placement
library(dplyr)    # For data manipulation
library(tidyr)    # For pivot_longer
library(ggpubr)
library(clusterProfiler)
library(enrichplot)
library(fgsea)
library(GSEABase)
library(msigdbr)
library(patchwork)
library(rlang)
library(org.Hs.eg.db)

# ------------------------------
# Set a seed for reproducibility
# ------------------------------
set.seed(44)

# ---------------------
# Set working directory
# ---------------------
setwd("/Volumes/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/scripts")


# *************************
# Functions
# *************************

# ------------------------------------------------------
# visualize C1 vs C2 expression of certain gene Function
# ------------------------------------------------------

visualize_C1_vs_C2_PROX1_expression <- function(c2_cleaned_names_value, rna_seq_data_value, pluvicto_master_sheet_value) {
  c2_matches <- grep(
    paste(c2_cleaned_names_value, collapse = "|"), 
    colnames(rna_seq_data_value),
    value = TRUE
  )
  
  c2_matches   # all matching column names
  df_new <- rna_seq_data_value[c2_matches]   # subset dataframe
  only_prox1 <- df_new["PROX1", ]
  
  sort_keys <- sub(".*?(\\d+)_E_C(\\d)", "\\1.\\2", colnames(only_prox1))
  sort_keys_numeric <- as.numeric(sort_keys)
  
  # Get sorted column names
  sorted_cols <- colnames(only_prox1)[order(sort_keys_numeric)]
  
  # Reorder the data frame
  only_prox1_sorted <- t(only_prox1[ , sorted_cols])
  only_prox1_sorted <- as.data.frame(only_prox1_sorted)
  only_prox1_sorted$patient_id <- sub("_E_C[12]$", "", rownames(only_prox1_sorted))
  
  pluvicto_master_sheet_value$patient_id <- rownames(pluvicto_master_sheet_value)
  
  df_merged <- merge(only_prox1_sorted, pluvicto_master_sheet_value, by.x = "patient_id")
  
  rownames(df_merged) <- paste0(df_merged$patient_id, "_", sub(".*_E_C", "E_C", rownames(only_prox1_sorted)))
  df_merged <- df_merged[, c("PROX1", "response")]
  
  return(df_merged)
}


# ---------------------------------------------------------------------------------------
# Scatter Plot of Expression vs NE-10, NEPC Score (ULP and Deep), CCP, and Basal Function
# ---------------------------------------------------------------------------------------

# USED use_all_cycles check at the beginning 

expression_vs_metadata_scatter_plot <- function() {
  list_of_columns = c("NEPC.Score_ULP_C1","NEPC.Fraction_deep_C1","NE10_C1","CCP_C1","Basal_C1")
  prox1_vals <- data.frame(
    patient_id = colnames(rna_seq_data),
    PROX1 = as.numeric(rna_seq_data["PROX1", ])
  )
  
  for (i in list_of_columns) {
    df_temp <- data.frame(
      patient_id = rownames(pluvicto_master_sheet),
      new_column = pluvicto_master_sheet[[i]],
      responder_group = pluvicto_master_sheet$response
    )
    
    df_merged = merge(prox1_vals, df_temp, by="patient_id")
    
    scatter_plot <- ggplot(df_merged, aes(x = PROX1, y = new_column, color = responder_group)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = FALSE, col = "grey") +
      stat_cor(aes(x = PROX1, y = new_column), method = "pearson", inherit.aes = FALSE) +
      scale_color_manual(values = c(
        "non-responder" = "red",
        "responder" = "blue")) +
      labs(x = "PROX1 expression", y = i, color = "Responder group") +
      ggtitle(paste0("Scatter Plot of ", i, " Score vs PROX1 Expression"))
    
    time <- format(Sys.time(), "%Y%m%d_%H%M%S")
    
    ggsave(file.path("./outputs", paste0(i, "_scatter_plot_", time, ".png")), plot = scatter_plot, width = 10, height = 8, dpi = 300)
    
  }
}

# -------------------------------------
# Scatter plot for AR vs PROX1 Function
# -------------------------------------

# DO IT WITH ALL PATIENT 
# MAKE SURE PATIENTS IN BOTH DATAFRAMES HAVE SAME NOMENCLATURE FOR PATIENT ID
AR_vs_PROX1_scatter_plot <- function(rna_seq_data_value, progression_data_value, all_patients) {
  prox1_vals <- data.frame(
    patient_id = colnames(rna_seq_data_value),
    PROX1 = as.numeric(rna_seq_data_value["PROX1", ])
  )

    ar_vals <- data.frame(
    patient_id = colnames(rna_seq_data_value),
    AR = as.numeric(rna_seq_data_value["AR", ])
  )

    df_temp <- data.frame(
    patient_id = rownames(progression_data_value),
    responder_group = progression_data_value$progression_group
  )
  
    df_merged = merge(prox1_vals, ar_vals, by="patient_id")
  if (all_patients) {
    df_merged = merge(df_merged, df_temp, by="patient_id", all.x = TRUE)
    all_patients_title_addition <- "all_patients_"
  } else {
    df_merged = merge(df_merged, df_temp, by="patient_id")
    all_patients_title_addition <- "responders_vs_nonresponders"
  }
  
  scatter_plot <- ggplot(df_merged, aes(x = PROX1, y = AR, color = responder_group)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, col = "grey") +
    stat_cor(aes(x = PROX1, y = AR), method = "pearson", inherit.aes = FALSE) +
    scale_color_manual(values = c(
      "non-responder" = "red",
      "responder" = "blue"
      ), na.value = "grey") +
    labs(x = "PROX1 expression", y = "AR", color = "Responder group") +
    ggtitle(paste0("Scatter Plot of AR vs PROX1 Expression"))
  
  time <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  ggsave(file.path("./outputs", paste0("AR_vs_PROX1_scatter_plot_", all_patients_title_addition, time, ".png")), plot = scatter_plot, width = 10, height = 8, dpi = 300)
  
  return(df_merged)
}

# -----------------------------
# Removing C2 patients Function
# -----------------------------

remove_c2_patients_function <- function(rna_seq_data) {
  # Getting C2 patients
  c2_columns <- grep("C2", colnames(rna_seq_data))
  c2_column_names <- colnames(rna_seq_data)[c2_columns]
  
  # Removing C2 patients
  if(length(c2_columns) > 0) {
    cat("Removing", length(c2_columns), "C2 columns\n")
    print(c2_column_names)
    rna_seq_data <- rna_seq_data[, -c2_columns]
  }
  # Need to return c2_columns eventually
  return(list(rna_seq_data = rna_seq_data, c2_column_names = c2_column_names))
}

# -----------------------------------------------------------
# Paired Box-plot between c1 and c2 prox1 expression Function
# -----------------------------------------------------------

paired_box_plot_c1_vs_c2_function <- function(rna_seq_data_value, c2_column_names_value) {
  c2_column_names_value <- c2_column_names
  all_rows <- rownames(rna_seq_data_value)
  rna_seq_data_value$timepoint <- ifelse(all_rows %in% c2_column_names_value, "C2", "C1")
  rna_seq_data_value$timepoint <- factor(rna_seq_data_value$timepoint, levels = c("C1", "C2"))
  rna_seq_data_value <- rna_seq_data_value %>% mutate(name_stripped = gsub("_[A-Z]_C[0-9]+$", "", rownames(rna_seq_data_value)))
  
  rna_seq_data_value <- rna_seq_data_value %>% mutate(PROX1_AR = (PROX1 / AR))
  
  column_values_of_interest <- c("PROX1", "AR", "PROX1_AR")
  
  for (column in column_values_of_interest) {
    box_plot = ggplot(rna_seq_data_value, aes(x = timepoint, y = .data[[column]], group = name_stripped)) +
      geom_boxplot(aes(group = timepoint), outlier.shape = NA) +
      geom_line(color = "grey", alpha = 0.6) +
      geom_point(aes(color = response)) + 
      scale_color_manual(values = c(
        "non-responder" = "red",
        "responder" = "blue",
        "undefined" = "grey")) + 
      labs(title = paste0("Paired Box Plot (C1 and C2 patients) for ", column),
           x = "Timepoint",
           y = paste0(column, " Expression"),
           color = "Response")
      
    time <- format(Sys.time(), "%Y%m%d_%H%M%S")
    
    ggsave(file.path("./outputs", paste0(column, "_box_plot_C1_vs_C2_", time, ".png")), plot = box_plot, width = 10, height = 8, dpi = 300)
    
  }

  return(rna_seq_data_value)
}

# --------------------------------
# add_genes_to_data_frame_function
# --------------------------------

add_genes_to_data_frame_function <- function(df_merged_value, original_rna_seq_data_value, list_of_genes) {
  df_merged_value$sample_id <- rownames(df_merged_value)
  list_of_genes <- c("AR")
  for (gene in list_of_genes) {
    df_of_AR <- t(original_rna_seq_data_value[gene, ])
    df_of_AR <- data.frame(df_of_AR)
    df_of_AR$sample_id <- rownames(df_of_AR)
    df_merged_value <- merge(df_merged_value, df_of_AR, by = "sample_id")
  }
  rownames(df_merged_value) <- df_merged_value$sample_id
  df_merged_value$sample_id <- NULL
  return(df_merged_value)
}

# ------------------------------------------------------------
# High vs low AR expression bins for PROX1 expression Function
# ------------------------------------------------------------
high_vs_low_ar_for_prox1_function <- function(df_merged_value) {
  df_merged_value <- df_merged
  df_merged_value <- df_merged_value %>% mutate(PROX1_value = ifelse(PROX1 > 1, "high", "low")) 
  df_merged_value$PROX1_value <- factor(df_merged_value$PROX1_value, levels = c("high", "low"))
  
  box_plot = ggplot(df_merged_value, aes(x = PROX1_value, y = AR)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgrey", alpha = 0.05) +
    # geom_line(color = "grey", alpha = 0.6) +
    geom_point(aes(color = responder_group)) +
    scale_color_manual(values = c(
      "non-responder" = "red",
      "responder" = "blue",
      "undefined" = "grey")) +
    labs(title = paste0("Box Plot (high PROX1 vs low PROX1) for AR expression"),
         x = "PROX1 Expression (low = < 1 CPM, high = > 1 CPM)",
         y = "AR Expression",
         color = "Response") +
    stat_compare_means(
      comparisons = list(c("low", "high")), # groups to compare
      method = "t.test",                    # or "wilcox.test"
      label = "p.format"                    # shows nicely formatted p-value
    )

  time <- format(Sys.time(), "%Y%m%d_%H%M%S")

  ggsave(file.path("./outputs", paste0("PROX1_box_plot_high_vs_low_PROX1_", time, ".png")), plot = box_plot, width = 10, height = 8, dpi = 300)

  return(df_merged_value)
}

# ----------------
# Zscore dataframe
# ----------------

zscore_df <- function (df, normalize_across) {

    if (normalize_across == "rows") {
    normalize_across = 1
  } else {
    normalize_across = 2
  }
  df_zscore <- t(apply(df, normalize_across, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
  
  df_zscore <- as.data.frame(df_zscore)
  rownames(df_zscore) <- rownames(df)
  colnames(df_zscore) <- colnames(df)
  
  return(df_zscore)
}

# *************************

# **********************
# LOAD DATA
# **********************

rna_seq_data <- read.table("", 
                         sep="\t", header = TRUE, row.names=1, check.names=FALSE)

progression_data <- read.table("", 
                             sep="\t", header = TRUE, row.names=1, check.names=FALSE)

pluvicto_master_sheet <- read.table("", 
                                    sep=",", header = TRUE, check.names=FALSE)

gene_universal_reference_data <- read.table("",
                                            sep="\t", header = TRUE, row.names=1, check.names=FALSE)

tfbs_data <- read.table("",
                        sep="\t", header = TRUE, check.names = FALSE)

currated_gene_list <-read.csv("")
gene_ids <- currated_gene_list$gene_id

original_rna_seq_data <- rna_seq_data
origional_progression_data <- progression_data
origional_pluvicto_master_sheet <- pluvicto_master_sheet

# ********************
# Preliminary Analysis
# ********************

# -------------------------
# Call AR vs PROX1 function
# -------------------------

# Creates a scatter plot of AR vs PROX1 and box plot of high and low PROX1 expression (removes C2 and cleans patient's names)
call_AR_vs_PROX1_function <- FALSE
all_patients <- TRUE
if (call_AR_vs_PROX1_function) {
  
  # Remove and clean patients
  remove_c2_patients_list <- remove_c2_patients_function(rna_seq_data)
  rna_seq_data_only_c1 <- remove_c2_patients_list$rna_seq_data
  original_names <- colnames(rna_seq_data_only_c1)
  cleaned_names <- gsub("_[A-Z]_C[0-9]+$", "", original_names)
  
  # Update column names
  colnames(rna_seq_data_only_c1) <- cleaned_names
  rna_seq_data_only_c1_zscore <- zscore_df(rna_seq_data_only_c1)
  
  # Calls scatter plot
  df_merged <- AR_vs_PROX1_scatter_plot(rna_seq_data_only_c1_zscore, progression_data, all_patients)
  
  # Calls box plot
  high_vs_low_ar_for_prox1_function_output <- high_vs_low_ar_for_prox1_function(df_merged)
}

# ***************************
# Use TFBS data if wanted
# ***************************
volcano_plot_file_title <- ""

# Sets up Volcano plot's title and file name as well as replaces rna_seq_data with tfbs data
use_tfbs <- FALSE
if (use_tfbs) {
  # Set RNA-seq data to tfbs data (removes first six columns of metadata)
  rna_seq_data <- progression_data[, -(1:6)]
  rna_seq_data_checkpoint <- data.frame(t(rna_seq_data)) # Transpose to tfbs x patient
  volcano_plot_file_title <- paste(volcano_plot_file_title, "TFBS", sep = "")
  volcano_plot_file_title_checkpoint <- volcano_plot_file_title
  
} else {
  rna_seq_data_checkpoint <- rna_seq_data
  volcano_plot_file_title <- paste(volcano_plot_file_title, "genes", sep = "")
  volcano_plot_file_title_checkpoint <- volcano_plot_file_title
}

# *****************************************
# Set responder vs non-responder categories
# *****************************************

columns <- c("survival_days", "PSA_prog_days", "PSA_Progression", "tfx_prog_days", "T_cycles")

# Make rownames ids
rownames(pluvicto_master_sheet) <- pluvicto_master_sheet$Sample_ID

for (column in columns) {
  # For binary indicators 
  if (column == "PSA_Progression") {
    pluvicto_master_sheet <- pluvicto_master_sheet %>% mutate(!!paste0("progression_group_", column) :=ifelse(.data[[column]] == 1 , "non-responder", "responder"))
    next
  }
  
  if (column == "T_cycles") {
    column_temp <- paste(column, "1_vs_6", sep = "_")
    pluvicto_master_sheet <- pluvicto_master_sheet %>% mutate(!!paste0("progression_group_", column_temp) := case_when(
                                                                                            T_cycles == 6 ~ "responder",
                                                                                            T_cycles == 5 ~ NA_character_,
                                                                                            T_cycles == 4 ~ NA_character_,
                                                                                            T_cycles == 3 ~ NA_character_,
                                                                                            T_cycles == 2 ~ NA_character_,
                                                                                            T_cycles == 1 ~ "non-responder",
                                                                                            TRUE ~ NA_character_))
    column_temp <- paste(column, "2_vs_6", sep = "_")
    pluvicto_master_sheet <- pluvicto_master_sheet %>% mutate(!!paste0("progression_group_", column_temp) := case_when(
                                                                                            T_cycles == 6 ~ "responder",
                                                                                            T_cycles == 5 ~ NA_character_,
                                                                                            T_cycles == 4 ~ NA_character_,
                                                                                            T_cycles == 3 ~ NA_character_,
                                                                                            T_cycles == 2 ~ "non-responder",
                                                                                            T_cycles == 1 ~ NA_character_,
                                                                                            TRUE ~ NA_character_))
    column_temp <- paste(column, "1_2_vs_5_6", sep = "_")
    pluvicto_master_sheet <- pluvicto_master_sheet %>% mutate(!!paste0("progression_group_", column_temp) := case_when(
                                                                                            T_cycles == 6 ~ "responder",
                                                                                            T_cycles == 5 ~ "responder",
                                                                                            T_cycles == 4 ~ NA_character_,
                                                                                            T_cycles == 3 ~ NA_character_,
                                                                                            T_cycles == 2 ~ "non-responder",
                                                                                            T_cycles == 1 ~ "non-responder",
                                                                                            TRUE ~ NA_character_))
    column_temp <- paste(column, "1-5_vs_6", sep = "_")
    pluvicto_master_sheet <- pluvicto_master_sheet %>% mutate(!!paste0("progression_group_", column_temp) := case_when(
                                                                                            T_cycles == 6 ~ "responder",
                                                                                            T_cycles == 5 ~ "non-responder",
                                                                                            T_cycles == 4 ~ "non-responder",
                                                                                            T_cycles == 3 ~ "non-responder",
                                                                                            T_cycles == 2 ~ "non-responder",
                                                                                            T_cycles == 1 ~ "non-responder",
                                                                                            TRUE ~ NA_character_))
    next
  }
  
  if (column == "survival_days") {
    pluvicto_master_sheet <- pluvicto_master_sheet %>% mutate(!!paste0("progression_group_", column, "_252_cutoff") := ifelse(.data[[column]] < 252, "non-responder", "responder"))
  }
  
  quantile_df <- pluvicto_master_sheet[!is.na(pluvicto_master_sheet[[column]]), ]
  quantile_vector <- quantile(quantile_df[[column]], probs = c(0.25, 0.5, 0.75))
  print(paste(column, quantile_vector, sep = " "))
  pluvicto_master_sheet <- pluvicto_master_sheet %>% mutate(!!paste0("progression_group_", column, "_median") := ifelse(.data[[column]] < quantile_vector[2], "non-responder", "responder"))
  pluvicto_master_sheet <- pluvicto_master_sheet %>% mutate(!!paste0("progression_group_", column, "_quartile") := case_when(.data[[column]] <= quantile_vector[1] ~ "non-responder",
                                                                                            .data[[column]] >= quantile_vector[3] ~ "responder",
                                                                                            TRUE ~ NA_character_))
}

# **************
# START ANALYSIS
# **************

# List of progression groups
progression_vars <- c(
  "progression_group_survival_days_252_cutoff", # 1
  "progression_group_survival_days_median",     # 2
  "progression_group_survival_days_quartile",   # 3
  "progression_group_PSA_prog_days_median",     # 4
  "progression_group_PSA_prog_days_quartile",   # 5
  "progression_group_PSA_Progression",          # 6
  "progression_group_tfx_prog_days_median",     # 7
  "progression_group_tfx_prog_days_quartile",   # 8
  "progression_group_T_cycles_1_vs_6",          # 9
  "progression_group_T_cycles_2_vs_6",          # 10
  "progression_group_T_cycles_1_2_vs_5_6",      # 11
  "progression_group_T_cycles_1-5_vs_6"         # 12
)

progression_group_id <- progression_vars[9]

volcano_plot_file_title <- volcano_plot_file_title_checkpoint
  
# Set Meta data to progression_data variable 
progression_data <- pluvicto_master_sheet
names(progression_data)[names(progression_data) == progression_group_id] <- "progression_group"
progression_data <- progression_data[!is.na(progression_data$progression_group), ]

# User Needs to Change each time program runs
volcano_plot_file_title <- paste(volcano_plot_file_title, progression_group_id, sep = "_")

cat("non-responders: ", sum(progression_data$progression_group == "non-responder"))
cat("responders: ",sum(progression_data$progression_group == "responder"))

# --------------------------------------------
# Filter out all samples based on TFx
# --------------------------------------------
if (use_tfbs) {
  progression_data <- progression_data %>% filter(TFx_C1 >= 0.03)
  
  volcano_plot_file_title <- paste(volcano_plot_file_title, "TFx_10_perc_cutoff", sep = "_")
} else {
  progression_data <- progression_data %>% filter(TFx_C1 >= 0.10)

  volcano_plot_file_title <- paste(volcano_plot_file_title, "TFx_10_perc_cutoff", sep = "_")
}


cat("non-responders: ", sum(progression_data$progression_group == "non-responder"))
cat("responders: ",sum(progression_data$progression_group == "responder"))

# ---------------------------------------
# Setting data frames equal to each other
# ---------------------------------------

rna_seq_data <- rna_seq_data_checkpoint

cat("Original number of samples:", ncol(rna_seq_data), "\n")
cat("Original column names (first 10):", head(colnames(rna_seq_data), 10), "\n")

# Getting C2 patients and Removing C2 patients
remove_c2_patients_list <- remove_c2_patients_function(rna_seq_data)

rna_seq_data <- remove_c2_patients_list$rna_seq_data
c2_column_names <- remove_c2_patients_list$c2_column_names

# Clean column names by removing suffixes like _E_C1, _S_C1, etc.
# This removes everything after the last underscore that matches the pattern _[A-Z]_C[0-9]
original_names <- colnames(rna_seq_data)
cleaned_names <- gsub("_[A-Z]_C[0-9]+$", "", original_names)
c2_cleaned_names <- gsub("_[A-Z]_C[0-9]+$", "", c2_column_names)

# This might be in the wrong place and might not work (changed from changing data frames to original data frames in hope of helping)
visualize_C1_vs_C2_PROX1_expression_decision <- FALSE
if (visualize_C1_vs_C2_PROX1_expression_decision) {
  df_merged <- visualize_C1_vs_C2_PROX1_expression(c2_cleaned_names, original_rna_seq_data, origional_pluvicto_master_sheet)
  df_merged <- add_genes_to_data_frame_function(df_merged, original_rna_seq_data, c("AR"))
  df_merged <- paired_box_plot_c1_vs_c2_function_output <- paired_box_plot_c1_vs_c2_function(df_merged, c2_column_names)
  df_merged <- high_vs_low_ar_for_prox1_function(df_merged)
}


cat("Sample of name cleaning:\n")
for(i in 1:min(10, length(original_names))) {
  cat(original_names[i], "->", cleaned_names[i], "\n")
}

# Update column names
colnames(rna_seq_data) <- cleaned_names

# Check data dimensions
cat("Expression data dimensions:", dim(rna_seq_data), "\n")
cat("Clinical data dimensions:", dim(progression_data), "\n")
 
clean_colnames <- trimws(colnames(rna_seq_data))
clean_rownames <- trimws(rownames(progression_data))

common_patients <- intersect(clean_colnames, clean_rownames)
different <- setdiff(clean_colnames, clean_rownames)
cat("Number of common patients:", length(common_patients), "\n")
cat("Common patients: ", common_patients, "\n")

if (length(common_patients) == 0){
  stop()
}

# Subset data to common patients
rna_seq_data <- rna_seq_data[, common_patients]
progression_data <- progression_data[common_patients, , drop = FALSE]

# Aligning progression_data to rna_seq_data
progression_data <- progression_data[colnames(rna_seq_data), , drop = FALSE]

# Check data dimensions
cat("Expression data dimensions:", dim(rna_seq_data), "\n")
cat("Clinical data dimensions:", dim(progression_data), "\n")

pre_filtered_rna_seq_data <- rna_seq_data
pre_filtered_progression_data <- progression_data


cat("non-responders: ", sum(progression_data$progression_group == "non-responder"))
cat("responders: ",sum(progression_data$progression_group == "responder"))

# ----------------------
# Filtering RNA-seq or TFBS data
# ----------------------

# Set a check point so you dont need to run above code
rna_seq_data <- pre_filtered_rna_seq_data

# if data is tfbs -> normalize, convert to fold change, and find and correct pvalue
if (use_tfbs) {
  
  # Create group factor
  progression_factor <- factor(progression_data[["progression_group"]], levels = c("non-responder" , "responder"))
  tumor_fraction <- progression_data$TFx_C1
  
  tfbs_regression_estimate_df <- data.frame(
                                  tfbs = rownames(rna_seq_data),
                                  estimate = NA_real_,
                                  pval = NA_real_,
                                  pval_adj = NA_real_
                                )
  
  # Regress out TFx
  for (i in seq_len(nrow(rna_seq_data))) {
    tfbs <- rownames(rna_seq_data)[i]
    
    temp_df <- data.frame(
                y = progression_factor,
                tfbs = as.numeric(rna_seq_data[tfbs, ]),
                tumor_fraction = tumor_fraction
              )
    
    model <- glm(y ~ tfbs + tumor_fraction, data = temp_df, family = binomial)
    coef_summary <- summary(model)$coefficients
    
    tfbs_regression_estimate_df$estimate[i] <- coef_summary["tfbs", "Estimate"]
    tfbs_regression_estimate_df$pval[i] <- coef_summary["tfbs", "Pr(>|z|)"]
  }
  
  tfbs_regression_estimate_df <- tfbs_regression_estimate_df %>% arrange(pval)
  
  print(tfbs_regression_estimate_df)
  
  tfbs_regression_estimate_df$pval_adj <- p.adjust(tfbs_regression_estimate_df$pval, method = "fdr")
  
  tfbs_regression_estimate_df_filtered <- tfbs_regression_estimate_df %>% filter(pval_adj <= 0.2)
  
  if (nrow(tfbs_regression_estimate_df_filtered) < 1) {
    print(paste("No tfbs's after regressing out TFx from ", progression_group_id))
    next
  } else {
    print(paste(progression_group_id, " has tfbs's"))
  }
  
  rna_seq_data_logFC_pvalue_df <- data.frame(gene = rownames(rna_seq_data),
                                        logFC = NA, pvalue = NA)
  
  contrast_position <- "Responder_vs_Nonresponder"
  
  for (i in seq_len(nrow(rna_seq_data))) {
    vals <- as.numeric(rna_seq_data[i, ])
    g1 <- vals[progression_factor == levels(progression_factor)[1]]
    g2 <- vals[progression_factor == levels(progression_factor)[2]]
    
    if (contrast_position == "Responder_vs_Nonresponder") {
      eps <- (1 * (10^-7))
      rna_seq_data_logFC_pvalue_df$logFC[i] <- (mean(g2) + eps) / (mean(g1) + eps)
    } else {
      rna_seq_data_logFC_pvalue_df$logFC[i] <- (mean(g1) + eps) / (mean(g2) + eps)
    }
    
    if (all(!is.na(g2)) | all(!is.na(g1))) {
      mann_u <- wilcox.test(g2, g1)
      rna_seq_data_logFC_pvalue_df$pvalue[i] <- mann_u$p.value
    }
  }
  
  rna_seq_data_logFC_pvalue_df$adj.P.Val <- p.adjust(rna_seq_data_logFC_pvalue_df$pvalue, method = "BH")
  
  rna_seq_data_logFC_pvalue_df <- na.omit(rna_seq_data_logFC_pvalue_df)
  
  deg_results <- rna_seq_data_logFC_pvalue_df

} else {
  if ("PROX1" %in% row.names(rna_seq_data)) {
    print("PROX1 in rna_seq_data")
  } else {
    print("PROX1 not in rna_seq_data")
  }
  
  # Filter lowly expressed genes (optional but recommended)
  # Keep genes with CPM > 1 in at least half of the samples
  min_samples <- ceiling(ncol(rna_seq_data) * .1)
  keep <- rowSums(rna_seq_data > 1) >= min_samples
  cat("Genes before filtering:", nrow(rna_seq_data), "\n")
  rna_seq_data <- rna_seq_data[keep, ]
  cat("Genes after filtering:", nrow(rna_seq_data), "\n")
  
  if ("PROX1" %in% row.names(rna_seq_data)) {
    print("PROX1 in rna_seq_data")
  } else {
    print("PROX1 not in rna_seq_data")
  }
  
  # Remove genes with >= 50% NaN values
  nan_pct <- rowMeans(is.na(rna_seq_data))
  genes_to_keep <- nan_pct < 0.5
  genes_to_remove <- nan_pct > 0.5
  cat("Removing genes: ", genes_to_remove, "\n")
  cat("Before gene filtering:", nrow(rna_seq_data), "genes\n")
  rna_seq_data <- rna_seq_data[genes_to_keep, , drop = FALSE]
  cat("After removing genes with ≥50% NaN:", nrow(rna_seq_data), "genes\n")
  
  if ("PROX1" %in% row.names(rna_seq_data)) {
    print("PROX1 in rna_seq_data")
  } else {
    print("PROX1 not in rna_seq_data")
  }
  
  # Variance filtering
  var <- sort(apply(rna_seq_data, 1, var), decreasing = TRUE)
  var <- names(var)[1:floor(length(var)*0.20)]
  cat("Before gene filtering:", nrow(rna_seq_data), "genes\n")
  rna_seq_data <- rna_seq_data[var, ]
  cat("After removing genes with ≥50% NaN:", nrow(rna_seq_data), "genes\n")
  
  if ("PROX1" %in% row.names(rna_seq_data)) {
    print("PROX1 in rna_seq_data")
  } else {
    print("PROX1 not in rna_seq_data")
  }
  
  # ----------------------------------------
  # Control for TFx on a gene by gene basis
  # ----------------------------------------
  
  corr_list <- list()
  for (row in rownames(rna_seq_data)) {
    gene <- unlist(rna_seq_data[row, ])
    tfx <- progression_data$TFx_C1
    spearman_corr <- cor.test(x = gene, y = tfx, method = "spearman")
    pearson_corr <- cor.test(x = gene, y = tfx, method = "pearson")
    corr_list[[row]] <- list(spearman = spearman_corr, pearson = pearson_corr)
    
    temp_df <- data.frame(gene = gene, tfx = tfx)
  }
  
  corr_list <- corr_list[order(sapply(corr_list, function(x) x$spearman$p.value))]
  head(corr_list)
  
  temp_df <- bind_rows(lapply(names(corr_list), function(name) {
    spearman <- corr_list[[name]]$spearman
    data.frame(gene = name, spearman_rho = spearman$estimate)
  }))
  
  ggplot(temp_df, aes(x = spearman_rho)) +
    geom_density(alpha = 0.5)
  
  # ----------------------
  # Creating design
  # ----------------------
  
  new_df <- progression_data %>% filter(is.na(TFx_C1))

  setdiff(rownames(progression_data), colnames(rna_seq_data))
  dim(progression_data)
  dim(rna_seq_data)
  
  # Creating a factor group
  progression_data$design_group <- factor(progression_data$progression_group)
  
  design <- model.matrix(~0 + design_group + TFx_C1, data = progression_data)
  colnames(design) <- make.names(c(levels(progression_data$design_group), "TFx_C1"))
  
  # ----------------------------------------
  # Perform differential expression analysis
  # ----------------------------------------
  
  # Expression data is log2(CPM+1), suitable for limma-trend
  
  fit_limma <- lmFit(rna_seq_data, design)
  
  contrast_position <- "Responder_vs_Nonresponder"
  contrast_matrix <- makeContrasts(Responder_vs_Nonresponder = responder - non.responder, levels = design)
  
  fit_contrast <- contrasts.fit(fit_limma, contrast_matrix)
  
  fit_ebayes <- eBayes(fit_contrast, trend = TRUE)
  
  deg_results_tt <- topTable(fit_ebayes, coef = "Responder_vs_Nonresponder", n = Inf, adjust.method = "BH", sort.by = "P") 
  deg_results <- as.data.frame(deg_results_tt)
  deg_results$gene <- rownames(deg_results)

}

# ***********************
# Save DEG results as tsv
# ***********************

# Edit deg_results
# save_table <- deg_results
# rownames(save_table) <- NULL
# 
# file_name <- paste0("./outputs/data-tables/", volcano_plot_file_title, "deg_results.tsv")
# 
# write.table(save_table, 
#             file_name, 
#             sep = "\t", 
#             row.names = FALSE, 
#             col.names = !file.exists(file_name),
#             quote = FALSE,
#             append = file.exists(file_name)
#             )



# *********************
# Volcano Plot
# *********************

if (use_tfbs) {
  cutoff <- 1
} else {
  cutoff <- 0
}
deg_results$significant <- ifelse(
  deg_results$adj.P.Val < 0.05 & !is.na(deg_results$adj.P.Val),
  ifelse(
    deg_results$logFC > cutoff, "Significantly up-regulated",
    ifelse(deg_results$logFC < cutoff,
           "Significantly down-regulated",
           "Significantly neutral"
    )
  ),
  "Not Significant"
)

sig_genes_volcano <- deg_results %>% dplyr::filter(adj.P.Val < 0.05 & !is.na(adj.P.Val))

# total_gene_dict[[volcano_plot_file_title]] <- save_table$gene

sig_genes_up_regulated <- sig_genes_volcano %>% dplyr::filter(logFC > cutoff)
sig_genes_down_regulated <- sig_genes_volcano %>% dplyr::filter(logFC < cutoff)

high_adj_p_val_genes <- head(sig_genes_volcano)

sig_genes_up_regulated$gene
sig_genes_down_regulated$gene

# Select top 15 up-regulated and top 15 down-regulated significant genes by fold change for labeling
top_up_regulated_genes_fc <- sig_genes_volcano %>%
  dplyr::filter(logFC > cutoff) %>%
  dplyr::arrange(dplyr::desc(logFC))

top_up_regulated_genes_fc <- head(top_up_regulated_genes_fc, n = 10)

top_down_regulated_genes_fc <- sig_genes_volcano %>%
  dplyr::filter(logFC < cutoff) %>%
  dplyr::arrange(logFC) 

top_down_regulated_genes_fc <- head(top_down_regulated_genes_fc, n = 10)

top_genes_volcano <- unique(rbind(top_up_regulated_genes_fc, top_down_regulated_genes_fc))

prox1_gene_list <- c("PROX1")

volcano_plot <- ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
  geom_point(alpha = 0.6, size = 2.25) + # Increased from 1.5 to 2.25 (50% increase)
  scale_color_manual(values = c(
    "Significantly down-regulated" = "red",
    "Significantly neutral" = "black",
    "Not Significant" = "gray",
    "Significantly up-regulated" = "blue")) +
  theme_bw(base_size = 12) +
  # labs(title = paste(paste0("Volcano Plot: ", contrast_position, " | ", volcano_plot_file_title)),
  #      x = "Log2 Fold Change (Negative = Down-regulated, Positive = Up-regulated)",
  #      y = "-Log10 Adjusted P-value", color = "significant") +
  labs(title = paste(paste0("Volcano Plot: Pluvicto_DEG (Completed Trial vs Single Cycle)")),
       subtitle = "Controlled for TFx",
       x = "Log2 Fold Change (Negative = Down-regulated, Positive = Up-regulated)",
       y = "-Log10 Adjusted P-value", color = "significant") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  # geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  # geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  ggrepel::geom_text_repel(data = top_genes_volcano, aes(label = gene),
                           size = 5.25, # Increased from 3.5 to 5.25 (50% increase)
                           fontface = "plain", # Added bold font
                           box.padding = 0.5,
                           max.overlaps = Inf,
                           segment.color = "black",
                           min.segment.length = 0.1) +
  ggrepel::geom_text_repel(data = deg_results %>% dplyr::filter(gene %in% gene_ids),
                           aes(label = gene),
                           color = "gold",
                           size = 5.25,
                           fontface = "bold",
                           box.padding = 0.5,
                           max.overlaps = Inf,
                           segment.color = "gold",
                           min.segment.length = 0.1) +
  ggrepel::geom_text_repel(data = deg_results %>% dplyr::filter(gene %in% prox1_gene_list),
                           aes(label = gene),
                           color = "green",
                           size = 5.25,
                           fontface = "bold",
                           box.padding = 0.5,
                           max.overlaps = Inf,
                           segment.color = "black",
                           min.segment.length = 0.1)

volcano_plot

time <- format(Sys.time(), "%Y%m%d_%H%M%S")

volcano_plot_file_title <- paste(volcano_plot_file_title, "ctrl_for_TFx", sep = "_")

ggsave(file.path("./outputs", paste0(volcano_plot_file_title, "_DGEA_volcano_plot_", time, ".png")), plot = volcano_plot, width = 10, height = 8, dpi = 300)

}

# -----------------------------
# Density Plot of certain genes
# -----------------------------
gene_of_intereset <- "PROX1"

expression_list <- as.list(pre_filtered_rna_seq_data[gene_of_intereset, ])
numeric_values <- unlist(expression_list)
numeric_values <- as.numeric(numeric_values)

# Calculate density
density_estimate <- density(numeric_values)

# Plot density curve
plot(density_estimate, main = "Density Plot of Data", xlab = "Value", ylab = "Density")

# Convert data to data frame
denstiy_plot_df <- data.frame(expression = numeric_values)
group <- ifelse(design[,"responder"] == 1, "responder", "non.responder")
denstiy_plot_df$group <- group

# For one line
density_plot <- ggplot(denstiy_plot_df, aes(x = expression, fill = group)) +
                        geom_density(fill = "steelblue", alpha = 0.5) +
                        geom_rug(aes(color = group), sides = "b") +  # adds tick marks for each patient value
                        labs(title = paste0("Density of expression values for ", gene_of_intereset," across patients"),
                             x = "Expression", y = "Density")

# Split lines into responder vs non responder
density_plot <- ggplot(denstiy_plot_df, aes(x = expression, fill = group)) +
                        geom_density(alpha = 0.4) +
                        geom_rug(aes(color = group), sides = "b") +
                        labs(
                          title = paste0("Density of expression values for ", gene_of_intereset, " across patients"),
                          x = "Expression",
                          y = "Density"
                        )
  
time <- format(Sys.time(), "%Y%m%d_%H%M%S")

ggsave(file.path("./outputs", paste0(gene_of_intereset, "_density_plot_", time, ".png")), plot = density_plot, width = 10, height = 8, dpi = 300)


# -------------------------
# GSEA on significant genes
# -------------------------

# Viewing gene lists
sig_genes_up_regulated 
sig_genes_down_regulated 

create_gsea_list_entrezid <- function(df) {
  # Filtering down the original gene list to gene symbol and logFC
  if (use_all_features) {
    original_gene_list <- df$logFCxtvalue
  } else {
    original_gene_list <- df$logFC
  }
  
  names(original_gene_list) <- df$gene
  
  # Omit NAs and sort by logFC
  filtered_gene_list <- na.omit(original_gene_list)
  filtered_gene_list = sort(filtered_gene_list, decreasing = TRUE)
  
  # See how many genes were omitted
  difs <- length(filtered_gene_list) - length(original_gene_list)
  cat(abs(difs), " genes were omitted due to NAs")
  
  # Convert gene symbols with entrezid
  gene_symbols <- names(filtered_gene_list)
  entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Replace original gene list symbols with entrezids
  gene_list_entrez <- filtered_gene_list[entrez_ids$SYMBOL]
  names(gene_list_entrez) <- entrez_ids$ENTREZID
  
  # Remove any NA values
  gene_list_entrez <- gene_list_entrez[!is.na(names(gene_list_entrez))]
  
  # See difference between original gene list and gene list with entrezid
  difs <- length(filtered_gene_list) - length(gene_list_entrez)
  cat(abs(difs), " genes were lost in converting to EntrezID")
  
  return(gene_list_entrez)
}

get_ref_pathways <- function(category = "H", subcategory = NULL) {
  # Available categories:
  # H = Hallmark gene sets
  # C1 = Positional gene sets  
  # C2 = Curated gene sets (includes KEGG, REACTOME, BIOCARTA)
  # C3 = Regulatory target gene sets
  # C4 = Computational gene sets
  # C5 = Ontology gene sets (GO)
  # C6 = Oncogenic signature gene sets
  # C7 = Immunologic signature gene sets
  # C8 = Cell type signature gene sets
  
  if (!is.null(subcategory)) {
    # For C2, you can specify subcategories like "KEGG", "REACTOME", "BIOCARTA"
    gene_sets <- msigdbr(species = "Homo sapiens", 
                         category = category, 
                         subcategory = subcategory)
  } else {
    gene_sets <- msigdbr(species = "Homo sapiens", category = category)
  }
  
  # Convert to list format for fgsea
  pathways <- split(gene_sets$entrez_gene, gene_sets$gs_name)
  
  cat("=== Gene Set Information ===\n")
  cat("Category:", category, "\n")
  if (!is.null(subcategory)) cat("Subcategory:", subcategory, "\n")
  cat("Number of gene sets:", length(pathways), "\n")
  cat("Gene set size range:", range(lengths(pathways)), "\n")
  
  return(pathways)
}

run_gsea_go <- function(gene_list, pathways)  {
  gsea_multilevel <- fgseaMultilevel(
                      pathways = pathways,
                      stats = gene_list,
                      minSize = 3,
                      maxSize = 600
                    )
  
  return(gsea_multilevel)
}

plot_fgsea_results <- function(fgsea_results, pathway, pvalue_cutoff) {
  
  # Get top pathways
  top_pathways <- fgsea_results %>%
    filter(padj < pvalue_cutoff)

  # Dot plot
  p1 <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway, NES))) +
    geom_point(aes(size = size, color = padj)) +
    scale_color_gradient(low = "red", high = "blue", name = "Adj. P-value") +
    scale_size_continuous(name = "Gene Set Size") +
    theme_minimal() +
    labs(title = paste0("Pathway Enrichment Analysis: ", pathway),
      subtitle = paste0("Top enriched pathways in responders (adj. p < ", pvalue_cutoff, ")"),
         x = "Normalized Enrichment Score (NES)",
         y = NULL) +
    theme(
          axis.text.y = element_text(size = 8, face = "bold"),
          axis.text.x = element_text(size = 8, face = "bold")
          )

  print(p1)
  
  return(p1)
}

# Set up params for GSEA

pathway_list <- list(
  H = list(category = "H", subcategory = NULL),
  C2_REACTOME = list(category = "C2", subcategory = "REACTOME"),
  C2_KEGG = list(category = "C2", subcategory = "KEGG_MEDICUS"),
  C2_PID = list(category = "C2", subcategory = "PID"),
  C5_BP = list(category = "C5", subcategory = "GO:BP"),
  C5_MF = list(category = "C5", subcategory = "GO:MF"),
  C5_HPO = list(category = "C5", subcategory = "HPO"),
  C3_TFT = list(category = "C3", subcategory = "TFT:GTRD"),
  C6_ONCO = list(category = "C6", subcategory = NULL)
)

pvalue_cutoff <- 0.10

use_all_features <- TRUE
if (use_all_features) {
  deg_results <- deg_results %>% mutate(logFCxtvalue = logFC * t)
  all_features_list <- create_gsea_list_entrezid(deg_results)
  sig_regulated_gene_list <- list(all_features = all_features_list)
} else {
  up_regulated_gene_list <- create_gsea_list_entrezid(sig_genes_up_regulated)
  down_regulated_gene_list <- create_gsea_list_entrezid(sig_genes_down_regulated)
  
  up_regulated_gene_list <- head(up_regulated_gene_list, n = 200)
  down_regulated_gene_list <- head(down_regulated_gene_list, n = 200)
  
  sig_regulated_gene_list <- list(up_reg = up_regulated_gene_list, down_reg = down_regulated_gene_list)
  
}

# Sort lists
for (list in sig_regulated_gene_list) {
  sig_regulated_gene_list$list <- sort(sig_regulated_gene_list$list)
}

all_plots <- list()
parent_dir <- "./outputs"

for (gene_list in names(sig_regulated_gene_list)) {
  output_dir <- file.path(parent_dir, gene_list, volcano_plot_file_title)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  for (pathway in names(pathway_list)) {
    ref_pathways <- get_ref_pathways(pathway_list[[pathway]]$category, pathway_list[[pathway]]$subcategory)
    gsea_multilevel <- run_gsea_go(sig_regulated_gene_list[[gene_list]], ref_pathways)
    head(gsea_multilevel)
    plot_fgsea_result <- plot_fgsea_results(gsea_multilevel, pathway, pvalue_cutoff)
    all_plots[[pathway]] <- plot_fgsea_result
    
    time <- format(Sys.time(), "%Y%m%d_%H%M%S")
    file_name <- paste0(gene_list, "_in_resp_", pathway, "_", pvalue_cutoff, "_gsea_dot_plot_", time, ".png")
    ggsave(file.path(output_dir, paste0(volcano_plot_file_title, file_name)), plot = plot_fgsea_result, width = 10, height = 8, dpi = 300)
  }
  
  combined_plot <- wrap_plots(all_plots, ncol = 3) + 
    plot_annotation(title = paste0("GSEA Results - All Collections - ", gene_list, " in responders"),
                    theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
  
  time <- format(Sys.time(), "%Y%m%d_%H%M%S")
  file_name <- paste0(gene_list, "_in_resp_all_gsea_dot_plots_", time, ".png")
  ggsave(file.path(parent_dir, paste0(volcano_plot_file_title, file_name)), plot = combined_plot, width = 25, height = 18, dpi = 300)
  
}


# *******************************
# Venn diagram of DEG data tables
# *******************************
install.packages("ggVennDiagram")
library(ggVennDiagram)
names(total_gene_dict)
venn_diagram <- ggVennDiagram(total_gene_dict[c(1, 4, 6)])

time <- format(Sys.time(), "%Y%m%d_%H%M%S")
file_name <- paste0("venn_diagram", time, ".png")
ggsave(file.path("./outputs", paste0(volcano_plot_file_title, file_name)), plot = venn_diagram, width = 25, height = 18, dpi = 300)











