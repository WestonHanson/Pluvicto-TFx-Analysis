install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("impute")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggpubr")
install.packages("reshape2")
install.packages("conflicted")

library("limma")
library("edgeR")
library("impute")
library(ggplot2)
library(ggrepel)  # For better text label placement
library(dplyr)    # For data manipulation
library(tidyr)    # For pivot_longer
library(ggpubr)
library(reshape2)
library(conflicted)

# ------------------------------
# Set a seed for reporducibility
# ------------------------------
set.seed(44)

# ---------------------
# Set working directory
# ---------------------
setwd("/Volumes/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/scripts")


# ---------------------------------------------
# Set source (to use functions from other file)
# ---------------------------------------------
# source("/Volumes/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/scripts/RNA_seq_analysis_v3.R")

# *************************
# Functions
# *************************

# -----------------------
# convert to CPM function
# -----------------------
convert_to_CPM_function <- function (df) {
  for (column in colnames(df)) {
    column_sum <- sum(df[[column]], na.rm = TRUE)
    df[[column]] <- (log2(((df[[column]] / column_sum) * 1e6)+ 1) ) 
  }
  return(df)
}

# --------------------------------
# Create a new data frame function
# --------------------------------
create_new_data_frame <- function (df, rows, columns, column_of_intrest = NULL) {
  matched_cols <- colnames(df)[Reduce(`|`, lapply(columns, function(p) grepl(p, colnames(df))))]
  if (!is.null(column_of_intrest)) {
    cols_to_keep <- unique(c(column_of_intrest, matched_cols))
    matched_rows <- df[[column_of_intrest]] %in% rows
  } else {
    cols_to_keep <- matched_cols
    matched_rows <- rownames(df) %in% rows
  }
  
  new_df <- df[matched_rows, cols_to_keep, drop = FALSE]
  return(new_df)
}

# --------------------------------
# Convert NCBI gene to gene symbol
# --------------------------------
ncbi_to_gene_symbol <- function (gene_reference, rna_seq_df) {
  gene_reference$ncbi_number <- gsub("^.*?:", "", rownames(gene_reference))
  rna_seq_df$ncbi_number <- rownames(rna_seq_df)
  
  colnames(gene_reference)
  colnames(rna_seq_df)
  
  merged_df <- rna_seq_df %>% 
    left_join (
      gene_reference %>% dplyr::select(ncbi_number, Gene_symbol),
      by = "ncbi_number"
    )
  
  head(merged_df$Gene_symbol)
  
  return(merged_df)
}

# ------------------
# Combine Dataframes
# ------------------

combine_data_frames <- function(df_1, df_2) {
  combined_rows <- base::union(rownames(df_1), rownames(df_2))
  
  df_1_aligned <- df_1[combined_rows, , drop = TRUE]
  rownames(df_1_aligned) <- combined_rows
  
  missing_rows <- base::setdiff(rownames(df_2), rownames(df_1))
  for (r in missing_rows) {
    df_1_aligned[r, ] <- rep(df_2[r, , drop = TRUE], ncol(df_1))
  }
  
  return(df_1_aligned)
}

# ---------------------------------------------------
# Create dictionary of genes and metadata for heatmap
# ---------------------------------------------------
create_df_for_heatmap <- function(
    TAN_predicted_rna_seq_data, 
    pluvicto_predicted_rna_seq_data,
    TAN_predicted_AR_10,
    pluvicto_predicted_AR_10,
    patinets_of_interest, 
    genes_of_interest, 
    metadata_of_interest
  ) {
  heatmap_df <- data.frame(gene = genes_of_interest, stringsAsFactors = FALSE)
  
  for (patient_name in names(patinets_of_interest)) {
    patient <- patinets_of_interest[[patient_name]]
    
    # For Pluvicto name
    pluvicto_id <- patient$pluvicto
    filtered_cols <- pluvicto_predicted_rna_seq_data[, grepl(pluvicto_id, colnames(pluvicto_predicted_rna_seq_data))]
    filtered_rows <- filtered_cols[rownames(filtered_cols) %in% genes_of_interest, ]
    filtered_rows$gene <- rownames(filtered_rows)
    
    heatmap_df <- merge(heatmap_df, filtered_rows, by = "gene", all.x = TRUE)

    #For TAN name
    TAN_id <- patient$TAN
    filtered_cols <- TAN_predicted_rna_seq_data[, grepl(TAN_id, colnames(TAN_predicted_rna_seq_data)), drop = FALSE]
    filtered_rows <- filtered_cols[rownames(filtered_cols) %in% genes_of_interest, , drop = FALSE]
    filtered_rows$gene <- rownames(filtered_rows)

    heatmap_df <- merge(heatmap_df, filtered_rows, by = "gene", all.x = TRUE)
    
    # For Pluvicto AR-10
    filtered_cols <- pluvicto_predicted_AR_10[, grepl(pluvicto_id, colnames(pluvicto_predicted_AR_10)), drop = FALSE]
    filtered_rows <- filtered_cols[rownames(filtered_cols) %in% metadata_of_interest, , drop = FALSE]
    filtered_rows$gene <- rownames(filtered_rows)
    
    heatmap_df_t <- data.frame(t(heatmap_df), colname = rownames(t(heatmap_df)), row.names = NULL)
    filtered_rows_t <- data.frame(t(filtered_rows), colname = rownames(t(filtered_rows)), row.names = NULL)
    merge_temp <- merge(heatmap_df_t, filtered_rows_t, by = "colname", all = TRUE)
    rownames(merge_temp) <- merge_temp$colname
    merge_temp$colname <- NULL
    heatmap_df <- as.data.frame(t(merge_temp))
    rownames(heatmap_df) <- NULL
    
    # For TAN AR-10
    filtered_cols <- TAN_predicted_AR_10[, grepl(TAN_id, colnames(TAN_predicted_AR_10)), drop = FALSE]
    filtered_rows <- filtered_cols[rownames(filtered_cols) %in% metadata_of_interest, , drop = FALSE]
    filtered_rows$gene <- rownames(filtered_rows)
    
    heatmap_df[nrow(heatmap_df)-1, paste0("PC-TAN_", TAN_id)] <- filtered_rows["ARG.10", paste0("PC-TAN_", TAN_id)]
    heatmap_df[nrow(heatmap_df), paste0("PC-TAN_", TAN_id)] <- filtered_rows["NE.10", paste0("PC-TAN_", TAN_id)]
    
    rownames(heatmap_df) <- heatmap_df$gene
    heatmap_df$gene <- NULL
  }
  
  return(heatmap_df)
}

# ----------------
# Zscore dataframe
# ----------------

zscore_df <- function (df) {
  df_zscore <- t(apply(df, 1, function(x) {
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

TAN_raw_rna_seq_data <-read.table("../data-files/UW_TAN_CRPC_208_hg38_raw_counts.txt", 
                                  sep="\t", header = TRUE, row.names=1, check.names=FALSE)

TAN_true_AR_10 <- read.table("../data-files/TAN_RNA_seq_truth_gsva.txt", 
                                       header = TRUE, row.names=1, check.names=FALSE)

TAN_predicted_rna_seq_data <- read.table("../data-files/predictions_tumor_TAN.tsv", 
                               sep="\t", header = TRUE, row.names=1, check.names=FALSE)

pluvicto_predicted_rna_seq_data <- read.table("../data-files/predictions_tumor_pluvicto.tsv", 
                           sep="\t", header = TRUE, row.names=1, check.names=FALSE)

pluvicto_predicted_AR_10 <- read.table("../data-files/Pluvicto_predicted_gsva.txt", 
                                              header = TRUE, row.names=1, check.names=FALSE)

TAN_predicted_AR_10 <- read.table("../data-files/PC-TAN_predicted_gsva.txt", 
                                       header = TRUE, row.names=1, check.names=FALSE)

gene_universal_reference_data <- read.table("../data-files/MANE.GRCh38.v1.3.UniversalReference_2.tsv", 
                                         sep="\t", header = TRUE, row.names=1, check.names=FALSE)

pluvicto_master_sheet <- read.table("../data-files/pluvicto_C1_master.csv", 
                                    sep=",", header = TRUE, row.names=1, check.names=FALSE)

# **************
# START ANALYSIS
# **************

# ------------------
# Select patient ids
# ------------------

pluvicto_id <- "patient_id"
TAN_id <- "patient_id"

# -----------------------------------------------------------------------
# Converting raw RNA and comparing genes and patients to predicted values
# -----------------------------------------------------------------------

# Convert TAN RNA-seq to log(CPM + 1)
TAN_raw_rna_seq_data_cleaned <- na.omit(TAN_raw_rna_seq_data)
TAN_raw_rna_seq_data_converted <- convert_to_CPM_function(TAN_raw_rna_seq_data_cleaned)

# Zscore TAN RNA-seq data
TAN_raw_rna_seq_data_zscore <- zscore_df(TAN_raw_rna_seq_data_converted)

# Test to see if zscore worked
mean_column <- colMeans(TAN_raw_rna_seq_data_zscore, na.rm = TRUE)
std_column <- apply(TAN_raw_rna_seq_data_zscore, 2, sd, na.rm = TRUE)
hist(mean_column, main = "Column Means", xlab = "Mean")
hist(std_column, main = "Column SDs", xlab = "SD")

# Convert RNA-seq gene ids to common name
TAN_raw_rna_seq_data_converted_new_gene_ids <- ncbi_to_gene_symbol(gene_universal_reference_data, TAN_raw_rna_seq_data_zscore)

# Filter data down to only specific genes and patients 
rows_to_convert <- c("PROX1", "AR", "PDPN", "LYVE1", "CD34", "hnRNPK", "HIF1A", "HIF1AN", "MMP14", "PPARG")
columns_to_convert <- c(TAN_id)
column_of_intrest <- "Gene_symbol"
TAN_raw_rna_seq_data_converted_filtered <- create_new_data_frame(TAN_raw_rna_seq_data_converted_new_gene_ids, rows_to_convert, columns_to_convert, column_of_intrest)

# Zscore TAN true AR.10 data
TAN_true_AR_10_t <- scale(t(TAN_true_AR_10))
TAN_true_AR_10_zscore <- as.data.frame(t(TAN_true_AR_10_t))

# Test to see if zscore worked
mean_column <- colMeans(TAN_true_AR_10_zscore, na.rm = TRUE)
std_column <- apply(TAN_true_AR_10_zscore, 2, sd, na.rm = TRUE)
hist(mean_column, main = "Column Means", xlab = "Mean")
hist(std_column, main = "Column SDs", xlab = "SD")


# Collect NE.10 and ARG.10
rows_to_convert <- c("ARG.10", "NE.10")
columns_to_convert <- c(TAN_id)
TAN_true_AR_10_converted_filtered <- create_new_data_frame(TAN_true_AR_10_zscore, rows_to_convert, columns_to_convert)
colnames(TAN_true_AR_10_converted_filtered) <- paste(colnames(TAN_true_AR_10_converted_filtered), "_avg", sep ="")


# Zscore TAN predicted RNA seq data
TAN_predicted_rna_seq_data_t <- scale(t(TAN_predicted_rna_seq_data))
TAN_predicted_rna_seq_data_zscore <- as.data.frame(t(TAN_predicted_rna_seq_data_t))

# Test to see if zscore worked
mean_column <- colMeans(TAN_predicted_rna_seq_data_zscore, na.rm = TRUE)
std_column <- apply(TAN_predicted_rna_seq_data_zscore, 2, sd, na.rm = TRUE)
hist(mean_column, main = "Column Means", xlab = "Mean")
hist(std_column, main = "Column SDs", xlab = "SD")

# Filter predicted values to only specific genes and patients
rows_to_convert <- c("PROX1", "AR")
columns_to_convert <- c(TAN_id)
TAN_predicted_rna_seq_data_zscore$Gene_symbol <- rownames(TAN_predicted_rna_seq_data_zscore)
TAN_predicted_rna_seq_data_filtered <- create_new_data_frame(TAN_predicted_rna_seq_data_zscore, rows_to_convert, columns_to_convert)

# -------------------------------------------------------------------------------------
# Making heatmap for certain patients and certain genes and metadata across time points
# -------------------------------------------------------------------------------------

patients_of_interest <- list(patient_1 = list(pluvicto = pluvicto_id, TAN = TAN_id))
down_regulated <- c("CD34", "MMP14", "PRRX1", "PPARG")
up_regulated <- c("PDPN", "LYVE1")
interacting <- c("hnRNPK", "HIF1A", "HIF1AN")
main <- c("PROX1", "AR")
genes_of_interest <- c(down_regulated, up_regulated, interacting, main)
metadata_of_interest <- c("ARG.10", "NE.10")

# Zscore all dataframes
TAN_predicted_rna_seq_data_filtered_zscore <- zscore_df(TAN_predicted_rna_seq_data)
pluvicto_predicted_rna_seq_data_zscore <- zscore_df(pluvicto_predicted_rna_seq_data)
TAN_predicted_AR_10_zscore <- zscore_df(TAN_predicted_AR_10)
pluvicto_predicted_AR_10_zscore <- zscore_df(pluvicto_predicted_AR_10)

combined_pluvicto_df <- create_df_for_heatmap(
  TAN_predicted_rna_seq_data_filtered_zscore, 
  pluvicto_predicted_rna_seq_data_zscore,
  TAN_predicted_AR_10_zscore,
  pluvicto_predicted_AR_10_zscore,
  patients_of_interest, 
  genes_of_interest, 
  metadata_of_interest
)

combined_pluvicto_df <- na.omit(combined_pluvicto_df)

# Make dataframe for combined TAN data
TAN_raw_rna_seq_data_converted_filtered_copy <- TAN_raw_rna_seq_data_converted_filtered
rownames(TAN_raw_rna_seq_data_converted_filtered_copy) <- TAN_raw_rna_seq_data_converted_filtered_copy$Gene_symbol
TAN_raw_rna_seq_data_converted_filtered_copy$Gene_symbol <- NULL
combined_TAN_data <- combine_data_frames(TAN_raw_rna_seq_data_converted_filtered_copy, TAN_true_AR_10_converted_filtered)

# Merge TAN and pluvicto dataframes
heatmap_df <- merge(combined_pluvicto_df, combined_TAN_data, by = "row.names", all = TRUE)
rownames(heatmap_df) <- heatmap_df$Row.names
heatmap_df$Row.names <- NULL

# Create heatmap
areas_of_interest <- c(genes_of_interest, metadata_of_interest)
patient_list <- unlist(patients_of_interest, use.names = FALSE)

# Update dataframe to numeric
heatmap_df_numeric <- heatmap_df %>% mutate(across(where(is.character), as.numeric))
heatmap_df_numeric$data <- rownames(heatmap_df_numeric)

# Transform dataframe to a matrix - long style
heatmap_df_long <- melt(heatmap_df_numeric, id.vars = "data")

# Setting suffixes to genes
label_map <- setNames(genes_of_interest, genes_of_interest)
label_map[genes_of_interest %in% down_regulated] <- paste0(genes_of_interest[genes_of_interest %in% down_regulated], "_down")
label_map[genes_of_interest %in% up_regulated] <- paste0(genes_of_interest[genes_of_interest %in% up_regulated], "_up")
label_map[genes_of_interest %in% interacting] <- paste0(genes_of_interest[genes_of_interest %in% interacting], "_int")
label_map[genes_of_interest %in% main] <- paste0(genes_of_interest[genes_of_interest %in% main], "_main")
label_map[metadata_of_interest] <- metadata_of_interest

# Keep only genes that are actually in the dataframe
present_genes <- intersect(rownames(heatmap_df_numeric), genes_of_interest)
present_genes <- c(present_genes, metadata_of_interest)

# Extract their suffix labels
final_labels <- label_map[present_genes]

heatmap_plot <- ggplot(heatmap_df_long, aes(variable, data, fill = value)) +
                      geom_tile() +
                      scale_fill_gradient2(
                        low = "#4575B4",
                        mid = "white",
                        high = "#D73027",
                        midpoint = 0,
                        space = "Lab"
                      ) + 
                      scale_y_discrete(limits = present_genes, labels = final_labels) +
                      labs(title = paste0(pluvicto_id, " Heatmap Across Timepoints"),
                           x = "Patient and Timepoint",
                           y = "Area of interest")

time <- format(Sys.time(), "%Y%m%d_%H%M%S")

patient_list_collapsed <- paste(unlist(patient_list), collapse = "_")
areas_of_interest_collapsed <- paste(unlist(areas_of_interest), collapse = "_")

ggsave(file.path("./outputs", paste0(patient_list_collapsed, "_", areas_of_interest_collapsed, "_zscore_heatmap_", time, ".png")), plot = heatmap_plot, width = 10, height = 8, dpi = 300)




