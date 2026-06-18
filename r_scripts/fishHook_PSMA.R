# Author: Weston Hanson
# Place: Fred Hutchinson Cancer Center, Seattle, WA
# Date Created: 02/12/26
# Purpose: To use the R package fishHook to find fragile sites in the Pluvicto data.

# !!!!!!!!!!!!!!!!!!
# Load Packages
# !!!!!!!!!!!!!!!!!!

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(ask = FALSE, update = TRUE)

BiocManager::version()

############################################################
## Install required Bioconductor packages
############################################################

BiocManager::install(c(
  "fishHook",
  "GenomicRanges",
  "IRanges",
  "rtracklayer",
  "gUtils",
  "gTrack",
  "BSgenome.Hsapiens.UCSC.hg38",
  "TxDb.Hsapiens.UCSC.hg38.knownGene"
), ask = FALSE)

############################################################
## Install CRAN dependencies
############################################################

install.packages(c(
  "data.table",
  "devtools"
))

############################################################
## Install skitools from GitHub
############################################################

devtools::install_github("mskilab-org/gChain")

devtools::install_github("mskilab-org/skitools")

############################################################
## Validate installation
############################################################

BiocManager::valid()

############################################################
## Load libraries
############################################################

library(fishHook)
library(gTrack)
library(rtracklayer)
library(skitools)
library(data.table)
library(GenomeInfoDb)
library(GenomicRanges)
require(tidyverse)
require(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(org.Hs.eg.db)
library(GenomicRanges)
library(Biostrings)
library(rtracklayer)

# ---------------------------------------------------------------------------------------------------------------------------------------

fishHook_function <- function(sv_file_path, output_file_path, covariates_granges) {
  ####################################
  # Load Files
  ####################################
  sv_dt = fread(sv_file_path)
  
  ####################################
  # Define Breakpoints
  ####################################
  bp1 = GRanges(seqnames = sv_dt$chromosome_1, ranges = IRanges(start = sv_dt$start_1, end = sv_dt$start_1))
  bp2 = GRanges(seqnames = sv_dt$chromosome_2, ranges = IRanges(start = sv_dt$start_2, end = sv_dt$start_2))
  
  sv_breakpoints = c(bp1, bp2)
  seqlevelsStyle(sv_breakpoints) = "UCSC"
  sv_breakpoints = keepStandardChromosomes(sv_breakpoints, pruning.mode = "coarse")
  
  # Explicitly reorder seqlevels to standard order
  standard_order = paste0("chr", c(1:22, "X"))
  sv_breakpoints = keepSeqlevels(sv_breakpoints, intersect(standard_order, seqlevels(sv_breakpoints)), pruning.mode = "coarse")
  seqlevels(sv_breakpoints) = intersect(standard_order, seqlevels(sv_breakpoints))
  
  ####################################
  # hg38 model
  ####################################
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  genes_gr = genes(txdb)
  
  seqlevelsStyle(genes_gr) = "UCSC"
  genes_gr = keepStandardChromosomes(genes_gr, pruning.mode = "coarse")
  
  gene_symbols = mapIds(org.Hs.eg.db, keys = names(genes_gr), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  
  mcols(genes_gr)$symbol = gene_symbols
  
  ####################################
  # Harmonize chromosomes (added on)
  ####################################
  std_chr = paste0("chr", c(1:22, "X"))
  valid_chr_genes = intersect(seqlevels(genes_gr), std_chr)
  valid_chr_sv = intersect(seqlevels(sv_breakpoints), std_chr)
  
  # Use only chromosomes that exist in BOTH objects
  common_chr = intersect(valid_chr_genes, valid_chr_sv)
  
  genes_gr = keepSeqlevels(genes_gr, common_chr, pruning.mode = "coarse")
  sv_breakpoints = keepSeqlevels(sv_breakpoints, common_chr, pruning.mode = "coarse")
  
  genome_uscs = BSgenome.Hsapiens.UCSC.hg38
  # Create a proper Seqinfo object from UCSC genome for common chromosomes only
  si = seqinfo(genome_uscs)[common_chr]
  seqinfo(sv_breakpoints) = si
  seqinfo(genes_gr) = si
  
  ####################################
  # Define eligible territory
  ####################################
  eligible = genes_gr
  
  ####################################
  # Checking inputs
  ####################################
  cat("=== genes_gr ===\n")
  show(genes_gr)
  show(seqlevels(genes_gr))
  show(seqinfo(genes_gr))
  
  cat("=== sv_breakpoints ===\n")
  show(sv_breakpoints)
  show(seqlevels(sv_breakpoints))
  show(seqinfo(sv_breakpoints))
  
  cat("=== eligible ===\n")
  show(eligible)
  show(seqlevels(eligible))
  show(seqinfo(eligible))
  
  cat("=== hypotheses ===\n")
  show(seqinfo(genes_gr))
  
  cat("=== events ===\n")
  show(seqinfo(sv_breakpoints))
  
  cat("=== eligible ===\n")
  show(seqinfo(eligible))

  ####################################
  # Create covariate BEFORE Fish object
  ####################################
  cov_obj = NULL
  if(!is.null(covariates_granges)){
    # Filter covariate to match the harmonized chromosomes
    cov_features = keepSeqlevels(covariates_granges, common_chr, pruning.mode = "coarse")
    # Set matching seqinfo
    seqinfo(cov_features) = si
    
    # Show filtered covariates
    cat("=== covariates (filtered) ===\n")
    print(cov_features)
    show(seqlevels(cov_features))
    show(seqinfo(eligible))
    
    # Create Cov object from filtered GRanges
    cov_obj = suppressWarnings(Cov(cov_features, name="GC_content", field="GC_content", type="numeric"))
    cat("=== Cov object created successfully ===\n")
    print(cov_obj)
    cat("=== Cov content ===\n")
    print(cov_obj$data)
    print(cov_obj$data[[1]])
    print(cov_obj$data[[1]][is.na(cov_obj$data[[1]]$A)])
    print(cov_obj$data[[1]][is.na(cov_obj$data[[1]]$C)])
    print(cov_obj$data[[1]][is.na(cov_obj$data[[1]]$G)])
    print(cov_obj$data[[1]][is.na(cov_obj$data[[1]]$T)])
    print(cov_obj$data[[1]][is.na(cov_obj$data[[1]]$GC_content)])
  }

  ####################################
  # Run fishHook
  ####################################
  # We have loaded our hypothesis (genes), events, and eligible - now we run fishHook
  # Note: Object creation triggers counting of how many events are in the eligible portion of each hypthesis interval so the 'idcol' parameter makes sure each counts of each interval 
  #       gets at most one event. 
  
  fish = Fish(hypotheses = genes_gr, events = sv_breakpoints, eligible = eligible)
  
  # Assign covariate after Fish creation
  if(!is.null(cov_obj)){
    tryCatch({
      suppressWarnings({
        # Get the seqinfo directly from the fish object's hypotheses to ensure perfect match
        fish_seqinfo = seqinfo(fish@hypotheses)
        # Extract the underlying GRanges from the Cov object and update its seqinfo
        # Note: we can't directly modify Cov internals, so we'll try assignment as-is
        fish$covariates = cov_obj
        cat("=== Covariate assigned to fish object ===\n")
      })
    }, error = function(e) {
      cat("Error assigning covariate:", conditionMessage(e), "\n")
      cat("Continuing without covariates...\n")
      # Don't stop - let the model run without covariates
    })
  }
  
  # Calculate and display
  cat("=== Running fish$score() ===\n")
  tryCatch({
    suppressWarnings({
      fish$score()
      head(fish$res %Q% order(p))
      fish$qqp(plotly = FALSE)
    })
  }, error = function(e) {
    cat("Error during fish$score():", conditionMessage(e), "\n")
    cat("Full error:\n")
    print(e)
    # Continue anyway to see if we can still save results
    cat("Attempting to continue...\n")
  }, warning = function(w) {
    cat("Warning during fish$score():", conditionMessage(w), "\n")
  })
  
  ####################################
  # Save results
  ####################################
  res = as.data.table(fish$res)
  
  res_ordered = res[order(res$fdr)]
  
  fwrite(res_ordered, file = output_file_path, sep = ",")
  
  return("fishHook completed.")
}

####################################
# Create list of files
####################################

files_list = list.files(path = "/Volumes/ha_g/user/rpatton/PSMA_data/annotatedSVs_2025/annotatedCalls_filtered/", pattern = "*.txt", full.names = TRUE)

# format list into dataframe
files_dt = data.table(
  path        = files_list,
  patient_id  = sub("\\..*$", "", basename(files_list))
)

# Add output path for no covariates
files_dt[, no_covariates_out := file.path(
  "/Volumes/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/scripts/outputs/fishHook-analysis/no-covariates",
  paste0(patient_id, "_fishHook_sv_results_no_covariates", ".csv")
)]

# Add output path for covariates
files_dt[, covariates_out := file.path(
  "/Volumes/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/scripts/outputs/fishHook-analysis/covariates",
  paste0(patient_id, "_fishHook_sv_results_covariates", ".csv")
)]

####################################
# Run fishHook with no covariates
####################################
# call fishHook_no_covariates with no covariates
out_list = mapply(fishHook_function, sv_file_path = files_dt$path, output_file_path = files_dt$no_covariates_out, SIMPLIFY = FALSE, MoreArgs = list(covariates_granges = NULL))

####################################
# Create covariates
####################################

# ----------------------------------------------------------
# Create nt_composition.hg38.rds (for GC content covariate)
# ----------------------------------------------------------

genome = BSgenome.Hsapiens.NCBI.GRCh38
keep_chr = c(as.character(1:22), "X")
genome_info = seqinfo(genome)[keep_chr]

bin_size = 1e7
bins = tileGenome(genome_info, tilewidth = bin_size, cut.last.tile.in.chrom = TRUE)
seqs = getSeq(genome, bins)

nt_freq = letterFrequency(seqs, letters = c("A", "C", "G", "T"), as.prob = TRUE)

# Turn chromosomes into UCSC style
seqlevelsStyle(bins) = "UCSC"
bins = keepStandardChromosomes(bins, pruning.mode = "coarse")

# Update seqlengths from UCSC genome
genome_ucsc = BSgenome.Hsapiens.UCSC.hg38
seqlengths(bins) = seqlengths(genome_ucsc)[seqlevels(bins)]

# Add frequencies for each base
bins$A = nt_freq[, "A"]
bins$C = nt_freq[, "C"]
bins$G = nt_freq[, "G"]
bins$T = nt_freq[, "T"]

# Save bins as RDS file
file_nt_comp = "/Volumes/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/data-files/nt_composition_v2.hg38.rds"
saveRDS(bins, file_nt_comp)

# Read in RDS file
data_nt_comp = readRDS(file_nt_comp)

# Calculate GC frequency
data_nt_comp$GC_content = data_nt_comp$G + data_nt_comp$C
data_nt_comp
seqinfo(data_nt_comp)

# Try different covariate
reptimedata_hg38 = gr.sub(import.bedGraph("/Volumes/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/data-files/RT_NHEK_Keratinocytes_Int92817591_hg38.bedgraph", 'chr', ''))
seqlevelsStyle(reptimedata_hg38) = "UCSC"
reptime = Cov(reptimedata_hg38, field = 'score', name = 'ReplicationTiming')

####################################
# Run fishHook with covariates
####################################
# call fishHook_function with covariates (pass the GRanges so it can be filtered per file)
out_list = mapply(fishHook_function, sv_file_path = files_dt$path, output_file_path = files_dt$covariates_out, SIMPLIFY = FALSE, MoreArgs = list(covariates_granges = data_nt_comp))





