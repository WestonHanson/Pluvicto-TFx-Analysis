# Author: Weston Hanson
# Date: 04/29/26
# Place: Fred Hutchinson Cancer Center - Seattle, WA
# Purpose: To visualize the Pluvicto structual variants. 

options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages('devtools')
install.packages('testthat')
install.packages("data.table")
install.packages("gGnome")

## allows dependencies that throw warnings to install
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)

devtools::install_github('mskilab/gGnome')

suppressPackageStartupMessages({
  library(data.table)
  library(gGnome)
  library(GenomicRanges)
})

# -----------------------------
# Parse arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("
Usage:
  Rscript run_gGnome.R <sample_name> <jabba_rds> <cov_rds> <outdir>

Example:
  Rscript run_gGnome.R FH0084_E_C1 patient.jabba.simple.rds cov.raw.rds output_dir
")
}

sample_name <- args[1]
jabba_file  <- args[2]
cov_file    <- args[3]
cov_fixed_file <- args[4]
outdir      <- args[5]

cat("[INFO] Sample:", sample_name)
cat("[INFO] Graph:", graph_file)
cat("[INFO] Coverage:", cov_file)
cat("[INFO] Output dir:", outdir)

# -----------------------------
# Fix coverage if needed
# -----------------------------
cat("[INFO] Reading coverage RDS...")
cov <- readRDS(cov_file)

if (!inherits(cov, "GRanges")) {
  stop(\"Coverage file must be a GRanges object\")
}

cols <- colnames(mcols(cov))

if (!("coverage" %in% cols)) {
  if ("count" %in% cols) {
    cat("[INFO] Converting 'count' -> 'coverage')
    mcols(cov)$coverage <- mcols(cov)$count
    mcols(cov)$count <- NULL

    fixed_cov_file <- sub("\\.rds$", ".fixed.rds", cov_file)
    saveRDS(cov, fixed_cov_file)
    cov_file <- fixed_cov_file

    cat("[INFO] Saved fixed coverage to:", cov_file)
  } else {
    stop("Coverage GRanges must contain 'coverage' or 'count' column")
  }
} else {
  cat("[INFO] Coverage already has 'coverage' column")
}

# -----------------------------
# Build data.table
# -----------------------------
jsdt <- data.table(
  sample   = sample_name,
  graph    = graph_file,
  coverage = cov_file
)

# -----------------------------
# Run gGnome
# -----------------------------
cat("[INFO] Running gGnome.js...")

gGnome.js(
  data   = jsdt,
  outdir = outdir,
  cov.col = \"coverage\"
)

cat("[INFO] Done! Output written to:", outdir)

# cov = readRDS(cov_file)
# # keep only standard chromosomes
# cov <- keepStandardChromosomes(cov, pruning.mode = "coarse")
# # optional: enforce naming style (very important)
# seqlevelsStyle(cov) <- "NCBI"   # gives "1,2,3"
# # OR
# # seqlevelsStyle(cov) <- "UCSC" # gives "chr1,chr2"
# # fix coverage column
# mcols(cov)$coverage <- mcols(cov)$count
# mcols(cov)$count <- NULL
# # save a fixed version
# saveRDS(cov, cov_fixed_file)
# 
# # generate Table
# jsdt = data.table(
#   sample = 'FH0084_E_C1',
#   graph = FH0084_E_C1_RDS,
#   coverage = cov_fixed_file
# )
# 
# # Generate gGnome.js project
# gGnome.js(
#   data = jsdt,
#   outdir = out_dir,
#   cov.col = "coverage",
#   append = TRUE
# )
