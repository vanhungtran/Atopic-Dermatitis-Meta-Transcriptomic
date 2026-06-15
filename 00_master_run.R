#!/usr/bin/env Rscript
# =============================================================================
# AD Meta-Analysis Pipeline --- Master Runner
# Cross-Platform Meta-Analysis of the Atopic Dermatitis Transcriptomic
# Disease Signature Across 39 Datasets (2,138 Individuals)
# =============================================================================
# Author: Van Hung Huynh-Tran
# Contact: luctran.biostats@gmail.com
# R version: 4.5.0
# =============================================================================
#
# This script runs the full analysis in order. Each step sources its own
# sub-script, which loads only the packages it needs.
#
# Before running:
#   1. Place all stage1_GSE*.csv and meta_*_effects.csv files in the
#      analysis/ directory (curated from GEO; see METHODS.md for details).
#   2. Install required packages (see DEPENDENCIES section below).
#   3. Set your working directory to analysis/, then source this file.
#
# Usage from PowerShell:
#   Set-Location ".../disease_signature/analysis"
#   & "C:/Program Files/R/R-4.5.0/bin/Rscript.exe" scripts/00_master_run.R
#
# =============================================================================

# ---- Dependencies -----------------------------------------------------------
# Install once:
#   install.packages(c("metafor","limma","DESeq2","GEOquery","gprofiler2",
#                      "igraph","ggplot2","ggrepel","pheatmap","dplyr",
#                      "data.table","org.Hs.eg.db"))
#
# Bioconductor packages:
#   BiocManager::install(c("limma","DESeq2","GEOquery","gprofiler2",
#                          "org.Hs.eg.db","hgu133plus2.db","hgu133a.db"))

SCRIPTS <- "."                         # all scripts are in this directory
DATA_DIR <- ".."                       # curated data files are in analysis/

# NOTE: Individual scripts contain hardcoded paths to the analysis/
# directory. If porting to a different machine, update the `out <-` 
# assignments in each script, or set your working directory to match.
# The scripts are provided as-executed for reproducibility.

# ---- Step 0: Normalise gene names (run after any data updates) --------------
cat("\n=== Step 0: Gene-name normalisation ===\n")
source(file.path(SCRIPTS, "00_normalize_gene_names.R"))

# ---- Step 7: Trait vs State classification ----------------------------------
cat("\n=== Step 7: Trait/State classification ===\n")
source(file.path(SCRIPTS, "07_causal_trait_state.R"))

# ---- Step 8: GWAS literature annotation (supporting) ------------------------
cat("\n=== Step 8: GWAS literature annotation ===\n")
source(file.path(SCRIPTS, "08_causal_coloc.R"))

# ---- Step 10: Treatment-response directionality -----------------------------
cat("\n=== Step 10: Treatment response ===\n")
source(file.path(SCRIPTS, "10_causal_treatment.R"))

# ---- Step 11: Causal evidence integration -----------------------------------
cat("\n=== Step 11: Causal summary ===\n")
source(file.path(SCRIPTS, "11_causal_summary.R"))

# ---- Step 12: Disease biomarker signature -----------------------------------
cat("\n=== Step 12: Disease biomarker signature ===\n")
source(file.path(SCRIPTS, "12_disease_biomarker_signature.R"))

# ---- Step 13: Leave-k-out sensitivity ---------------------------------------
cat("\n=== Step 13: k-Sensitivity analysis ===\n")
source(file.path(SCRIPTS, "13_k_sensitivity.R"))

# ---- Step 14: AD/LS decomposition (three-group split) -----------------------
cat("\n=== Step 14: AD/LS decomposition ===\n")
source(file.path(SCRIPTS, "14_ad_ls_decomposition.R"))

# ---- Step 15: Drug repurposing (druggability + signature reversal) ----------
cat("\n=== Step 15: Drug repurposing ===\n")
if (!nzchar(Sys.getenv("RUN_FULL_SIGNATURESEARCH"))) {
  Sys.setenv(RUN_FULL_SIGNATURESEARCH = "1")
}
source(file.path(SCRIPTS, "15_combined_drug_repurposing.R"))

# ---- Step 16: Regenerate all figures ----------------------------------------
cat("\n=== Step 16: Regenerate all figures ===\n")
source(file.path(SCRIPTS, "16_regenerate_all_figures.R"))

cat("\n=== Pipeline complete. ===\n")
cat("Output CSV and PDF files are in the analysis/ directory.\n")
cat("Compile the manuscript: pdflatex manuscript_disease_signature.tex\n")

# =============================================================================
# INPUT FILES (place in analysis/ before running):
#   stage1_disease_GSE*.csv         per-dataset disease-contrast DE results
#   stage1_GSE*.csv                 per-dataset lesion-contrast DE results
#   meta_disease_effects.csv        pooled disease-axis random-effects MA
#   meta_lesion_effects.csv         pooled lesion-axis random-effects MA
#   consensus_disease_signature.csv consensus vote-count table
#   dataset_contrast_annotation.csv study-level metadata
#   logFC_disease_matrix.csv        disease-axis logFC matrix
#   logFC_lesion_matrix.csv         lesion-axis logFC matrix
#   SE_disease_matrix.csv           disease-axis SE matrix
#   SE_lesion_matrix.csv            lesion-axis SE matrix
#
# KEY OUTPUTS:
#   disease_biomarker_signature_final.csv  135-gene field-defect signature
#   core_disease_signature_final.csv       consensus-ranked disease genes
#   causal_summary.csv                     causal evidence integration
#   ad_ls_decomposition_summary.csv        three-group partition counts
#   plot_*.pdf                             all manuscript figures
# =============================================================================
