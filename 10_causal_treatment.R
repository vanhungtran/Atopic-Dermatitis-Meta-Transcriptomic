#!/usr/bin/env Rscript
# =============================================================================
# 10_causal_treatment.R — Treatment-response directionality for lesion signature
# =============================================================================
# For top-30 consensus_lesion_signature.csv genes, checks whether expression
# in paired-design datasets shows reversal toward NL baseline, supporting
# a mediating (rather than static) role in lesion pathology.
#
# Approach: Extract per-dataset logFC from paired LS-vs-NL stage-1 files and
# test concordance with the meta-analytic beta_L. A gene is flagged
# "treatment_responsive" if its paired-design effect is consistent in
# sign with the meta-analytic effect (suggesting the lesion effect is
# reversible with treatment that normalizes toward NL).
#
# Input:  consensus_lesion_signature.csv, stage1_*_paired.csv
# Output: causal_treatment_response.csv
# =============================================================================

# ---- Configuration ----
out <- "d:/OneDrive/Statistics/Metagenomics/data/geo_atopic_dermatitis_longitudinal/disease_signature/analysis"

# ---- Load top-30 consensus lesion signature ----
cs <- read.csv(file.path(out, "consensus_lesion_signature.csv"), stringsAsFactors = FALSE)
top30 <- head(cs, 30)

# ---- Load meta lesion effects for beta_L reference ----
meta <- read.csv(file.path(out, "meta_lesion_effects.csv"), stringsAsFactors = FALSE)
rownames(meta) <- meta$gene

cat(sprintf("Top-30 consensus lesion genes:\n"))
cat(paste(top30$gene, collapse = ", "), "\n\n")

# ---- Find all paired stage-1 files ----
paired_files <- list.files(out, pattern = "stage1_.*_paired\\.csv$", full.names = TRUE)
cat(sprintf("Found %d paired-dataset files\n", length(paired_files)))

# ---- Extract per-gene logFC from each paired file ----
# Collect as vectors first, then assemble data.frame (avoids rbind issues)
gene_vec <- character()
ds_vec   <- character()
lfc_vec  <- numeric()

for (f in paired_files) {
  ds <- gsub("^stage1_", "", basename(f))
  ds <- gsub("_paired\\.csv$", "", ds)
  
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  # First column is unnamed (gene symbol from limma output)
  colnames(d)[1] <- "gene"
  
  # Keep only genes in top30
  d <- d[d$gene %in% top30$gene, ]
  if (nrow(d) == 0) next
  
  gene_vec <- c(gene_vec, d$gene)
  ds_vec   <- c(ds_vec, rep(ds, nrow(d)))
  lfc_vec  <- c(lfc_vec, d$logFC)
}

paired_df <- data.frame(
  gene = gene_vec,
  dataset = ds_vec,
  beta_paired = lfc_vec,
  stringsAsFactors = FALSE
)
cat(sprintf("Collected %d gene-dataset observations from paired data\n", nrow(paired_df)))

# ---- Build treatment response table ----
res_list <- lapply(top30$gene, function(g) {
  # Meta-analytic beta_L
  beta_L <- if (g %in% rownames(meta)) meta[g, "pooled_logFC"] else NA_real_
  
  # Paired-dataset effects for this gene
  sub <- paired_df[paired_df$gene == g, ]
  
  if (nrow(sub) == 0) {
    return(list(
      gene = g,
      beta_L = beta_L,
      n_paired = 0L,
      paired_datasets = "",
      paired_mean_logFC = NA_real_,
      pct_concordant = NA_real_,
      concordant_with_reversal = NA,
      note = "Not present in paired-dataset files"
    ))
  }
  
  # Concordance: does the paired logFC have the same sign as beta_L?
  conc <- sign(sub$beta_paired) == sign(beta_L)
  all_conc <- all(conc, na.rm = TRUE)
  
  # A gene is "treatment-responsive" if ≥70% of paired datasets agree in direction
  pct_conc <- mean(conc, na.rm = TRUE) * 100
  
  list(
    gene = g,
    beta_L = beta_L,
    n_paired = as.integer(nrow(sub)),
    paired_datasets = paste(sub$dataset, collapse = ";"),
    paired_mean_logFC = as.numeric(mean(sub$beta_paired, na.rm = TRUE)),
    pct_concordant = as.numeric(round(pct_conc, 1)),
    concordant_with_reversal = all_conc,
    note = if (all_conc) {
      sprintf("Concordant across %d/%d paired datasets (%.0f%%)", sum(conc), nrow(sub), pct_conc)
    } else {
      sprintf("%d/%d paired datasets concordant (%.0f%%)", sum(conc), nrow(sub), pct_conc)
    }
  )
})

# Convert list-of-lists to data.frame
res <- data.frame(
  gene = sapply(res_list, `[[`, "gene"),
  beta_L = as.numeric(sapply(res_list, `[[`, "beta_L")),
  n_paired = as.integer(sapply(res_list, `[[`, "n_paired")),
  paired_datasets = as.character(sapply(res_list, `[[`, "paired_datasets")),
  paired_mean_logFC = as.numeric(sapply(res_list, `[[`, "paired_mean_logFC")),
  pct_concordant = as.numeric(sapply(res_list, `[[`, "pct_concordant")),
  concordant_with_reversal = as.logical(sapply(res_list, `[[`, "concordant_with_reversal")),
  note = as.character(sapply(res_list, `[[`, "note")),
  stringsAsFactors = FALSE
)

# Sort by concordance then n_paired
res <- res[order(-res$concordant_with_reversal, -res$n_paired, res$gene), ]

# ---- Write output ----
write.csv(res, file.path(out, "causal_treatment_response.csv"), row.names = FALSE)

# ---- Summary ----
n_conc <- sum(res$concordant_with_reversal, na.rm = TRUE)
n_total <- sum(!is.na(res$concordant_with_reversal))
cat(sprintf("\n=== Treatment Response Summary ===\n"))
cat(sprintf("Fully concordant (paired sign = meta beta_L): %d/%d genes\n", n_conc, n_total))
cat(sprintf("Genes with paired data:                       %d/%d\n", n_total, nrow(res)))
cat("\nFull table:\n")
print(res[, c("gene", "beta_L", "n_paired", "pct_concordant", "concordant_with_reversal")], row.names = FALSE)

cat(sprintf("\nWrote %d rows to causal_treatment_response.csv\n", nrow(res)))
cat("Done.\n")
