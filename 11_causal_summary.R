#!/usr/bin/env Rscript
# =============================================================================
# 11_causal_summary.R — Integrate all causal evidence into a single summary
# =============================================================================
# Joins results from Steps 1, 2 and 4 into one unified table per candidate gene:
#   (1) trait/state classification,
#   (2) GWAS literature support  — a single annotation column. NOTE: the former
#       "coloc" (Step 2) and "MR" (Step 3) outputs are the SAME literature
#       look-up (is the gene in >=2 curated AD GWAS catalogs); no coloc.abf() or
#       TwoSampleMR was run (sumstats unavailable). They are therefore collapsed
#       into ONE support flag here, not counted twice.
#   (4) treatment-response directionality — the primary causal-direction evidence.
#
# Input:  causal_trait_state_classification.csv, causal_coloc_results.csv,
#         causal_treatment_response.csv   (causal_mr_results.csv NOT used —
#         redundant with coloc; see note above)
# Output: causal_summary.csv
# =============================================================================

# ---- Configuration ----
out <- "d:/OneDrive/Statistics/Metagenomics/data/geo_atopic_dermatitis_longitudinal/disease_signature/analysis"

# ---- Load step outputs ----
ts  <- read.csv(file.path(out, "causal_trait_state_classification.csv"), stringsAsFactors = FALSE)
col <- read.csv(file.path(out, "causal_coloc_results.csv"), stringsAsFactors = FALSE)
trt <- read.csv(file.path(out, "causal_treatment_response.csv"), stringsAsFactors = FALSE)

cat(sprintf("Loaded: trait/state=%d, gwas-literature=%d, treatment=%d genes\n",
            nrow(ts), nrow(col), nrow(trt)))

# ---- Build per-gene summary from candidate list ----
cand <- readLines(file.path(out, "causal_candidate_genes.txt"))

res <- do.call(rbind, lapply(cand, function(g) {
  # Trait/state info
  ts_row <- ts[ts$gene == g, ]
  category <- if (nrow(ts_row) > 0) ts_row$category[1] else NA
  trait_like <- category %in% c("trait_only", "trait_and_state_concordant")
  
  # GWAS literature support (single annotation; NOT a coloc/MR analysis).
  # The former Step-2 "coloc" and Step-3 "MR" outputs are the same curated-
  # catalogue look-up, so they contribute ONE support flag, not two.
  col_row <- col[col$gene == g, ]
  gwas_literature_support <- if (nrow(col_row) > 0) col_row$gwas_support[1] else FALSE
  n_gwas <- if (nrow(col_row) > 0) col_row$n_gwas_catalogs[1] else 0
  eqtl_skin <- if (nrow(col_row) > 0) col_row$eqtl_skin[1] else FALSE

  # Treatment response (only defined for top-30 lesion consensus genes)
  trt_row <- trt[trt$gene == g, ]
  treatment_responsive <- if (nrow(trt_row) > 0) {
    if (!is.na(trt_row$concordant_with_reversal[1])) trt_row$concordant_with_reversal[1] else NA
  } else NA
  
  # Count supports (max 3): trait-like, GWAS literature support (one flag),
  # and treatment-response directionality.
  supports <- sum(c(
    trait_like,
    gwas_literature_support,
    if (!is.na(treatment_responsive)) treatment_responsive else FALSE
  ))
  
  # Tier assignment (3 = all supports, 2, 1, 0)
  tier <- if (supports >= 3) "Tier 1"
          else if (supports >= 2) "Tier 2"
          else if (supports >= 1) "Tier 3"
          else "Associational only"

  data.frame(
    gene = g,
    category = category,
    trait_like = trait_like,
    gwas_literature_support = gwas_literature_support,
    n_gwas_catalogs = n_gwas,
    eqtl_skin = eqtl_skin,
    treatment_responsive = treatment_responsive,
    n_supports = supports,
    causal_evidence_tier = tier,
    stringsAsFactors = FALSE
  )
}))

# Sort by tier then gene
res <- res[order(factor(res$causal_evidence_tier,
                 levels = c("Tier 1", "Tier 2", "Tier 3", "Associational only")),
                 res$gene), ]

# ---- Write output ----
write.csv(res, file.path(out, "causal_summary.csv"), row.names = FALSE)

# ---- Summary report ----
cat("\n========================================\n")
cat("   CAUSAL EVIDENCE SUMMARY\n")
cat("========================================\n\n")

tier_counts <- table(res$causal_evidence_tier)
for (tn in names(tier_counts)) {
  cat(sprintf("%s: %d genes\n", tn, tier_counts[tn]))
}
cat(sprintf("\nTotal candidate genes: %d\n", nrow(res)))

# Print Tier 1 and Tier 2 genes
cat("\n--- Tier 1 genes (≥3 supports) ---\n")
t1 <- res[res$causal_evidence_tier == "Tier 1", ]
if (nrow(t1) > 0) {
  print(t1[, c("gene", "category", "trait_like", "gwas_literature_support", "treatment_responsive")])
} else {
  cat("None reached Tier 1.\n")
  cat("NOTE: none of the candidate genes appear in the curated AD GWAS\n")
  cat("catalogues used here. This is NOT interpreted as 'novelty': the genetic\n")
  cat("reference is incomplete (no live GWAS/eQTL sumstats were available), so an\n")
  cat("absence of overlap against a partial loci list is uninterpretable as a\n")
  cat("test of genetic support. GWAS is reported only as a supporting annotation.\n")
}

cat("\n--- Tier 2 genes (2 supports) ---\n")
t2 <- res[res$causal_evidence_tier == "Tier 2", ]
if (nrow(t2) > 0) {
  print(t2[, c("gene", "category", "trait_like", "gwas_literature_support", "treatment_responsive")])
} else {
  cat("None.\n")
}

cat(sprintf("\nWrote %d rows to causal_summary.csv\n", nrow(res)))
cat("Done.\n")
