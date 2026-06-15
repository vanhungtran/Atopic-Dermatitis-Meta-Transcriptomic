#!/usr/bin/env Rscript
# =============================================================================
# 07_causal_trait_state.R — Trait vs. State classification of AD genes
# =============================================================================
# Splits genes into trait-like (already dysregulated in clinically-normal AD
# skin -> candidate driver/predisposing) vs state-like (only changes in active
# lesions -> reactive) by joining the two pooled meta-analysis contrasts.
#
#   beta_D : AD non-lesional vs healthy control   (disease / trait axis)
#   beta_L : AD lesional vs AD non-lesional        (lesion / state axis)
#
# Input:  meta_disease_effects.csv, meta_lesion_effects.csv
# Output: causal_trait_state_classification.csv
#         plot_trait_state_quadrant.png / .pdf
# =============================================================================

suppressMessages({ library(ggplot2) })

# ---- Configuration ----
out <- "d:/OneDrive/Statistics/Metagenomics/data/geo_atopic_dermatitis_longitudinal/disease_signature/analysis"
source(file.path(out, "scripts", "plotting_utils.R"))
alpha    <- 0.05   # padj significance threshold
fc_min   <- 0.5    # min |pooled log2FC| on the significant axis (drop trivial effects)
# Housekeeping / structural genes whose tiny, ultra-precise effects are almost
# certainly normalization artifacts rather than AD biology. Dropped before
# classification. Pattern matches ribosomal protein (RPS*/RPL*), elongation/
# initiation factors, and a few canonical housekeeping symbols.
hk_pattern <- "^(RPS|RPL|RPLP|MRPS|MRPL|EEF|EIF|ACTB|ACTG1|GAPDH|TUBB|TUBA|B2M|HPRT1|PPIA|UBC|YWHAZ)[0-9A-Z]*$"

# ---- Load the two pooled contrasts ----
D <- read.csv(file.path(out, "meta_disease_effects.csv"), stringsAsFactors = FALSE)
L <- read.csv(file.path(out, "meta_lesion_effects.csv"),  stringsAsFactors = FALSE)

# Keep only the columns we need, renamed per-axis (gene names already HGNC-
# normalized upstream; this join assumes that — see METHODS Effect Size section)
D <- data.frame(gene = D$gene, beta_D = D$pooled_logFC, padj_D = D$padj,
                I2_D = D$I2, k_D = D$k, stringsAsFactors = FALSE)
L <- data.frame(gene = L$gene, beta_L = L$pooled_logFC, padj_L = L$padj,
                I2_L = L$I2, k_L = L$k, stringsAsFactors = FALSE)

# Full outer join: a gene may be measured in only one axis
m <- merge(D, L, by = "gene", all = TRUE)

cat(sprintf("Genes: %d in disease, %d in lesion, %d in union\n",
            nrow(D), nrow(L), nrow(m)))

# ---- Drop housekeeping / structural genes (normalization artifacts) ----
hk <- grepl(hk_pattern, m$gene)
cat(sprintf("Removed %d housekeeping/structural genes (e.g. %s)\n",
            sum(hk), paste(head(m$gene[hk], 6), collapse = ", ")))
m <- m[!hk, ]

# ---- Classify ----
# Significance now also requires a non-trivial fold-change on the relevant axis
sigD <- !is.na(m$padj_D) & m$padj_D < alpha & abs(m$beta_D) >= fc_min
sigL <- !is.na(m$padj_L) & m$padj_L < alpha & abs(m$beta_L) >= fc_min
same_sign <- sign(m$beta_D) == sign(m$beta_L)

m$category <- "ns"                                              # neither sig
m$category[ sigD & !sigL] <- "trait_only"                      # driver candidate
m$category[!sigD &  sigL] <- "state_only"                      # reactive
m$category[ sigD &  sigL &  same_sign] <- "trait_and_state_concordant"
m$category[ sigD &  sigL & !same_sign] <- "trait_and_state_discordant"

# Drop genes significant in neither contrast from the saved table
res <- m[m$category != "ns", ]
res <- res[order(factor(res$category,
                 levels = c("trait_and_state_concordant", "trait_only",
                            "trait_and_state_discordant", "state_only")),
                 -pmax(abs(res$beta_D), abs(res$beta_L), na.rm = TRUE)), ]

write.csv(res, file.path(out, "causal_trait_state_classification.csv"),
          row.names = FALSE)

# Candidate gene list for downstream causal steps (coloc / MR):
# union of trait_only + concordant genes — the driver candidates.
cand <- sort(unique(res$gene[res$category %in%
              c("trait_only", "trait_and_state_concordant")]))
writeLines(cand, file.path(out, "causal_candidate_genes.txt"))
cat(sprintf("\nWrote %d candidate genes to causal_candidate_genes.txt\n",
            length(cand)))

cat("\n=== Category counts (padj < 0.05) ===\n")
print(table(res$category))

# ---- Quadrant plot (only genes scored in BOTH axes) ----
pdat <- m[!is.na(m$beta_D) & !is.na(m$beta_L) & m$category != "ns", ]
cols <- c(trait_only = ad_colors$trait,
          state_only = ad_colors$lesion,
          trait_and_state_concordant = ad_colors$both,
          trait_and_state_discordant = ad_colors$discordant)

# Label the most extreme genes per category for readability
lab <- do.call(rbind, lapply(split(pdat, pdat$category), function(d) {
  d[order(-(abs(d$beta_D) + abs(d$beta_L)))[seq_len(min(8, nrow(d)))], ]
}))

p <- ggplot(pdat, aes(beta_D, beta_L, color = category)) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = ad_colors$muted) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = ad_colors$muted) +
  geom_point(alpha = 0.68, size = 1.6) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene),
                           size = 3, max.overlaps = 20, show.legend = FALSE,
                           min.segment.length = 0, segment.color = ad_colors$muted) +
  scale_color_manual(values = cols) +
  labs(x = expression(beta[D]~": AD non-lesional vs healthy control (trait axis)"),
       y = expression(beta[L]~": AD lesional vs non-lesional (state axis)"),
       title = "Trait vs state classification of AD-associated genes",
       subtitle = sprintf("%d classified genes after housekeeping and |log2FC| filters",
                          nrow(res)),
       color = NULL) +
  theme_ad_pub(base_size = 11)

save_pub_pair("plot_trait_state_quadrant", p, width = 8.5, height = 7.2,
              dpi = 600, out_dir = out)

cat("\nWrote: causal_trait_state_classification.csv,",
    "plot_trait_state_quadrant.{png,pdf}\n")
cat("Done.\n")
