#!/usr/bin/env Rscript
# =============================================================================
# 12_disease_biomarker_signature.R — Final AD Disease Biomarker Signature
# =============================================================================
# Produces one canonical, artifact-free disease-signature table by:
#   Stage 3: recomputing the composite biomarker score on the pooled
#            disease-axis effects (beta_D: AD non-lesional vs healthy control),
#   Stage 4: joining the per-dataset consensus vote count (n_sig),
#   then applying the same housekeeping + |log2FC| filters used in step 07 so
#   the headline signature is not contaminated by normalization artifacts.
#
# Input:  meta_disease_effects.csv, consensus_disease_signature.csv
# Output: disease_biomarker_signature_final.csv
#         plot_disease_signature_top30.png / .pdf
# =============================================================================

suppressMessages({ library(ggplot2) })

# ---- Configuration ----
out <- "d:/OneDrive/Statistics/Metagenomics/data/geo_atopic_dermatitis_longitudinal/disease_signature/analysis"
source(file.path(out, "scripts", "plotting_utils.R"))
alpha   <- 0.05    # padj cutoff
fc_min  <- 0.5     # min |pooled log2FC|
top_n   <- 30      # genes to plot
hk_pattern <- "^(RPS|RPL|RPLP|MRPS|MRPL|EEF|EIF|ACTB|ACTG1|GAPDH|TUBB|TUBA|B2M|HPRT1|PPIA|UBC|YWHAZ)[0-9A-Z]*$"

# ---- Load pooled disease effects ----
D <- read.csv(file.path(out, "meta_disease_effects.csv"), stringsAsFactors = FALSE)
D_all <- D
cat(sprintf("Pooled disease genes: %d\n", nrow(D)))

# ---- Load consensus vote counts (how many datasets each gene is sig in) ----
C <- read.csv(file.path(out, "consensus_disease_signature.csv"), stringsAsFactors = FALSE)
C <- C[, c("gene", "n_sig")]

# ---- Filter: housekeeping/structural + trivial effect sizes ----
hk <- grepl(hk_pattern, D$gene)
cat(sprintf("Removing %d housekeeping/structural genes\n", sum(hk)))
D <- D[!hk, ]
D <- D[!is.na(D$padj) & D$padj < alpha & abs(D$pooled_logFC) >= fc_min, ]
cat(sprintf("After padj<%.2f & |log2FC|>=%.1f: %d genes\n", alpha, fc_min, nrow(D)))

# ---- Stage 3: (re)compute composite biomarker score ----
# score = -log10(p) * |log2FC| * (1 - I2) * sqrt(k)
# I2 in the file is on a 0-100 scale -> convert to 0-1.
I2_frac <- pmin(pmax(D$I2 / 100, 0), 1)
D$score <- -log10(D$pval) * abs(D$pooled_logFC) * (1 - I2_frac) * sqrt(D$k)

# ---- Stage 4: attach consensus vote count ----
nbefore <- nrow(D)
D <- merge(D, C, by = "gene", all.x = TRUE)
D$n_sig[is.na(D$n_sig)] <- 0
stopifnot(nrow(D) == nbefore)   # left join must not change row count

# ---- Assemble final, ordered signature table ----
keep_cols <- c("gene", "pooled_logFC", "se", "ci_lower", "ci_upper",
               "pval", "padj", "I2", "k", "n_sig", "score")
sig <- D[, keep_cols]
sig$direction <- ifelse(sig$pooled_logFC > 0, "up", "down")
sig <- sig[order(-sig$score), ]

write.csv(sig, file.path(out, "disease_biomarker_signature_final.csv"),
          row.names = FALSE)

cat(sprintf("\nFinal disease biomarker signature: %d genes (%d up, %d down)\n",
            nrow(sig), sum(sig$direction == "up"), sum(sig$direction == "down")))
cat("\nTop 20 by composite score:\n")
print(head(sig[, c("gene", "pooled_logFC", "padj", "I2", "k", "n_sig", "score")], 20))

# ---- Plot top-N as a ranked fold-change bar (colored by consensus) ----
pd <- head(sig, top_n)
pd$gene <- factor(pd$gene, levels = rev(pd$gene))
p <- ggplot(pd, aes(pooled_logFC, gene, fill = n_sig)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3,
                 color = ad_colors$ink, linewidth = 0.35) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = ad_colors$muted) +
  scale_fill_gradientn(colors = c("#D9EAF7", ad_colors$disease, ad_colors$both),
                       name = "# datasets\nsignificant") +
  labs(x = expression("Pooled "*log[2]*"FC  (AD non-lesional vs healthy control)"),
       y = NULL,
       title = sprintf("Top %d AD disease biomarkers", top_n),
       subtitle = "Random-effects pooled effect; whiskers show 95% CI") +
  theme_ad_pub(base_size = 11) +
  theme(axis.text.y = element_text(face = "bold"))

save_pub_pair("plot_disease_signature_top30", p, width = 8.5, height = 7.2,
              dpi = 600, out_dir = out)

cat("\nWrote: disease_biomarker_signature_final.csv,",
    "plot_disease_signature_top30.{png,pdf}\nDone.\n")
