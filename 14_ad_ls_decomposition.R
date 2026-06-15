#!/usr/bin/env Rscript
# =============================================================================
# 14_ad_ls_decomposition.R — Decompose the AD transcriptome into three groups
# =============================================================================
# Using the three sample types in this dataset (AD_LS, AD_NL, HC) and the two
# resulting axes, partition genes into the three groups the analysis can resolve
# WITHOUT any external disease:
#
#   AD axis  (beta_D) = AD_NL vs HC      -> "what AD does to skin at baseline"
#   LS axis  (beta_L) = AD_LS vs AD_NL   -> "what forming a lesion adds"
#
#   group "AD_related"  : sig on AD axis only          (trait of having AD)
#   group "LS_related"  : sig on LS axis only          (lesion reaction)
#   group "AD_and_LS"   : sig on both, same direction  (core AD genes)
#   (genes sig on both with OPPOSITE direction are reported as "discordant")
#
# Input:  causal_trait_state_classification.csv  (from 07_causal_trait_state.R)
# Output: ad_ls_gene_decomposition.csv
#         ad_ls_decomposition_summary.csv
#         plot_ad_ls_decomposition.png / .pdf
# =============================================================================

suppressMessages({ library(ggplot2) })

out <- "d:/OneDrive/Statistics/Metagenomics/data/geo_atopic_dermatitis_longitudinal/disease_signature/analysis"
source(file.path(out, "scripts", "plotting_utils.R"))

cls <- read.csv(file.path(out, "causal_trait_state_classification.csv"),
                stringsAsFactors = FALSE)
cat(sprintf("Loaded %d classified genes\n", nrow(cls)))

# ---- Relabel internal categories -> user-facing AD / LS / both naming ----
remap <- c(
  trait_only                 = "AD_related",
  state_only                 = "LS_related",
  trait_and_state_concordant = "AD_and_LS",
  trait_and_state_discordant = "discordant"
)
cls$group <- remap[cls$category]

# Order columns for a clean, self-explanatory table
keep <- c("gene", "group",
          "beta_D", "padj_D", "beta_L", "padj_L",
          "I2_D", "k_D", "I2_L", "k_L")
keep <- keep[keep %in% names(cls)]
res <- cls[, keep]

# Sort: AD_and_LS first (most interpretable), then AD_related, then the rest,
# each by combined effect magnitude.
grp_order <- c("AD_and_LS", "AD_related", "discordant", "LS_related")
res <- res[order(factor(res$group, levels = grp_order),
                 -pmax(abs(res$beta_D), abs(res$beta_L), na.rm = TRUE)), ]

write.csv(res, file.path(out, "ad_ls_gene_decomposition.csv"), row.names = FALSE)

# ---- Summary counts ----
summ <- as.data.frame(table(res$group), stringsAsFactors = FALSE)
names(summ) <- c("group", "n_genes")
summ <- summ[match(grp_order, summ$group), ]
summ <- summ[!is.na(summ$group), ]
write.csv(summ, file.path(out, "ad_ls_decomposition_summary.csv"), row.names = FALSE)

cat("\n=== AD / LS decomposition (HK + |log2FC|>=0.5 filtered) ===\n")
print(summ, row.names = FALSE)
cat("\nAD_and_LS genes (core AD):\n")
print(res$gene[res$group == "AD_and_LS"])
cat("\nAD_related genes (AD trait, baseline):\n")
print(res$gene[res$group == "AD_related"])

# ---- Quadrant figure: density-style overview with focused labels ----
pd <- cls[!is.na(cls$beta_D) & !is.na(cls$beta_L), ]
pd$group <- factor(remap[pd$category],
                   levels = c("AD_and_LS", "AD_related",
                              "LS_related", "discordant"))
pd$beta_D_plot <- pmax(pmin(pd$beta_D, 4), -4)
pd$beta_L_plot <- pmax(pmin(pd$beta_L, 3), -3)
cols <- c(AD_and_LS = ad_colors$both, AD_related = ad_colors$trait,
          LS_related = "#AEB8C2", discordant = ad_colors$discordant)

# Label half of AD_and_LS + AD_related genes (stratified sampling)
lab <- pd[pd$group %in% c("AD_and_LS", "AD_related"), ]
lab$label_score <- abs(lab$beta_D) + 0.5 * abs(lab$beta_L)
lab <- lab[order(-lab$label_score), ]
# Sample ~half: keep every other by score rank within each group
set.seed(42)
lab <- do.call(rbind, lapply(split(lab, lab$group), function(x) {
  x[seq_len(ceiling(nrow(x) / 2)), , drop = FALSE]
}))
rownames(lab) <- NULL
lab <- lab[order(-lab$label_score), ]
count_labels <- summ
names(count_labels)[names(count_labels) == "n_genes"] <- "n"
count_labels$label <- paste0(count_labels$group, "\n", count_labels$n, " genes")
count_labels$x <- c(3.15, -3.2, 3.15, -3.2)[match(count_labels$group,
                                                c("AD_and_LS", "AD_related",
                                                  "LS_related", "discordant"))]
count_labels$y <- c(1.15, 2.10, -2.35, -2.10)[match(count_labels$group,
                                                c("AD_and_LS", "AD_related",
                                                  "LS_related", "discordant"))]

p <- ggplot(pd, aes(beta_D_plot, beta_L_plot, color = group)) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = ad_colors$muted) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = ad_colors$muted) +
  geom_point(data = pd[pd$group == "LS_related", ], color = "#CBD5E1",
             alpha = 0.28, size = 1.0) +
  geom_point(data = pd[pd$group != "LS_related", ],
             aes(size = abs(beta_D) + abs(beta_L)), alpha = 0.82) +
  geom_label(data = count_labels, aes(x = x, y = y, label = label, fill = group),
             color = "white", fontface = "bold", size = 3.1,
             linewidth = 0, show.legend = FALSE, inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 2.2,
                           max.overlaps = 80, show.legend = FALSE,
                           min.segment.length = Inf, segment.color = NA,
                           box.padding = 0.3, force = 0.5) +
  scale_color_manual(values = cols, name = NULL,
                     labels = c("AD + LS (core)", "AD-related (trait)",
                                "LS-related (lesion)", "discordant")) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_size_continuous(range = c(1.7, 4.8), guide = "none") +
  labs(x = expression(beta[D]*": AD non-lesional vs healthy control  (AD axis)"),
       y = expression(beta[L]*": AD lesional vs non-lesional  (LS axis)"),
       title = "Decomposition of the AD transcriptome",
       subtitle = expression(beta[D]*" and "*beta[L]*" are random-effects meta-analytic log"[2]*" fold changes; axes capped for readability")) +
  theme_ad_pub(base_size = 11) +
  coord_cartesian(xlim = c(-4.2, 4.2), ylim = c(-3.1, 3.1))

save_pub_pair("plot_ad_ls_decomposition", p, width = 8.5, height = 7.2,
              dpi = 600, out_dir = out)

cat("\nWrote: ad_ls_gene_decomposition.csv, ad_ls_decomposition_summary.csv,",
    "plot_ad_ls_decomposition.{png,pdf}\nDone.\n")
