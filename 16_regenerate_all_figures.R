#!/usr/bin/env Rscript
# =============================================================================
# regenerate_all_figures.R — Master figure generation for AD meta-analysis manuscript
# =============================================================================
# Reads normalized meta-analysis outputs and regenerates all manuscript figures.
# Run after normalize_gene_names.R and downstream scripts (07-14).
# =============================================================================

suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(data.table)
  library(dplyr)
  library(pheatmap)
  library(igraph)
})

# ---- Configuration ----
out <- "d:/OneDrive/Statistics/Metagenomics/data/geo_atopic_dermatitis_longitudinal/disease_signature/analysis"
setwd(out)
source(file.path(out, "scripts", "plotting_utils.R"))

# Housekeeping pattern
hk_pattern <- "^(ACTB|ACTG1|B2M|GAPDH|HPRT1|PPIA|RPL[0-9]|RPLP[0-9]|RPS[0-9]|RPSA|TBP|TUBB|TUBA|UBC|YWHAZ|EEF|EIF|MRPS|MRPL|TFRC|SDHA|PGK1|VIM|HMBS)"
alpha <- 0.05
fc_min <- 0.5
theme_pub <- theme_ad_pub(base_size = 11)

# ---- Load data ----
cat("Loading normalized data...\n")
L <- read.csv("meta_lesion_effects.csv", stringsAsFactors = FALSE)
D <- read.csv("meta_disease_effects.csv", stringsAsFactors = FALSE)
ts <- read.csv("causal_trait_state_classification.csv", stringsAsFactors = FALSE)
db <- read.csv("disease_biomarker_signature_final.csv", stringsAsFactors = FALSE)
dec <- read.csv("ad_ls_gene_decomposition.csv", stringsAsFactors = FALSE)
lodo <- read.csv("sensitivity_lodo.csv", stringsAsFactors = FALSE)
cs <- read.csv("consensus_disease_signature.csv", stringsAsFactors = FALSE)
lb <- read.csv("lesion_biomarkers.csv", stringsAsFactors = FALSE)

if (file.exists("logFC_lesion_matrix.csv")) {
  lfc_mat <- as.matrix(read.csv("logFC_lesion_matrix.csv", row.names = 1))
}
if (file.exists("coexpr_gene_modules.csv")) {
  mods <- read.csv("coexpr_gene_modules.csv", stringsAsFactors = FALSE)
}
if (file.exists("coexpr_module_summary.csv")) {
  mod_summary <- read.csv("coexpr_module_summary.csv", stringsAsFactors = FALSE)
}
if (file.exists("coexpr_module_trait.csv")) {
  mod_trait <- read.csv("coexpr_module_trait.csv", stringsAsFactors = FALSE)
}
if (file.exists("gene_hub_scores_v2.csv")) {
  hubs <- read.csv("gene_hub_scores_v2.csv", stringsAsFactors = FALSE)
} else if (file.exists("hub_gene_scores.csv")) {
  hubs <- read.csv("hub_gene_scores.csv", stringsAsFactors = FALSE)
}

module_palette <- setNames(
  c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#B07AA1",
    "#EDC948", "#7B3294", "#9C755F", "#BAB0AC", "#2A9D8F", "#D55E00"),
  paste0("M", seq_len(12))
)

# Pick the most *specific* informative enrichment term for a module, instead of
# the raw smallest-p term (which is almost always a giant root-level GO term
# such as "multicellular organismal process" or "cytoplasm" — statistically
# top but biologically empty, and prone to duplicate labels across modules).
# Rule: among significant terms, prefer those whose term_size falls in an
# interpretable range; try progressively looser caps so a label is always
# returned, and within the chosen pool prefer specific pathway sources
# (KEGG/Reactome/WikiPathways) before broad GO, breaking ties by p_value.
# Return the FULL ranked candidate-term list for a module (best first), instead
# of a single term. Ranking removes giant root-level GO terms via a size cap,
# prefers readable GO:BP over idiosyncratic pathway names, then breaks ties by
# p_value. The ranked list lets us de-duplicate labels across modules: if two
# modules share a top term, the weaker one falls back to its next distinct term.
rank_terms <- function(e) {
  if (is.null(e) || nrow(e) == 0) return(NULL)
  ts <- suppressWarnings(as.numeric(e$term_size))
  e$.ts <- ifelse(is.na(ts), Inf, ts)
  # term_size band: smaller (more specific) is better, but keep a floor of 5.
  e$.band <- ifelse(e$.ts < 5, 9,
             ifelse(e$.ts <= 500, 1,
             ifelse(e$.ts <= 800, 2,
             ifelse(e$.ts <= 1500, 3, 4))))
  # Prefer GO:BP (readable biological processes); pathway DBs (KEGG/Reactome/WP)
  # and GO:MF/CC are fallbacks, as their names are often idiosyncratic or
  # off-context for skin (e.g. "Cytoskeleton in muscle cells", virus WikiPaths).
  src_rank <- c("GO:BP" = 1,
                "KEGG" = 2, "REAC" = 2, "WP" = 2,
                "GO:MF" = 3, "GO:CC" = 3)
  e$.sr <- ifelse(e$source %in% names(src_rank), src_rank[e$source], 4)
  e[order(e$.band, e$.sr, e$p_value), , drop = FALSE]
}

# Top terms with p above this are too weak to assert a clean function; flagged.
weak_p_thresh <- 1e-3

module_annotation <- NULL
if (exists("mod_summary") && nrow(mod_summary) > 0) {
  module_annotation <- mod_summary
  names(module_annotation) <- sub("^n_genes$", "n_genes", names(module_annotation))
  module_annotation$top_term <- NA_character_
  module_annotation$source <- NA_character_
  module_annotation$term_p <- NA_real_

  # Load ranked candidates per module.
  ranked <- list()
  for (m in module_annotation$module) {
    f <- file.path(out, paste0("enrichment_module_", m, ".csv"))
    if (file.exists(f)) {
      e <- read.csv(f, stringsAsFactors = FALSE)
      if (nrow(e) > 0) ranked[[m]] <- rank_terms(e)
    }
  }

  # Greedy de-duplication: assign modules in order of their strongest term's
  # p_value (most confident first), each taking the best candidate term not yet
  # claimed by another module so labels stay distinct.
  order_modules <- names(ranked)
  if (length(order_modules) > 0) {
    best_p <- sapply(order_modules, function(m) ranked[[m]]$p_value[1])
    order_modules <- order_modules[order(best_p)]
  }
  used_terms <- character(0)
  for (m in order_modules) {
    cand <- ranked[[m]]
    pick <- cand[!(cand$term_name %in% used_terms), , drop = FALSE]
    if (nrow(pick) == 0) pick <- cand            # all taken -> allow a repeat
    chosen <- pick[1, ]
    used_terms <- c(used_terms, chosen$term_name)
    label <- chosen$term_name
    if (is.na(chosen$p_value) || chosen$p_value > weak_p_thresh) {
      label <- paste0(label, " (weak)")
    }
    sel <- module_annotation$module == m
    module_annotation$top_term[sel] <- label
    module_annotation$source[sel]   <- chosen$source
    module_annotation$term_p[sel]   <- chosen$p_value
  }
  if (exists("mod_trait") && nrow(mod_trait) > 0) {
    module_annotation <- merge(module_annotation, mod_trait, by = "module", all.x = TRUE)
  }
  module_annotation$top_term[is.na(module_annotation$top_term)] <- "No enrichment annotation"
  module_annotation$module_label <- sprintf("%s  |  %s", module_annotation$module,
                                            module_annotation$top_term)
  write.csv(module_annotation, "coexpr_module_annotation.csv", row.names = FALSE)
}

# =============================================================================
# FIGURE 1: Lesion Volcano Plot
# =============================================================================
cat("Figure 1: Lesion Volcano...\n")
Lf <- L[L$k >= 5 & !grepl(hk_pattern, L$gene), ]
Lf$sig <- ifelse(Lf$padj < alpha & abs(Lf$pooled_logFC) >= fc_min, "FDR<0.05 & |log2FC|≥0.5",
          ifelse(Lf$padj < alpha, "FDR<0.05", "NS"))
Lf$sig <- factor(Lf$sig, levels = c("NS", "FDR<0.05", "FDR<0.05 & |log2FC|≥0.5"))
Lf$nlogP <- -log10(Lf$padj)

top_lab <- head(Lf[order(Lf$padj), ], 20)
g <- ggplot(Lf, aes(x = pooled_logFC, y = nlogP, color = sig)) +
  geom_point(alpha = 0.42, size = 0.75) +
  scale_color_manual(values = c("#B8C2CC", ad_colors$disease, ad_colors$up), name = "") +
  geom_vline(xintercept = c(-fc_min, fc_min), linetype = "dashed",
             color = ad_colors$muted, alpha = 0.55) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed",
             color = ad_colors$muted, alpha = 0.55) +
  geom_text_repel(data = top_lab, aes(label = gene), size = 2.9, max.overlaps = 25,
                  min.segment.length = 0, segment.color = ad_colors$muted) +
  labs(x = expression(log[2]~"Fold Change (LS vs NL)"),
       y = expression(-log[10](P[adj])),
       title = "AD lesional vs non-lesional meta-analysis",
       subtitle = sprintf("%s genes tested; red marks FDR < 0.05 and |log2FC| >= %.1f",
                          comma(nrow(Lf)), fc_min)) +
  theme_pub
save_pub_pair("plot_lesion_volcano", g, width = 8.5, height = 7.2, dpi = 500)

# =============================================================================
# FIGURE 2: Heterogeneity vs Effect Size
# =============================================================================
cat("Figure 2: Heterogeneity vs Effect Size...\n")
Lf$pooled_logFC_plot <- pmax(pmin(Lf$pooled_logFC, 4), -4)
Lf$outlier_fc <- abs(Lf$pooled_logFC) > 4
L_sig <- Lf[Lf$sig == "FDR<0.05 & |log2FC|≥0.5", ]
g <- ggplot(Lf, aes(x = pooled_logFC_plot, y = I2)) +
  stat_bin_2d(bins = 70, aes(fill = after_stat(log10(count))),
              alpha = 0.88, show.legend = TRUE) +
  geom_point(data = L_sig, aes(color = ifelse(pooled_logFC > 0, "Up", "Down"),
                               shape = outlier_fc),
             alpha = 0.72, size = 0.75) +
  scale_color_manual(values = direction_palette[c("Up", "Down")], name = "") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "Within display range",
                                `TRUE` = "|log2FC| > 4"),
                     name = "") +
  scale_fill_gradientn(colors = c("#F8FAFC", "#9CC9DF", ad_colors$disease, ad_colors$both),
                       name = expression(log[10]("gene count"))) +
  geom_hline(yintercept = 50, linetype = "dashed", color = ad_colors$up, alpha = 0.55) +
  labs(x = expression(log[2]~"Fold Change"), y = expression(I^2~"(%)"),
       title = "Effect size and heterogeneity landscape",
       subtitle = "Display capped at +/-4 log2FC; colored points meet FDR and effect-size criteria") +
  coord_cartesian(xlim = c(-4.15, 4.15), ylim = c(0, 100)) +
  guides(shape = "none") +
  theme_ad_pub(base_size = 11, legend = "right")
save_pub_pair("plot_lesion_i2_vs_fc", g, width = 8.2, height = 6.2, dpi = 500)

# =============================================================================
# FIGURE 3: I² Histogram
# =============================================================================
cat("Figure 3: I² Histogram...\n")
g <- ggplot(Lf, aes(x = I2)) +
  geom_histogram(bins = 60, fill = ad_colors$disease, alpha = 0.82, color = "white") +
  geom_vline(xintercept = 50, linetype = "dashed", color = ad_colors$up) +
  labs(x = expression(I^2~"(%)"), y = "Number of Genes",
       title = "Distribution of between-study heterogeneity",
       subtitle = "Lesion contrast after k and housekeeping filters") +
  theme_pub
save_pub_pair("plot_lesion_i2_hist", g, width = 7.2, height = 5.2, dpi = 500)

# =============================================================================
# FIGURE 4: Forest Plot — Top 16 Lesion Biomarkers (strict)
# =============================================================================
cat("Figure 4: Forest Plot...\n")
# Strict: |logFC|>1, I2<50%
L_strict <- L[L$k >= 5 & L$padj < 0.05 & abs(L$pooled_logFC) > 1 & L$I2 < 50, ]
L_strict <- L_strict[!grepl(hk_pattern, L_strict$gene), ]
L_strict <- L_strict[order(L_strict$pooled_logFC), ]
L_strict$gene <- factor(L_strict$gene, levels = L_strict$gene)
L_strict$direction <- ifelse(L_strict$pooled_logFC > 0, "Up", "Down")

g <- ggplot(L_strict, aes(x = pooled_logFC, y = gene, color = direction)) +
  geom_point(aes(size = k), alpha = 0.85) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25, linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = ad_colors$muted, alpha = 0.65) +
  scale_color_manual(values = direction_palette[c("Up", "Down")], guide = "none") +
  scale_size_continuous(range = c(2.5, 6), name = "k (datasets)") +
  labs(x = expression(log[2]~"Fold Change (95% CI)"),
       y = "", title = "Lesional biomarker forest plot") +
  theme_pub +
  theme(axis.text.y = element_text(size = 9, face = "bold"),
        legend.position = "bottom")
save_pub_pair("plot_lesion_forest", g, width = 9, height = 7, dpi = 500)

# =============================================================================
# FIGURE 5: Lesion vs Disease Contrast Comparison
# =============================================================================
cat("Figure 5: Lesion vs Disease...\n")
Df <- D[D$k >= 5 & !grepl(hk_pattern, D$gene), ]
common <- merge(
  Lf[, c("gene", "pooled_logFC", "padj")],
  Df[, c("gene", "pooled_logFC", "padj")],
  by = "gene", suffixes = c("_L", "_D")
)
common$sig_type <- "Neither"
common$sig_type[common$padj_L < alpha & abs(common$pooled_logFC_L) >= fc_min &
                common$padj_D >= alpha] <- "Lesion only"
common$sig_type[common$padj_D < alpha & abs(common$pooled_logFC_D) >= fc_min &
                common$padj_L >= alpha] <- "Disease only"
common$sig_type[common$padj_L < alpha & abs(common$pooled_logFC_L) >= fc_min &
                common$padj_D < alpha & abs(common$pooled_logFC_D) >= fc_min] <- "Both"

g <- ggplot(common, aes(x = pooled_logFC_D, y = pooled_logFC_L, color = sig_type)) +
  geom_point(alpha = 0.58, size = 0.85) +
  scale_color_manual(values = c("Neither" = "#D7DCE2", "Lesion only" = ad_colors$lesion,
                                 "Disease only" = ad_colors$disease, "Both" = ad_colors$both),
                     name = "") +
  geom_hline(yintercept = 0, color = ad_colors$muted, alpha = 0.45) +
  geom_vline(xintercept = 0, color = ad_colors$muted, alpha = 0.45) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = ad_colors$muted, alpha = 0.4) +
  labs(x = expression(log[2]~"FC (Disease: NL vs HC)"),
       y = expression(log[2]~"FC (Lesion: LS vs NL)"),
       title = "Lesion and disease effect-size decomposition",
       subtitle = "Each point is a gene measured on both transcriptomic axes") +
  theme_pub
save_pub_pair("plot_lesion_vs_disease", g, width = 8.2, height = 7.2, dpi = 500)

# =============================================================================
# FIGURE 6: Lesional Directional Balance (combined strict + permissive)
# =============================================================================
cat("Figure 6: Directional Balance...\n")
# Strict: |logFC|>1, I2<50%
L_strict <- L[L$k >= 5 & L$padj < 0.05 & abs(L$pooled_logFC) > 1 & L$I2 < 50, ]
L_strict <- L_strict[!grepl(hk_pattern, L_strict$gene), ]
L_strict$direction <- ifelse(L_strict$pooled_logFC > 0, "Up", "Down")
dc_s <- as.data.frame(table(L_strict$direction))
names(dc_s) <- c("Direction", "Count")
dc_s$Threshold <- sprintf("Strict\n(%d genes)", nrow(L_strict))

# Permissive: |logFC|>0.5, I2<50%
L_perm <- L[L$k >= 5 & L$padj < 0.05 & abs(L$pooled_logFC) > 0.5 & L$I2 < 50, ]
L_perm <- L_perm[!grepl(hk_pattern, L_perm$gene), ]
L_perm$direction <- ifelse(L_perm$pooled_logFC > 0, "Up", "Down")
dc_p <- as.data.frame(table(L_perm$direction))
names(dc_p) <- c("Direction", "Count")
dc_p$Threshold <- sprintf("Permissive\n(%d genes)", nrow(L_perm))

dir_count <- rbind(dc_s, dc_p)
dir_count$Threshold <- factor(dir_count$Threshold,
  levels = c(sprintf("Strict\n(%d genes)", nrow(L_strict)),
             sprintf("Permissive\n(%d genes)", nrow(L_perm))))

g <- ggplot(dir_count, aes(x = Direction, y = Count, fill = Direction)) +
  geom_col(width = 0.54, color = "white", linewidth = 0.35) +
  scale_fill_manual(values = direction_palette[c("Up", "Down")], guide = "none") +
  geom_text(aes(label = Count), vjust = -0.55, size = 4.2, fontface = "bold",
            color = ad_colors$ink) +
  facet_wrap(~ Threshold, scales = "free_y") +
  labs(y = "Number of Significant Genes",
       x = NULL,
       title = "Directional balance of lesional expression changes",
       subtitle = "FDR<0.05, k>=5, no HK. Strict: |log2FC|>1, I2<50%. Permissive: |log2FC|>0.5, I2<50%.") +
  theme_pub +
  theme(strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "grey95"),
        plot.subtitle = element_text(size = 8, color = "grey40"))
save_pub_pair("plot_lesion_direction", g, width = 8, height = 5.5, dpi = 500)

# =============================================================================
# FIGURE 7: Consensus Distribution
# =============================================================================
cat("Figure 7: Consensus Distribution...\n")
if (exists("cs") && nrow(cs) > 0 && "n_sig" %in% names(cs)) {
  g <- ggplot(cs, aes(x = n_sig)) +
    geom_histogram(bins = 30, fill = ad_colors$disease, alpha = 0.88, color = NA) +
    labs(x = "Number of Significant Datasets", y = "Number of Genes",
         title = "Consensus disease-signature distribution",
         subtitle = "Higher values indicate replication across more datasets") +
    theme_pub
  save_pub_pair("plot_consensus_distribution", g, width = 7.2, height = 5.2, dpi = 500)
}

# =============================================================================
# FIGURE 8: LODO Sensitivity
# =============================================================================
cat("Figure 8: LODO Sensitivity...\n")
if (exists("lodo") && nrow(lodo) > 0) {
  g <- ggplot(lodo, aes(x = reorder(dataset, rank_corr),
                         y = rank_corr)) +
    geom_col(fill = ad_colors$disease, alpha = 0.88, color = "white", linewidth = 0.2) +
    geom_hline(yintercept = mean(lodo$rank_corr, na.rm = TRUE),
               linetype = "dashed", color = ad_colors$up) +
    coord_flip() +
    labs(x = "Dataset Removed", y = "Rank Correlation with Full Signature",
         title = "Leave-one-dataset-out sensitivity",
         subtitle = "Dashed line marks the mean rank correlation") +
    theme_pub + theme(axis.text.y = element_text(size = 7))
  save_pub_pair("plot_lodo_sensitivity", g, width = 9.2, height = 7.2, dpi = 500)
}

# =============================================================================
# FIGURE 9: Biomarker Heatmap
# =============================================================================
cat("Figure 9: Biomarker Heatmap...\n")
if (exists("lfc_mat")) {
  top_bm <- head(db[order(db$score, decreasing = TRUE), ], 50)
  bm_genes <- intersect(top_bm$gene, rownames(lfc_mat))
  if (length(bm_genes) >= 10) {
    bm_mat <- lfc_mat[bm_genes, , drop = FALSE]
    bm_mat <- bm_mat[, colSums(is.na(bm_mat)) < nrow(bm_mat) * 0.5, drop = FALSE]
    bm_mat[is.na(bm_mat)] <- 0
    pdf("plot_lesion_heatmap.pdf", width = 12, height = 10)
    pheatmap(bm_mat, scale = "none", cluster_rows = TRUE, cluster_cols = TRUE,
             color = heat_colors(100),
             main = "Top Disease Biomarkers Across Lesion Datasets",
             fontsize_row = 6, fontsize_col = 6,
             breaks = seq(-3, 3, length.out = 101))
    dev.off()
    png("plot_lesion_heatmap.png", width = 12, height = 10, units = "in", res = 300)
    pheatmap(bm_mat, scale = "none", cluster_rows = TRUE, cluster_cols = TRUE,
             color = heat_colors(100),
             main = "Top Disease Biomarkers Across Lesion Datasets",
             fontsize_row = 6, fontsize_col = 6,
             breaks = seq(-3, 3, length.out = 101))
    dev.off()
  }
}

# =============================================================================
# FIGURE 10: Co-expression Module Heatmap (if module data available)
# =============================================================================
cat("Figure 10: Co-expression Modules...\n")
if (!is.null(module_annotation) && nrow(module_annotation) > 0) {
  ma <- module_annotation
  ma$n_genes <- as.numeric(ma$n_genes)
  ma$module <- factor(ma$module, levels = ma$module[order(ma$n_genes)])
  ma$term_short <- scales::wrap_format(34)(ma$top_term)
  ma$module_term <- paste0(ma$module, "\n", ma$term_short)
  ma$module_term <- factor(ma$module_term, levels = ma$module_term[order(ma$n_genes)])
  ma$neglog_p <- pmin(-log10(pmax(ma$term_p, .Machine$double.xmin)), 80)
  g <- ggplot(ma, aes(x = n_genes, y = module_term, fill = module)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.25) +
    geom_text(aes(label = scales::comma(n_genes)), hjust = -0.12,
              size = 3, color = ad_colors$ink) +
    scale_fill_manual(values = module_palette, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    labs(x = "Genes assigned to module", y = NULL,
         title = "Annotated co-expression modules",
         subtitle = "Each module is labelled with its top enrichment term") +
    theme_ad_pub(base_size = 10) +
    theme(axis.text.y = element_text(face = "bold", size = 8.2))
  save_pub_pair("plot_coexpr_modules", g, width = 9.2, height = 7.6, dpi = 500)

  if (exists("mod_trait") && nrow(mod_trait) > 0) {
    mt <- merge(mod_trait, ma[, c("module", "top_term", "n_genes")],
                by = "module", all.x = TRUE)
    mt$module <- factor(mt$module, levels = mt$module[order(mt$cor)])
    mt$sig <- ifelse(mt$pval < 0.001, "FDR/P < 0.001",
              ifelse(mt$pval < 0.01, "P < 0.01",
              ifelse(mt$pval < 0.05, "P < 0.05", "NS")))
    mt$label <- paste0(mt$module, " (", scales::comma(mt$n_genes), ")")
    mt$label <- factor(mt$label, levels = mt$label[order(mt$cor)])
    g <- ggplot(mt, aes(x = cor, y = label, color = sig)) +
      geom_vline(xintercept = 0, color = ad_colors$muted, linewidth = 0.35) +
      geom_segment(aes(x = 0, xend = cor, yend = label),
                   linewidth = 1.0, alpha = 0.85) +
      geom_point(aes(size = -log10(pmax(pval, .Machine$double.xmin))),
                 alpha = 0.95) +
      scale_color_manual(values = c("FDR/P < 0.001" = ad_colors$both,
                                    "P < 0.01" = ad_colors$up,
                                    "P < 0.05" = ad_colors$gold,
                                    "NS" = "#AEB8C2"),
                         name = NULL) +
      scale_size_continuous(name = expression(-log[10](P)), range = c(2.5, 6.5)) +
      labs(x = "Correlation with lesion-effect PC1", y = NULL,
           title = "Module-trait association",
           subtitle = "Point size reflects association strength; labels include module size") +
      theme_ad_pub(base_size = 10, legend = "right") +
      theme(axis.text.y = element_text(face = "bold"))
    save_pub_pair("plot_coexpr_module_trait", g, width = 8.8, height = 6.8, dpi = 500)
  }
}

if (exists("mods") && nrow(mods) > 0 && exists("lfc_mat")) {
  mod_genes <- intersect(mods$gene, rownames(lfc_mat))
  if (length(mod_genes) > 100) {
    set.seed(42)
    mod_genes <- sample(mod_genes, min(500, length(mod_genes)))
  }
  if (length(mod_genes) >= 20) {
    mod_mat <- lfc_mat[mod_genes, , drop = FALSE]
    mod_mat <- mod_mat[, colSums(is.na(mod_mat)) < nrow(mod_mat) * 0.5, drop = FALSE]
    mod_mat[is.na(mod_mat)] <- 0
    annotation_row <- data.frame(Module = factor(mods$module[match(mod_genes, mods$gene)]),
                                 row.names = mod_genes)
    ann_colors <- list(Module = module_palette[intersect(names(module_palette),
                                                         levels(annotation_row$Module))])
    pdf("plot_coexpr_modules_heatmap.pdf", width = 14, height = 10)
    pheatmap(mod_mat, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
             annotation_row = annotation_row,
             annotation_colors = ann_colors,
             color = heat_colors(100),
             main = "Co-expression Module Profiles Across Datasets",
             fontsize_row = 3, fontsize_col = 6, show_rownames = FALSE)
    dev.off()
    png("plot_coexpr_modules_heatmap.png", width = 14, height = 10, units = "in", res = 300)
    pheatmap(mod_mat, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
             annotation_row = annotation_row,
             annotation_colors = ann_colors,
             color = heat_colors(100),
             main = "Co-expression Module Profiles Across Datasets",
             fontsize_row = 3, fontsize_col = 6, show_rownames = FALSE)
    dev.off()
  }
}

# =============================================================================
# FIGURE 10b: Hub Centrality and Hub Heatmap
# =============================================================================
cat("Figure 10b: Hub Gene Figures...\n")
if (exists("hubs") && nrow(hubs) > 0) {
  hub_dt <- hubs
  if (!"hub_score" %in% names(hub_dt)) {
    hub_dt$hub_score <- rank(-hub_dt$degree, ties.method = "min")
  }
  if (!"community" %in% names(hub_dt)) {
    hub_dt$community <- 1
  }
  hub_dt <- merge(hub_dt, cs[, c("gene", "n_sig")], by = "gene", all.x = TRUE)
  hub_dt <- merge(hub_dt, L[, c("gene", "pooled_logFC", "padj")], by = "gene", all.x = TRUE)
  hub_dt$n_sig[is.na(hub_dt$n_sig)] <- 0
  hub_dt$direction <- ifelse(hub_dt$pooled_logFC >= 0, "Up", "Down")
  hub_dt$effect_size <- abs(hub_dt$pooled_logFC)
  hub_top <- head(hub_dt[order(-hub_dt$n_sig, -hub_dt$effect_size, hub_dt$hub_score), ], 30)
  hub_top$gene <- factor(hub_top$gene, levels = rev(hub_top$gene))
  g <- ggplot(hub_top, aes(x = n_sig, y = gene)) +
    geom_segment(aes(x = 0, xend = n_sig, yend = gene),
                 color = "#CBD5E1", linewidth = 0.9) +
    geom_point(aes(size = effect_size, fill = direction),
               shape = 21, color = "white", stroke = 0.35, alpha = 0.95) +
    scale_fill_manual(values = direction_palette[c("Up", "Down")], name = "Lesion\ndirection") +
    scale_size_continuous(name = expression("|pooled log"[2]*"FC|"), range = c(2.8, 7.8)) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = "Consensus datasets with nominal significance",
         y = NULL,
         title = "Top hub genes in the lesion consensus network",
         subtitle = "Selected from the hub network; ranked by consensus support and lesion effect magnitude") +
    theme_ad_pub(base_size = 10, legend = "right") +
    theme(axis.text.y = element_text(face = "bold", size = 8))
  save_pub_pair("plot_hub_centrality", g, width = 8.8, height = 7.4, dpi = 500)

  if (exists("lfc_mat")) {
    hub_genes <- intersect(as.character(hub_top$gene), rownames(lfc_mat))
    if (length(hub_genes) >= 8) {
      hub_mat <- lfc_mat[hub_genes, , drop = FALSE]
      hub_mat <- hub_mat[, colSums(is.na(hub_mat)) < nrow(hub_mat) * 0.5, drop = FALSE]
      hub_mat[is.na(hub_mat)] <- 0
      row_ann <- data.frame(
        Community = factor(hub_top$community[match(hub_genes, as.character(hub_top$gene))]),
        n_sig = hub_top$n_sig[match(hub_genes, as.character(hub_top$gene))],
        row.names = hub_genes
      )
      pdf("plot_hub_heatmap.pdf", width = 10.5, height = 8.8)
      pheatmap(hub_mat, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
               annotation_row = row_ann,
               color = heat_colors(101),
               breaks = seq(-2.5, 2.5, length.out = 102),
               main = "Top Hub Gene Profiles Across Lesion Datasets",
               fontsize_row = 7.5, fontsize_col = 6)
      dev.off()
      png("plot_hub_heatmap.png", width = 10.5, height = 8.8, units = "in", res = 450)
      pheatmap(hub_mat, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
               annotation_row = row_ann,
               color = heat_colors(101),
               breaks = seq(-2.5, 2.5, length.out = 102),
               main = "Top Hub Gene Profiles Across Lesion Datasets",
               fontsize_row = 7.5, fontsize_col = 6)
      dev.off()
    }
  }
}

# =============================================================================
# FIGURE 11: Dataset Contributions
# =============================================================================
cat("Figure 11: Dataset Contributions...\n")
if (exists("lfc_mat")) {
  n_genes_per_ds <- colSums(!is.na(lfc_mat))
  ds_df <- data.frame(dataset = names(n_genes_per_ds), n_genes = n_genes_per_ds)
  ds_df <- ds_df[order(-ds_df$n_genes), ]
  ds_df$dataset <- factor(ds_df$dataset, levels = ds_df$dataset)
  g <- ggplot(ds_df, aes(x = dataset, y = n_genes / 1000)) +
    geom_col(fill = ad_colors$disease, alpha = 0.88, color = "white", linewidth = 0.2) +
    coord_flip() +
    labs(x = "", y = "Genes Tested (thousands)",
         title = "Gene coverage per dataset",
         subtitle = "Lesion log2FC matrix used for co-expression and heatmaps") +
    theme_pub + theme(axis.text.y = element_text(size = 6))
  save_pub_pair("plot_dataset_contributions", g, width = 8.2, height = 7.2, dpi = 500)
}

# =============================================================================
# FIGURE 12: Top Consensus Genes
# =============================================================================
cat("Figure 12: Top Consensus Genes...\n")
if (exists("cs") && nrow(cs) > 0) {
  top_cs <- head(cs[order(-cs$n_sig), ], 30)
  top_cs$gene <- factor(top_cs$gene, levels = rev(top_cs$gene))
  g <- ggplot(top_cs, aes(x = n_sig, y = gene)) +
    geom_col(fill = ad_colors$trait, alpha = 0.88, color = "white", linewidth = 0.2) +
    labs(x = "Number of Significant Datasets (out of 30)",
         y = "", title = "Top consensus disease-signature genes",
         subtitle = "Ranked by number of significant contributing datasets") +
    theme_pub + theme(axis.text.y = element_text(size = 8))
  save_pub_pair("plot_top30_consensus", g, width = 8.2, height = 7.2, dpi = 500)
}

# Legacy manuscript filenames retained as aliases to the refreshed figures.
alias_pairs <- list(
  c("plot_lesion_volcano", "plot_volcano"),
  c("plot_lesion_forest", "plot_forest"),
  c("plot_lesion_direction", "plot_directional_balance"),
  c("plot_lesion_i2_vs_fc", "plot_heterogeneity")
)
for (pair in alias_pairs) {
  for (ext in c("pdf", "png")) {
    src <- paste0(pair[1], ".", ext)
    dst <- paste0(pair[2], ".", ext)
    if (file.exists(src)) file.copy(src, dst, overwrite = TRUE)
  }
}

cat("\n=== All figures regenerated ===\n")
cat("Output directory:", out, "\n")
