#!/usr/bin/env Rscript
# =============================================================================
# k_sensitivity.R — Sensitivity of results to k threshold (k≥5 vs k≥4)
# =============================================================================
out <- "d:/OneDrive/Statistics/Metagenomics/data/geo_atopic_dermatitis_longitudinal/disease_signature/analysis"

hk_pattern <- "^(RPS|RPL|RPLP|MRPS|MRPL|EEF|EIF|ACTB|ACTG1|GAPDH|TUBB|TUBA|B2M|HPRT1|PPIA|UBC|YWHAZ)[0-9A-Z]*$"
alpha <- 0.05
fc_min <- 0.5

# ---- Load data ----
D <- read.csv(file.path(out, "meta_disease_effects.csv"), stringsAsFactors=FALSE)
L <- read.csv(file.path(out, "meta_lesion_effects.csv"), stringsAsFactors=FALSE)

# ---- Filter function ----
run_at_k <- function(D, L, kmin) {
  Df <- D[!grepl(hk_pattern, D$gene) & D$k >= kmin & !is.na(D$padj) & D$padj < alpha & abs(D$pooled_logFC) >= fc_min, ]
  Lf <- L[!grepl(hk_pattern, L$gene) & L$k >= kmin & !is.na(L$padj) & L$padj < alpha & abs(L$pooled_logFC) >= fc_min, ]
  list(D = Df, L = Lf, D_genes = unique(Df$gene), L_genes = unique(Lf$gene))
}

r5 <- run_at_k(D, L, 5)
r4 <- run_at_k(D, L, 4)

# ---- GWAS genes of interest ----
gwas <- c("IL13","IL4R","IL4","IL5","FLG","TSLP","STAT3","GATA3",
          "IL1RL1","IL18R1","IL18RAP","RAD50","KIF3A","OVOL1",
          "CLEC16A","LRRC32","IL2RA","IL7R","IL6R","CARD11","IKZF3")

# ---- Build comparison table ----
cat("=== GWAS genes: k≥5 vs k≥4 ===\n")
cat(sprintf("%-10s %6s %6s %10s %10s %10s %10s\n",
            "Gene","k_D5","k_D4","padj_D5","padj_D4","log2FC_D","I2_D"))
for(g in gwas) {
  d5 <- D[D$gene == g & D$k >= 5, ]
  d4 <- D[D$gene == g & D$k >= 4, ]
  sig5 <- if(nrow(d5)>0 && d5$padj < alpha && abs(d5$pooled_logFC) >= fc_min) "*" else ""
  sig4 <- if(nrow(d4)>0 && d4$padj < alpha && abs(d4$pooled_logFC) >= fc_min) "*" else ""
  k5 <- if(nrow(d5)>0) d5$k[1] else NA
  k4 <- if(nrow(d4)>0) d4$k[1] else NA
  p5 <- if(nrow(d5)>0) sprintf("%.2e",d5$padj[1]) else "NA"
  p4 <- if(nrow(d4)>0) sprintf("%.2e",d4$padj[1]) else "NA"
  fc <- if(nrow(d4)>0) sprintf("%.3f",d4$pooled_logFC[1]) else "NA"
  i2 <- if(nrow(d4)>0) sprintf("%.1f%%",d4$I2[1]) else "NA"
  status <- ""
  if(is.na(k5) && !is.na(k4)) status <- "NEW at k≥4"
  else if(!is.na(k5) && sig5=="" && sig4=="*") status <- "SIG at k≥4"
  else if(sig5=="*") status <- "OK at k≥5"
  cat(sprintf("%-10s %5s%s %5s%s %10s %10s %10s %10s  %s\n",
              g, ifelse(is.na(k5),"-",k5), sig5,
              ifelse(is.na(k4),"-",k4), sig4, p5, p4, fc, i2, status))
}

# ---- Disease-axis genes gained at k≥4 ----
cat("\n=== Disease-axis genes gained at k≥4 ===\n")
d4_new <- setdiff(r4$D_genes, r5$D_genes)
cat(sprintf("k≥5 disease genes: %d  |  k≥4 disease genes: %d  |  Gained: %d\n",
            length(r5$D_genes), length(r4$D_genes), length(d4_new)))
if(length(d4_new) > 0) {
  new_info <- D[D$gene %in% d4_new, ]
  new_info <- new_info[order(new_info$padj), ]
  cat("\nTop gained disease-axis genes (by padj):\n")
  print(head(new_info[, c("gene","pooled_logFC","padj","I2","k")], 20))
}

# ---- Trait/state changes ----
cat("\n=== Trait/state classification: genes that change at k≥4 ===\n")
# Which genes that are ns at k≥5 become classified at k≥4?
ts5 <- read.csv(file.path(out, "causal_trait_state_classification.csv"), stringsAsFactors=FALSE)
Df5 <- D[!grepl(hk_pattern, D$gene) & D$k >= 5, c("gene","pooled_logFC","padj","k")]
Lf5 <- L[!grepl(hk_pattern, L$gene) & L$k >= 5, c("gene","pooled_logFC","padj","k")]
Df4 <- D[!grepl(hk_pattern, D$gene) & D$k >= 4, c("gene","pooled_logFC","padj","k")]
Lf4 <- L[!grepl(hk_pattern, L$gene) & L$k >= 4, c("gene","pooled_logFC","padj","k")]

# Genes that were missing from D-axis at k≥5 but present at k≥4
D_new <- setdiff(Df4$gene[Df4$padj < alpha & abs(Df4$pooled_logFC) >= fc_min],
                 Df5$gene[Df5$padj < alpha & abs(Df5$pooled_logFC) >= fc_min])
cat(sprintf("D-axis sig genes gained at k≥4: %d\n", length(D_new)))
if(length(D_new) > 0) {
  gwas_in_new <- intersect(D_new, gwas)
  if(length(gwas_in_new) > 0) cat("GWAS genes among them:", paste(gwas_in_new, collapse=", "), "\n")
  other_new <- setdiff(D_new, gwas)
  cat("Top other new D-axis genes:\n")
  D_ranks <- D[D$gene %in% other_new & D$padj < alpha & abs(D$pooled_logFC) >= fc_min, ]
  D_ranks <- D_ranks[order(D_ranks$padj), ]
  print(head(D_ranks[, c("gene","pooled_logFC","padj","I2","k")], 15))
}

cat("\nDone.\n")
