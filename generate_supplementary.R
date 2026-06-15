#!/usr/bin/env Rscript
# =============================================================================
# generate_supplementary.R — Generate Supplementary Tables S1 and S2
# =============================================================================

suppressMessages({
  library(data.table)
})

out <- "d:/OneDrive/Statistics/Metagenomics/data/geo_atopic_dermatitis_longitudinal/disease_signature/analysis"

L <- fread(file.path(out, "meta_lesion_effects.csv"))
D <- fread(file.path(out, "meta_disease_effects.csv"))
db <- fread(file.path(out, "disease_biomarker_signature_final.csv"))
dec <- fread(file.path(out, "ad_ls_gene_decomposition.csv"))
ds_ann <- fread(file.path(out, "dataset_contrast_annotation.csv"))

hk_pattern <- "^(RPS|RPL|RPLP|MRPS|MRPL|EEF|EIF|ACTB|ACTG1|GAPDH|TUBB|TUBA|B2M|HPRT1|PPIA|UBC|YWHAZ)"

# ===========================================================================
# SUPPLEMENTARY TABLE S1: Full lesion biomarker results
# ===========================================================================
cat("Generating Supplementary Table S1...\n")

# Strict biomarkers: |logFC|>1, I2<50%, k>=5, padj<0.05, no HK
lb_strict <- L[L$k >= 5 & L$padj < 0.05 & abs(L$pooled_logFC) > 1 & L$I2 < 50]
lb_strict <- lb_strict[!grepl(hk_pattern, lb_strict$gene)]
lb_strict$direction <- ifelse(lb_strict$pooled_logFC > 0, "Up", "Down")
lb_strict <- lb_strict[order(-abs(lb_strict$pooled_logFC))]

# Permissive: |logFC|>0.5, I2<50%, k>=5, padj<0.05, no HK
lb_perm <- L[L$k >= 5 & L$padj < 0.05 & abs(L$pooled_logFC) > 0.5 & L$I2 < 50]
lb_perm <- lb_perm[!grepl(hk_pattern, lb_perm$gene)]
lb_perm$direction <- ifelse(lb_perm$pooled_logFC > 0, "Up", "Down")
lb_perm <- lb_perm[order(-abs(lb_perm$pooled_logFC))]

# Disease-significant genes (no HK)
dsig <- D[D$padj < 0.05 & D$k >= 5]
dsig <- dsig[!grepl(hk_pattern, dsig$gene)]
dsig$direction <- ifelse(dsig$pooled_logFC > 0, "Up", "Down")
dsig <- dsig[order(dsig$padj)]

# Combine into one supplementary table
s1_strict <- lb_strict[, .(
  gene, pooled_logFC, se, ci_lower, ci_upper, pval, padj, I2, tau2, k, direction,
  rank_strict = .I
)]

s1_perm <- lb_perm[, .(
  gene, pooled_logFC, se, ci_lower, ci_upper, pval, padj, I2, tau2, k, direction,
  rank_permissive = .I
)]

s1_disease <- dsig[, .(
  gene, logFC_disease = pooled_logFC, se_disease = se, padj_disease = padj,
  I2_disease = I2, k_disease = k, direction_disease = direction
)]

# Merge permissive ranks with disease info
s1 <- merge(s1_perm, s1_disease, by = "gene", all.x = TRUE, suffixes = c("", ".d"))
# Add strict biomarker flag
s1$strict_biomarker <- s1$gene %in% lb_strict$gene
s1 <- s1[order(-abs(pooled_logFC))]

fwrite(s1, file.path(out, "Supplementary_Table_S1.csv"))
cat(sprintf("  S1: %d rows (perm) / %d strict biomarkers / %d disease-sig\n",
            nrow(s1_perm), nrow(lb_strict), nrow(dsig)))

# ===========================================================================
# SUPPLEMENTARY TABLE S2: Decomposition, Drug Repurposing, k-Sensitivity
# ===========================================================================
cat("Generating Supplementary Table S2...\n")

# Sheet A: AD/LS decomposition gene lists
dec_out <- dec[, .(gene, group, beta_D, beta_L, padj_D, padj_L, I2_D, I2_L, k_D, k_L)]
setnames(dec_out, c("gene","decomposition_group","logFC_disease","logFC_lesion",
                     "padj_disease","padj_lesion","I2_disease","I2_lesion",
                     "k_disease","k_lesion"))

# Sheet B: Drug repurposing results
dr <- fread(file.path(out, "drug_repurposing_results.csv"))

# Sheet C: Disease biomarker signature
db_out <- db[, .(gene, pooled_logFC, padj, I2, k, direction, score, rank = .I)]
setnames(db_out, c("gene","logFC","padj","I2","k","direction","biomarker_score","biomarker_rank"))

# Write as multi-sheet CSV (use a marker row)
fwrite(dec_out, file.path(out, "Supplementary_Table_S2_decomposition.csv"))
fwrite(dr, file.path(out, "Supplementary_Table_S2_drug_repurposing.csv"))
fwrite(db_out, file.path(out, "Supplementary_Table_S2_biomarkers.csv"))

cat(sprintf("  S2: %d decomposition + %d drug pairs + %d disease biomarkers\n",
            nrow(dec_out), nrow(dr), nrow(db_out)))

# ===========================================================================
# Also generate a combined Supplementary PDF with key tables
# ===========================================================================
cat("Generating supplementary PDF placeholder...\n")
# Write a simple combined CSV that can be the "Supplementary File 2"
s2_combined <- rbind(
  data.table(section = "=== AD/LS DECOMPOSITION ===", content = sprintf("%d genes total (%d AD-related, %d AD-and-LS, %d LS-related)",
              nrow(dec), sum(dec$group=="AD_related"), sum(dec$group=="AD_and_LS"), sum(dec$group=="LS_related"))),
  data.table(section = "AD-related genes", content = paste(sort(dec$gene[dec$group=="AD_related"]), collapse=", ")),
  data.table(section = "AD-and-LS genes", content = paste(sort(dec$gene[dec$group=="AD_and_LS"]), collapse=", ")),
  data.table(section = "", content = ""),
  data.table(section = "=== DISEASE BIOMARKER SIGNATURE ===", content = sprintf("%d genes (%d up, %d down)", nrow(db), sum(db$direction=="up"), sum(db$direction=="down"))),
  data.table(section = "Top 20 biomarkers", content = paste(head(db$gene, 20), collapse=", ")),
  data.table(section = "", content = ""),
  data.table(section = "=== DRUG REPURPOSING ===", content = sprintf("%d drug-gene pairs across %d unique genes", nrow(dr), length(unique(dr$gene)))),
  data.table(section = "Drug categories", content = paste(unique(dr$category), collapse=", "))
)
fwrite(s2_combined, file.path(out, "Supplementary_File_2_Summary.csv"))

cat("\n=== Supplementary files generated ===\n")
cat("Supplementary_Table_S1.csv — Full lesion biomarker results with disease-axis annotations\n")
cat("Supplementary_Table_S2_decomposition.csv — AD/LS decomposition gene lists\n")
cat("Supplementary_Table_S2_drug_repurposing.csv — Drug repurposing results\n")
cat("Supplementary_Table_S2_biomarkers.csv — Disease biomarker signature\n")
cat("Supplementary_File_2_Summary.csv — Combined summary\n")
cat("\nDone.\n")
