#!/usr/bin/env Rscript
# Reviewer-requested checks: disease contrast gene coverage, MSH5-SAPCD1 RNA-seq breakdown
setwd("d:/OneDrive/Statistics/Metagenomics/data/geo_atopic_dermatitis_longitudinal/disease_signature/analysis")

cat("=== 1. Disease contrast k distribution ===\n")
meta_dz <- read.csv("meta_disease_effects.csv")
kt <- aggregate(list(n = meta_dz$gene), by = list(k = meta_dz$k), FUN = length)
kt <- kt[order(-kt$n), ]
print(kt[1:15, ])
cat(sprintf("\nGenes with k>=5: %d / %d total\n", sum(meta_dz$k >= 5, na.rm = TRUE), nrow(meta_dz)))

cat("\n=== 2. Disease contrast per-dataset coverage ===\n")
logFC_dz <- read.csv("logFC_disease_matrix.csv", row.names = 1)
SE_dz    <- read.csv("SE_disease_matrix.csv", row.names = 1)
coverage <- data.frame(
  dataset = colnames(logFC_dz),
  genes_measured = colSums(!is.na(logFC_dz)),
  stringsAsFactors = FALSE
)
coverage <- coverage[order(-coverage$genes_measured), ]
print(coverage)

cat("\n=== 3. Disease contrast gene overlap ===\n")
# genes present in >=5 disease datasets
gene_presence <- rowSums(!is.na(logFC_dz))
genes_k5 <- names(gene_presence[gene_presence >= 5])
cat(sprintf("Genes in >=5 disease datasets: %d\n", length(genes_k5)))

cat("\n=== 4. MSH5-SAPCD1 RNA-seq vs microarray ===\n")
meta_les <- read.csv("meta_lesion_effects.csv")
msh5 <- meta_les[meta_les$gene == "MSH5-SAPCD1", ]
if (nrow(msh5) > 0) {
  cat(sprintf("MSH5-SAPCD1 lesion: log2FC=%.3f, p=%.2e, padj=%.2e, I2=%.1f%%, k=%d\n",
              msh5$pooled_logFC[1], msh5$pval[1], msh5$padj[1], msh5$I2[1], msh5$k[1]))
}

# Check RNA-seq datasets
ds_annot <- read.csv("dataset_contrast_annotation.csv")
rna_seq_ds <- ds_annot$GEO[grepl("RNA-seq", ds_annot$Platform, ignore.case = TRUE)]
cat(sprintf("\nRNA-seq datasets (n=%d): %s\n", length(rna_seq_ds), paste(rna_seq_ds, collapse = ", ")))

# Per-dataset logFC for MSH5-SAPCD1 in lesion contrast
logFC_les <- read.csv("logFC_lesion_matrix.csv", row.names = 1)
if ("MSH5-SAPCD1" %in% rownames(logFC_les)) {
  msh5_fc <- as.numeric(logFC_les["MSH5-SAPCD1", ])
  names(msh5_fc) <- colnames(logFC_les)
  msh5_fc <- msh5_fc[!is.na(msh5_fc)]
  cat(sprintf("\nMSH5-SAPCD1 per-dataset log2FC (n=%d datasets):\n", length(msh5_fc)))
  for (nm in names(msh5_fc)) {
    is_rna <- grepl("RNA-seq", ds_annot$Platform[ds_annot$GEO == nm], ignore.case = TRUE)
    cat(sprintf("  %-20s  log2FC=%+.3f  %s\n", nm, msh5_fc[nm], ifelse(length(is_rna)>0 && is_rna, "[RNA-seq]", "[microarray]")))
  }
}

cat("\nDone.\n")
