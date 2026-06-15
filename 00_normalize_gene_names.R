#!/usr/bin/env Rscript
# =============================================================================
# normalize_gene_names.R — Retroactive probe→HGNC mapping for meta-analysis files
# =============================================================================
# Maps Affymetrix/Illumina probe IDs to HGNC gene symbols in meta_*_effects.csv
# and consensus_disease_signature.csv. For multi-probe→same-gene collisions,
# keeps the row with the largest |pooled_logFC|.
# =============================================================================

suppressMessages({
  library(hgu133plus2.db)
  library(hgu133a.db)
  library(org.Hs.eg.db)
})

out <- "d:/OneDrive/Statistics/Metagenomics/data/geo_atopic_dermatitis_longitudinal/disease_signature/analysis"

# ---- Build probe→symbol lookup table ----
# Use select() because as.list(hgu133plus2SYMBOL) keys by ENTREZID, not PROBEID
cat("Building probe-to-symbol lookup...\n")

build_map <- function(db, db_name) {
  tryCatch({
    # Get all PROBEID→SYMBOL mappings via select
    mapped <- select(db, keys=keys(db, "PROBEID"), columns="SYMBOL", keytype="PROBEID")
    # For multi-mapping probes, take the first (most common) symbol
    mapped <- mapped[!is.na(mapped$SYMBOL), ]
    mapped <- mapped[!duplicated(mapped$PROBEID), ]
    m <- setNames(mapped$SYMBOL, mapped$PROBEID)
    cat(sprintf("  %s: %d probes mapped\n", db_name, length(m)))
    return(m)
  }, error = function(e) {
    cat(sprintf("  %s: ERROR %s\n", db_name, e$message))
    return(c())
  })
}

hgu133p2_map <- build_map(hgu133plus2.db, "hgu133plus2")
hgu133a_map  <- build_map(hgu133a.db, "hgu133a")

# Combine: hgu133p2 first, hgu133a fallback
probe_to_symbol <- hgu133p2_map
for (probe in names(hgu133a_map)) {
  if (is.null(probe_to_symbol[[probe]]) || is.na(probe_to_symbol[[probe]])) {
    probe_to_symbol[[probe]] <- hgu133a_map[[probe]]
  }
}

cat(sprintf("  Combined lookup: %d probes\n", length(probe_to_symbol)))

# ---- Probe detection pattern ----
is_probe <- function(x) grepl("^[0-9]+(_[a-z])?_at$", x) | grepl("^ILMN_", x)

# ---- Normalize a single data frame ----
normalize_genes <- function(df, gene_col = "gene", effect_col = "pooled_logFC") {
  genes <- df[[gene_col]]
  probe_idx <- which(is_probe(genes))
  
  if (length(probe_idx) == 0) {
    cat(sprintf("  No probes found in %d rows\n", nrow(df)))
    return(df)
  }
  
  cat(sprintf("  Found %d probe IDs in %d rows (%.1f%%)\n",
              length(unique(genes[probe_idx])), nrow(df),
              100*length(probe_idx)/nrow(df)))
  
  # Map probes
  new_genes <- genes
  n_mapped <- 0L
  for (i in seq_along(probe_idx)) {
    p <- genes[probe_idx[i]]
    sym <- probe_to_symbol[p]  # single-bracket returns NA for missing names
    if (!is.na(sym)) {
      new_genes[probe_idx[i]] <- sym
      n_mapped <- n_mapped + 1L
    }
  }
  cat(sprintf("  Mapped %d/%d probes to gene symbols\n", n_mapped, length(probe_idx)))
  
  # Remove rows with unmapped probes (probes that stayed as probe IDs)
  still_probes <- which(is_probe(new_genes))
  if (length(still_probes) > 0) {
    cat(sprintf("  Removing %d rows with unmapped probes\n", length(still_probes)))
    df <- df[-still_probes, ]
    new_genes <- new_genes[-still_probes]
  }
  
  df$gene_original <- genes
  df$gene <- new_genes
  
  # Collapse multi-probe→same-gene: keep largest |effect|
  n_before <- nrow(df)
  df$abs_eff <- abs(df[[effect_col]])
  df <- df[order(-df$abs_eff), ]
  df <- df[!duplicated(df$gene), ]
  df$abs_eff <- NULL
  n_after <- nrow(df)
  cat(sprintf("  Collapsed %d -> %d rows (removed %d duplicate genes)\n",
              n_before, n_after, n_before - n_after))
  
  return(df)
}

# ---- Process meta_lesion_effects.csv ----
cat("\n=== Normalizing meta_lesion_effects.csv ===\n")
L <- read.csv(file.path(out, "meta_lesion_effects.csv"), stringsAsFactors=FALSE)
L_norm <- normalize_genes(L, "gene", "pooled_logFC")
write.csv(L_norm[, setdiff(names(L_norm), "gene_original")],
          file.path(out, "meta_lesion_effects.csv"), row.names=FALSE)
cat(sprintf("  Wrote %d rows\n", nrow(L_norm)))

# ---- Process meta_disease_effects.csv ----
cat("\n=== Normalizing meta_disease_effects.csv ===\n")
D <- read.csv(file.path(out, "meta_disease_effects.csv"), stringsAsFactors=FALSE)
D_norm <- normalize_genes(D, "gene", "pooled_logFC")
write.csv(D_norm[, setdiff(names(D_norm), "gene_original")],
          file.path(out, "meta_disease_effects.csv"), row.names=FALSE)
cat(sprintf("  Wrote %d rows\n", nrow(D_norm)))

# ---- Process consensus_disease_signature.csv ----
cat("\n=== Normalizing consensus_disease_signature.csv ===\n")
CS <- read.csv(file.path(out, "consensus_disease_signature.csv"), stringsAsFactors=FALSE)
CS_norm <- normalize_genes(CS, "gene", "n_sig")
write.csv(CS_norm[, setdiff(names(CS_norm), "gene_original")],
          file.path(out, "consensus_disease_signature.csv"), row.names=FALSE)
cat(sprintf("  Wrote %d rows\n", nrow(CS_norm)))

# ---- Process disease_biomarkers.csv ----
cat("\n=== Normalizing disease_biomarkers.csv ===\n")
DB <- read.csv(file.path(out, "disease_biomarkers.csv"), stringsAsFactors=FALSE)
DB_norm <- normalize_genes(DB, "gene", "pooled_logFC")
write.csv(DB_norm[, setdiff(names(DB_norm), "gene_original")],
          file.path(out, "disease_biomarkers.csv"), row.names=FALSE)
cat(sprintf("  Wrote %d rows\n", nrow(DB_norm)))

# ---- Quick verification ----
cat("\n=== Post-normalization verification ===\n")
Lv <- read.csv(file.path(out, "meta_lesion_effects.csv"), stringsAsFactors=FALSE)
Dv <- read.csv(file.path(out, "meta_disease_effects.csv"), stringsAsFactors=FALSE)
probes_L <- sum(is_probe(Lv$gene))
probes_D <- sum(is_probe(Dv$gene))
cat(sprintf("Lesion probes remaining: %d / %d (%.1f%%)\n",
            probes_L, nrow(Lv), 100*probes_L/nrow(Lv)))
cat(sprintf("Disease probes remaining: %d / %d (%.1f%%)\n",
            probes_D, nrow(Dv), 100*probes_D/nrow(Dv)))

# Check specific examples
cat("\nExample mappings:\n")
examples <- c("221619_s_at","204580_at","209125_at","203476_at","214580_x_at")
for (ex in examples) {
  sym <- probe_to_symbol[ex]
  if (is.na(sym)) sym <- "NOT MAPPED"
  cat(sprintf("  %s -> %s\n", ex, sym))
}

cat("\nDone. Original gene names preserved in gene_original column (removed from output).\n")
