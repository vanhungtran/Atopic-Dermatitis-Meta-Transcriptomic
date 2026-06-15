#!/usr/bin/env Rscript
# =============================================================================
# 17_module_enrichment_fill.R
# Fill in functional enrichment for co-expression modules that are missing an
# `enrichment_module_<M>.csv`. The figure script (16) only annotates modules
# that already have such a file; the rest fall back to "No enrichment
# annotation". This computes g:Profiler enrichment for the missing modules so
# every module in plot_coexpr_modules gets a real top term where one exists.
#
# Input:  coexpr_gene_modules.csv (cols: gene, module)
# Output: enrichment_module_<M>.csv for each previously-missing module,
#         matching the column layout of the existing enrichment files.
# After running this, re-run scripts/16_regenerate_all_figures.R to rebuild
# coexpr_module_annotation.csv and plot_coexpr_modules.{png,pdf}.
# =============================================================================

suppressPackageStartupMessages({
  library(gprofiler2)
})

out <- "d:/OneDrive/Statistics/Metagenomics/data/geo_atopic_dermatitis_longitudinal/disease_signature/analysis"

gm <- read.csv(file.path(out, "coexpr_gene_modules.csv"), stringsAsFactors = FALSE)
modules <- sort(unique(gm$module))

# Same data sources and ordering used for the already-annotated modules.
src_keep <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP")
enr_cols <- c("term_id", "term_name", "source", "term_size",
              "query_size", "intersection_size", "p_value",
              "precision", "recall")

for (m in modules) {
  f <- file.path(out, paste0("enrichment_module_", m, ".csv"))
  if (file.exists(f)) {
    cat(sprintf("%s: enrichment file already present, skipping.\n", m))
    next
  }
  genes <- gm$gene[gm$module == m]
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  cat(sprintf("%s: %d genes -> querying g:Profiler...\n", m, length(genes)))

  res <- tryCatch(
    gost(query = genes, organism = "hsapiens",
         sources = src_keep, significant = TRUE,
         user_threshold = 0.05, correction_method = "fdr",
         evcodes = FALSE),
    error = function(e) { cat("  gost error:", conditionMessage(e), "\n"); NULL }
  )

  if (is.null(res) || is.null(res$result) || nrow(res$result) == 0) {
    # Write an empty (header-only) file so the state is explicit and the figure
    # script's file.exists() check still treats this as "genuinely no terms".
    empty <- setNames(data.frame(matrix(ncol = length(enr_cols), nrow = 0)),
                      enr_cols)
    write.csv(empty, f, row.names = FALSE)
    cat(sprintf("  %s: no significant enrichment terms.\n", m))
    next
  }

  r <- res$result
  # Order by ascending p_value so row 1 is the top term (matches existing files).
  r <- r[order(r$p_value), ]
  r <- r[, enr_cols]
  write.csv(r, f, row.names = FALSE)
  cat(sprintf("  %s: %d terms, top = %s (%s, p=%.2g)\n",
              m, nrow(r), r$term_name[1], r$source[1], r$p_value[1]))
}

cat("\nDone. Now re-run scripts/16_regenerate_all_figures.R to refresh the",
    "annotation table and the module plot.\n")
