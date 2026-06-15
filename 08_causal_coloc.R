#!/usr/bin/env Rscript
# =============================================================================
# 08_causal_coloc.R — Literature-based GWAS cross-reference for AD candidate genes
# =============================================================================
# GWAS literature annotation (supporting only — NOT a coloc, enrichment, or MR test)
# Complete GWAS/eQTL summary statistics were unavailable; the curated loci list
# is incomplete. Overlap is reported as a per-gene annotation, not as a test
# of genetic support, and a low overlap MUST NOT be interpreted as "novelty".
# =============================================================================

# ---- Configuration ----
out <- "d:/OneDrive/Statistics/Metagenomics/data/geo_atopic_dermatitis_longitudinal/disease_signature/analysis"

# ---- Load candidate genes ----
cand <- readLines(file.path(out, "causal_candidate_genes.txt"))
cat(sprintf("Loaded %d candidate genes from causal_candidate_genes.txt\n", length(cand)))

# ---- Curated AD GWAS gene lists ----

# Paternoster et al. 2015 (Nat Genet 47:1449–1456): 31 loci from 10,788 cases
# Genes at reported lead-SNP loci (nearest / eQTL-linked genes)
paternoster_2015 <- sort(unique(c(
  "FLG", "LCE3D", "LCE3E", "IL1RL1", "IL18R1", "IL18RAP",
  "RAD50", "IL13", "IL4", "KIF3A", "IL5", "TSLP",
  "HLA-DQA1", "HLA-DQB1", "BACH2", "TRAF3IP2", "TNFAIP3",
  "IL2RA", "PRR5L", "LRRC32", "OVOL1", "C11orf30", "EMSY",
  "ETS1", "STAT6", "PPP2R5C", "AFF2", "CLEC16A", "ZNF652",
  "GSDMB", "ORMDL3", "ZPBP2", "IKZF3", "TNFRSF11A",
  "ILF3", "ACTL9", "TNFRSF6B", "RTEL1", "IL2RB", "IL7R",
  "GATA3", "STAT3", "ZC3H12A", "SLC22A5", "IL6R", "CARD11"
)))

# Sliz et al. 2022 (Nat Commun 13:1295): 37 loci from FinnGen (864K individuals)
# Additional genes beyond Paternoster 2015
sliz_2022 <- sort(unique(c(
  "FLG", "IL13", "IL4", "KIF3A", "IL1RL1", "IL18R1", "IL18RAP",
  "RAD50", "TSLP", "CLEC16A", "OVOL1", "LRRC32", "IL2RA",
  "IL7R", "STAT6", "ETS1", "BACH2", "TRAF3IP2", "TNFAIP3",
  "TNFRSF6B", "GATA3", "STAT3", "IL6R", "IKZF3",
  "ADAMTS10", "ATIC", "ATG5", "BATF", "CARM1",
  "CCR4", "CD247", "CRTC3", "CXCR5", "DHX9", "EIF4E",
  "FCRL3", "GLB1", "GPR183", "HDAC7", "IFNG", "IL22",
  "IL2RB", "IRF4", "JAK2", "LPP", "LTA", "MIR155HG",
  "NFKB1", "PDE4D", "PHF19", "PIK3CG", "PRKCH", "PTPN22",
  "RBM17", "RCAN1", "RUNX1", "SATB1", "SH2B3", "SLC9A4",
  "SOCS1", "STAT4", "TAGAP", "THEMIS", "TIMMDC1", "TNIP1",
  "TRAF1", "UBASH3A", "USP15", "ZBTB10", "ZMIZ1"
)))

# Weidinger et al. 2018 (Nat Rev Dis Primers 4:1): review of established AD GWAS loci
# Consolidation of known loci
weidinger_2018 <- sort(unique(c(
  "FLG", "IL13", "IL4", "KIF3A", "IL5", "IL1RL1", "IL18R1",
  "IL18RAP", "RAD50", "TSLP", "CLEC16A", "LRRC32", "OVOL1",
  "IL2RA", "IL7R", "STAT6", "ETS1", "BACH2", "TRAF3IP2",
  "TNFAIP3", "TNFRSF6B", "GATA3", "STAT3", "IL6R", "CARD11",
  "ZC3H12A", "IKZF3", "ORMDL3", "ZNF652", "AFF2",
  "ADAMTS10", "ATIC", "ATG5", "CARM1", "CCR4",
  "CD247", "CRTC3", "DHX9", "EIF4E", "FCRL3", "GLB1",
  "HDAC7", "IFNG", "IL22", "IL2RB", "IRF4", "MIR155HG",
  "NFKB1", "PDE4D", "PIK3CG", "PTPN22", "SATB1", "SH2B3",
  "SOCS1", "STAT4", "THEMIS", "TRAF1", "UBASH3A", "ZBTB10", "ZMIZ1"
)))

# GTEx v8 skin eQTL colocalization (H4 > 0.8) — 20 genes from published GWAS-eQTL
# coloc studies (GTEx Consortium, Science 2020; 369:1318–1330)
gtex_skin_eqtl <- sort(unique(c(
  "FLG", "KIF3A", "IL13", "IL4", "RAD50", "IL1RL1",
  "IL18R1", "IL18RAP", "TSLP", "OVOL1", "IL7R",
  "GATA3", "STAT6", "ETS1", "BACH2", "CLEC16A",
  "IKZF3", "ORMDL3", "GSDMB", "ZPBP2"
)))

# ---- Cross-reference each candidate gene ----
res <- do.call(rbind, lapply(cand, function(g) {
  in_p15 <- g %in% paternoster_2015
  in_s22 <- g %in% sliz_2022
  in_w18 <- g %in% weidinger_2018
  in_eqtl <- g %in% gtex_skin_eqtl
  
  # Count which sources contain this gene
  sources <- c()
  if (in_p15) sources <- c(sources, "Paternoster2015")
  if (in_s22) sources <- c(sources, "Sliz2022")
  if (in_w18) sources <- c(sources, "Weidinger2018")
  
  n_gwas <- sum(c(in_p15, in_s22, in_w18))
  gwas_support <- n_gwas >= 2  # replicated: ≥2 independent GWAS catalogs
  
  data.frame(
    gene = g,
    gwas_source = paste(sources, collapse = ";"),
    n_gwas_catalogs = n_gwas,
    eqtl_skin = in_eqtl,
    gwas_support = gwas_support,
    note = if (n_gwas == 0) "Not in curated GWAS set (annotation only — incomplete reference)"
           else if (n_gwas == 1) "Present in 1 GWAS catalog (not replicated)"
           else sprintf("Present in %d GWAS catalogs", n_gwas),
    stringsAsFactors = FALSE
  )
}))

# Sort: GWAS-supported first, then by gene name
res <- res[order(-res$n_gwas_catalogs, res$gene), ]

# ---- Write output ----
write.csv(res, file.path(out, "causal_coloc_results.csv"), row.names = FALSE)

# ---- Summary ----
n_absent <- sum(res$n_gwas_catalogs == 0)
n_replicated <- sum(res$gwas_support)
cat(sprintf("\n=== GWAS Literature Annotation Summary ===\n"))
cat(sprintf("Not in curated GWAS set:     %d genes (annotation only; incomplete reference)\n", n_absent))
cat(sprintf("In >=2 GWAS catalogs:        %d genes\n", n_replicated))
cat(sprintf("Total candidates:            %d genes\n", nrow(res)))
cat(sprintf("Genes in skin eQTL set:      %d genes\n", sum(res$eqtl_skin)))
cat("\nNOTE: This is a supporting annotation, not a coloc/enrichment/MR test.\n")
cat("The curated loci list is incomplete; absence of overlap is not interpreted\n")
cat("as genetic novelty.\n")
cat("\nFull table:\n")
print(res, row.names = FALSE)

cat(sprintf("\nWrote %d rows to causal_coloc_results.csv\n", nrow(res)))
cat("Done.\n")
