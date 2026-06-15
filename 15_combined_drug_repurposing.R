#!/usr/bin/env Rscript
# =============================================================================
# 15_combined_drug_repurposing.R -- two-tier drug repurposing analysis
# =============================================================================
# Tier 1: Druggability annotation of AD/lesion/both biomarker groups using
# curated interactions plus live Open Targets tractability and DGIdb GraphQL
# annotation where available.
#
# Tier 2: LINCS/CMap signature-reversal analysis through signatureSearch. The
# script runs against an ExperimentHub reference database when available and
# writes explicit method status rows otherwise.
# =============================================================================

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(httr2)
  library(jsonlite)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

out <- "d:/OneDrive/Statistics/Metagenomics/data/geo_atopic_dermatitis_longitudinal/disease_signature/analysis"
source(file.path(out, "scripts", "plotting_utils.R"))

alpha <- 0.05
fc_min <- 0.5
known_ad_drugs <- c("dupilumab", "tralokinumab", "lebrikizumab",
                    "abrocitinib", "upadacitinib", "baricitinib",
                    "ruxolitinib", "crisaborole", "tacrolimus",
                    "pimecrolimus", "delgocitinib", "hydrocortisone",
                    "dexamethasone", "prednisolone", "methotrexate",
                    "ciclosporin", "cyclosporine", "mycophenolic acid",
                    "azathioprine", "sirolimus")

normalize_perturbagen <- function(x) {
  tolower(sub("__.*$", "", as.character(x)))
}

safe_chr <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_character_)
  paste(sort(unique(na.omit(as.character(x)))), collapse = "; ")
}

stage_is_clinical <- function(stage) {
  grepl("Approved|Phase", stage, ignore.case = TRUE)
}

stage_is_approved <- function(stage) {
  grepl("Approved", stage, ignore.case = TRUE)
}

read_axis_table <- function() {
  dec <- fread(file.path(out, "ad_ls_gene_decomposition.csv"))
  dec <- dec[, .(gene, target_axis = group, beta_D, padj_D, beta_L, padj_L)]
  dec_core <- dec[target_axis %in% c("AD_related", "AD_and_LS")]

  db <- fread(file.path(out, "disease_biomarker_signature_final.csv"))
  db <- db[, .(gene, target_axis = "AD_related",
               beta_D = pooled_logFC, padj_D = padj,
               beta_L = NA_real_, padj_L = NA_real_)]

  lb <- fread(file.path(out, "lesion_biomarkers.csv"))
  lb <- lb[, .(gene, target_axis = "LS_related",
               beta_D = NA_real_, padj_D = NA_real_,
               beta_L = pooled_logFC, padj_L = padj)]

  axes <- rbindlist(list(dec_core, db, lb), fill = TRUE)
  axes[, target_axis := fifelse(gene %in% dec_core[target_axis == "AD_and_LS", gene],
                                "AD_and_LS", target_axis)]
  axes <- axes[, .(
    target_axis = if ("AD_and_LS" %in% target_axis) "AD_and_LS" else target_axis[1],
    beta_D = beta_D[which.max(!is.na(beta_D))],
    padj_D = padj_D[which.max(!is.na(padj_D))],
    beta_L = beta_L[which.max(!is.na(beta_L))],
    padj_L = padj_L[which.max(!is.na(padj_L))]
  ), by = gene]
  axes[]
}

axis_map <- read_axis_table()

# ---- Curated drug-target interactions retained from the previous analysis ----
drug_map <- rbindlist(list(
  data.table(gene = "CCL17", axis = "Disease", logFC = 1.434, padj = 6.10e-06,
             drug = "Mogamulizumab (anti-CCR4)", category = "Immunotherapy", stage = "Phase III"),
  data.table(gene = "GATA3", axis = "Disease", logFC = -0.167, padj = 4.83e-02,
             drug = "Corticosteroids; JAK inhibitors", category = "Anti-inflammatory", stage = "Approved"),
  data.table(gene = "GJB2", axis = "Disease", logFC = 0.578, padj = 4.78e-11,
             drug = "Tonabersat (connexin blocker)", category = "Ion channel", stage = "Phase II"),
  data.table(gene = "CLDN1", axis = "Disease", logFC = -0.489, padj = 2.04e-07,
             drug = "Corticosteroids; Retinoids", category = "Barrier repair", stage = "Approved"),
  data.table(gene = "TNFSF10", axis = "Disease", logFC = 0.554, padj = 1.65e-02,
             drug = "TRAIL-R agonists", category = "Immunotherapy", stage = "Phase I/II"),
  data.table(gene = "TNC", axis = "Disease", logFC = 0.851, padj = 9.50e-04,
             drug = "Anti-tenascin antibodies", category = "Anti-fibrotic", stage = "Preclinical"),
  data.table(gene = "BGN", axis = "Disease", logFC = -0.443, padj = 1.90e-04,
             drug = "Decorin mimetics; TGF-beta modulators", category = "Anti-fibrotic", stage = "Preclinical"),
  data.table(gene = "RHOA", axis = "Disease", logFC = 0.121, padj = 3.87e-02,
             drug = "Fasudil (ROCK inhibitor); Statins", category = "Anti-inflammatory", stage = "Approved"),
  data.table(gene = "CIDEA", axis = "Disease", logFC = -0.755, padj = 1.60e-04,
             drug = "PPAR-gamma agonists", category = "Metabolic", stage = "Approved"),
  data.table(gene = "KLF6", axis = "Disease", logFC = 0.180, padj = 2.13e-02,
             drug = "Anti-fibrotic agents", category = "Anti-fibrotic", stage = "Preclinical"),
  data.table(gene = "MSN", axis = "Disease", logFC = 0.448, padj = 4.57e-02,
             drug = "ERM inhibitors (investigational)", category = "Investigational", stage = "Preclinical"),
  data.table(gene = "RAD54L", axis = "Disease", logFC = 0.229, padj = 2.41e-02,
             drug = "PARP inhibitors", category = "Oncology", stage = "Approved"),
  data.table(gene = "FKBP14", axis = "Disease", logFC = -0.186, padj = 1.05e-03,
             drug = "FKBP ligands; mTOR pathway", category = "Immunosuppressant", stage = "Approved"),
  data.table(gene = "S100A7", axis = "Lesion", logFC = 1.86, padj = 1.45e-05,
             drug = "S100A7/psoriasin inhibitor", category = "Anti-inflammatory", stage = "Investigational"),
  data.table(gene = "S100A8", axis = "Lesion", logFC = 2.58, padj = 2.93e-07,
             drug = "Tasquinimod (calprotectin inhibitor)", category = "Anti-inflammatory", stage = "Phase II"),
  data.table(gene = "S100A9", axis = "Lesion", logFC = 3.43, padj = 7.90e-17,
             drug = "Tasquinimod (calprotectin inhibitor)", category = "Anti-inflammatory", stage = "Phase II"),
  data.table(gene = "IL36G", axis = "Lesion", logFC = 1.53, padj = 1.98e-07,
             drug = "Spesolimab (anti-IL36R)", category = "Biologic", stage = "Approved"),
  data.table(gene = "SERPINB3", axis = "Lesion", logFC = 2.68, padj = 1.48e-22,
             drug = "Serpin/SCCA modulators", category = "Protease inhibitor", stage = "Investigational"),
  data.table(gene = "AKR1B10", axis = "Lesion", logFC = 2.61, padj = 2.43e-23,
             drug = "Epalrestat (aldose reductase inh.)", category = "Metabolic", stage = "Approved"),
  data.table(gene = "KRT16", axis = "Lesion", logFC = 1.09, padj = 2.00e-05,
             drug = "Keratin modulators; wound healing", category = "Dermatologic", stage = "Investigational"),
  data.table(gene = "KRT6A", axis = "Lesion", logFC = 1.11, padj = 8.42e-06,
             drug = "Keratin 6A modulators", category = "Dermatologic", stage = "Investigational"),
  data.table(gene = "MMP12", axis = "Lesion", logFC = 2.63, padj = 1.10e-08,
             drug = "Doxycycline (MMP inhibitor)", category = "Anti-remodeling", stage = "Approved"),
  data.table(gene = "DEFB4A", axis = "Lesion", logFC = 1.40, padj = 5.42e-08,
             drug = "Beta-defensin modulators", category = "Antimicrobial", stage = "Investigational"),
  data.table(gene = "CXCL9", axis = "Lesion", logFC = 0.96, padj = 5.31e-04,
             drug = "CXCR3 antagonists", category = "Anti-inflammatory", stage = "Investigational"),
  data.table(gene = "CXCL10", axis = "Lesion", logFC = 1.20, padj = 6.71e-06,
             drug = "CXCR3 antagonists", category = "Anti-inflammatory", stage = "Investigational"),
  data.table(gene = "CXCL11", axis = "Lesion", logFC = 1.05, padj = 8.39e-05,
             drug = "CXCR3 antagonists", category = "Anti-inflammatory", stage = "Investigational"),
  data.table(gene = "PI3", axis = "Lesion", logFC = 2.67, padj = 1.35e-17,
             drug = "Elafin/SKALP peptide therapy", category = "Anti-inflammatory", stage = "Preclinical"),
  data.table(gene = "CHI3L1", axis = "Lesion", logFC = 0.95, padj = 3.12e-03,
             drug = "Anti-CHI3L1 antibodies", category = "Anti-fibrotic", stage = "Preclinical"),
  data.table(gene = "POSTN", axis = "Lesion", logFC = 1.15, padj = 6.48e-05,
             drug = "Anti-periostin antibodies", category = "Anti-fibrotic", stage = "Investigational"),
  data.table(gene = "ERBB4", axis = "Lesion", logFC = -2.34, padj = 1.10e-10,
             drug = "Neratinib (pan-HER TKI)", category = "Oncology/TKI", stage = "Approved"),
  data.table(gene = "BTC", axis = "Lesion", logFC = -1.58, padj = 8.48e-06,
             drug = "EGF-family targeted therapy", category = "Growth factor", stage = "Investigational"),
  data.table(gene = "IL4R", axis = "Lesion", logFC = 1.054, padj = 3.41e-19,
             drug = "Dupilumab (anti-IL4Ralpha)", category = "AD Biologic Target", stage = "Approved"),
  data.table(gene = "IL13", axis = "Lesion", logFC = 0.695, padj = 2.00e-01,
             drug = "Tralokinumab; Lebrikizumab (anti-IL13)", category = "AD Biologic Target", stage = "Approved"),
  data.table(gene = "IL13RA1", axis = "Lesion", logFC = 0.357, padj = 1.26e-05,
             drug = "IL-13 receptor complex", category = "AD Biologic Target", stage = "Investigational"),
  data.table(gene = "IL13RA2", axis = "Lesion", logFC = 1.198, padj = 6.81e-06,
             drug = "IL-13 decoy receptor", category = "AD Biologic Target", stage = "Investigational")
))

drug_map <- merge(drug_map, axis_map[, .(gene, target_axis)], by = "gene", all.x = TRUE)
drug_map[is.na(target_axis) & axis == "Disease", target_axis := "AD_related"]
drug_map[is.na(target_axis) & axis == "Lesion", target_axis := "LS_related"]

query_dgidb_gene <- function(gene) {
  q <- paste0(
    "query q($genes:[String!]){ interactions(geneNames:$genes, first:100){ ",
    "totalCount nodes { drug { name approved } gene { name } interactionTypes { type directionality } ",
    "sources { sourceDbName } } } }"
  )
  tryCatch({
    resp <- request("https://dgidb.org/api/graphql") |>
      req_body_json(list(query = q, variables = list(genes = list(gene)))) |>
      req_timeout(20) |>
      req_error(is_error = function(resp) FALSE) |>
      req_perform()
    if (resp_status(resp) >= 300) stop(resp_body_string(resp))
    x <- resp_body_json(resp, simplifyVector = FALSE)
    nodes <- x$data$interactions$nodes
    sources <- unlist(lapply(nodes, function(n) {
      vapply(n$sources, function(s) s$sourceDbName %||% NA_character_, character(1))
    }), use.names = FALSE)
    approved <- unlist(lapply(nodes, function(n) isTRUE(n$drug$approved)), use.names = FALSE)
    data.table(
      gene = gene,
      dgidb_interaction_count = x$data$interactions$totalCount %||% length(nodes),
      dgidb_approved_interaction = any(approved, na.rm = TRUE),
      dgidb_sources = safe_chr(sources),
      dgidb_status = "queried"
    )
  }, error = function(e) {
    data.table(gene = gene, dgidb_interaction_count = NA_integer_,
               dgidb_approved_interaction = NA,
               dgidb_sources = NA_character_,
               dgidb_status = paste("query_failed:", conditionMessage(e)))
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a

query_dgidb_batch <- function(genes) {
  q <- paste0(
    "query q($genes:[String!]){ interactions(geneNames:$genes, first:1000){ ",
    "totalCount nodes { drug { name approved } gene { name } interactionTypes { type directionality } ",
    "sources { sourceDbName } } } }"
  )
  tryCatch({
    resp <- request("https://dgidb.org/api/graphql") |>
      req_body_json(list(query = q, variables = list(genes = as.list(genes)))) |>
      req_timeout(60) |>
      req_error(is_error = function(resp) FALSE) |>
      req_perform()
    if (resp_status(resp) >= 300) stop(resp_body_string(resp))
    x <- resp_body_json(resp, simplifyVector = FALSE)
    nodes <- x$data$interactions$nodes
    if (length(nodes) == 0) {
      return(data.table(gene = genes, dgidb_interaction_count = 0L,
                        dgidb_approved_interaction = FALSE,
                        dgidb_sources = NA_character_,
                        dgidb_status = "queried_batch"))
    }
    flat <- rbindlist(lapply(nodes, function(n) {
      data.table(
        gene = n$gene$name %||% NA_character_,
        drug = n$drug$name %||% NA_character_,
        approved = isTRUE(n$drug$approved),
        source = safe_chr(vapply(n$sources, function(s) s$sourceDbName %||% NA_character_,
                                 character(1)))
      )
    }), fill = TRUE)
    summ <- flat[!is.na(gene), .(
      dgidb_interaction_count = uniqueN(drug),
      dgidb_approved_interaction = any(approved, na.rm = TRUE),
      dgidb_sources = safe_chr(unlist(strsplit(source, "; ", fixed = TRUE))),
      dgidb_status = "queried_batch"
    ), by = gene]
    merge(data.table(gene = genes), summ, by = "gene", all.x = TRUE)[
      is.na(dgidb_interaction_count),
      `:=`(dgidb_interaction_count = 0L,
           dgidb_approved_interaction = FALSE,
           dgidb_status = "queried_batch")]
  }, error = function(e) {
    rbindlist(lapply(genes, query_dgidb_gene), fill = TRUE)
  })
}

query_open_targets <- function(symbol) {
  ens <- tryCatch(
    AnnotationDbi::mapIds(org.Hs.eg.db, keys = symbol, keytype = "SYMBOL",
                          column = "ENSEMBL", multiVals = "first"),
    error = function(e) NA_character_
  )
  if (is.na(ens) || !nzchar(ens)) {
    return(data.table(gene = symbol, ensembl_id = NA_character_,
                      open_targets_tractability = NA_character_,
                      open_targets_known_drug_count = NA_integer_,
                      open_targets_status = "no_ensembl_mapping"))
  }

  q <- paste0(
    "query target($ensemblId:String!){ target(ensemblId:$ensemblId){ ",
    "id approvedSymbol tractability { label modality value } } }"
  )
  tryCatch({
    resp <- request("https://api.platform.opentargets.org/api/v4/graphql") |>
      req_body_json(list(query = q, variables = list(ensemblId = unname(ens)))) |>
      req_timeout(20) |>
      req_error(is_error = function(resp) FALSE) |>
      req_perform()
    if (resp_status(resp) >= 300) stop(resp_body_string(resp))
    x <- resp_body_json(resp, simplifyVector = FALSE)
    tr <- x$data$target$tractability
    true_tr <- Filter(function(z) isTRUE(z$value), tr)
    labels <- vapply(true_tr, function(z) paste(z$modality, z$label, sep = ":"),
                     character(1))
    approved_count <- sum(vapply(true_tr, function(z) z$label == "Approved Drug",
                                  logical(1)))
    data.table(gene = symbol, ensembl_id = unname(ens),
               open_targets_tractability = safe_chr(labels),
               open_targets_known_drug_count = approved_count,
               open_targets_status = "queried")
  }, error = function(e) {
    data.table(gene = symbol, ensembl_id = unname(ens),
               open_targets_tractability = NA_character_,
               open_targets_known_drug_count = NA_integer_,
               open_targets_status = paste("query_failed:", conditionMessage(e)))
  })
}

genes_to_annotate <- sort(unique(c(axis_map$gene, drug_map$gene)))
cat(sprintf("Annotating %d biomarker/drug-target genes with Open Targets + DGIdb\n",
            length(genes_to_annotate)))
ot <- rbindlist(lapply(genes_to_annotate, query_open_targets), fill = TRUE)
dg <- query_dgidb_batch(genes_to_annotate)

curated_counts <- drug_map[, .(
  curated_drug_count = uniqueN(drug),
  curated_approved = any(stage_is_approved(stage)),
  curated_clinical = any(stage_is_clinical(stage)),
  curated_drugs = safe_chr(drug)
), by = gene]

anno <- merge(axis_map, curated_counts, by = "gene", all.x = TRUE)
anno <- merge(anno, ot, by = "gene", all.x = TRUE)
anno <- merge(anno, dg, by = "gene", all.x = TRUE)
anno[is.na(curated_drug_count), curated_drug_count := 0L]
anno[is.na(curated_approved), curated_approved := FALSE]
anno[is.na(curated_clinical), curated_clinical := FALSE]
anno[is.na(dgidb_interaction_count), dgidb_interaction_count := 0L]
anno[is.na(dgidb_approved_interaction), dgidb_approved_interaction := FALSE]
anno[is.na(open_targets_known_drug_count), open_targets_known_drug_count := 0L]

anno[, druggability_tier := fifelse(
  curated_approved | dgidb_approved_interaction | open_targets_known_drug_count > 0,
  "Tier 1: approved/clinical target",
  fifelse(curated_clinical | dgidb_interaction_count > 0,
          "Tier 2: drug-gene evidence",
          fifelse(!is.na(open_targets_tractability) & nzchar(open_targets_tractability),
                  "Tier 3: tractability only",
                  "Unannotated")))
]

fwrite(anno[order(target_axis, druggability_tier, gene)],
       file.path(out, "drug_repurposing_druggability_annotation.csv"))

drug_map <- merge(drug_map, anno[, .(gene, druggability_tier,
                                     open_targets_tractability,
                                     open_targets_known_drug_count,
                                     dgidb_interaction_count,
                                     dgidb_sources,
                                     target_axis)],
                  by = c("gene", "target_axis"), all.x = TRUE)

# ---- Signature reversal with signatureSearch ----
symbol_to_entrez <- function(symbols) {
  mapped <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = unique(symbols), keytype = "SYMBOL",
    column = "ENTREZID", multiVals = "first"))
  na.omit(unname(mapped[symbols]))
}

build_query <- function(effects, axis_label, top_n = 150) {
  effects <- effects[!is.na(padj) & padj < alpha & abs(logFC) >= fc_min]
  effects <- effects[order(-abs(logFC))]
  effects <- head(effects, top_n)
  up <- symbol_to_entrez(effects[logFC > 0, gene])
  down <- symbol_to_entrez(effects[logFC < 0, gene])
  list(axis = axis_label, upset = unique(up), downset = unique(down),
       n_up = length(unique(up)), n_down = length(unique(down)))
}

run_signature_search <- function(query_obj, ref_preference = c("cmap", "lincs")) {
  ref_preference <- match.arg(ref_preference)
  if (!identical(Sys.getenv("RUN_FULL_SIGNATURESEARCH"), "1")) {
    return(data.table(query_axis = query_obj$axis, perturbagen = known_ad_drugs,
                      score = NA_real_, rank = seq_along(known_ad_drugs),
                      known_ad_drug = TRUE, validation_hit = FALSE,
                      reference_db = ref_preference,
                      method_status = paste0(
                        "not_run_set_RUN_FULL_SIGNATURESEARCH_1_for_full_",
                        ref_preference, "_ExperimentHub_search")))
  }
  if (!requireNamespace("signatureSearch", quietly = TRUE) ||
      !requireNamespace("ExperimentHub", quietly = TRUE)) {
    return(data.table(query_axis = query_obj$axis, perturbagen = known_ad_drugs,
                      score = NA_real_, rank = seq_along(known_ad_drugs),
                      known_ad_drug = TRUE, validation_hit = FALSE,
                      reference_db = ref_preference,
                      method_status = "not_run_signatureSearch_or_ExperimentHub_missing"))
  }
  if (query_obj$n_up < 5 || query_obj$n_down < 5) {
    return(data.table(query_axis = query_obj$axis, perturbagen = known_ad_drugs,
                      score = NA_real_, rank = seq_along(known_ad_drugs),
                      known_ad_drug = TRUE, validation_hit = FALSE,
                      reference_db = ref_preference,
                      method_status = sprintf("not_run_insufficient_query_genes_up_%d_down_%d",
                                              query_obj$n_up, query_obj$n_down)))
  }

  tryCatch({
    suppressMessages({
      eh <- ExperimentHub::ExperimentHub()
      ref_id <- if (ref_preference == "lincs") "EH3226" else "EH3223"
      refdb <- eh[[ref_id]]
      qsig <- signatureSearch::qSig(
        query = list(upset = query_obj$upset, downset = query_obj$downset),
        gess_method = if (ref_preference == "lincs") "LINCS" else "CMAP",
        refdb = refdb
      )
      gess <- if (ref_preference == "lincs") {
        signatureSearch::gess_lincs(qsig, sortby = "NCS", tau = FALSE,
                                    chunk_size = 5000, workers = 1,
                                    addAnnotations = FALSE)
      } else {
        signatureSearch::gess_cmap(qsig, chunk_size = 5000, workers = 1,
                                   addAnnotations = FALSE)
      }
      res <- as.data.table(signatureSearch::result(gess))
    })
    score_col <- intersect(c("NCSct", "NCS", "WTCS", "scaled_score", "effect"), names(res))[1]
    pert_col <- intersect(c("pert", "pert_iname", "compound", "drug_name"), names(res))[1]
    if (is.na(score_col) || is.na(pert_col)) stop("signatureSearch result lacks score/perturbagen columns")
    res <- res[order(get(score_col))]
    res <- head(res, 200)
    out_dt <- data.table(
      query_axis = query_obj$axis,
      perturbagen = as.character(res[[pert_col]]),
      score = as.numeric(res[[score_col]]),
      rank = seq_len(nrow(res)),
      known_ad_drug = normalize_perturbagen(res[[pert_col]]) %in% known_ad_drugs,
      validation_hit = normalize_perturbagen(res[[pert_col]]) %in% known_ad_drugs,
      reference_db = ref_preference,
      method_status = paste0("run_", ref_preference, "_", ref_id)
    )
    controls <- data.table(query_axis = query_obj$axis, perturbagen = known_ad_drugs,
                           score = NA_real_, rank = NA_integer_,
                           known_ad_drug = TRUE,
                           validation_hit = known_ad_drugs %in% normalize_perturbagen(out_dt$perturbagen),
                           reference_db = ref_preference,
                           method_status = paste0("validation_controls_", ref_preference))
    rbindlist(list(out_dt, controls), fill = TRUE)
  }, error = function(e) {
    data.table(query_axis = query_obj$axis, perturbagen = known_ad_drugs,
               score = NA_real_, rank = seq_along(known_ad_drugs),
               known_ad_drug = TRUE, validation_hit = FALSE,
               reference_db = ref_preference,
               method_status = paste("not_run_reference_unavailable:", conditionMessage(e)))
  })
}

D <- fread(file.path(out, "meta_disease_effects.csv"))[
  , .(gene, logFC = pooled_logFC, padj)]
L <- fread(file.path(out, "meta_lesion_effects.csv"))[
  , .(gene, logFC = pooled_logFC, padj)]
both_genes <- axis_map[target_axis == "AD_and_LS", gene]
both <- rbindlist(list(
  D[gene %in% both_genes],
  L[gene %in% both_genes]
))[, .(logFC = mean(logFC, na.rm = TRUE), padj = min(padj, na.rm = TRUE)), by = gene]

queries <- list(
  build_query(D, "AD_related"),
  build_query(L, "LS_related"),
  build_query(both, "AD_and_LS", top_n = 80)
)

sig_res <- rbindlist(c(
  lapply(queries, run_signature_search, ref_preference = "cmap"),
  lapply(queries, run_signature_search, ref_preference = "lincs")
), fill = TRUE)
fwrite(sig_res, file.path(out, "drug_repurposing_signature_reversal.csv"))

# ---- Figures ----
dm <- drug_map[order(target_axis, -abs(logFC))]
dm[, gene_label := factor(sprintf("%s (%s)", gene, signed_label(logFC)),
                          levels = rev(sprintf("%s (%s)", gene, signed_label(logFC))))]

category_palette <- c(
  "Immunotherapy" = "#4E79A7",
  "Anti-inflammatory" = "#D55E00",
  "Ion channel" = "#59A14F",
  "Barrier repair" = "#9C755F",
  "Anti-fibrotic" = "#B07AA1",
  "Metabolic" = "#F28E2B",
  "Investigational" = "#8D99AE",
  "Oncology" = "#CC79A7",
  "Protease inhibitor" = "#76B7B2",
  "Biologic" = "#E15759",
  "Dermatologic" = "#EDC948",
  "Anti-remodeling" = "#AF7AA1",
  "Antimicrobial" = "#86BCB6",
  "Oncology/TKI" = "#FF9DA7",
  "Immunosuppressant" = "#6D6875",
  "Growth factor" = "#9C6644",
  "AD Biologic Target" = "#111827"
)

g1 <- ggplot(dm, aes(x = logFC, y = gene_label, color = category, shape = stage)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = ad_colors$muted) +
  geom_point(aes(size = -log10(padj)), alpha = 0.9, stroke = 0.8) +
  facet_grid(target_axis ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = category_palette, name = "Drug category") +
  scale_shape_manual(values = c("Approved" = 16, "Phase II" = 17, "Phase III" = 15,
                                "Phase I/II" = 8, "Investigational" = 1,
                                "Preclinical" = 4),
                     name = "Clinical stage") +
  scale_size_continuous(name = expression(-log[10](P[adj])), range = c(2.2, 6.5)) +
  labs(x = expression(log[2]~"Fold Change on source axis"),
       y = NULL,
       title = "Drug repurposing targets across AD, lesion, and shared axes",
       subtitle = "Curated target overlap annotated with Open Targets tractability and DGIdb evidence") +
  theme_ad_pub(base_size = 10, legend = "right") +
  theme(axis.text.y = element_text(size = 7.5, face = "bold"),
        strip.text.y = element_text(angle = 0))
save_pub_pair("plot_drug_repurposing_coverage", g1, width = 11.5, height = 9.2,
              dpi = 500, out_dir = out)

anno_plot <- anno[target_axis %in% c("AD_related", "LS_related", "AD_and_LS"),
                  .N, by = .(target_axis, druggability_tier)]
anno_plot[, target_axis := factor(target_axis,
                                  levels = c("AD_related", "LS_related", "AD_and_LS"))]
g2 <- ggplot(anno_plot, aes(target_axis, N, fill = druggability_tier)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.35) +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.4) +
  scale_fill_manual(values = c(
    "Tier 1: approved/clinical target" = ad_colors$both,
    "Tier 2: drug-gene evidence" = ad_colors$disease,
    "Tier 3: tractability only" = ad_colors$gold,
    "Unannotated" = "#AEB8C2"
  ), name = NULL) +
  labs(x = NULL, y = "Number of biomarker genes",
       title = "Druggability annotation by transcriptomic axis",
       subtitle = "Open Targets tractability plus DGIdb and curated drug-gene evidence") +
  theme_ad_pub(base_size = 11, legend = "right")
save_pub_pair("plot_drug_repurposing_druggability", g2, width = 8.2, height = 5.4,
              dpi = 500, out_dir = out)

sig_plot <- sig_res[!is.na(rank) & !is.na(score)][order(rank)]
if (nrow(sig_plot) > 0) {
  sig_plot <- sig_plot[, head(.SD, 15), by = query_axis]
  sig_plot[, perturbagen := factor(perturbagen, levels = rev(unique(perturbagen)))]
  g3 <- ggplot(sig_plot, aes(score, perturbagen, color = query_axis)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = ad_colors$muted) +
    geom_point(size = 2.8, alpha = 0.9) +
    facet_grid(query_axis ~ ., scales = "free_y", space = "free_y") +
    scale_color_manual(values = axis_palette, guide = "none") +
    labs(x = "signatureSearch connectivity score (negative = reversal)",
         y = NULL,
         title = "CMap signature-reversal candidates",
         subtitle = "Top inverse perturbagens for AD, lesion, and shared signatures") +
    theme_ad_pub(base_size = 10)
  save_pub_pair("plot_drug_repurposing_signature_reversal", g3, width = 8.8, height = 9,
                dpi = 500, out_dir = out)
}

summary_dt <- rbindlist(list(
  data.table(analysis_tier = "Tier 1: Open Targets + DGIdb annotation",
             n_rows = nrow(anno),
             n_validation_hits = NA_integer_,
             method_status = safe_chr(unique(c(anno$open_targets_status, anno$dgidb_status)))),
  data.table(analysis_tier = "Tier 2: LINCS/CMap signature reversal",
             n_rows = nrow(sig_res),
             n_validation_hits = sum(sig_res$validation_hit, na.rm = TRUE),
             method_status = safe_chr(unique(sig_res$method_status)))
), fill = TRUE)
fwrite(summary_dt, file.path(out, "drug_repurposing_summary.csv"))
fwrite(dm, file.path(out, "drug_repurposing_results.csv"))

cat(sprintf("Total curated drug-gene pairs: %d\n", nrow(dm)))
cat(sprintf("Druggability annotation genes: %d\n", nrow(anno)))
cat(sprintf("Signature-reversal rows: %d\n", nrow(sig_res)))
cat("Known AD drug validation hits recovered: ",
    sum(sig_res$validation_hit, na.rm = TRUE), "\n", sep = "")
cat("Wrote drug_repurposing_results.csv, drug_repurposing_druggability_annotation.csv,\n",
    "drug_repurposing_signature_reversal.csv, drug_repurposing_summary.csv\n", sep = "")
