# AD Transcriptomic Disease Signature — Analysis Code

This repository contains the complete R analysis pipeline for:

**"A Cross-Platform Meta-Analysis of the Atopic Dermatitis Transcriptomic
Disease Signature Across 39 Datasets (2,138 Individuals)"**

Van Hung Huynh-Tran — luctran.biostats@gmail.com

---

## Requirements

- R ≥ 4.5.0
- R packages: `metafor`, `limma`, `DESeq2`, `GEOquery`, `gprofiler2`, `igraph`,
  `ggplot2`, `ggrepel`, `pheatmap`, `dplyr`, `data.table`, `org.Hs.eg.db`
- Bioconductor packages: `limma`, `DESeq2`, `GEOquery`, `gprofiler2`,
  `org.Hs.eg.db`, `hgu133plus2.db`, `hgu133a.db`

Install missing packages before running:

```r
install.packages(c("metafor", "igraph", "ggplot2", "ggrepel", "pheatmap",
                    "dplyr", "data.table"))
BiocManager::install(c("limma", "DESeq2", "GEOquery", "gprofiler2",
                        "org.Hs.eg.db", "hgu133plus2.db", "hgu133a.db"))
```

## Input Data

Before running, place the following curated files in the `analysis/` directory
(available from Zenodo, DOI: 20707787, or from the corresponding author upon
reasonable request):

- `stage1_disease_GSE*.csv` — per-dataset disease-contrast DE results
- `stage1_GSE*.csv` — per-dataset lesion-contrast DE results
- `meta_disease_effects.csv` — pooled disease-axis random-effects meta-analysis
- `meta_lesion_effects.csv` — pooled lesion-axis random-effects meta-analysis
- `consensus_disease_signature.csv` — consensus vote-count table
- `dataset_contrast_annotation.csv` — study-level metadata
- `logFC_disease_matrix.csv`, `logFC_lesion_matrix.csv` — logFC matrices
- `SE_disease_matrix.csv`, `SE_lesion_matrix.csv` — standard error matrices

## Running the Pipeline

From the `analysis/` directory:

```r
source("zenodo_code/00_master_run.R")
```

Or from the command line (Windows PowerShell):

```powershell
Set-Location ".../disease_signature/analysis"
& "C:/Program Files/R/R-4.5.0/bin/Rscript.exe" zenodo_code/00_master_run.R
```

## Script Inventory

| Script | Description |
|--------|-------------|
| `00_master_run.R` | Master runner — executes all steps in order |
| `00_normalize_gene_names.R` | Gene-name normalisation across platforms |
| `07_causal_trait_state.R` | Trait vs. state classification |
| `08_causal_coloc.R` | GWAS literature annotation |
| `10_causal_treatment.R` | Treatment-response directionality |
| `11_causal_summary.R` | Causal evidence integration |
| `12_disease_biomarker_signature.R` | 135-gene field-defect signature |
| `13_k_sensitivity.R` | Leave-k-out sensitivity analysis |
| `14_ad_ls_decomposition.R` | Three-group AD/LS partition |
| `15_combined_drug_repurposing.R` | Druggability + signature reversal |
| `16_regenerate_all_figures.R` | Regenerate all manuscript figures |
| `17_module_enrichment_fill.R` | Module enrichment gap-filling |
| `99_validate_drug_repurposing_outputs.R` | Drug repurposing validation |
| `generate_supplementary.R` | Supplementary analyses |
| `plotting_utils.R` | Shared plotting functions |
| `_gen_supp_table_s1.R` | Supplementary Table S1 generation |
| `_reviewer_checks.R` | Reviewer-requested diagnostics |

## Output

All outputs are written to the `analysis/` directory. Key files:

- `disease_biomarker_signature_final.csv` — 135-gene field-defect signature
- `core_disease_signature_final.csv` — consensus-ranked disease genes
- `causal_summary.csv` — causal evidence integration
- `ad_ls_decomposition_summary.csv` — three-group partition counts
- `plot_*.pdf` — all manuscript figures

## License

MIT. See `LICENSE` file in the Zenodo deposit.

## Citation

If you use this code, please cite the accompanying manuscript and the Zenodo
deposit (DOI: 20707787).
