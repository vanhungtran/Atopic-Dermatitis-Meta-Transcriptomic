# Statistical Methods — AD Disease Signature Meta-Analysis

## Study Design
We conducted a two-contrast random-effects meta-analysis of 39 transcriptomic
datasets from the Gene Expression Omnibus (GEO), each comparing atopic
dermatitis (AD) skin to non-lesional or healthy control skin. Datasets were
classified using a 2×2 design: AD status (AD, HC) × Lesion status (LS, NL).
Healthy controls have only NL samples; AD patients contribute LS and/or NL
samples. This yields two independent effect estimates per gene:

1. **Lesion effect** (`beta_L`): AD lesional vs AD non-lesional skin (22
   datasets, 12 with paired within-patient biopsies). Estimated via linear
   models blocking for patient when possible.
2. **Disease effect** (`beta_D`): AD non-lesional vs healthy control skin
   (15 datasets). Estimated via unpaired comparison.

## Effect Size Extraction
For each dataset, differential expression was computed using **limma**
(microarray) or **DESeq2** (RNA-seq) with the appropriate contrast. Effect
sizes are reported as log2 fold change (log2FC) with standard errors. For
paired-design datasets (same patient, LS and NL biopsies), patient identity
was included as a blocking factor: `~ patient + tissue`. All gene identifiers
were mapped to HGNC gene symbols.

**Gene-name normalization (critical for all joins/merges)**: Because each
dataset's platform annotation (Affymetrix probe IDs, Illumina probe IDs,
RNA-seq gene IDs/Ensembl, etc.) uses different identifier schemes and gene
symbols drift over time, every per-dataset table must be mapped to a single
current HGNC symbol set *before* any merge — including the `meta_*_effects`
pooling, the `consensus_*_signature` joins, the trait/state classification
join in the Causal Inference section, and any join against external GWAS/eQTL
gene lists. Use `limma::alias2SymbolTable()` or
`HGNChelper::checkGeneSymbols()` to resolve outdated aliases to current
symbols, and collapse duplicate probes/IDs mapping to the same symbol by
taking the value with the largest absolute effect (or average, if blocking
makes this inappropriate) before writing each `stage1_*.csv`. Verify gene
counts before/after normalization to catch silent symbol collisions.

## Meta-Analysis
Per-gene random-effects meta-analysis was performed using restricted maximum
likelihood (REML) estimation via the `metafor` package (Viechtbauer, 2010):

`theta_hat = sum(w_i * theta_i) / sum(w_i),  w_i = 1 / (v_i + tau^2)`

where `theta_hat_i` is the log2FC from dataset `i`, `v_i` its sampling
variance, and `tau^2` the between-study variance estimated by REML.

Genes present in fewer than 5 datasets were excluded from pooling.

### Heterogeneity
Between-study heterogeneity was quantified using Cochran's Q test, the I²
statistic with 95% confidence intervals, and the absolute between-study
variance `tau^2`. Prediction intervals (95% PI) were computed to estimate the
range of true effects in future studies (Harrer et al., 2022, Ch 4).

### Biomarker Ranking
Genes were ranked by a composite score:
`score = -log10(p) * |log2FC| * (1 - I^2) * sqrt(k)`
penalizing heterogeneity and rewarding replication breadth.

### Subgroup Analysis
Pre-specified subgroup analyses tested whether pooled effects differed by:
platform (microarray vs RNA-seq), study design (paired vs unpaired), and
tissue type (skin biopsy, tape-strip, PBMC, laser-capture microdissection).
Between-subgroup differences were tested using Cochran's Q.

### Meta-Regression
Continuous moderators (sample size, publication year) were tested via
mixed-effects meta-regression with permutation-based p-values when k < 20.

## Publication Bias
For genes present in ≥10 datasets, funnel plot asymmetry was assessed via
Egger's regression test. Contour-enhanced funnel plots distinguished
publication bias from heterogeneity. Trim-and-fill (Duval & Tweedie) provided
adjusted pooled estimates as sensitivity analysis.

## Causal Inference Extensions
The two-contrast design above is purely associational. The following analyses
add causal evidence to (a) the disease effect (`beta_D`: AD-NL vs HC) and
(b) the lesion effect (`beta_L`: AD-LS vs AD-NL). Each subsection specifies
exact inputs, methods, and output files so it can be run as a standalone
script (`07_causal_*.R`) added to the existing pipeline directory
(`analysis/github/scripts/`, alongside `02_meta_analysis.R` etc.).

### 1. Trait vs. State Classification (no external data required)
**Objective**: Separate genes that are already dysregulated in clinically
normal AD skin (trait-like, candidate causal/predisposing) from genes that
only change in active lesions (state-like, likely reactive to inflammation).

**Inputs**: `meta_disease_effects.csv` (`beta_D`, padj) and
`meta_lesion_effects.csv` (`beta_L`, padj), joined on `gene`.

**Filtering**: Before classification, (i) housekeeping/structural genes whose
tiny but ultra-precise pooled effects are normalization artifacts rather than
biology are removed by symbol pattern (ribosomal proteins `RPS*/RPL*/RPLP*/
MRPS*/MRPL*`, elongation/initiation factors `EEF*/EIF*`, and canonical
housekeeping genes `ACTB, ACTG1, GAPDH, TUBB*, TUBA*, B2M, HPRT1, PPIA, UBC,
YWHAZ`); and (ii) a gene is only counted as "significant" on an axis if it
additionally has |pooled log2FC| ≥ 0.5 on that axis, so trivial high-`k`
effects do not dominate. Both thresholds are exposed as `hk_pattern` and
`fc_min` at the top of the script.

**Method**: For each retained gene significant in at least one contrast,
classify into one of four categories:
- `trait_only`: significant `beta_D`, non-significant `beta_L`
- `state_only`: significant `beta_L`, non-significant `beta_D`
- `trait_and_state_concordant`: both significant, same sign
- `trait_and_state_discordant`: both significant, opposite sign

**Output**: `causal_trait_state_classification.csv` with columns `gene,
beta_D, padj_D, beta_L, padj_L, category`. Generate a 2x2 quadrant scatter
plot (`plot_trait_state_quadrant.png`) of `beta_D` vs `beta_L` colored by
category. `trait_only` and `trait_and_state_concordant` genes are the
strongest candidates for drivers of AD susceptibility rather than lesion
byproducts.

### 2. GWAS Literature Annotation (supporting only)
**Objective**: Annotate trait/state candidate genes with whether they appear
among curated AD GWAS loci, as supporting per-gene context. This is an
annotation, NOT a hypothesis test: it is not an enrichment, colocalization, or
Mendelian-randomization analysis. Because the curated loci list is incomplete
and no GWAS/eQTL summary statistics were available, a low overlap must NOT be
interpreted as genetic novelty or as evidence against a genetic role.

**Data source**: Curated gene lists from three published AD GWAS meta-analyses:
- Paternoster et al. 2015 (Nat Genet 47:1449–1456): 31 loci from 10,788 cases
  / 30,047 controls of European ancestry.
- Sliz et al. 2022 (Nat Commun 13:1295): 37 loci from FinnGen (864K
  individuals).
- Weidinger et al. 2018 (Nat Rev Dis Primers 4:1): review of established AD
  GWAS loci.

Additionally, 20 genes with known GTEx v8 skin eQTL colocalization (H4 > 0.8)
from published GWAS-eQTL coloc studies (GTEx Consortium, Science 2020;
369:1318–1330) are flagged.

**Method**: For each gene in the step-1 candidate list (`trait_only` +
`trait_and_state_concordant`), check membership in the union of the three
GWAS gene sets. A gene is flagged `gwas_support = TRUE` if it appears in ≥2
of the three catalogs (replicated). Skin eQTL colocalization is recorded as a
separate boolean flag.

**Output**: `causal_coloc_results.csv` with columns `gene, gwas_source,
eqtl_skin, gwas_support, note`. This replaces the formal `coloc::coloc.abf()`
analysis that requires complete GWAS + eQTL summary statistics downloads
(not feasible in the current environment due to GWAS Catalog REST API
endpoint changes).

### 3. Mendelian Randomization — REMOVED (folded into step 2)
A separate "literature-level MR" step previously flagged `mr_support = TRUE`
when a gene appeared in ≥2 GWAS catalogs. This is the SAME curated-catalogue
look-up as step 2, so counting it as a second, independent support double-counts
one piece of evidence and inflates the causal-evidence tiers. It has been
collapsed into the single `gwas_support` flag from step 2. No formal two-sample
MR (`TwoSampleMR`) was run — it requires complete GWAS + eQTL summary
statistics, which were unavailable. `causal_mr_results.csv` is no longer used.

### 4. Lesion-Effect Directionality via Treatment Response
**Objective**: For the `beta_L` (lesion) signature, test whether top genes
normalize toward NL/HC levels upon successful treatment — supporting a
causal/mediating role in lesion pathology rather than a static marker of
chronicity.

**Inputs**: Datasets in the inventory with paired pre/post-treatment biopsies
(check `dataset_contrast_annotation.csv` for designs such as
GSE283265 LS_baseline vs follow-up, GSE133477/GSE130588/GSE133385/GSE58558/
GSE99802/GSE60709/GSE111053/GSE111055/GSE59294 "Day1" vs later timepoints,
which already have paired LS-vs-NL contrasts available as
`stage1_*_paired.csv`).

**Method**: For each top-30 `consensus_lesion_signature.csv` gene, fit
`~ timepoint + (1|patient)` (or paired t-test where only two timepoints
exist) on log-expression, restricted to lesional samples. A gene is
"treatment-responsive" if its fold-change reverses sign relative to
`beta_L` with treatment (i.e., moves back toward the NL baseline).

**Output**: `causal_treatment_response.csv` with columns `gene, dataset,
beta_L, beta_treatment_response, concordant_with_reversal (TRUE/FALSE)`.
Summarize the proportion of top-30 genes that are treatment-responsive as
supporting evidence that the lesion signature reflects active, reversible
disease processes rather than fixed tissue remodeling.

### Reporting
Combine results from steps 1, 2 and 4 into a single `causal_summary.csv` (one
row per gene in the union of all candidate lists) with boolean columns
`trait_like, gwas_literature_support, treatment_responsive`, and a final
`causal_evidence_tier`. There are THREE possible supports (Tier 1: all 3;
Tier 2: 2; Tier 3: 1; none: associational only). GWAS contributes a single
support (the step-2 annotation); it is not double-counted. Report this tier alongside the existing
biomarker `score` in the manuscript's results table.

## Bayesian Sensitivity Analysis
For the top 100 biomarker genes, Bayesian random-effects meta-analysis was
performed with weakly informative priors (Half-Cauchy(0, 0.5) on $\tau$) via
the `brms` package. Posterior probabilities of |pooled log2FC| > 1 were
reported alongside frequentist results.

## Drug Repurposing
Drug repurposing was implemented as a two-tier analysis in
`scripts/15_combined_drug_repurposing.R`.

**Tier 1: Druggability annotation.** Disease-axis biomarkers, lesion-axis
biomarkers, AD-and-LS genes, and curated AD drug targets were annotated with
Open Targets tractability labels and DGIdb drug-gene interaction counts. Gene
symbols were mapped to Ensembl identifiers with `org.Hs.eg.db` before Open
Targets GraphQL queries. DGIdb GraphQL interactions were summarized per gene as
interaction counts, source databases, and whether any approved drug interaction
was present. Each gene was assigned a `druggability_tier`: Tier 1 for approved
or clinical target evidence, Tier 2 for drug-gene evidence, Tier 3 for
tractability-only evidence, or unannotated. The output table is
`drug_repurposing_druggability_annotation.csv`; the same annotations are joined
back into `drug_repurposing_results.csv`.

**Tier 2: LINCS/CMap signature reversal.** For the AD-related, LS-related, and
AD-and-LS signatures, up- and down-regulated gene sets were converted to Entrez
IDs and queried with `signatureSearch` against CMap (`EH3223`) and LINCS
(`EH3226`) ExperimentHub references. Negative connectivity scores were treated
as signature-reversal candidates. A validation panel of established AD
therapies (including topical/systemic corticosteroid, calcineurin-inhibitor,
immunosuppressant, biologic, PDE4, and JAK-inhibitor drugs) was carried through
the output to test whether known AD therapies were recovered among inverse
perturbagens. Results are written to
`drug_repurposing_signature_reversal.csv` and summarized in
`drug_repurposing_summary.csv`.

## Software
All analyses were performed in R 4.5.0 using packages: `metafor` (v4.x),
`meta`, `limma`, `DESeq2`, `ggplot2`, `gprofiler2`, `pheatmap`, `brms`,
`httr2`, `jsonlite`, `org.Hs.eg.db`, `ExperimentHub`, and `signatureSearch`.
The analysis pipeline and all figure-generating code are available at
[repository URL]. Reporting follows PRISMA 2020 guidelines.

## References
- Harrer M, Cuijpers P, Furukawa TA, Ebert DD. *Doing Meta-Analysis with R.* CRC Press, 2022.
- Viechtbauer W. Conducting meta-analyses in R with the metafor package. *J Stat Softw.* 2010;36(3):1-48.
