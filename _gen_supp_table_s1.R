#!/usr/bin/env Rscript
# Generate Supplementary Table S1: Full lesion biomarker meta-analysis results

x <- read.csv("disease_biomarkers.csv", stringsAsFactors = FALSE)
x <- x[order(-x$score), ]

# Format columns
x$gene      <- gsub("_", "\\_", x$gene, fixed = TRUE)
x$gene      <- gsub("%", "\\%", x$gene, fixed = TRUE)
x$gene      <- gsub("&", "\\&", x$gene, fixed = TRUE)
x$gene      <- gsub("#", "\\#", x$gene, fixed = TRUE)
x$logFC_fmt <- sprintf("%.2f", x$pooled_logFC)
x$ci_fmt    <- sprintf("[%.2f, %.2f]", x$ci_lower, x$ci_upper)
x$pval_fmt  <- ifelse(x$pval < 1e-100, "$<$1e-100", sprintf("%.1e", x$pval))
x$padj_fmt  <- ifelse(x$padj < 1e-100, "$<$1e-100", sprintf("%.1e", x$padj))
x$I2_fmt    <- sprintf("%.1f", x$I2)
x$score_fmt <- sprintf("%.1f", x$score)

out <- file("supp_table_s1.tex", "w")

cat("% Supplementary Table S1: Top 200 lesion biomarker meta-analysis results\n", file = out)
cat("% Auto-generated from disease_biomarkers.csv\n\n", file = out)
cat("\\subsection*{Supplementary Table S1: Full Lesion Biomarker Meta-Analysis Results}\n", file = out)
cat("\\noindent Top 200 genes ranked by composite score ($-\\log_{10}P \\times |\\log_2\\text{FC}| \\times (1-I^2/100) \\times \\sqrt{k}$).\n", file = out)
cat("\\vspace{4pt}\n\n", file = out)

cat("\\begin{longtable}{lrrrrrrr}\n", file = out)
cat("\\caption{Top 200 lesion biomarkers from random-effects meta-analysis of 31 datasets (LS vs NL/HC).}\\\\\n", file = out)
cat("\\toprule\n", file = out)
cat("\\textbf{Gene} & \\textbf{log$_2$FC} & \\textbf{95\\% CI} & \\textbf{$P$} & \\textbf{FDR} & \\textbf{$I^2$ (\\%)} & \\textbf{$k$} & \\textbf{Score} \\\\\n", file = out)
cat("\\midrule\n", file = out)
cat("\\endfirsthead\n\n", file = out)

cat("\\toprule\n", file = out)
cat("\\textbf{Gene} & \\textbf{log$_2$FC} & \\textbf{95\\% CI} & \\textbf{$P$} & \\textbf{FDR} & \\textbf{$I^2$ (\\%)} & \\textbf{$k$} & \\textbf{Score} \\\\\n", file = out)
cat("\\midrule\n", file = out)
cat("\\endhead\n\n", file = out)

cat("\\bottomrule\n", file = out)
cat("\\endfoot\n\n", file = out)

for (i in seq_len(nrow(x))) {
  line <- sprintf("%s & %s & %s & %s & %s & %s & %d & %s \\\\\n",
    x$gene[i], x$logFC_fmt[i], x$ci_fmt[i], x$pval_fmt[i], x$padj_fmt[i],
    x$I2_fmt[i], x$k[i], x$score_fmt[i])
  cat(line, file = out)
}

cat("\\end{longtable}\n", file = out)
close(out)

cat("Generated supp_table_s1.tex with", nrow(x), "rows\n")
