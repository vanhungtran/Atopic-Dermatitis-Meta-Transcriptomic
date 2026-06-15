#!/usr/bin/env Rscript
# Shared plotting helpers for the AD disease-signature pipeline.

suppressMessages({
  library(ggplot2)
  library(scales)
})

ad_colors <- list(
  ink = "#1F2933",
  muted = "#64748B",
  grid = "#E5E7EB",
  panel = "#FAFAF8",
  up = "#C0392B",
  down = "#2166AC",
  lesion = "#D55E00",
  disease = "#0072B2",
  both = "#7B3294",
  neutral = "#CBD5E1",
  trait = "#009E73",
  discordant = "#CC79A7",
  gold = "#E69F00",
  teal = "#008B8B"
)

axis_palette <- c(
  "Disease" = ad_colors$disease,
  "Lesion" = ad_colors$lesion,
  "AD_and_LS" = ad_colors$both,
  "AD_related" = ad_colors$trait,
  "LS_related" = ad_colors$lesion,
  "Both" = ad_colors$both,
  "Neither" = "#D7DCE2",
  "Disease only" = ad_colors$disease,
  "Lesion only" = ad_colors$lesion,
  "discordant" = ad_colors$discordant
)

direction_palette <- c(
  "Up" = ad_colors$up,
  "Down" = ad_colors$down,
  "up" = ad_colors$up,
  "down" = ad_colors$down
)

stage_palette <- c(
  "Approved" = "#124559",
  "Phase III" = "#2A9D8F",
  "Phase II" = "#E9C46A",
  "Phase I/II" = "#F4A261",
  "Investigational" = "#8D99AE",
  "Preclinical" = "#6D6875"
)

theme_ad_pub <- function(base_size = 11, legend = "bottom") {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(color = ad_colors$ink),
      plot.title = element_text(face = "bold", size = base_size + 3,
                                margin = margin(b = 5)),
      plot.subtitle = element_text(size = base_size, color = ad_colors$muted,
                                   margin = margin(b = 10)),
      plot.caption = element_text(size = base_size - 2, color = ad_colors$muted,
                                  hjust = 0),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = ad_colors$ink),
      panel.grid.major = element_line(color = ad_colors$grid, linewidth = 0.28),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = ad_colors$panel, color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = legend,
      legend.title = element_text(face = "bold"),
      legend.key = element_rect(fill = "white", color = NA),
      strip.text = element_text(face = "bold", color = ad_colors$ink),
      strip.background = element_rect(fill = "#EEF2F5", color = NA)
    )
}

format_p <- function(x) {
  ifelse(is.na(x), NA_character_,
         ifelse(x < 1e-3, scientific(x, digits = 2), number(x, accuracy = 0.001)))
}

save_pub_plot <- function(filename, plot, width, height, dpi = 450) {
  ext <- tools::file_ext(filename)
  if (tolower(ext) == "pdf") {
    ggsave(filename, plot, width = width, height = height, device = cairo_pdf)
  } else {
    ggsave(filename, plot, width = width, height = height, dpi = dpi,
           bg = "white")
  }
}

save_pub_pair <- function(stem, plot, width, height, dpi = 450, out_dir = ".") {
  save_pub_plot(file.path(out_dir, paste0(stem, ".pdf")), plot, width, height, dpi)
  save_pub_plot(file.path(out_dir, paste0(stem, ".png")), plot, width, height, dpi)
}

signed_label <- function(x, digits = 2) {
  sprintf(paste0("%+.", digits, "f"), x)
}

heat_colors <- function(n = 101) {
  colorRampPalette(c("#1E3A8A", "#F8FAFC", "#B91C1C"))(n)
}
