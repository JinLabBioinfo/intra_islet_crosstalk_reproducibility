#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig7
# Rscript Fig7f.VolcanoPlot.r
# Rscript Fig7f.VolcanoPlot.r label
# Rscript Fig7f.VolcanoPlot.r nolabel
genes_to_label <- c("CYP26B1", "SYT16", "DHRS3", "CYP26A1", "RARB", 
"MMP1", "POSTN", "FGG", "CXCL5")
optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(ggplot2)
library(ggrepel)
library(ggbreak)

args <- commandArgs(trailingOnly = TRUE)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

source_data_dir <- normalizePath(file.path(script_dir, "..", "..", "SourceData", "Fig7f.VolcanoPlot"), mustWork = FALSE)
input_file <- file.path(source_data_dir, "Sig.HumanIslets_Sorted_CtrlRetinol.DEseq2.all.noNA.csv")

setwd(script_dir)

res_all_no_na <- read.csv(input_file, stringsAsFactors = FALSE)

if (!"Geneid" %in% colnames(res_all_no_na)) {
  stop("Column 'Geneid' was not found in the input file.")
}

res_all_no_na$significance <- with(
  res_all_no_na,
  ifelse(padj < 0.05 & abs(log2FoldChange) > 0.5, "Significant", "Not Significant")
)

res_all_no_na$color <- with(
  res_all_no_na,
  ifelse(
    significance == "Not Significant",
    "grey",
    ifelse(log2FoldChange < -0.5, "#158BB8", ifelse(log2FoldChange > 0.5, "red", "grey"))
  )
)

res_all_no_na$neg_log10_padj <- -log10(res_all_no_na$padj)

ymax <- max(res_all_no_na$neg_log10_padj, na.rm = TRUE)
upper_gap <- max(20, floor(ymax) - 2)
sig_line <- -log10(0.05)
label_mode <- if (length(args) >= 1) tolower(args[1]) else "label"
add_labels <- !label_mode %in% c("nolabel", "no_label", "none", "false", "0")

p_ggbreak <- ggplot(
  res_all_no_na,
  aes(x = log2FoldChange, y = neg_log10_padj, color = color)
) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_identity() +
  theme_bw() +
  labs(
    title = "Volcano Plot (label selected genes)",
    x = "Log2 Fold Change",
    y = "-Log10(FDR)"
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = sig_line, linetype = "dashed", color = "black") +
  scale_y_break(c(20, upper_gap), space = 0.05, scales = 0.7) +
  theme(
    text = element_text(color = "black"),
    panel.grid = element_blank()
  )

if (add_labels) {
  p_ggbreak <- p_ggbreak +
    geom_text_repel(
      data = subset(res_all_no_na, Geneid %in% genes_to_label),
      aes(label = Geneid),
      size = 3,
      box.padding = 0.3,
      max.overlaps = 20,
      color = "black"
    )
}

ggsave(
  "Fig7f.VolcanoPlot.pdf",
  plot = p_ggbreak,
  device = "pdf",
  width = 8,
  height = 10
)

ggsave(
  "Fig7f.VolcanoPlot.png",
  plot = p_ggbreak,
  device = "png",
  width = 8,
  height = 10,
  dpi = 600
)
