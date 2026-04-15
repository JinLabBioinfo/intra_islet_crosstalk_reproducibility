#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig7
# Rscript Fig7g.GSEAlollipopPlot.r

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(ggplot2)
library(dplyr)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

source_data_dir <- normalizePath(file.path(script_dir, "..", "..", "SourceData", "Fig7g.GSEAlollipopPlot"), mustWork = FALSE)
input_file <- file.path(source_data_dir, "GSEA_Result.csv")

setwd(script_dir)

GSEA_Table <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

plot_df <- GSEA_Table %>%
  mutate(
    padj_cat = factor(
      case_when(
        padj > 0.05 ~ "> 0.05",
        padj <= 0.05 & padj > 0.01 ~ "0.05-0.01",
        padj <= 0.01 & padj > 0.0001 ~ "0.01-0.0001",
        padj <= 0.0001 ~ "< 0.0001"
      ),
      levels = c("> 0.05", "0.05-0.01", "0.01-0.0001", "< 0.0001")
    ),
    Handwriting = factor(Handwriting, levels = Handwriting[order(NES)])
  )

p <- ggplot(plot_df, aes(x = NES, y = Handwriting)) +
  geom_segment(
    aes(x = 0, xend = NES, y = Handwriting, yend = Handwriting, color = NES),
    size = 1
  ) +
  geom_point(
    aes(size = padj_cat, fill = NES),
    shape = 21, color = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_size_manual(
    name = "Adjusted p-value",
    values = c(6, 8, 10, 12),
    drop = FALSE
  ) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "NES"
  ) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    guide = "none"
  ) +
  theme_bw() +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    title = "GSEA Lollipop Plot"
  )

ggsave("Fig7g.GSEAlollipopPlot.pdf", plot = p, width = 10, height = 10)
