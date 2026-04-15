#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig2
# Rscript Fig2b.2.KS_Test_p_value_BarPlot.r

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(ggplot2)
library(dplyr)
library(tidyr)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

source_data_dir <- normalizePath(file.path(script_dir, "..", "..", "SourceData", "Fig2b.DonorRisk"), mustWork = FALSE)
summary_file <- file.path(source_data_dir, "All_CellTypes_Summary.csv")

setwd(script_dir)

CellTypeBarPlot <- read.csv(summary_file, stringsAsFactors = FALSE)

CellTypeBarPlot <- CellTypeBarPlot %>%
  filter(!CellType %in% c("ActivatePSC", "QuiescentPSC")) %>%
  mutate(
    neg_log10_KS = -log10(pmax(KS_Test_p_value, .Machine$double.xmin))
  ) %>%
  arrange(desc(neg_log10_KS)) %>%
  mutate(
    CellType = factor(CellType, levels = rev(CellType))
  )

KS_bar_plot <- ggplot(CellTypeBarPlot, aes(x = neg_log10_KS, y = CellType)) +
  geom_col(width = 0.7, fill = "#5E616D") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "#474747", linewidth = 1) +
  geom_vline(xintercept = -log10(0.01), linetype = "dashed", color = "black", linewidth = 1) +
  annotate("text", x = -log10(0.05) - 0.18, y = 0.55, label = "-log10(0.05)", color = "#474747", size = 4, hjust = -0.01, vjust = 0.5) +
  annotate("text", x = -log10(0.01), y = 0.55, label = "-log10(0.01)", color = "black", size = 4, hjust = -0.05, vjust = 0.5) +
  labs(
    title = expression(-log[10]("KS Test p-value by Cell Type")),
    x = expression(-log[10]("KS_Test_p_value")),
    y = "Cell Type"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "none"
  ) +
  expand_limits(x = max(CellTypeBarPlot$neg_log10_KS) * 1.12)

ggsave(
  filename = "KS_Test_p_value_BarPlot.pdf",
  plot = KS_bar_plot,
  width = 10,
  height = 8
)
