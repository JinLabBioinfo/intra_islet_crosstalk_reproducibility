#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig1
# Rscript Fig1g.DnFunctionalMap.r

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(dplyr)
library(tidyr)
library(ggplot2)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

source_data_dir <- normalizePath(file.path(script_dir, "..", "..", "SourceData"), mustWork = FALSE)
input_file <- file.path(source_data_dir, "T2D_DN_signatures_Functional.csv")

setwd(script_dir)

bin_fdr <- function(x) {
  if (is.na(x)) {
    return(NA_character_)
  }
  if (x > 0.05) {
    return(">0.05")
  }
  if (x > 1e-2) {
    return("0.05")
  }
  if (x > 1e-4) {
    return("0.01")
  }
  if (x > 1e-6) {
    return("1e-04")
  }
  if (x > 1e-8) {
    return("1e-06")
  }
  if (x > 1e-10) {
    return("1e-08")
  }
  if (x > 1e-12) {
    return("1e-10")
  }
  if (x > 1e-16) {
    return("1e-12")
  }
  if (x > 1e-25) {
    return("1e-16")
  }
  if (x > 1e-40) {
    return("1e-25")
  }
  return("1e-40")
}

plot_colors <- c(
  ">0.05" = "white",
  "0.05" = "#becedb",
  "0.01" = "#87abc6",
  "1e-04" = "#7aa0bd",
  "1e-06" = "#6d95b4",
  "1e-08" = "#608aab",
  "1e-10" = "#3A6A8F",
  "1e-12" = "#356386",
  "1e-16" = "#2D5675",
  "1e-25" = "#20425C",
  "1e-40" = "#0B2132"
)

plot_levels <- names(plot_colors)
celltype_order <- c("Beta", "Alpha", "Delta", "PP", "PSC", "Endothelial", "Duct", "Acinar")

functional_df <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

if (grepl("Names", names(functional_df)[1])) {
  names(functional_df)[1] <- "Names"
}

functional_df$TermLabel <- paste0(functional_df$Names, " (", functional_df$Gene_Count, ")")

fdr_cols <- paste0(celltype_order, "_FDR")
count_cols <- paste0(celltype_order, ".Counts")
fdr_long <- functional_df %>%
  select(TermLabel, all_of(fdr_cols)) %>%
  pivot_longer(
    cols = all_of(fdr_cols),
    names_to = "CellType",
    values_to = "FDR"
  ) %>%
  mutate(
    CellType = sub("_FDR$", "", CellType),
    FDR_Bin = vapply(FDR, bin_fdr, character(1)),
    FDR_Bin = factor(FDR_Bin, levels = plot_levels)
  )

count_long <- functional_df %>%
  select(TermLabel, all_of(count_cols)) %>%
  pivot_longer(
    cols = all_of(count_cols),
    names_to = "CellType",
    values_to = "Counts"
  ) %>%
  mutate(CellType = sub("\\.Counts$", "", CellType))

plot_df <- fdr_long %>%
  left_join(count_long, by = c("TermLabel", "CellType"))

plot_df$CellType <- factor(plot_df$CellType, levels = celltype_order)
plot_df$TermLabel <- factor(plot_df$TermLabel, levels = rev(functional_df$TermLabel))

p1 <- ggplot(plot_df) +
  aes(CellType, TermLabel, fill = FDR_Bin) +
  geom_tile(colour = "black", linewidth = 0.3) +
  geom_text(aes(label = Counts), size = 5, color = "black") +
  scale_fill_manual("FDR", values = plot_colors, drop = FALSE) +
  labs(x = "Cell Type", y = "Terms") +
  theme_classic() +
  theme(
    legend.key = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 18, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black")
  )

pdf("Fig1g.DnFunctionalMap.pdf", height = 13, width = 10)
print(p1)
dev.off()

ggsave("Fig1g.DnFunctionalMap.png", p1, width = 11, height = 13, dpi = 600)
