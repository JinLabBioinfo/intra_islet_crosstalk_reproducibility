#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig2
# Rscript Fig2b.1.BarplotForDonor.r

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
setwd(script_dir)

CellType <- c(
  "Beta",
  "PSC",
  "Alpha",
  "Acinar",
  "Delta",
  "PP",
  "Endothelial",
  "Duct"
)

readCSV <- c(
  file.path(source_data_dir, "Beta.Combine.T2D_Risk.csv"),
  file.path(source_data_dir, "PSC.Combine.T2D_Risk.csv"),
  file.path(source_data_dir, "Alpha.Combine.T2D_Risk.csv"),
  file.path(source_data_dir, "Acinar.Combine.T2D_Risk.csv"),
  file.path(source_data_dir, "Delta.Combine.T2D_Risk.csv"),
  file.path(source_data_dir, "PP.Combine.T2D_Risk.csv"),
  file.path(source_data_dir, "Endothelial.Combine.T2D_Risk.csv"),
  file.path(source_data_dir, "Duct.Combine.T2D_Risk.csv")
)

cell_data_list <- setNames(readCSV, CellType) |>
  lapply(read.csv, stringsAsFactors = FALSE, check.names = FALSE)

cell_data_list <- lapply(names(cell_data_list), function(celltype) {
  df <- cell_data_list[[celltype]]
  df$CellType <- celltype
  df
})

names(cell_data_list) <- CellType

cell_data_list_fixed <- lapply(names(cell_data_list), function(celltype) {
  df <- cell_data_list[[celltype]]
  if (names(df)[1] == "") {
    names(df)[1] <- "X"
  }
  colnames(df)[which(colnames(df) == celltype)] <- "Value"
  df[, c("X", "Donor", "Value", "Rank", "Group", "CellType")]
})

combined_df <- do.call(rbind, cell_data_list_fixed)

CellType_plot <- c(
  "Beta",
  "PSC",
  "Acinar",
  "Alpha",
  "Delta",
  "Endothelial",
  "PP",
  "Duct"
)

combined_df$CellType <- factor(combined_df$CellType, levels = CellType_plot)

pcombine <- ggplot(combined_df, aes(x = Rank, color = Group)) +
  geom_vline(aes(xintercept = Rank, color = Group), alpha = 0.7, linewidth = 1) +
  facet_grid(CellType ~ ., scales = "free_y", switch = "y") +
  scale_color_manual(values = c("H" = "blue", "T2D" = "red")) +
  theme_bw() +
  labs(title = NULL, x = "Rank (Sorted Donors)", y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    strip.text.y.left = element_text(angle = 0, size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray90", color = NA),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black")
  )

ggsave("Combined.DonorRisk.VerticalLines.pdf", pcombine, width = 5, height = 4)
ggsave("Combined.DonorRisk.VerticalLines.png", pcombine, width = 5, height = 4, dpi = 600)
