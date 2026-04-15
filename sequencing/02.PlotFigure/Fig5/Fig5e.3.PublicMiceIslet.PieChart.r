#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig5
# Rscript Fig5e.3.PublicMiceIslet.PieChart.r

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

require(RColorBrewer)
require(dplyr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(ggforce)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

source_data_dir <- normalizePath(file.path(script_dir, "..", "..", "SourceData", "Fig5e.3.PublicMiceIslet.PieChart"), mustWork = FALSE)
input_file <- file.path(source_data_dir, "NewDataset_CellTypeNew_MetaData.csv")

setwd(script_dir)

meta <- read.csv(input_file, stringsAsFactors = FALSE)

required_cols <- c("newCluster", "CellTypeNew")
if (!all(required_cols %in% colnames(meta))) {
  stop("Input metadata must contain: ", paste(required_cols, collapse = ", "))
}
meta <- meta[, required_cols]

celltype_counts <- table(meta$CellTypeNew)
celltype_df <- as.data.frame(celltype_counts)
colnames(celltype_df) <- c("CellType", "Count")
celltype_df$Percentage <- (celltype_df$Count / sum(celltype_df$Count)) * 100
celltype_df$Label <- paste0(celltype_df$CellType, " (", round(celltype_df$Percentage, 1), "%)")
celltype_df$CellType <- factor(
  celltype_df$CellType,
  levels = c("Beta", "Alpha", "Other Endocrine", "SC", "Other Non-Endocrine")
)
celltype_df <- celltype_df[order(celltype_df$CellType), ]

celltype_df$Fraction <- celltype_df$Count / sum(celltype_df$Count)
sc_target_pos <- "right"
target_angle <- c(right = 0, top = pi/2, left = pi, bottom = -pi/2)[sc_target_pos]

base_end <- 2 * pi * cumsum(celltype_df$Fraction)
base_start <- c(0, head(base_end, -1))
base_mid <- (base_start + base_end) / 2
sc_mid_base <- base_mid[celltype_df$CellType == "SC"][1]
angle_offset <- target_angle - sc_mid_base

celltype_df$end <- angle_offset + 2 * pi * cumsum(celltype_df$Fraction)
celltype_df$start <- c(angle_offset, head(celltype_df$end, -1))
celltype_df$mid <- (celltype_df$start + celltype_df$end) / 2

up_shift <- 0.08
celltype_df$x0 <- 0
celltype_df$y0 <- ifelse(celltype_df$CellType == "SC", up_shift, 0)

cell_colors <- c(
  "Beta" = "#1f78b4",
  "Alpha" = "#22A2C3",
  "Other Endocrine" = "#7da7d7",
  "SC" = "#806332",
  "Other Non-Endocrine" = "#E8B004"
)

p <- ggplot(celltype_df) +
  geom_arc_bar(
    aes(x0 = x0, y0 = y0, r0 = 0, r = 1, start = start, end = end, fill = CellType),
    color = "white"
  ) +
  coord_fixed() +
  scale_fill_manual(values = cell_colors, labels = celltype_df$Label) +
  labs(title = "Public Mice Islet Cell Type Composition", fill = "Cell Type") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )

ggsave("Fig5e.3.PublicMiceIslet.PieChart.pdf", p, width = 7, height = 6)
ggsave("Fig5e.3.PublicMiceIslet.PieChart.png", p, width = 7, height = 6, dpi = 300)
