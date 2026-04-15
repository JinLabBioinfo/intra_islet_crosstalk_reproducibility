#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig1
# Rscript Fig1c.InHouseCohort.UMAP.r

# Load environment script when available
optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

# Set working directory and source-data path
script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

source_data_dir <- normalizePath(file.path(script_dir, "..", "..", "SourceData"), mustWork = FALSE)
input_file <- file.path(source_data_dir, "InHouseCohort.UMAP.csv")

setwd(script_dir)

# Load necessary libraries
require(RColorBrewer)
require(dplyr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

# Define the order of cell types
order <- c(
  "Beta",
  "Alpha",
  "Proliferating_Alpha",
  "Delta",
  "PP",
  "Epsilon",
  "Acinar",
  "Duct",
  "TFF1_Duct",
  "Activated_PSC",
  "Quiescent_PSC",
  "Proliferating_PSC",
  "Endothelial",
  "Macrophage",
  "Dendritic",
  "Schwann"
)

# Define colors for each cell type
mycolors <- c(
  "#1F78B4", "#A6CEE3", "#8acef3", "#33A02C", "#B2DF8A", "#008B8B",
  "#E31A1C", "#FF7F00", "#fca956", "#FB9A99", "#efafae", "#af6a69",
  "#FDBF6F", "#FF69B4", "#CAB2D6", "#6A3D9A"
)

data.infoandUMAP <- read.csv(input_file, stringsAsFactors = FALSE)
str(data.infoandUMAP)
table(data.infoandUMAP$CellType)
unique(data.infoandUMAP$CellType)

# Set the CellType factor levels
data.infoandUMAP$CellType <- factor(data.infoandUMAP$CellType, levels = order)

p1 <- ggplot(data.infoandUMAP, aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
  geom_point(size = 0.06) +
  scale_color_manual(values = mycolors) +
  labs(color = "Cell type") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggsave("DropSeq.LV1.all_Ordered.png", p1, width = 6, height = 4.5, dpi = 600)

pdf("01.DropSeq.LV1.all_Ordered.pdf", height = 4.5, width = 6)
print(p1)
dev.off()

write.csv(table(data.infoandUMAP$CellType), "DropSeqCellType.csv")
