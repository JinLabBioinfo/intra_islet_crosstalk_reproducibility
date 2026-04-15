#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig1
# Rscript Fig1d.PANCDB_Cohort.UMAP.r

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
input_file <- file.path(source_data_dir, "PANCDB_Cohort.UMAP.csv")

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
  "CELA1_Acinar",
  "Duct",
  "TFF1_Duct",
  "SFRP5_Duct",
  "PSC_Activated",
  "PSC_Quiescent",
  "Endothelial",
  "Macrophage"
)

# Define colors for each cell type
mycolors <- c(
  "#1F78B4", "#A6CEE3", "#8acef3", "#33A02C", "#B2DF8A", "#008B8B",
  "#E31A1C", "#9B1E64", "#FF7F00", "#fca956", "#FBB612",
  "#FB9A99", "#efafae", "#FDBF6F", "#FF69B4"
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

ggsave("PANCDB_Cohort.UMAP.png", p1, width = 6, height = 4.5, dpi = 600)

pdf("Fig1d.PANCDB_Cohort.UMAP.pdf", height = 4.5, width = 6)
print(p1)
dev.off()

write.csv(table(data.infoandUMAP$CellType), "PANCDB_Cohort.CellType.csv")
