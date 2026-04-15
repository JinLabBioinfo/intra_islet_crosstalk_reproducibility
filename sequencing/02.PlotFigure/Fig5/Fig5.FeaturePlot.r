#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig5
# Rscript Fig5.FeaturePlot.r /path/to/object.rds Rbp4
# Rscript Fig5.FeaturePlot.r /path/to/object.rds Rbp4 Ttr Stra6

lib_seurat52 = "~/.Rlibs_seurat52"
optional_extra_lib = Sys.getenv("FIG5_FEATUREPLOT_EXTRA_LIB", "")
if (nzchar(optional_extra_lib)) {
  .libPaths(c(lib_seurat52, optional_extra_lib, .libPaths()))
} else {
  .libPaths(c(lib_seurat52, .libPaths()))
}

# load ggplot2 first from the right library
library(ggplot2, lib.loc = lib_seurat52)

# then load the rest
library(Seurat)
library(patchwork)
library(rtracklayer)
library(dplyr)
library(combinat)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript Fig5.FeaturePlot.r /path/to/object.rds GENE1 [GENE2 ...]")
}

raw_object_path <- normalizePath(args[1], mustWork = TRUE)
plot_genes <- args[-1]

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

setwd(script_dir)
pbmc_filtered <- readRDS(raw_object_path)

# Save helper
save_plot <- function(p, name, w = 5, h = 5, dpi = 600) {
  ggsave(paste0(name, ".pdf"), plot = p, width = w, height = h, dpi = dpi)
  ggsave(paste0(name, ".png"), plot = p, width = w, height = h, dpi = dpi)
}

for (gene_name in plot_genes) {
  p1 <- FeaturePlot(
    pbmc_filtered,
    features = gene_name,
    cols = c("gray", "red"),
    order = TRUE,
    combine = FALSE
  )[[1]]

  save_plot(p1, paste0(gene_name, "_FeaturePlot"))
}
