#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig1
# Rscript Fig1e.AllCellTypeCirclizeHeatMap.r

# Load environment script when available
optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(dplyr)
library(plyr)
library(circlize)
library(ComplexHeatmap)
library(grid)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

source_data_dir <- normalizePath(
  file.path(script_dir, "..", "..", "SourceData", "Fig1e.AllCellTypeCirclizeHeatMap"),
  mustWork = FALSE
)
raw_data_dir <- file.path(source_data_dir, "RawData")
gene_dir <- file.path(source_data_dir, "genes")

setwd(script_dir)

get_merged_data <- function(raw_data_dir) {
  df.m <- data.frame()

  file_list <- list.files(
    path = raw_data_dir,
    pattern = "^Merged\\.T2D\\..*\\.csv$",
    full.names = TRUE
  )

  for (file in file_list) {
    group <- sub("^Merged\\.T2D\\.(.*)\\.csv$", "\\1", basename(file))
    df <- read.csv(file, stringsAsFactors = FALSE, row.names = 1)
    colnames(df) <- c("DropSeq T2D status (All)", "PANCDB T2D status (All)")
    df$gene <- rownames(df)
    df$CellType <- group
    df.m <- rbind(df.m, df)
  }

  df.m
}

get_gene_lists <- function(gene_dir) {
  gene_files <- list.files(path = gene_dir, pattern = "\\.genes$", full.names = TRUE)
  gene_lists <- list()

  for (file in gene_files) {
    celltype <- sub("\\.genes$", "", basename(file))
    gene_lists[[celltype]] <- readLines(file)
  }

  gene_lists
}

order_by_gene_lists <- function(df.m, gene_lists) {
  df.m.ordered_list <- list()

  for (ct in names(gene_lists)) {
    gene_list <- gene_lists[[ct]]
    df_sub <- df.m[df.m$CellType == ct & df.m$gene %in% gene_list, ]
    df_sub$gene <- factor(df_sub$gene, levels = gene_list, ordered = TRUE)
    df_sub <- df_sub[order(df_sub$gene), ]
    df_sub$Rank_by_List <- seq_len(nrow(df_sub))
    df.m.ordered_list[[ct]] <- df_sub
  }

  do.call(rbind, df.m.ordered_list)
}

draw_circular_heatmap <- function(plot_data, include_labels = TRUE, file_type = "pdf", file_name) {
  if (file_type == "pdf") {
    pdf(file_name, width = 8, height = 8)
  } else {
    png(file_name, width = 7000, height = 7000, res = 1200)
  }

  col_fun <- colorRamp2(c(-1, 0, 1), c("#2166ac", "#f7f7f7", "#b2182b"))

  celltype_order <- c("Beta", "Alpha", "Delta", "PP", "PSC", "Duct", "Endoth", "Acinar")
  color_vector <- c("#1F78B4", "#A6CEE3", "#33A02C", "#B2DF8A", "#FB9A99", "#FF7F00", "#FDBF6F", "#E31A1C")
  color_map <- setNames(color_vector, celltype_order)

  plot_data$CellType <- mapvalues(plot_data$CellType, from = c("Endothelial"), to = c("Endoth"))
  plot_data$CellType <- factor(plot_data$CellType, levels = celltype_order)
  split.m <- plot_data$CellType

  circos.par(
    start.degree = 90,
    track.margin = c(0, 0.002),
    cell.padding = c(0, 0, 0, 0),
    gap.after = c(rep(4, length(celltype_order) - 1), 40)
  )

  part.row <- 1 / 28

  circos.heatmap(
    plot_data[, c("DropSeq T2D status (All)", "PANCDB T2D status (All)")],
    split = split.m,
    col = col_fun,
    cluster = FALSE,
    track.height = 15 * part.row,
    bg.border = "black",
    bg.lwd = 1,
    bg.lty = 1,
    show.sector.labels = FALSE
  )

  if (include_labels) {
    circos.track(
      track.index = get.current.track.index(),
      panel.fun = function(x, y) {
        if (CELL_META$sector.index == "Acinar") {
          labels <- c("PANCDB T2D status (All)", "DropSeq T2D status (All)")
          circos.text(
            CELL_META$cell.xlim[2] + 10,
            seq_along(labels) - 0.5,
            labels,
            cex = 0.6,
            adj = c(0, 0.8),
            facing = "inside"
          )
        }
      },
      bg.border = NA
    )
  }

  circos.heatmap(
    as.character(split.m),
    split = split.m,
    col = color_map,
    track.height = 0.02
  )

  if (include_labels) {
    circos.trackPlotRegion(
      ylim = c(0, 1),
      track.height = 0.03,
      bg.border = NA,
      panel.fun = function(x, y) {
        circos.text(
          x = mean(CELL_META$xlim),
          y = 0.5,
          labels = CELL_META$sector.index,
          facing = "inside",
          niceFacing = TRUE,
          cex = 0.6
        )
      }
    )
  }

  lgd <- Legend(
    at = c(-1, 0, 1),
    col_fun = col_fun,
    title = "Slope",
    title_position = "topcenter"
  )
  draw(lgd, x = unit(20, "mm"), y = unit(3, "mm"), just = c("right", "bottom"))

  circos.clear()
  dev.off()
}

df.m <- get_merged_data(raw_data_dir)
gene_lists <- get_gene_lists(gene_dir)
df.m.all_ordered_by_list <- order_by_gene_lists(df.m, gene_lists)

draw_circular_heatmap(
  plot_data = df.m.all_ordered_by_list,
  include_labels = TRUE,
  file_type = "pdf",
  file_name = "Fig1e.AllCellTypeCirclizeHeatMap.pdf"
)

draw_circular_heatmap(
  plot_data = df.m.all_ordered_by_list,
  include_labels = TRUE,
  file_type = "png",
  file_name = "Fig1e.AllCellTypeCirclizeHeatMap.png"
)

draw_circular_heatmap(
  plot_data = df.m.all_ordered_by_list,
  include_labels = FALSE,
  file_type = "png",
  file_name = "Fig1e.AllCellTypeCirclizeHeatMap.NoLabel.png"
)
