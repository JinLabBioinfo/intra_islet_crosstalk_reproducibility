#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig6
# Rscript Fig6c.1.RetinolDerivatives.r

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(tidyverse)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(usethis)
library(ggrepel)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(grid)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

source_data_dir <- normalizePath(file.path(script_dir, "..", "..", "SourceData", "Fig6c.1.RetinolDerivatives"), mustWork = FALSE)
setwd(script_dir)

ListData <- c(
  "RetinylEster" = file.path(source_data_dir, "01.RetinylEster.csv"),
  "11-cis-Retinol" = file.path(source_data_dir, "02.11-cis-Retinol.csv"),
  "11-cis-Retinal" = file.path(source_data_dir, "03.11-cis-Retinal.csv"),
  "all-trans-Retinal" = file.path(source_data_dir, "04.all-trans-Retinal.csv")
)

read_matrix <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  colnames(df)[1] <- "SampleID"

  df_clean <- df %>%
    mutate(across(-SampleID, ~ as.numeric(gsub(",", "", as.character(.x)))))

  mat <- df_clean %>%
    column_to_rownames("SampleID") %>%
    as.matrix()

  return(mat)
}

mat_list <- lapply(ListData, read_matrix)

OrderSampleType <- c("PSC_Medium", "PSC_Lysate", "EndoC_Medium", "EndoC_Lysate")
OrderCondition <- c("Ctrl", "RBP4", "Rol", "Combine")
OrderSampleName <- as.vector(t(outer(OrderSampleType, OrderCondition, paste, sep = "_")))
OrderCellType <- c("PSC", "EndoC")

my_cols <- c(
  "Ctrl" = "#4C566A",
  "RBP4" = "#1F78B4",
  "Rol" = "#E6A800",
  "Combine" = "#2E7D32"
)

calc_y_axis <- function(x) {
  y_max <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(y_max) || y_max <= 0) {
    y_max <- 1
  }

  y_exp <- floor(log10(y_max))
  y_unit <- 10^y_exp
  n_tick <- ceiling(y_max / y_unit)
  if (n_tick < 3 && y_exp > 0) {
    y_exp <- y_exp - 1
    y_unit <- 10^y_exp
    n_tick <- ceiling(y_max / y_unit)
  }
  y_breaks <- seq(0, n_tick * y_unit, by = y_unit)

  list(exp = y_exp, unit = y_unit, breaks = y_breaks)
}

global_y_cfg <- calc_y_axis(unlist(lapply(mat_list, as.numeric)))

for (met_name in names(mat_list)) {
  retinyl_mat <- mat_list[[met_name]]

  retinyl_df <- as.data.frame(retinyl_mat) %>%
    tibble::rownames_to_column("SampleType")

  retinyl_long <- retinyl_df %>%
    tidyr::pivot_longer(
      cols = -SampleType,
      names_to = "Condition",
      values_to = "Value"
    )

  retinyl_long <- retinyl_long %>%
    dplyr::mutate(
      Condition = dplyr::recode(
        Condition,
        "Ctrl_Med" = "Ctrl",
        "RBP4_Med" = "RBP4",
        "Rol_Med" = "Rol",
        "Combine_Med" = "Combine"
      ),
      Condition = factor(Condition, levels = OrderCondition),
      SampleType = factor(SampleType, levels = OrderSampleType),
      SampleName = paste(SampleType, Condition, sep = "_"),
      SampleName = factor(SampleName, levels = OrderSampleName),
      IsNA = is.na(Value),
      ValuePlot = dplyr::if_else(IsNA, 0, Value)
    )

  for (cell_type in OrderCellType) {
    group_sample_types <- paste0(cell_type, c("_Medium", "_Lysate"))

    plot_df <- retinyl_long %>%
      dplyr::filter(SampleType %in% group_sample_types) %>%
      dplyr::mutate(
        SampleType = factor(as.character(SampleType), levels = group_sample_types)
      )

    p <- ggplot(plot_df, aes(x = SampleType, y = ValuePlot, fill = Condition)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, na.rm = TRUE) +
      scale_fill_manual(values = my_cols) +
      scale_y_continuous(
        breaks = global_y_cfg$breaks,
        labels = function(x) as.character(round(x / global_y_cfg$unit)),
        limits = c(0, max(global_y_cfg$breaks))
      ) +
      labs(
        title = paste0(met_name, " - ", cell_type),
        x = "SampleType",
        y = paste0("Peak area (x10^", global_y_cfg$exp, ")"),
        fill = "Condition"
      ) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        strip.text = element_text(color = "black")
      )

    safe_name <- gsub("[^A-Za-z0-9]+", "_", met_name)
    ggsave(
      file.path(script_dir, paste0("Fig6c.1.", safe_name, "_", cell_type, "_barplot.pdf")),
      plot = p, width = 5, height = 3
    )
  }
}
