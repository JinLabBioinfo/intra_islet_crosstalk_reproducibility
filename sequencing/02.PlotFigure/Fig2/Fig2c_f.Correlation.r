#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig2
# Rscript Fig2c_f.Correlation.r Beta
# Rscript Fig2c_f.Correlation.r Alpha
# Rscript Fig2c_f.Correlation.r Delta
# Rscript Fig2c_f.Correlation.r PP
# Rscript Fig2c_f.Correlation.r PSC
# Rscript Fig2c_f.Correlation.r Endothelial
# Rscript Fig2c_f.Correlation.r Duct
# Rscript Fig2c_f.Correlation.r Acinar

args <- commandArgs(trailingOnly = TRUE)
CellType <- if (length(args) >= 1) args[1] else "Beta"

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(rlang)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

source_data_dir <- normalizePath(file.path(script_dir, "..", "..", "SourceData", "GetStatisticsDonorScore", "00.RawObj"), mustWork = FALSE)
input_file <- file.path(source_data_dir, paste0("Disease_Status_", CellType, "_All.rds"))
donor_info_file <- normalizePath(file.path(script_dir, "..", "..", "SourceData", "GetStatisticsDonorScore", "Combine_DonorID.csv"), mustWork = FALSE)

setwd(script_dir)

if (!file.exists(input_file)) {
  stop("Input RDS not found for cell type: ", CellType)
}

Obj <- readRDS(input_file)
DonorData <- Obj[["beta.RNA.PCA.20bin.ob"]][["data.info.withbin"]]
CombineDonorID <- read.csv(donor_info_file, stringsAsFactors = FALSE)
colnames(CombineDonorID)[1] <- "Sample"
CombineDonorID$Sample <- trimws(CombineDonorID$Sample)
CombineDonorID <- CombineDonorID[, c(1, 8)]

DonorData$pseudo.index.balanced.scaled <- (DonorData$pseudo.index.balanced - min(DonorData$pseudo.index.balanced)) /
  (max(DonorData$pseudo.index.balanced) - min(DonorData$pseudo.index.balanced))
DonorData <- merge(DonorData, CombineDonorID, by = "Sample", all.x = TRUE)

T2D_Risk <- DonorData %>%
  group_by(Donor) %>%
  summarize(Median_Pseudo_Index = median(pseudo.index.balanced.scaled, na.rm = TRUE), .groups = "drop")

metadata <- DonorData %>%
  select(Donor, Disease_Status, BMI, HbA1C, Age, Sex) %>%
  distinct()

indexTable <- T2D_Risk %>%
  left_join(metadata, by = "Donor")

plot_age_vs_risk_by_group <- function(data, x_col, y_col, group_col,
                                      save_path = NULL, save_format = c("pdf", "png"),
                                      width = 8, height = 6, dpi = 600,
                                      y_limits = c(0.34, 0.7),
                                      colors = c("Y" = "red", "N" = "blue"),
                                      shapes = c("Y" = 16, "N" = 15),
                                      dotsize = 3) {
  save_format <- match.arg(save_format)

  data_clean <- data %>%
    mutate(
      !!sym(x_col) := as.numeric(ifelse(!!sym(x_col) == "N/A", NA, !!sym(x_col))),
      !!sym(group_col) := factor(!!sym(group_col))
    ) %>%
    drop_na(!!sym(x_col), !!sym(y_col), !!sym(group_col))

  lm_model_y <- lm(as.formula(paste(y_col, "~", x_col)), data = data_clean %>% filter(!!sym(group_col) == "Y"))
  lm_model_n <- lm(as.formula(paste(y_col, "~", x_col)), data = data_clean %>% filter(!!sym(group_col) == "N"))

  format_slope <- function(intercept, slope) {
    slope_sign <- ifelse(slope >= 0, "+", "-")
    paste0("y = ", format(intercept, digits = 3), " ", slope_sign, " ", format(abs(slope), digits = 3), "x")
  }

  formula_y <- format_slope(lm_model_y$coefficients[1], lm_model_y$coefficients[2])
  formula_n <- format_slope(lm_model_n$coefficients[1], lm_model_n$coefficients[2])
  pvalue_y <- format(summary(lm_model_y)$coefficients[2, 4], digits = 3)
  pvalue_n <- format(summary(lm_model_n)$coefficients[2, 4], digits = 3)

  subtitle <- paste(
    "Y: ", formula_y, " | p =", pvalue_y, "\n",
    "N: ", formula_n, " | p =", pvalue_n
  )

  p <- ggplot(data_clean, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(aes(color = !!sym(group_col), shape = !!sym(group_col)), size = dotsize, alpha = 1) +
    geom_smooth(
      method = "lm",
      aes(color = !!sym(group_col), fill = !!sym(group_col)),
      linetype = "solid",
      se = TRUE,
      alpha = 0.1
    ) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    coord_cartesian(ylim = y_limits) +
    labs(
      x = x_col,
      y = y_col,
      title = paste(x_col, "vs", y_col, "by", group_col),
      subtitle = subtitle
    ) +
    facet_wrap(as.formula(paste("~", group_col)), scales = "fixed") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
      strip.text = element_text(color = "black"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(color = "black"),
      plot.subtitle = element_text(color = "black"),
      legend.title = element_text(color = "black"),
      legend.text = element_text(color = "black")
    )

  if (!is.null(save_path)) {
    out_path <- paste0(save_path, ".", save_format)
    ggsave(out_path, plot = p, width = width, height = height, dpi = dpi, units = "in")
    message("Plot saved to: ", out_path)
  }

  return(p)
}

p_age <- plot_age_vs_risk_by_group(
  data = indexTable,
  x_col = "Age",
  y_col = "Median_Pseudo_Index",
  group_col = "Disease_Status",
  save_path = paste0(CellType, "_Age_vs_T2D_Risk"),
  save_format = "pdf"
)

p_bmi <- plot_age_vs_risk_by_group(
  data = indexTable,
  x_col = "BMI",
  y_col = "Median_Pseudo_Index",
  group_col = "Disease_Status",
  save_path = paste0(CellType, "_BMI_vs_T2D_Risk"),
  save_format = "pdf"
)

print(p_age)
print(p_bmi)
