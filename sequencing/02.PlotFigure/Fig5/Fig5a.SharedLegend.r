#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig5
# Rscript Fig5a.SharedLegend.r

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

setwd(script_dir)

sc_to_beta_dir <- normalizePath(
  file.path(script_dir, "..", "..", "SourceData", "Fig5a.1.SCtoBeta"),
  mustWork = FALSE
)
beta_to_sc_dir <- normalizePath(
  file.path(script_dir, "..", "..", "SourceData", "Fig5a.2.BetaToSC"),
  mustWork = FALSE
)

get_log2_fc <- function(df) {
  df %>%
    select(interaction_name_2, dataset, prob) %>%
    pivot_wider(names_from = dataset, values_from = prob) %>%
    mutate(
      log2_fc = ifelse(!is.na(H) & !is.na(T2D), log2(T2D / H), NA_real_)
    ) %>%
    pull(log2_fc)
}

sc_to_beta_up <- read.csv(file.path(sc_to_beta_dir, "Selected.Up.csv"), stringsAsFactors = FALSE)
sc_to_beta_dn <- read.csv(file.path(sc_to_beta_dir, "Selected.DN.csv"), stringsAsFactors = FALSE)
beta_to_sc_up <- read.csv(file.path(beta_to_sc_dir, "Selected.Increased_in_T2D_Prob_pval.csv"), stringsAsFactors = FALSE)
beta_to_sc_dn <- read.csv(file.path(beta_to_sc_dir, "Selected.Decreased_in_T2D_Prob_pval.csv"), stringsAsFactors = FALSE)

all_fc <- c(
  get_log2_fc(sc_to_beta_up),
  get_log2_fc(sc_to_beta_dn),
  get_log2_fc(beta_to_sc_up),
  get_log2_fc(beta_to_sc_dn)
)

all_fc <- all_fc[!is.na(all_fc)]

if (length(all_fc) == 0) {
  stop("No valid log2 fold-change values were found.")
}

min_val <- min(all_fc)
max_val <- max(all_fc)

legend_plot <- ggplot(data.frame(x = 0, y = 0, fc = 0), aes(x = x, y = y, color = fc)) +
  geom_point(size = 4) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(min_val, max_val),
    name = "log2(T2D / H)"
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(color = "black", size = 11),
    legend.text = element_text(color = "black", size = 10)
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      barheight = unit(60, "mm"),
      barwidth = unit(5, "mm")
    )
  )

legend_grob <- get_legend(legend_plot)
legend_only <- ggdraw(legend_grob)

ggsave("Fig5a.SharedLegend.pdf", plot = legend_only, width = 2.0, height = 3.2, units = "in")
ggsave("Fig5a.SharedLegend.png", plot = legend_only, width = 2.0, height = 3.2, units = "in", dpi = 600)
