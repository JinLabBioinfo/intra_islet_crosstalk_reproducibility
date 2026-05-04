#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig5
# Rscript Fig5a.CellChat.DropSeq.Plot.r

library(ggplot2)
library(dplyr)
library(gridExtra)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

setwd(script_dir)

source_data_root <- normalizePath(
  file.path(script_dir, "..", "..", "SourceData"),
  mustWork = TRUE
)

psc_to_beta_file <- file.path(
  source_data_root,
  "Fig5a.1.SCtoBeta",
  "Fig5.PSC_To_Beta.DropSeq.diff_prob.T2D_minus_NonT2D.Final.csv"
)
beta_to_psc_file <- file.path(
  source_data_root,
  "Fig5a.2.BetaToSC",
  "Fig5.Beta_To_PSC.DropSeq.diff_prob.T2D_minus_NonT2D.Final.csv"
)

PSC_To_Beta_plot_df <- read.csv(psc_to_beta_file, stringsAsFactors = FALSE)
Beta_To_PSC_plot_df <- read.csv(beta_to_psc_file, stringsAsFactors = FALSE)

prepare_plot_df <- function(plot_df) {
  plot_df %>%
    mutate(
      interaction_name_2 = factor(
        interaction_name_2,
        levels = rev(interaction_name_2)
      ),
      diff_prob_plot = ifelse(
        is.na(diff_prob.T2D_minus_NonT2D),
        0,
        diff_prob.T2D_minus_NonT2D
      )
    )
}

PSC_To_Beta_plot_df <- prepare_plot_df(PSC_To_Beta_plot_df)
Beta_To_PSC_plot_df <- prepare_plot_df(Beta_To_PSC_plot_df)

all_diff_prob <- c(
  PSC_To_Beta_plot_df$diff_prob.T2D_minus_NonT2D,
  Beta_To_PSC_plot_df$diff_prob.T2D_minus_NonT2D
)

min_val <- min(all_diff_prob, na.rm = TRUE)
max_val <- max(all_diff_prob, na.rm = TRUE)

common_color_scale <- scale_color_gradient2(
  low = "blue",
  mid = "white",
  high = "red",
  midpoint = 0,
  limits = c(min_val, max_val),
  name = "T2D - NonT2D probability"
)

base_theme <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(color = "black", hjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.grid = element_blank()
  )

make_interaction_plot <- function(plot_df, title, label_side = c("left", "right")) {
  label_side <- match.arg(label_side)

  plot <- ggplot(
    plot_df,
    aes(x = interaction_name_2, y = 0, color = diff_prob_plot)
  ) +
    geom_point(size = 6) +
    common_color_scale +
    labs(title = title, x = NULL, y = NULL) +
    coord_flip() +
    base_theme +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )

  if (label_side == "right") {
    plot <- plot +
      scale_x_discrete(position = "top") +
      theme(axis.text.y = element_text(hjust = 0))
  }

  plot
}

Beta_To_PSC_plot <- make_interaction_plot(
  Beta_To_PSC_plot_df,
  title = "Beta to PSC",
  label_side = "left"
)

PSC_To_Beta_plot <- make_interaction_plot(
  PSC_To_Beta_plot_df,
  title = "PSC to Beta",
  label_side = "right"
)

legend_plot <- ggplot(
  Beta_To_PSC_plot_df,
  aes(x = interaction_name_2, y = 0, color = diff_prob_plot)
) +
  geom_point(size = 6) +
  common_color_scale +
  coord_flip() +
  base_theme +
  theme(legend.position = "right")

get_legend <- function(plot) {
  plot_grob <- ggplotGrob(plot)
  legend_index <- which(sapply(plot_grob$grobs, function(x) x$name) == "guide-box")

  if (length(legend_index) == 0) {
    stop("No legend was found in the plot.")
  }

  plot_grob$grobs[[legend_index[1]]]
}

shared_legend <- get_legend(legend_plot)

combined_plot <- arrangeGrob(
  Beta_To_PSC_plot,
  PSC_To_Beta_plot,
  shared_legend,
  ncol = 3,
  widths = c(1, 1, 0.35)
)

output_prefix <- "Fig5a.CellChat.DropSeq.diff_prob.T2D_minus_NonT2D.combined_legend"

ggsave(
  filename = paste0(output_prefix, ".pdf"),
  plot = combined_plot,
  device = "pdf",
  width = 10,
  height = 5,
  units = "in"
)

ggsave(
  filename = paste0(output_prefix, ".png"),
  plot = combined_plot,
  width = 10,
  height = 5,
  units = "in",
  dpi = 300
)
