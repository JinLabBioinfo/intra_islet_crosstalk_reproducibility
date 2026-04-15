#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig5
# Rscript Fig5a.2.BetaToSC.Plot.r

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(dplyr)
library(tidyr)
library(ggplot2)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

source_data_dir <- normalizePath(file.path(script_dir, "..", "..", "SourceData", "Fig5a.2.BetaToSC"), mustWork = FALSE)
selected_up_file <- file.path(source_data_dir, "Selected.Increased_in_T2D_Prob_pval.csv")
selected_dn_file <- file.path(source_data_dir, "Selected.Decreased_in_T2D_Prob_pval.csv")

setwd(script_dir)

selected.Increase <- read.csv(selected_up_file, stringsAsFactors = FALSE)
selected.Decreased <- read.csv(selected_dn_file, stringsAsFactors = FALSE)

log2fc_gg1r <- selected.Increase %>%
  select(interaction_name_2, dataset, prob) %>%
  pivot_wider(names_from = dataset, values_from = prob) %>%
  mutate(
    log2_fc = ifelse(!is.na(H) & !is.na(T2D), log2(T2D / H), NA_real_)
  ) %>%
  select(interaction_name_2, log2_fc)

log2fc_gg2r <- selected.Decreased %>%
  select(interaction_name_2, dataset, prob) %>%
  pivot_wider(names_from = dataset, values_from = prob) %>%
  mutate(
    log2_fc = ifelse(!is.na(H) & !is.na(T2D), log2(T2D / H), NA_real_)
  ) %>%
  select(interaction_name_2, log2_fc)

Increased.Signal <- c(
  "SPP1 - (ITGA5+ITGB1)",
  "SPP1 - (ITGAV+ITGB1)",
  "DLK1 - NOTCH2",
  "TGFB1 - (ACVR1B+TGFBR2)",
  "CNTN1 - NOTCH2",
  "BMP5 - (ACVR1+BMPR2)"
)

Decreased.Signal <- c(
  "RBP4 - STRA6",
  "SEMA4D - PLXNB2",
  "EDN3 - EDNRA",
  "EDN3 - EDNRB"
)

pick_fc <- function(df, sels, label) {
  df %>%
    filter(interaction_name_2 %in% sels) %>%
    mutate(direction = label)
}

plot_df <- bind_rows(
  pick_fc(log2fc_gg1r, Increased.Signal, "Increased in T2D"),
  pick_fc(log2fc_gg2r, Decreased.Signal, "Decreased in T2D")
) %>%
  mutate(
    interaction_name_2 = factor(
      interaction_name_2,
      levels = rev(c(Decreased.Signal, rev(Increased.Signal)))
    ),
    had_both = !is.na(log2_fc),
    log2_fc_plot = ifelse(is.na(log2_fc), 0, log2_fc)
  )

min_val <- min(plot_df$log2_fc, na.rm = TRUE)
max_val <- max(plot_df$log2_fc, na.rm = TRUE)

p <- ggplot(plot_df, aes(x = interaction_name_2, y = 0, color = log2_fc_plot)) +
  geom_point(size = 6) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limits = c(min_val, max_val),
    name = "log2(T2D / H)"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  coord_flip() +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave(
  filename = "Fig5a.2.BetaToSC.Plot.pdf",
  plot = p +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ),
  device = "pdf",
  width = 5,
  height = 5,
  units = "in"
)

p_right <- ggplot(plot_df, aes(x = interaction_name_2, y = 0, color = log2_fc_plot)) +
  geom_point(size = 6) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limits = c(min_val, max_val),
    name = "log2(T2D / H)"
  ) +
  labs(x = NULL, y = NULL) +
  coord_flip() +
  scale_x_discrete(position = "top") +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave(
  filename = "Fig5a.2.BetaToSC.Plot.right.pdf",
  plot = p_right +
    theme(
      axis.text.y = element_text(hjust = 0),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ),
  device = "pdf",
  width = 5,
  height = 5,
  units = "in"
)
