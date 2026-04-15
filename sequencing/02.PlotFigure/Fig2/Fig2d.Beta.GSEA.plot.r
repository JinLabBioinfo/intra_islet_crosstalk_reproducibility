#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig2
# Rscript Fig2d.Beta.GSEA.plot.r

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(ggplot2)
library(dplyr)
library(tidyr)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

source_data_dir <- normalizePath(file.path(script_dir, "..", "..", "SourceData", "Fig2d_g_GSEA"), mustWork = FALSE)
input_file <- file.path(source_data_dir, "Beta.combined_T2D_Aging_BMI.csv")

setwd(script_dir)

combined_table <- read.csv(input_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
colnames(combined_table) <- c(
  "pval_T2D",
  "NES_T2D",
  "pval_Aging",
  "NES_Aging",
  "pval_BMI",
  "NES_BMI"
)

df_long <- data.frame(
  Term = rownames(combined_table),
  NES_T2D = combined_table$NES_T2D,
  NES_Aging_Healthy = combined_table$NES_Aging,
  NES_BMI_Healthy = combined_table$NES_BMI,
  logp_T2D = -log10(combined_table$pval_T2D + 0.001),
  logp_Aging_Healthy = -log10(combined_table$pval_Aging + 0.001),
  logp_BMI_Healthy = -log10(combined_table$pval_BMI + 0.001)
)

set.seed(123)
num_clusters <- 2
nes_matrix <- combined_table[, c("NES_T2D", "NES_Aging")]
kmeans_result <- kmeans(nes_matrix, centers = num_clusters)
df_long$Cluster <- as.factor(kmeans_result$cluster)

df_plot_nes <- df_long %>%
  pivot_longer(
    cols = c("NES_T2D", "NES_Aging_Healthy", "NES_BMI_Healthy"),
    names_to = "Condition",
    values_to = "NES"
  )

df_plot_logp <- df_long %>%
  pivot_longer(
    cols = c("logp_T2D", "logp_Aging_Healthy", "logp_BMI_Healthy"),
    names_to = "Condition",
    values_to = "logp"
  )

df_plot_nes$Condition <- gsub("NES_", "", df_plot_nes$Condition)
df_plot_logp$Condition <- gsub("logp_", "", df_plot_logp$Condition)

df_plot <- merge(df_plot_nes, df_plot_logp, by = c("Term", "Condition", "Cluster"))

df_plot$Condition <- factor(df_plot$Condition, levels = c("BMI_Healthy", "Aging_Healthy", "T2D"))
df_plot$Shape <- ifelse(df_plot$NES > 0, "Up", "Down")
df_plot$Shape <- factor(df_plot$Shape, levels = c("Up", "Down"))

manual_term_order <- list(
  "2" = c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "WONG_ADULT_TISSUE_STEM_MODULE",
    "GOBP_POSITIVE_REGULATION_OF_RNA_METABOLIC_PROCESS",
    "GOBP_POSITIVE_REGULATION_OF_PROGRAMMED_CELL_DEATH",
    "GOMF_TRANSCRIPTION_REGULATOR_ACTIVITY",
    "ELVIDGE_HYPOXIA_BY_DMOG_UP",
    "GOBP_INTRACELLULAR_SIGNALING_CASSETTE",
    "GOBP_APOPTOTIC_PROCESS",
    "GOBP_REGULATION_OF_CELL_DIFFERENTIATION",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "REACTOME_SIGNALING_BY_INTERLEUKINS",
    "KEGG_TGF_BETA_SIGNALING_PATHWAY",
    "GOBP_REGULATION_OF_CELL_CYCLE",
    "GOBP_CELLULAR_RESPONSE_TO_STRESS",
    "WP_DNA_DAMAGE_RESPONSE"
  ),
  "1" = c(
    "SERVITJA_ISLET_HNF1A_TARGETS_DN",
    "GOCC_INNER_MITOCHONDRIAL_MEMBRANE_PROTEIN_COMPLEX",
    "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",
    "GOBP_OXIDATIVE_PHOSPHORYLATION",
    "GOMF_OXIDOREDUCTASE_ACTIVITY_ACTING_ON_NAD_P_H",
    "REACTOME_REGULATION_OF_BETA_CELL_DEVELOPMENT",
    "REACTOME_REGULATION_OF_GENE_EXPRESSION_IN_BETA_CELLS",
    "GOMF_PHOSPHATASE_INHIBITOR_ACTIVITY",
    "GOBP_RESPONSE_TO_RETINOIC_ACID",
    "GOBP_INSULIN_SECRETION",
    "HALLMARK_PANCREAS_BETA_CELLS",
    "HP_ABNORMAL_CIRCULATING_INSULIN_CONCENTRATION",
    "GOBP_POTASSIUM_ION_HOMEOSTASIS",
    "REACTOME_RETINOID_CYCLE_DISEASE_EVENTS",
    "HALLMARK_FATTY_ACID_METABOLISM"
  )
)

df_plot_ordered <- df_plot %>%
  group_split(Cluster) %>%
  lapply(function(df) {
    clust <- as.character(unique(df$Cluster))
    custom_levels <- manual_term_order[[clust]]
    df$Term <- factor(as.character(df$Term), levels = custom_levels)
    df
  }) %>%
  bind_rows()

df_plot_ordered$Cluster <- factor(df_plot_ordered$Cluster, levels = c("2", "1"))
df_plot_ordered$ClusterTerm <- with(df_plot_ordered, interaction(Cluster, Term, sep = "_"))
df_plot_ordered$ClusterTerm <- factor(
  df_plot_ordered$ClusterTerm,
  levels = unique(df_plot_ordered$ClusterTerm[order(df_plot_ordered$Cluster, df_plot_ordered$Term)])
)

arrow_up <- "⬆"
arrow_down <- "⬇"

p <- ggplot(df_plot_ordered, aes(y = Condition, x = ClusterTerm)) +
  geom_text(
    aes(label = ifelse(Shape == "Up", arrow_up, arrow_down), color = NES),
    size = 15
  ) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = range(df_plot_ordered$NES, na.rm = TRUE)
  ) +
  scale_x_discrete(labels = function(x) sub("^[^_]*_", "", x)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(color = "black", face = "bold"),
    plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(color = "black", size = 12),
    axis.title.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 10),
    legend.title = element_text(color = "black", size = 10),
    legend.text = element_text(color = "black", size = 9),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  labs(
    y = "Condition",
    x = "Hallmark Term",
    color = "NES"
  )

p <- p + theme(text = element_text(family = "DejaVu Sans"))

ggsave("Fig2d.Beta.GSEA.plot.png", p, width = 12, height = 6)
ggsave("Fig2d.Beta.GSEA.plot.pdf", p, width = 12, height = 6, device = cairo_pdf)
