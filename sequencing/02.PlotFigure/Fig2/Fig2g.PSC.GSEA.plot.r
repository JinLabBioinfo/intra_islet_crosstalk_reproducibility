#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig2
# Rscript Fig2g.PSC.GSEA.plot.r

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
input_file <- file.path(source_data_dir, "PSC.00.combined_T2D_Aging_BMI.csv")

setwd(script_dir)

combined_table <- read.csv(input_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
colnames(combined_table) <- c(
  "pval_T2D",
  "NES_T2D",
  "pval_Aging_Healthy",
  "NES_Aging_Healthy",
  "pval_BMI_Healthy",
  "NES_BMI_Healthy"
)

df_long <- data.frame(
  Term = rownames(combined_table),
  NES_T2D = combined_table$NES_T2D,
  NES_Aging_Healthy = combined_table$NES_Aging_Healthy,
  NES_BMI_Healthy = combined_table$NES_BMI_Healthy,
  logp_T2D = -log10(combined_table$pval_T2D + 0.001),
  logp_Aging_Healthy = -log10(combined_table$pval_Aging_Healthy + 0.001),
  logp_BMI_Healthy = -log10(combined_table$pval_BMI_Healthy + 0.001)
)

set.seed(123)
num_clusters <- 2
nes_matrix <- combined_table[, c("NES_T2D", "NES_Aging_Healthy")]
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
df_plot$Shape <- ifelse(df_plot$NES > 0, "Up", "Down")
df_plot$Shape <- factor(df_plot$Shape, levels = c("Up", "Down"))

manual_term_order <- list(
  "1" = c(
    "GOBP_COLLAGEN_METABOLIC_PROCESS",
    "GAVISH_3CA_METAPROGRAM_FIBROBLASTS_CAF_4",
    "GOBP_CELL_PROJECTION_MORPHOGENESIS",
    "GOBP_MUSCLE_TISSUE_DEVELOPMENT",
    "GOMF_PROTEIN_KINASE_ACTIVITY",
    "GOBP_MICROTUBULE_CYTOSKELETON_ORGANIZATION_INVOLVED_IN_MITOSIS",
    "GOBP_MITOTIC_CELL_CYCLE_PROCESS",
    "GOBP_SPINDLE_LOCALIZATION",
    "MURARO_PANCREAS_MESENCHYMAL_STROMAL_CELL",
    "REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS",
    "GOBP_CELL_SUBSTRATE_JUNCTION_ORGANIZATION",
    "GOBP_POSITIVE_REGULATION_OF_RNA_METABOLIC_PROCESS",
    "GOBP_REGULATION_OF_CELL_SUBSTRATE_ADHESION",
    "GOBP_CELL_MATRIX_ADHESION",
    "GOMF_FIBRONECTIN_BINDING"
  ),
  "2" = c(
    "GAVISH_3CA_METAPROGRAM_FIBROBLASTS_LIPID_METABOLISM",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "GOBP_PROTEIN_REFOLDING",
    "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",
    "GOBP_CELLULAR_RESPIRATION",
    "GOCC_LATE_ENDOSOME",
    "GOMF_OXIDOREDUCTASE_ACTIVITY",
    "GOMF_SIGNALING_RECEPTOR_REGULATOR_ACTIVITY",
    "HALLMARK_ADIPOGENESIS",
    "GOBP_CELLULAR_RESPONSE_TO_TOXIC_SUBSTANCE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "APPIERTO_RESPONSE_TO_FENRETINIDE_DN",
    "GOBP_REGULATION_OF_LIPID_STORAGE",
    "GOMF_LIPID_BINDING",
    "GOMF_RETINAL_DEHYDROGENASE_ACTIVITY"
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

df_plot_ordered$Cluster <- factor(df_plot_ordered$Cluster, levels = c("1", "2"))
df_plot_ordered$ClusterTerm <- with(df_plot_ordered, interaction(Cluster, Term, sep = "_"))
df_plot_ordered$ClusterTerm <- factor(
  df_plot_ordered$ClusterTerm,
  levels = unique(df_plot_ordered$ClusterTerm[order(df_plot_ordered$Cluster, df_plot_ordered$Term)])
)
df_plot_ordered$Condition <- factor(
  df_plot_ordered$Condition,
  levels = rev(c("T2D", "Aging_Healthy", "BMI_Healthy"))
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

ggsave("Fig2g.PSC.GSEA.plot.png", p, width = 12, height = 6)
ggsave("Fig2g.PSC.GSEA.plot.pdf", p, width = 12, height = 6, device = cairo_pdf)
