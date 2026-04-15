# Author: Jiachen Sun

source("/home/jxs2269/islet_cosmx/for_jasmine/src/pipeline_reproducibility_common/_utility_funs.r")

# Calculate geo-adjacency score
geo_adj <- return_geoadj(readRDS("for_jasmine/c2_list_AL189T.rds"))
geo_adj$cell_type <- recode(geo_adj$cell_type, "aStellate" = "Stellate", "qStellate" = "Stellate")
write.csv(
  geo_adj %>%
    filter(geo_adj$cell_type %in% c("Alpha", "Beta", "Stellate", "Acinar", "Ductal")),
  "for_jasmine/geo_adjacency_AL189T.csv",
  quote = FALSE
)

# Load annotated and raw Seurat objects for AL189T
obj_merged <- LoadSeuratRds("for_jasmine/obj_merged_AL189T_annotated.rds")
obj_raw <- LoadSeuratRds("for_jasmine/obj_raw_AL189T.rds")

# UMAP of all cells; endocrine in gold, nonendocrine in royalblue (simplified palette)
ggsave("for_jasmine/umap_all_AL189T.png",
  DimPlot(obj_merged,
    dims = c(2, 1),
    reduction = "umap",
    shuffle = TRUE,
    raster = FALSE,
    cols = c(
      "#FFD700", # Alpha
      "#FFD700", # Beta
      "#FFD700", # Delta
      "#FFD700", # Gamma
      "#4169E1", # aStellate
      "#4169E1", # qStellate
      "#4169E1", # Ductal
      "#4169E1", # Acinar
      "#4169E1", # Endothelial
      "#4169E1" # Immune
    )
  ) & theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) &
    NoLegend() &
    coord_equal(),
  width = 4, height = 4, bg = "transparent"
)

# UMAP on nonendocrine only; re-SCT, PCA, UMAP to visualize finer separation
ggsave("for_jasmine/umap_nonendocrine_AL189T.png", DimPlot(
  quick_reanalyze_subset(obj_merged, c("aStellate", "qStellate", "Ductal", "Acinar", "Endothelial", "Immune")),
  dims = c(2, 1),
  reduction = "umap",
  shuffle = TRUE,
  raster = FALSE,
  cols = c(
    "#FF4500", # aStellate
    "#FF4500", # qStellate
    "#006400", # Ductal
    "#A020F0", # Acinar
    "#00868B", # Endothelial
    "#FF00FF" # Immune
  )
) & theme(
  axis.line = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank()
) &
  NoLegend() &
  coord_equal(), width = 4, height = 4, bg = "transparent")

ggsave("for_jasmine/dotplot_marker_AL189T.pdf",
  DotPlot(obj_merged,
    features = c(Alpha, Beta, Delta, Gamma, aStellate, qStellate, Ductal, Acinar, Endothelial, Immune),
    group.by = "cell_type",
    cols = c("#D3D3D3", "#FF0000")
  ) +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  width = 6.5, height = 7.5, bg = "transparent"
)

# Define identities and split into endocrine vs exocrine for downstream analyses
Idents(obj_merged) <- "cell_type"
obj_menv <- subset(obj_merged, subset = region == "Endocrine")
obj_exo <- subset(obj_merged, subset = region == "Exocrine")

# Whole-slide spatial map (simplified endocrine vs nonendocrine palette)
Idents(obj_raw) <- "Unclassified"
Idents(obj_raw) <- factor(Idents(obj_merged), levels = c(levels(Idents(obj_merged)), "Unclassified", "QC-Dropped"))
obj_raw <- annotate_qc_dropped(obj_raw, cell_type_levels = levels(Idents(obj_raw)))

ggsave("for_jasmine/image_full_AL189T.png",
  ImageDimPlot(obj_raw, boundaries = "segmentation", border.size = 0.1, border.color = "#000000", col = c(
    "Alpha" = "#FFD700",
    "Beta" = "#FFD700",
    "Delta" = "#FFD700",
    "Gamma" = "#FFD700",
    "aStellate" = "#4169E1",
    "qStellate" = "#4169E1",
    "Ductal" = "#4169E1",
    "Acinar" = "#4169E1",
    "Endothelial" = "#4169E1",
    "Immune" = "#4169E1",
    "Unclassified" = "#BEBEBE",
    "QC-Dropped" = "#BEBEBE"
  )),
  width = 16000, height = 12800, units = "px", limitsize = FALSE
)

# Visualize one representative region
ggsave("for_jasmine/image_fov_AL189T.png",
  ImageDimPlot(
    subset(obj_raw,
      subset = fov %in% c(609, 610, 611, 634, 635, 636, 655, 656, 657)
    ),
    boundaries = "segmentation",
    border.size = 0.1,
    border.color = "#000000",
    col = c(
      "Alpha" = "#EEEE00",
      "Beta" = "#00EEEE",
      "Delta" = "#EEA9B8",
      "Gamma" = "#66CD00",
      "aStellate" = "#FF4500",
      "qStellate" = "#FF4500",
      "Ductal" = "#006400",
      "Acinar" = "#A020F0",
      "Endothelial" = "#00868B",
      "Immune" = "#FF00FF",
      "Unclassified" = "#BEBEBE",
      "QC-Dropped" = "#BEBEBE"
    )
  ),
  width = 4000, height = 3200, units = "px", limitsize = FALSE
)

DefaultAssay(obj_menv) <- "RNA"
DefaultAssay(obj_exo) <- "RNA"

# Find signature genes for stellate cells (CosMx vs Drop-seq reference)
obj_stellate <- process_stellate(obj_menv, obj_exo, extract_barcodes(readRDS("for_jasmine/c1_list_AL189T.rds")))
sig_dropseq <- read.csv("/mnt/jinstore/JinLab02/sxz694/01.drop_seq/05.Islet/08.new/04.CellTypeCompare/06.subcelltype.compare/02.PSC.marker/acti_quei.spec.genes.csv", row.names = 1)[azimuth_stellate_sig_list, ]
sig_cosmx <- FindMarkers(obj_stellate, group.by = "region", ident.1 = "Exocrine", ident.2 = "Endocrine", min.pct = 0, logfc.threshold = 0)[azimuth_stellate_sig_list, ]
dual_logfc_barplot(sig_dropseq, sig_cosmx, azimuth_stellate_sig_list, "for_jasmine/barplot_log2fc_AL189T.pdf")
