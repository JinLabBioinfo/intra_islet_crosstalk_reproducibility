# Author: Jiachen Sun

source("/home/jxs2269/islet_cosmx/for_jasmine/src/pipeline_reproducibility_common/_utility_funs.r")

# Run subpopulation correction for AL189T
obj_cosmx <- process_subpopulation(LoadSeuratRds("/mnt/jinstore/Archive01/pxg353/05.SpacialData/DropSeq.Endo.5kGenes.rds"), LoadSeuratRds("for_jasmine/obj_merged_AL189T_annotated.rds"))

# Correct cell type assignment (manual map by cluster index)
obj_cosmx <- RenameIdents(obj_cosmx, c(
    "0" = "Beta",
    "1" = "Alpha",
    "2" = "Alpha",
    "3" = "Beta",
    "4" = "Delta",
    "5" = "Gamma"
))

Idents(obj_cosmx) <- factor(Idents(obj_cosmx), levels = c(
    "Alpha",
    "Beta",
    "Delta",
    "Gamma"
))

# UMAP plot of corrected endocrine labels
ggsave("for_jasmine/umap_endocrine_AL189T.png",
    DimPlot(obj_cosmx,
        dims = c(2, 1),
        reduction = "umap",
        shuffle = TRUE,
        raster = FALSE,
        cols = c(
            "#EEEE00", # Alpha
            "#00EEEE", # Beta
            "#EEA9B8", # Delta
            "#66CD00" # Gamma
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

# Apply corrected endocrine labels back to full object
obj_cosmx_full <- LoadSeuratRds("for_jasmine/obj_merged_AL189T_annotated.rds")
Idents(obj_cosmx_full) <- Idents(obj_cosmx)
Idents(obj_cosmx_full) <- factor(Idents(obj_cosmx_full), levels = c(
    "Alpha",
    "Beta",
    "Delta",
    "Gamma",
    "aStellate",
    "qStellate",
    "Ductal",
    "Acinar",
    "Endothelial",
    "Immune"
))
obj_cosmx_full$cell_type <- Idents(obj_cosmx_full)

SaveSeuratRds(obj_cosmx_full, "for_jasmine/obj_merged_AL189T_annotated.rds")
