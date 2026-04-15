# Author: Jiachen Sun

source("/home/jxs2269/islet_cosmx/for_jasmine/src/pipeline_reproducibility_common/_utility_funs.r")

# Load annotated merged objects
obj_merged_ak090t <- LoadSeuratRds("for_jasmine/obj_merged_AK090T_annotated.rds")
obj_merged_al189t <- LoadSeuratRds("for_jasmine/obj_merged_AL189T_annotated.rds")
obj_merged_am010 <- LoadSeuratRds("for_jasmine/obj_merged_AM010_annotated.rds")

# Endocrine/Exocrine subsets and set assay
obj_menv_ak090t <- subset(obj_merged_ak090t, subset = region == "Endocrine")
obj_menv_al189t <- subset(obj_merged_al189t, subset = region == "Endocrine")
obj_menv_am010 <- subset(obj_merged_am010, subset = region == "Endocrine")

obj_exo_ak090t <- subset(obj_merged_ak090t, subset = region == "Exocrine")
obj_exo_al189t <- subset(obj_merged_al189t, subset = region == "Exocrine")
obj_exo_am010 <- subset(obj_merged_am010, subset = region == "Exocrine")

# Tag geolocation for abundance summaries
obj_exo_ak090t$geoloc <- "ak090t_exo"
obj_exo_al189t$geoloc <- "al189t_exo"
obj_exo_am010$geoloc <- "am010_exo"

# Coarse mapping from detailed cell_type to macro groups for abundance summaries
cell_type_mappings <- c(
    "Alpha" = "Endocrine",
    "Beta" = "Endocrine",
    "Delta" = "Endocrine",
    "Gamma" = "Endocrine",
    "aStellate" = "aStellate",
    "qStellate" = "qStellate",
    "Ductal" = "Nonendocrine",
    "Acinar" = "Nonendocrine",
    "Endothelial" = "Nonendocrine",
    "Immune" = "Nonendocrine",
    "Inconclusive" = "Nonendocrine",
    "Unclassified" = "Unclassified"
)

# Calculate stellate abundance within endocrine islets: c2 (islet + coat) vs c0 (core)
obj_c2_ak090t <- subset(obj_menv_ak090t, cells = extract_barcodes(readRDS("for_jasmine/c2_list_AK090T.rds")))
obj_c2_ak090t$geoloc <- "ak090t_c12"
obj_c2_ak090t$geoloc[colnames(obj_c2_ak090t) %fin% extract_barcodes(readRDS("for_jasmine/c0_list_AK090T.rds"))] <- "ak090t_c0"

obj_c2_al189t <- subset(obj_menv_al189t, cells = extract_barcodes(readRDS("for_jasmine/c2_list_AL189T.rds")))
obj_c2_al189t$geoloc <- "al189t_c12"
obj_c2_al189t$geoloc[colnames(obj_c2_al189t) %fin% extract_barcodes(readRDS("for_jasmine/c0_list_AL189T.rds"))] <- "al189t_c0"

obj_c2_am010 <- subset(obj_menv_am010, cells = extract_barcodes(readRDS("for_jasmine/c2_list_AM010.rds")))
obj_c2_am010$geoloc <- "am010_c12"
obj_c2_am010$geoloc[colnames(obj_c2_am010) %fin% extract_barcodes(readRDS("for_jasmine/c0_list_AM010.rds"))] <- "am010_c0"

# Compute macro-group abundances by region/geoloc
abd_df <- cbind(
    return_cell_abd(obj_exo_ak090t, cell_type_mappings),
    return_cell_abd(obj_exo_al189t, cell_type_mappings)[, -1],
    return_cell_abd(obj_exo_am010, cell_type_mappings)[, -1],
    return_cell_abd(obj_c2_ak090t, cell_type_mappings)[, -1],
    return_cell_abd(obj_c2_al189t, cell_type_mappings)[, -1],
    return_cell_abd(obj_c2_am010, cell_type_mappings)[, -1]
)
rownames(abd_df) <- abd_df$cell_type
abd_df <- abd_df[, -1][c("aStellate", "qStellate"), c("ak090t_exo", "al189t_exo", "am010_exo", "ak090t_c12", "al189t_c12", "am010_c12", "ak090t_c0", "al189t_c0", "am010_c0")]
write.csv(round(sweep(abd_df, 2, colSums(abd_df), FUN = "/"), 4), file = "for_jasmine/stellate_abundance.csv", quote = FALSE)
