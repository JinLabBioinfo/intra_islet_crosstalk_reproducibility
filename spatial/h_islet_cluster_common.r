# Author: Jiachen Sun

source("/home/jxs2269/islet_cosmx/for_jasmine/src/pipeline_reproducibility_common/_utility_funs.r")

# Load endocrine-only subsets
obj_menv_ak090t <- subset(LoadSeuratRds("for_jasmine/obj_merged_AK090T_annotated.rds"), subset = region == "Endocrine")
obj_menv_al189t <- subset(LoadSeuratRds("for_jasmine/obj_merged_AL189T_annotated.rds"), subset = region == "Endocrine")
obj_menv_am010 <- subset(LoadSeuratRds("for_jasmine/obj_merged_AM010_annotated.rds"), subset = region == "Endocrine")

# Map detailed cell types to macro groups for clustering/plots
obj_menv_ak090t <- RenameIdents(obj_menv_ak090t, c(
    "Alpha" = "Endocrine",
    "Beta" = "Endocrine",
    "Delta" = "Endocrine",
    "Gamma" = "Endocrine",
    "aStellate" = "Stellate",
    "qStellate" = "Stellate",
    "Ductal" = "Nonendocrine",
    "Acinar" = "Nonendocrine",
    "Inconclusive" = "Nonendocrine"
))
obj_menv_al189t <- RenameIdents(obj_menv_al189t, c(
    "Alpha" = "Endocrine",
    "Beta" = "Endocrine",
    "Delta" = "Endocrine",
    "Gamma" = "Endocrine",
    "aStellate" = "Stellate",
    "qStellate" = "Stellate",
    "Endothelial" = "Nonendocrine",
    "Ductal" = "Nonendocrine",
    "Acinar" = "Nonendocrine",
    "Immune" = "Nonendocrine"
))
obj_menv_am010 <- RenameIdents(obj_menv_am010, c(
    "Alpha" = "Endocrine",
    "Beta" = "Endocrine",
    "Delta" = "Endocrine",
    "Gamma" = "Endocrine",
    "aStellate" = "Stellate",
    "qStellate" = "Stellate",
    "Ductal" = "Nonendocrine",
    "Acinar" = "Nonendocrine",
    "Immune" = "Nonendocrine",
    "Inconclusive" = "Nonendocrine"
))

palette <- c(
    "Endocrine" = "#FFD700",
    "Stellate" = "#FF4500",
    "Nonendocrine" = "#4169E1",
    "Unclassified" = "#BEBEBE"
)

# Order identities consistently for both samples
Idents(obj_menv_ak090t) <- factor(Idents(obj_menv_ak090t), levels = names(palette))
Idents(obj_menv_al189t) <- factor(Idents(obj_menv_al189t), levels = names(palette))
Idents(obj_menv_am010) <- factor(Idents(obj_menv_am010), levels = names(palette))

obj_menv_ak090t$cell_type <- Idents(obj_menv_ak090t)
obj_menv_al189t$cell_type <- Idents(obj_menv_al189t)
obj_menv_am010$cell_type <- Idents(obj_menv_am010)

obj_raw_ak090t <- LoadSeuratRds("for_jasmine/obj_raw_AK090T.rds")
obj_raw_al189t <- LoadSeuratRds("for_jasmine/obj_raw_AL189T.rds")
obj_raw_am010 <- LoadSeuratRds("for_jasmine/obj_raw_AM010.rds")

# Initialize raw Idents and then sync from endocrine-only objects
Idents(obj_raw_ak090t) <- "Unclassified"
Idents(obj_raw_al189t) <- "Unclassified"
Idents(obj_raw_am010) <- "Unclassified"

Idents(obj_raw_ak090t) <- Idents(obj_menv_ak090t)
Idents(obj_raw_al189t) <- Idents(obj_menv_al189t)
Idents(obj_raw_am010) <- Idents(obj_menv_am010)

# Build islet cluster lists from raw + endocrine with cached graph/threshold
cluster_list_ak090t <- get_islet_cluster_list(
    obj_raw_ak090t,
    obj_menv_ak090t,
    "for_jasmine/graph_raw_AK090T.rds",
    "for_jasmine/c2_list_AK090T.rds",
    "for_jasmine/threshold_AK090T.rds",
    "AK090T",
    names(palette)
)
cluster_list_al189t <- get_islet_cluster_list(
    obj_raw_al189t,
    obj_menv_al189t,
    "for_jasmine/graph_raw_AL189T.rds",
    "for_jasmine/c2_list_AL189T.rds",
    "for_jasmine/threshold_AL189T.rds",
    "AL189T",
    names(palette)
)
cluster_list_am010 <- get_islet_cluster_list(
    obj_raw_am010,
    obj_menv_am010,
    "for_jasmine/graph_raw_AM010.rds",
    "for_jasmine/c2_list_AM010.rds",
    "for_jasmine/threshold_AM010.rds",
    "AM010",
    names(palette)
)

# Calculate shape statistics and export combined table
shape_stat_cluster_ak090t <- return_shape_stat(cluster_list_ak090t, "Endocrine")
shape_stat_cluster_al189t <- return_shape_stat(cluster_list_al189t, "Endocrine")
shape_stat_cluster_am010 <- return_shape_stat(cluster_list_am010, "Endocrine")
shape_stat_cluster_ak090t$sample <- "AK090T"
shape_stat_cluster_al189t$sample <- "AL189T"
shape_stat_cluster_am010$sample <- "AM010"
shape_stat_cluster <- rbind(shape_stat_cluster_ak090t, shape_stat_cluster_al189t, shape_stat_cluster_am010) %>%
    select(-ends_with("_circularity"))
write.csv(shape_stat_cluster, "for_jasmine/shape_stat_ic.csv", quote = FALSE)

# Plot individual islet clusters, sorted by size
dir.create("for_jasmine/image_ic_AK090T", recursive = TRUE)
dir.create("for_jasmine/image_ic_AL189T", recursive = TRUE)
dir.create("for_jasmine/image_ic_AM010", recursive = TRUE)

ic_index_ak090t <- c(7, 18, 22, 28, 34, 36, 42, 44, 59, 73, 91, 93)
ic_index_al189t <- c(17, 40, 51, 67, 70, 86, 114, 130, 148, 188, 189, 216)
ic_index_am010 <- c(45, 55, 57, 89, 106, 119, 155, 161, 168, 179, 207, 248)

max_size <- max(
    c(
        calculate_max_size(
            ic_index_ak090t,
            obj_raw_ak090t,
            cluster_list_ak090t
        ), calculate_max_size(
            ic_index_al189t,
            obj_raw_al189t,
            cluster_list_al189t
        ), calculate_max_size(
            ic_index_am010,
            obj_raw_am010,
            cluster_list_am010
        )
    )
)

individual_islet_plot_wrapper(
    ic_index_ak090t,
    annotate_qc_dropped(obj_raw_ak090t, cell_type_levels = names(palette)),
    cluster_list_ak090t,
    c(palette, "QC-Dropped" = "#7F7F7F"),
    "for_jasmine/image_ic_AK090T/image_ic_AK090T",
    max_size
)
individual_islet_plot_wrapper(
    ic_index_al189t,
    annotate_qc_dropped(obj_raw_al189t, cell_type_levels = names(palette)),
    cluster_list_al189t,
    c(palette, "QC-Dropped" = "#7F7F7F"),
    "for_jasmine/image_ic_AL189T/image_ic_AL189T",
    max_size
)
individual_islet_plot_wrapper(
    ic_index_am010,
    annotate_qc_dropped(obj_raw_am010, cell_type_levels = names(palette)),
    cluster_list_am010,
    c(palette, "QC-Dropped" = "#7F7F7F"),
    "for_jasmine/image_ic_AM010/image_ic_AM010",
    max_size
)
