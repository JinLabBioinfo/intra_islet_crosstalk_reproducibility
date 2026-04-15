# Author: Jiachen Sun

source("/home/jxs2269/islet_cosmx/for_jasmine/src/pipeline_reproducibility_common/_utility_funs.r")

obj <- merge_object_pair(LoadSeuratRds("for_jasmine/obj_menv_AL189T.rds"), LoadSeuratRds("for_jasmine/obj_exo_AL189T.rds"))

# Global clustering on merged object (SCT -> PCA -> graph -> SLM -> UMAP)
obj <- SCTransform(obj, clip.range = c(-10, 10), variable.features.n = 2000) %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30, reduction = "pca") %>%
    FindClusters(resolution = 0.8, algorithm = 3) %>%
    RunUMAP(dims = 1:30, reduction = "pca")

# Remove tiny clusters (<10 cells) for stability
obj <- subset(obj, idents = names(table(obj$seurat_clusters))[table(obj$seurat_clusters) >= 10])

# Coarse label mapping by cluster ID
obj <- RenameIdents(obj,
    "0" = "Ductal",
    "1" = "Acinar",
    "2" = "Endocrine",
    "3" = "Endocrine",
    "4" = "Stellate",
    "5" = "Acinar",
    "6" = "Acinar",
    "7" = "Endocrine",
    "8" = "Immune",
    "9" = "Endocrine",
    "10" = "Endocrine",
    "11" = "Stellate"
)

# Stellate subclustering to resolve activated/quiescent and endothelial subtypes
obj_subset_1 <- subset(obj, idents = c("Stellate")) %>%
    SCTransform(clip.range = c(-5, 5), variable.features.n = 2000) %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30, reduction = "pca") %>%
    FindClusters(resolution = 0.5, algorithm = 3) %>%
    RunUMAP(dims = 1:30, reduction = "pca")

# Rename and propagate refined identities back to the merged object
Idents(obj) <- Idents(RenameIdents(obj_subset_1,
    "0" = "aStellate",
    "1" = "aStellate",
    "2" = "qStellate",
    "3" = "aStellate",
    "4" = "Endothelial"
))

# Save annotated merged object for AL189T
SaveSeuratRds(obj, "for_jasmine/obj_merged_AL189T_annotated.rds")
