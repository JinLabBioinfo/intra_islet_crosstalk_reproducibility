# Author: Jiachen Sun

source("/home/jxs2269/islet_cosmx/for_jasmine/src/pipeline_reproducibility_common/_utility_funs.r")

obj <- merge_object_pair(LoadSeuratRds("for_jasmine/obj_menv_AM010.rds"), LoadSeuratRds("for_jasmine/obj_exo_AM010.rds"))

# Global clustering on merged object (SCT -> PCA -> graph -> SLM -> UMAP)
obj <- SCTransform(obj, clip.range = c(-15, 15), variable.features.n = 2000) %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:25, reduction = "pca") %>%
    FindClusters(resolution = 0.8, algorithm = 3) %>%
    RunUMAP(dims = 1:25, reduction = "pca")

# Remove tiny clusters (<10 cells) for stability
obj <- subset(obj, idents = names(table(obj$seurat_clusters))[table(obj$seurat_clusters) >= 10])

# Coarse label mapping by cluster ID
obj <- RenameIdents(obj,
    "0" = "Nonendocrine",
    "1" = "Acinar",
    "2" = "Acinar",
    "3" = "Ductal",
    "4" = "Endocrine",
    "5" = "Endocrine",
    "6" = "Inconclusive",
    "7" = "Endocrine",
    "8" = "Endocrine",
    "9" = "Endocrine",
    "10" = "Endocrine",
    "11" = "Stellate"
)

# Nonendocrine subclustering level 1
obj_subset_1 <- subset(obj, idents = c("Nonendocrine")) %>%
    SCTransform(clip.range = c(-10, 10), variable.features.n = 2000) %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:15, reduction = "pca") %>%
    FindClusters(resolution = 0.4, algorithm = 3) %>%
    RunUMAP(dims = 1:15, reduction = "pca")

# Rename and propagate refined identities back to the merged object
Idents(obj) <- Idents(RenameIdents(obj_subset_1,
    "0" = "Inconclusive",
    "1" = "Acinar",
    "2" = "Inconclusive",
    "3" = "Endocrine",
    "4" = "Immune",
    "5" = "Endocrine"
))

# Stellate subclustering
obj_subset_2 <- subset(obj, idents = c("Stellate")) %>%
    SCTransform(clip.range = c(-5, 5), variable.features.n = 2000) %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:15, reduction = "pca") %>%
    FindClusters(resolution = 0.3, algorithm = 3) %>%
    RunUMAP(dims = 1:15, reduction = "pca")

# Rename and propagate refined identities back to the merged object
Idents(obj) <- Idents(RenameIdents(obj_subset_2,
    "0" = "aStellate",
    "1" = "aStellate",
    "2" = "qStellate"
))

# Save annotated merged object for AM010
SaveSeuratRds(obj, "for_jasmine/obj_merged_AM010_annotated.rds")
