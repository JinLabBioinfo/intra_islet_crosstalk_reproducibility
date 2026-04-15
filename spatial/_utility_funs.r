# Author: Jiachen Sun

### Initialization Begins ###
# Load common settings
set.seed(42)
options(future.globals.maxSize = 2E6 * 1024^2)
num_cores <- 8

# Load common libraries
# This special loading sequence avoids error in VlnPlot.
library(ggplot2)
library(Seurat)
library(patchwork)
library(ggrepel)

library(igraph)

library(parallel)
library(fastmatch)
library(data.table)
library(Matrix)

library(dplyr)
library(tidyr)

# Set variables for paths
data_dir_var <- "for_jasmine"
magick_path_var <- "~/lib/magick-7.1.1/magick"
pigz_path_var <- "~/mambaforge/bin/pigz"

# Define marker genes for cell types
Alpha <- c("GCG", "TTR", "GC")
Beta <- c("INS", "IAPP", "DLK1")
Delta <- c("SST", "RBP4")
Gamma <- c()
aStellate <- c("COL1A1", "COL3A1", "COL6A3")
qStellate <- c("RGS5", "C11orf96", "ADIRF")
Ductal <- c("ANXA4", "SPP1", "KRT19", "KRT7")
Acinar <- c("PLA2G1B", "CPB1", "SPINK1")
Immune <- c("CD74", "HLA-DRB", "APOE")
Endothelial <- c("PLVAP", "FLT1")

# Azimuth website: https://azimuth.hubmapconsortium.org/references/#Human%20-%20Pancreas
azimuth_stellate_sig_list <- c("COL1A1", "COL1A2", "COL6A3", "COL3A1", "TIMP3", "TIMP1", "RGS5", "C11orf96", "FABP4", "ADIRF", "NDUFA4L2", "ESAM")
### Initialization Ends ###

# Read CosMx expression matrix (tx) .csv file in chunks.
create_sparse_matrix_chunked <- function(file_path, chunk_size = 5E5) {
    require(R.utils)

    total_lines <- countLines(file_path)[1]
    header_names <- colnames(fread(file_path, nrows = 0))
    gene_names <- setdiff(header_names, c("fov", "cell_ID"))

    i_list <- list()
    j_list <- list()
    x_list <- list()
    final_cell_names <- list()

    file_row_index <- 1
    cells_processed_count <- 0

    while (file_row_index < total_lines) {
        # Stream a CSV chunk without header then assign colnames from the first row.
        # Using skip + nrows keeps memory bounded, and we avoid loading the whole file.
        # Note: file_row_index starts at 1 to skip the header line in the CSV.
        chunk <- fread(
            file = file_path,
            data.table = FALSE,
            nrows = chunk_size,
            skip = file_row_index,
            header = FALSE
        )

        if (nrow(chunk) == 0) {
            break
        }

        colnames(chunk) <- header_names
        # Cells with cell_ID == 0 are background; remove them early to reduce work.
        chunk <- chunk[chunk$cell_ID != 0, ]

        if (nrow(chunk) == 0) {
            file_row_index <- file_row_index + chunk_size
            next
        }

        # Compose stable cell barcodes for the chunk.
        chunk_cell_ids <- paste0(as.character(chunk$cell_ID), "_", chunk$fov)
        # Restrict to gene columns only; order preserved from header.
        chunk_data <- subset(chunk, select = gene_names)
        # Sparse extraction: arr.ind returns (row=cell index within chunk, col=gene index).
        ind <- which(chunk_data != 0, arr.ind = TRUE)

        if (nrow(ind) > 0) {
            # Values at non-zero positions
            x_list[[length(x_list) + 1]] <- chunk_data[ind]
            # Row indices in the sparse matrix correspond to genes
            i_list[[length(i_list) + 1]] <- ind[, "col"]
            # Column indices correspond to global cell position = local row + offset
            j_list[[length(j_list) + 1]] <- ind[, "row"] + cells_processed_count
            # Append cell names for this chunk (aligned with columns)
            final_cell_names[[length(final_cell_names) + 1]] <- chunk_cell_ids
            # Advance global cell counter by the number of rows in this chunk
            cells_processed_count <- cells_processed_count + length(chunk_cell_ids)
        }

        file_row_index <- file_row_index + chunk_size
        rm(chunk, chunk_data, ind)
    }

    sparseMatrix(
        i = unlist(i_list, use.names = FALSE),
        j = unlist(j_list, use.names = FALSE),
        x = unlist(x_list, use.names = FALSE),
        dims = c(length(gene_names), cells_processed_count),
        dimnames = list(gene_names, unlist(final_cell_names, use.names = FALSE))
    )
}

# Optimized from Seurat::ReadNanostring, which causes memory overflow.
read_nanostring <- function(data_dir = data_dir_var,
                            mtx_file = NULL,
                            metadata_file = NULL,
                            segmentations_file = NULL,
                            genes_filter = NA_character_) {
    files <- c(
        matrix = mtx_file %||% "[_a-zA-Z0-9]*_exprMat_file.csv",
        metadata_file = metadata_file %||% "[_a-zA-Z0-9]*_metadata_file.csv",
        segmentations_file = segmentations_file %||% "[_a-zA-Z0-9]*-polygons.csv"
    )

    files <- vapply(X = files, FUN = function(x) {
        x <- as.character(x = x)
        if (isTRUE(x = dirname(path = x) == ".")) {
            sort(x = list.files(
                path = data_dir,
                pattern = x,
                recursive = FALSE,
                full.names = TRUE
            ), decreasing = TRUE)[1L]
        } else {
            x
        }
    }, FUN.VALUE = character(length = 1L), USE.NAMES = TRUE)

    files[!file.exists(files)] <- NA_character_

    # Read metadata once; used for both centroids and adding to Seurat object.
    md <- fread(
        file = files[["metadata_file"]],
        sep = ",", data.table = FALSE, verbose = FALSE
    )

    # Drop any unresolved parts
    files <- files[!is.na(x = files)]

    outs <- list(matrix = NULL, centroids = NULL, segmentations = NULL)
    for (otype in names(x = outs)) {
        outs[[otype]] <- switch(EXPR = otype,
            matrix = {
                # Stream the expression CSV as a sparse matrix, then filter by gene regex patterns.
                tx <- create_sparse_matrix_chunked(file_path = files[[otype]])
                tx <- tx[!rowSums(sapply(genes_filter, function(filter) grepl(filter, rownames(tx)))) > 0, , drop = FALSE]
                # Drop empty columns (cells with zero total counts)
                as.sparse(x = tx[, which(Matrix::colSums(tx) != 0)])
            },
            centroids = {
                data.frame(
                    x = md$CenterX_global_px,
                    y = md$CenterY_global_px,
                    cell = paste0(
                        as.character(md$cell_ID),
                        "_",
                        md$fov
                    ), stringsAsFactors = FALSE
                )
            },
            segmentations = {
                segs <- fread(
                    file = files[["segmentations_file"]],
                    sep = ",",
                    data.table = FALSE,
                    verbose = FALSE
                )
                data.frame(
                    x = segs$x_global_px, y = segs$y_global_px,
                    cell = paste0(
                        as.character(segs$cellID),
                        "_",
                        segs$fov
                    ), stringsAsFactors = FALSE
                )
            },
        )
    }

    outs[["metadata"]] <- md
    outs
}

# Adapted from Seurat::LoadNanostring, with transcript coordinates removed.
load_nanostring <- function(data_dir = data_dir_var,
                            fov = "image",
                            assay = "RNA",
                            prefix = "RNA") {
    data <- read_nanostring(
        data_dir = data_dir_var,
        type = c(
            "segmentations",
            "centroids"
        ),
        genes_filter = c(
            "Negative",
            "SystemControl"
        )
    )

    obj <- CreateSeuratObject(counts = data$matrix, assay = assay)

    # Create a Seurat image object from centroids/segmentation boundaries.
    coords <- CreateFOV(
        coords = list(
            centroids = CreateCentroids(data$centroids),
            segmentation = CreateSegmentation(data$segmentations)
        ),
        type = c(
            "segmentation",
            "centroids"
        ), assay = assay
    )

    # Keep cells present in both the image and the count matrix to avoid dangling entries.
    bcs <- intersect(Cells(obj), intersect(
        Cells(x = coords, boundary = "segmentation"),
        Cells(x = coords, boundary = "centroids")
    ))
    obj <- subset(x = obj, cells = bcs) # remove cells without corresponding polygon
    obj[[fov]] <- subset(x = coords, cells = bcs)

    rownames(data$metadata) <- paste0(as.character(data$metadata$cell_ID), "_", data$metadata$fov)
    obj <- AddMetaData(object = obj, metadata = data$metadata)
    obj@meta.data$fov <- as.factor(obj@meta.data$fov)
    obj@meta.data$orig.ident <- prefix

    # Prefix cell barcodes to make sample identity explicit across merges.
    RenameCells(obj, new.names = paste0(prefix, "_", colnames(obj)))
}

# Construct a adjacency graph with delaunay triangulation and k-nearest neighbors, then filter edges with a threshold.
# knn_k: neighborhood size for initial candidates; cutoff: quantile used for the distance threshold [0,1].
adj_graph_construction <- function(obj,
                                   threshold_path,
                                   add_cell_type = FALSE,
                                   knn_k = 10,
                                   cutoff = 0.9) {
    require(FNN)
    require(geometry)

    cell_coords <- as.data.frame(obj@images[[1]]$centroids@coords)
    barcodes <- colnames(obj)

    # Delaunay triangulation captures a robust proximity structure without a fixed radius.
    triangles <- delaunayn(cell_coords)
    # Assumes >= 3 points; for degenerate cases (n<3), triangulation is not defined.
    # Convert triangles to unique undirected edge list (i,j) pairs across all faces.
    tng_edges <- unlist(mclapply(seq_len(nrow(triangles)), function(i) {
        list(
            sort(c(triangles[i, 1], triangles[i, 2])),
            sort(c(triangles[i, 1], triangles[i, 3])),
            sort(c(triangles[i, 2], triangles[i, 3]))
        )
    }, mc.cores = num_cores), recursive = FALSE) %>%
        do.call(rbind, .) %>%
        unique() %>%
        as.data.table()
    colnames(tng_edges) <- c("i", "j")

    # kNN graph provides local neighborhood candidates for distance thresholding.
    knn_index <- get.knn(cell_coords, k = knn_k)$nn.index

    # If knn.k exceeds available neighbors, FNN returns up to n-1 neighbors.
    if (file.exists(threshold_path)) {
        threshold <- readRDS(threshold_path)
    } else {
        # Randomize iteration order for slight robustness; same set of cells covered.
        sample_size <- nrow(cell_coords)
        index_rand <- sample(nrow(cell_coords), sample_size)
        # cutoff is a quantile in [0,1]; higher => looser threshold (more edges).
        threshold <- quantile(as.numeric(unlist(mclapply(1:sample_size, FUN = function(i) {
            unlist(apply(cell_coords[knn_index[index_rand[i], ], ], MARGIN = 1, FUN = function(x) {
                sqrt((cell_coords[index_rand[i], ][1] - x[1])^2 +
                    (cell_coords[index_rand[i], ][2] - x[2])^2)
            }))
        }, mc.cores = num_cores))), cutoff)
        saveRDS(threshold, threshold_path)
    }

    # For each cell i, keep only neighbor indices within the threshold distance.
    # We replicate the i-th coordinate to match the number of neighbors and compute
    # Euclidean distances in a vectorized way.
    dist_results <- mclapply(seq_len(nrow(cell_coords)), FUN = function(i) {
        list(indices = knn_index[i, ][sqrt(
            rowSums(
                (cell_coords[knn_index[i, ], ] - (cell_coords[i, ] %>%
                    slice(rep(seq_len(n()), each = nrow(cell_coords[knn_index[i, ], ])))))^2
            )
        ) < threshold])
    }, mc.cores = num_cores)

    # Build a symmetric adjacency from the filtered neighbor sets.
    knn_edges <- as.data.table(
        as_edgelist(
            graph_from_adjacency_matrix(
                sparseMatrix(
                    i = rep(
                        seq_len(nrow(cell_coords)),
                        unlist(mclapply(dist_results,
                            function(x) length(x$indices),
                            mc.cores = num_cores
                        ))
                    ),
                    j = unlist(mclapply(dist_results,
                        function(x) x$indices,
                        mc.cores = num_cores
                    )),
                    x = 1,
                    dims = c(
                        nrow(cell_coords),
                        nrow(cell_coords)
                    )
                ),
                mode = "max"
            )
        )
    )
    colnames(knn_edges) <- c("i", "j")

    # Intersection of Delaunay and distance-filtered kNN edges for robust locality.
    merged_edges <- merge(tng_edges, knn_edges, by = c("i", "j"), all = FALSE)
    g <- graph_from_adjacency_matrix(
        sparseMatrix(
            i = merged_edges$i,
            j = merged_edges$j,
            x = 1,
            dims = c(
                nrow(cell_coords),
                nrow(cell_coords)
            )
        ),
        mode = "max" # ensure undirected adjacency (i<->j) from the symmetric matrix
    )

    V(g)$barcode <- barcodes
    V(g)$coord_x <- cell_coords[, 1]
    V(g)$coord_y <- cell_coords[, 2]

    if (add_cell_type) {
        V(g)$cell_type <- obj$cell_type
    } else {
        V(g)$cell_type <- "Unclassified"
    }

    g
}

# Iteratively expand to neighbors of existing vertices on a graph list.
expand_neighbors <- function(input_list, graph_raw, iter) {
    for (i in 1:iter) {
        # For each subgraph, gather 1-coated neighbors from the full graph by barcode.
        barcode_list <- mclapply(input_list, FUN = function(islet) {
            unlist(mclapply(V(islet)$barcode, FUN = function(barcode) {
                V(graph_raw)$barcode[neighbors(graph_raw, which(V(graph_raw)$barcode == barcode))]
            }, mc.cores = ceiling(num_cores / 4)))
        }, mc.cores = ceiling(num_cores / 2))

        # Expand by inducing a subgraph on the union of current vertices and neighbors.
        output_list <- mclapply(seq_along(input_list), FUN = function(i) {
            induced_subgraph(
                graph_raw,
                which(V(graph_raw)$barcode %fin%
                    c(
                        V(input_list[[i]])$barcode,
                        barcode_list[[i]]
                    ))
            )
        }, mc.cores = ceiling(num_cores / 2))
        input_list <- output_list
    }

    output_list
}

# Get a corase list of "islet candidates" Only used in preprocessing.
get_endo_list <- function(graph_qc) {
    # Restrict to endocrine cells and split into connected components.
    endo_graph <- induced_subgraph(graph_qc, which(V(graph_qc)$cell_type %fin% c("Endocrine")))
    mclapply(1:components(endo_graph)$no, FUN = function(i) {
        subgraph <- induced_subgraph(endo_graph, which(components(endo_graph)$membership == i))
    }, mc.cores = num_cores)
}

# Get a list of islet candidates from a unified adjacency graph. Sort and filter functionalities are included.
get_islet_list <- function(c3_graph, endo_cell_list, cutoff = 30) {
    # Seed with components composed of endocrine (or provided) cell types.
    islet_graph <- induced_subgraph(c3_graph, which(V(c3_graph)$cell_type %fin% endo_cell_list))
    islet_list <- mclapply(1:components(islet_graph)$no, FUN = function(i) {
        subgraph <- induced_subgraph(islet_graph, which(components(islet_graph)$membership == i))
    }, mc.cores = num_cores)

    # Sort by size (descending) and drop small fragments.
    islet_list <- islet_list[order(unlist(mclapply(islet_list, function(islet) {
        length(V(islet))
    }, mc.cores = num_cores)), decreasing = TRUE)]
    islet_list <- islet_list[sapply(islet_list, function(islet) {
        length(V(islet))
    }) > cutoff]

    # Precompute Euclidean edge weights (used by certain layouts/plots/statistics).
    islet_list <- mclapply(islet_list, function(g) {
        E(g)$weight <- sapply(E(g), function(edge) {
            nodes <- ends(g, edge)
            sqrt((V(g)$coord_x[nodes[, 1]] - V(g)$coord_x[nodes[, 2]])^2 +
                (V(g)$coord_y[nodes[, 1]] - V(g)$coord_y[nodes[, 2]])^2)
        })
        g
    }, mc.cores = num_cores)
}

# Get a list of 2-coated islets with neighbor expansion (wrapped).
get_c2_list <- function(islet_list, graph_raw, sample) {
    # Two iterations of neighbor expansion around endocrine seeds.
    c2_list <- expand_neighbors(islet_list, graph_raw, 2)
    c2_list <- c2_list[order(unlist(mclapply(c2_list, function(islet) {
        length(V(islet))
    }, mc.cores = num_cores)), decreasing = TRUE)]

    # Assign stable IDs per islet within this sample.
    c2_list <- mclapply(seq_along(c2_list), function(i) {
        islet <- c2_list[[i]]
        V(islet)$islet.id <- paste0(sample, "_", i)
        islet
    }, mc.cores = num_cores)

    # Edge weights for downstream plotting/metrics.
    c2_list <- mclapply(c2_list, function(g) {
        E(g)$weight <- sapply(E(g), function(edge) {
            nodes <- ends(g, edge)
            sqrt((V(g)$coord_x[nodes[, 1]] - V(g)$coord_x[nodes[, 2]])^2 +
                (V(g)$coord_y[nodes[, 1]] - V(g)$coord_y[nodes[, 2]])^2)
        })
        g
    }, mc.cores = num_cores)
}

# Use concave hull to identify islet boundaries and filter out cells outside the hull polygon.
get_c0_list <- function(c2_list, endo_cell_list) {
    require(concaveman)
    require(sp)

    mclapply(lapply(c2_list, FUN = function(islet) {
        # Compute concave hull using endocrine-only core to define boundary.
        endo_islet <- induced_subgraph(islet, which(V(islet)$cell_type %fin% endo_cell_list))
        endo_islet <- induced_subgraph(endo_islet, which(degree(endo_islet) > 0))
        hull <- concaveman(as.matrix(data.frame(V(endo_islet)$coord_x, V(endo_islet)$coord_y)), concavity = 2)
        hull <- hull[-nrow(hull), ]

        # Tag endocrine vertices that lie exactly on hull vertices (for boundary retention).
        V(endo_islet)$hull.no <- lapply(seq_along(V(endo_islet)), function(i) {
            bool <- paste0(V(endo_islet)$coord_x[i], V(endo_islet)$coord_y[i]) == paste0(hull[, 1], hull[, 2])
            if (any(bool)) {
                which(bool)
            } else {
                0
            }
        }) %>% unlist()

        # Propagate hull vertex indices back to the full islet vertex set (others = 0).
        V(islet)$hull.no <- sapply(seq_along(V(islet)), function(j) {
            if (V(islet)$barcode[j] %fin% V(endo_islet)$barcode) {
                V(endo_islet)$hull.no[which(V(endo_islet)$barcode == V(islet)$barcode[j])]
            } else {
                0
            }
        })

        # Ordered hull coordinates for polygon construction (close ring by repeating first point).
        hull_coords <- data.frame(
            x = V(islet)$coord_x[V(islet)$hull.no > 0],
            y = V(islet)$coord_y[V(islet)$hull.no > 0]
        )[order(V(islet)$hull.no[V(islet)$hull.no > 0]), ]
        hull_coords <- rbind(hull_coords, hull_coords[1, ])

        # Inside test: keep cells inside the concave polygon or on the hull itself.
        V(islet)$inside <- sapply(seq_along(V(islet)), function(i) {
            ifelse(!is.na(over(
                SpatialPoints(matrix(c(V(islet)$coord_x[i], V(islet)$coord_y[i]), ncol = 2)),
                SpatialPolygons(list(Polygons(list(Polygon(hull_coords)), "poly1")))
            )) || V(islet)$hull.no[i] > 0, TRUE, FALSE)
        })

        # Retain only inside vertices and keep the largest connected component (robustness).
        islet <- induced_subgraph(islet, which(V(islet)$inside))
        induced_subgraph(islet, which(components(islet)$membership == which.max(components(islet)$csize)))
    }), function(g) {
        # Refresh edge weights post-trimming.
        E(g)$weight <- sapply(E(g), function(edge) {
            nodes <- ends(g, edge)
            dx <- V(g)$coord_x[nodes[, 1]] - V(g)$coord_x[nodes[, 2]]
            dy <- V(g)$coord_y[nodes[, 1]] - V(g)$coord_y[nodes[, 2]]
            sqrt(dx^2 + dy^2)
        })
        g
    }, mc.cores = num_cores)
}

# Return a list of islet clusters.
get_islet_cluster_list <- function(obj_raw,
                                   obj_menv,
                                   graph_raw_path,
                                   graph_list_c2_path,
                                   threshold_path,
                                   sample_name,
                                   cell_type_levels,
                                   core_cell_type = "Endocrine",
                                   cutoff = 30) {
    # Pipeline overview for obj_ic (islet-cluster object):
    # 1) Build a spatial graph on obj_menv, get endocrine seeds and expand.
    # 2) Transfer cell_type labels from obj_menv to each expanded graph's vertices.
    # 3) Extract all unique barcodes from these graphs.
    # 4) Subset obj_raw to these barcodes to form obj_ic (use all cells for clustering).
    # Note: filtering by endocrine abundance is deferred until after component isolation.
    obj_ic <- subset(
        obj_raw,
        cells = colnames(
            subset(
                obj_raw,
                cells = unique(
                    unlist(
                        extract_barcodes(
                            lapply(
                                get_c2_list(
                                    get_endo_list(
                                        adj_graph_construction(
                                            obj_menv,
                                            threshold_path,
                                            TRUE
                                        )
                                    ),
                                    readRDS(graph_raw_path),
                                    sample_name
                                ), function(islet) {
                                    # Transfer cell_type labels from obj_menv to graph vertices.
                                    V(islet)$cell_type <- obj_menv$cell_type[V(islet)$barcode]
                                    V(islet)$cell_type <- factor(
                                        V(islet)$cell_type,
                                        levels = cell_type_levels
                                    )
                                    V(islet)$cell_type[is.na(V(islet)$cell_type)] <- "Unclassified"
                                    islet
                                }
                            )
                        )
                    )
                )
            )
        )
    )

    obj_ic$cell_type <- Idents(obj_ic)

    # The key difference is that we are using all cells to perform component isolation.
    ic_list <- get_islet_list(
        adj_graph_construction(
            obj_ic,
            threshold_path,
            TRUE
        ),
        cell_type_levels
    )

    # We do filtering here instead, using all endocrine cells in that islet cluster, annotate QC-dropped cells,
    # and remove Unclassified cells that are not part of the microenvironment.
    ic_list <- mclapply(
        annotate_qc_dropped(
            obj_raw,
            ic_list[
                unlist(
                    mclapply(
                        ic_list,
                        function(x) {
                            sum(V(x)$cell_type == core_cell_type)
                        },
                        mc.cores = num_cores
                    )
                ) > cutoff
            ],
            cell_type_levels
        ),
        function(islet) {
            induced_subgraph(
                islet,
                which(V(islet)$cell_type != "Unclassified")
            )
        },
        mc.cores = num_cores
    )

    # Finally, we sort by size again and assign stable IDs.
    ic_list <- ic_list[order(unlist(lapply(ic_list, function(islet) {
        length(V(islet))
    })), decreasing = TRUE)]

    lapply(seq_along(ic_list), function(i) {
        V(ic_list[[i]])$islet.id <- paste0(sample_name, "_", i)
        ic_list[[i]]
    })
}

# Wrapper function to read and prepare G_raw for postprocessing.
get_graph_raw_postprocess <- function(c3_graph, graph_raw_path) {
    c3_labels <- as.character(V(c3_graph)$cell_type)
    names(c3_labels) <- V(c3_graph)$barcode

    graph_raw <- readRDS(graph_raw_path)
    V(graph_raw)$cell_type <- "Unclassified"
    full_labels <- as.character(V(graph_raw)$cell_type)
    names(full_labels) <- V(graph_raw)$barcode

    # Overwrite graph_raw labels wherever c3_graph provides a label for that barcode.
    full_labels[names(full_labels) %fin% names(c3_labels)] <- c3_labels[names(c3_labels) %fin% names(full_labels)]
    V(graph_raw)$cell_type <- full_labels

    graph_raw
}

# A wrapper function to post-process islet data to generate c0, c1, c2 lists.
post_process_islet <- function(sample, data_dir = data_dir_var, endo_cell_list = c("Alpha", "Beta", "Delta", "Gamma")) {
    # Regenerate 3-coated graph on Endocrine subset only
    c3_graph <- adj_graph_construction(
        subset(LoadSeuratRds(paste0(data_dir, "/obj_merged_", sample, "_annotated.rds")), subset = region == "Endocrine"),
        paste0(data_dir, "/threshold_", sample, ".rds"), TRUE
    )
    graph_raw <- get_graph_raw_postprocess(c3_graph, paste0(data_dir, "/graph_raw_", sample, ".rds"))
    c2_list <- get_c2_list(get_islet_list(c3_graph, endo_cell_list), graph_raw, sample) # 2-coated expansion around endocrine cores
    c0_list <- get_c0_list(c2_list, endo_cell_list) # concave-hull trimming to core islets

    saveRDS(c0_list, paste0(data_dir, "/c0_list_", sample, ".rds"))
    saveRDS(c2_list, paste0(data_dir, "/c2_list_", sample, ".rds"))
    saveRDS(expand_neighbors(c0_list, graph_raw, 1), paste0(data_dir, "/c1_list_", sample, ".rds")) # 1-coated coat around c0
}

# This function handles data loading, quality control, sketching, clustering, and integration projection.
pre_process_islet_part_1 <- function(sample,
                                     data_dir = data_dir_var,
                                     qc_nFeature = 50,
                                     qc_nCount = 100,
                                     qc_negprobes = 5,
                                     sketch_ncells = 50000L,
                                     n_features = 1000,
                                     pca_dims = 1:30,
                                     n_trees = 500,
                                     cluster_res = 0.5,
                                     cluster_algo = 3,
                                     cluster_n_start = 50,
                                     cluster_n_iter = 50,
                                     umap_epochs = 500) {
    # Construct file paths using data_dir
    obj <- LoadSeuratRds(file.path(data_dir, paste0("obj_raw_", sample, ".rds")))

    # Quick FOV map
    make_annotated_fov_plot(obj, file.path(data_dir, paste0("fov_map_raw_", sample, ".png")), group_by = "fov")

    # G_raw construction (raw spatial adjacency graph; caches threshold)
    saveRDS(
        adj_graph_construction(
            obj,
            file.path(data_dir, paste0("threshold_", sample, ".rds"))
        ),
        file.path(data_dir, paste0("graph_raw_", sample, ".rds"))
    )

    # Quality control, sketching and primary cell typing (sketch to speed up clustering)
    obj <- subset(obj, subset = nFeature_RNA > qc_nFeature &
        nCount_RNA > qc_nCount &
        nCount_negprobes < qc_negprobes &
        qcCellsFlagged == FALSE &
        qcFlagsFOV == "Pass") %>%
        NormalizeData() %>%
        FindVariableFeatures(nfeatures = n_features)

    obj <- SketchData(obj, ncells = sketch_ncells, features = VariableFeatures(obj))
    DefaultAssay(obj) <- "sketch"

    obj <- obj %>%
        FindVariableFeatures(nfeatures = n_features) %>%
        ScaleData(vars.to.regress = "nCount_RNA") %>% # regress nCount_RNA for cleaner PCA
        RunPCA() %>%
        FindNeighbors(reduction = "pca", dims = pca_dims, n.trees = n_trees) %>%
        FindClusters(resolution = cluster_res, algorithm = cluster_algo, n.start = cluster_n_start, n.iter = cluster_n_iter) %>%
        RunUMAP(reduction = "pca", dims = pca_dims, n.epochs = umap_epochs)

    # Build full-data PCA
    ProjectIntegration(obj, assay = "RNA", sketched.assay = "sketch", reduction = "pca", reduction.name = "pca.full")
}

# This function identifies the endocrine cluster, constructs the c3 microenvironment, and saves the intermediate objects.
pre_process_islet_part_2 <- function(obj,
                                     sample,
                                     data_dir = data_dir_var,
                                     cluster_index = -1,
                                     cluster_min_size = 10,
                                     pca_dims = 1:30,
                                     coat_layers = 3) {
    # Assign "Endocrine" to a cluster based on cluster_index.
    cluster_counts <- table(Idents(obj))
    if (cluster_index >= 0) {
        selected_cluster <- as.character(cluster_index)
    } else {
        eligible <- cluster_counts[cluster_counts > cluster_min_size]
        selected_cluster <- names(if (cluster_index > 0) sort(eligible, decreasing = TRUE) else sort(eligible, decreasing = FALSE))[min(abs(cluster_index), length(eligible))]
    }

    new_ident <- setNames(rep("Nonendocrine", length(levels(obj))), levels(obj))
    new_ident[selected_cluster] <- "Endocrine"

    # Assign coarse labels
    obj <- RenameIdents(obj, new_ident)
    obj$cell_type <- Idents(obj)

    # Project coarse labels to full data
    obj <- ProjectData(obj, sketched.reduction = "pca", full.reduction = "pca.full", dims = pca_dims, refdata = list(cell_type = "cell_type"))

    # G_qc construction (graph after QC) and endocrine seeds
    graph_qc <- adj_graph_construction(
        obj,
        file.path(data_dir, paste0("threshold_", sample, ".rds")),
        TRUE
    )
    endo_list <- get_endo_list(graph_qc)

    # C3 dataset subset (3-coated microenvironment around endocrine seeds)
    DefaultAssay(obj) <- "RNA"
    obj_menv <- subset(obj, cells = unique(unlist(lapply(expand_neighbors(endo_list, graph_qc, coat_layers), function(islet) {
        V(islet)$barcode
    }))))

    # Save C3 dataset and QC object
    SaveSeuratRds(obj_menv, file.path(data_dir, paste0("obj_menv_", sample, ".rds")))
    SaveSeuratRds(obj, file.path(data_dir, paste0("obj_qc_", sample, ".rds")))

    # Mark microenvironment cells to guide manual exocrine region selection
    obj$menv_region <- colnames(obj) %fin% colnames(obj_menv)
    make_annotated_fov_plot(obj, file.path(data_dir, paste0("fov_map_menv_", sample, ".png")))

    obj
}

# This function handles the manual FOV selection for the exocrine region and saves the final subset.
pre_process_islet_part_3 <- function(obj,
                                     sample,
                                     exo_fov,
                                     data_dir = data_dir_var) {
    # Select exocrine region (manual FOV lists by sample)
    obj$target_fov <- obj$fov %in% exo_fov

    # Save exocrine region subset
    SaveSeuratRds(subset(obj, subset = target_fov == TRUE), file.path(data_dir, paste0("obj_exo_", sample, ".rds")))
}

# Merge Endocrine and Exocrine region Seurat objects with downsampling functionality.
merge_object_pair <- function(obj_menv,
                              obj_exo,
                              downsample_ratio = 0.5) {
    # Label regions and clean up assays
    obj_menv$region <- "Endocrine"
    obj_exo$region <- "Exocrine"

    obj_menv[["sketch"]] <- NULL
    obj_exo[["sketch"]] <- NULL

    # Downsample exocrine cells (default 50%)
    obj_exo <- subset(obj_exo, cells = sample(colnames(obj_exo), round(ncol(obj_exo) * downsample_ratio)))

    # Remove cells that overlap with the Endocrine object to avoid duplication
    obj_exo <- subset(obj_exo, cells = setdiff(colnames(obj_exo), colnames(obj_menv)))

    # Merge and unify layers
    JoinLayers(merge(obj_menv, obj_exo))
}

# Identify subpopulations that lacks definitive markers in CosMx data using a reference scRNA-seq dataset.
process_subpopulation <- function(obj_ref,
                                  obj_cosmx,
                                  ref_celltype_col = "CellType",
                                  query_subset_ident = "Endocrine",
                                  clusters_to_score = c("Alpha", "Proliferating_Alpha", "Beta", "Delta", "PP", "Epsilon"),
                                  ctm_genes = c("INS", "SST", "GCG", "TTR"),
                                  pca_dims = 1:10,
                                  cluster_res = 0.6,
                                  sct_clip_range = c(-10, 10),
                                  min_pct = 0.05,
                                  logfc_threshold = 0.25,
                                  p_val_adj = 0.01) {
    # Reduce to common features
    features <- intersect(rownames(obj_ref), rownames(obj_cosmx))

    # Subset objects and clean assays
    # Process Reference (Drop-seq)
    obj_ref <- DietSeurat(obj_ref, features = features, misc = FALSE)
    DefaultAssay(obj_ref) <- "RNA"
    obj_ref[["RNA"]] <- CreateAssay5Object(obj_ref[["RNA"]]$counts)

    # Process Query (CosMx)
    obj_cosmx <- DietSeurat(subset(obj_cosmx, idents = query_subset_ident), features = features, misc = FALSE)
    DefaultAssay(obj_cosmx) <- "RNA"
    obj_cosmx[["RNA"]] <- CreateAssay5Object(obj_cosmx[["RNA"]]$counts)

    # Run SCTransform workflow on Reference
    obj_ref <- obj_ref %>%
        SCTransform(verbose = FALSE) %>%
        RunPCA(verbose = FALSE) %>%
        RunUMAP(dims = pca_dims, verbose = FALSE)

    # Find signature genes from Reference
    # Set identity to the specified cell type column
    Idents(obj_ref) <- ref_celltype_col

    top_features <- FindAllMarkers(obj_ref, min.pct = min_pct, logfc.threshold = logfc_threshold, only.pos = TRUE, assay = "SCT", verbose = FALSE)
    top_features <- top_features[top_features$p_val_adj < p_val_adj, ]

    # Add module scores and re-cluster
    add_module_scores(
        obj = DietSeurat(obj_cosmx, features = top_features$gene, misc = FALSE) %>% # Filter CosMx to only the marker genes found
            SCTransform(variable.features.n = nrow(.), clip.range = sct_clip_range, verbose = FALSE) %>% # Run SCTransform workflow on Query using specific markers
            RunPCA(features = top_features$gene, verbose = FALSE) %>%
            RunUMAP(dims = pca_dims, verbose = FALSE),
        clusters = clusters_to_score,
        top_features = top_features,
        ctm_genes = ctm_genes
    ) %>%
        FindNeighbors(dims = pca_dims, verbose = FALSE) %>%
        FindClusters(resolution = cluster_res, algorithm = 3, verbose = FALSE)
}

# A quick SCTransform pipeline for certain identities for visualization purposes.
quick_reanalyze_subset <- function(obj,
                                   idents_to_keep = NULL,
                                   sct_clip_range = c(-5, 5),
                                   umap_dims = 1:40) {
    # Subset and run the standard SCTransform pipeline
    subset(obj, idents = idents_to_keep) %>%
        SCTransform(clip.range = sct_clip_range, verbose = FALSE) %>%
        RunPCA(verbose = FALSE) %>%
        RunUMAP(dims = umap_dims, verbose = FALSE)
}

# Return geo-adjacency scores for cells in a list of islets.
return_geoadj <- function(input_list) {
    do.call(rbind, lapply(input_list, function(x) {
        # Center coordinates within the islet then compute radial distance.
        coords <- as.matrix(data.frame(V(x)$coord_x, V(x)$coord_y))
        radius <- sqrt((coords[, 1] - mean(coords[, 1]))^2 + (coords[, 2] - mean(coords[, 2]))^2)

        data.frame(
            geo_adj_score = radius / mean(radius),
            cell_type = V(x)$cell_type
        )
    }))
}

# Return cell statistics for a list of islets.
return_cell_stat <- function(input_list, obj) {
    # Compute counts per cell_type for each islet, aligning to the full set of types in obj.
    stat <- as.data.frame(do.call(rbind, mclapply(mclapply(input_list, FUN = function(islet) {
        table(V(islet)$cell_type)
    }, mc.cores = num_cores), FUN = function(islet) {
        cell_types <- c(names(table(obj$cell_type)), "Unclassified")
        t(matrix <- matrix(sapply(cell_types, function(ct) {
            if (ct %fin% names(islet)) {
                islet[ct]
            } else {
                0
            }
        })))
    }, mc.cores = num_cores)))

    colnames(stat) <- c(names(table(obj$cell_type)), "Unclassified")
    stat$total <- rowSums(stat)

    stat
}

# Return shape statistics for a list of islets.
return_shape_stat <- function(input_list,
                              cell_type = "All",
                              concavity = 2,
                              log10 = FALSE) {
    require(concaveman)
    require(sf)

    if (cell_type != "All") {
        input_list <- lapply(input_list, function(islet) {
            induced_subgraph(islet, V(islet)$cell_type %in% cell_type)
        })
    }

    data.frame(
        row.names = mclapply(input_list, function(islet) {
            V(islet)$islet.id[1]
        }, mc.cores = num_cores) %>%
            unlist(),
        cell_number = mclapply(input_list, function(islet) {
            length(V(islet))
        }, mc.cores = num_cores) %>%
            unlist() %>%
            {
                if (log10) log10(.) else .
            },
        # Regularity: min/max radius from centroid of concave hull vertices.
        concave_regularity = lapply(input_list, function(islet) {
            coords <- as.matrix(data.frame(V(islet)$coord_x, V(islet)$coord_y))
            hull_points <- concaveman(coords, concavity = concavity)
            hull_points <- hull_points[-nrow(hull_points), ]
            delta_x <- hull_points[, 1] - mean(hull_points[, 1])
            delta_y <- hull_points[, 2] - mean(hull_points[, 2])
            radius <- sqrt(delta_x^2 + delta_y^2)
            min(radius) / max(radius)
        }) %>% unlist(),
        # Circularity: Polsby–Popper score on concave hull (4πA / P^2).
        concave_circularity = lapply(input_list, function(islet) {
            coords <- as.matrix(data.frame(V(islet)$coord_x, V(islet)$coord_y))
            hull_points <- concaveman(coords, concavity = concavity)
            concave_poly <- st_sfc(st_polygon(list(hull_points)))
            4 * pi * st_area(concave_poly) / st_perimeter(concave_poly)^2
        }) %>% unlist(),
        # Brittleness: compares concave vs convex boundary complexity and area.
        brittleness = lapply(input_list, function(islet) {
            coords <- as.matrix(data.frame(V(islet)$coord_x, V(islet)$coord_y))
            concave_poly <- st_sfc(st_polygon(list(concaveman(coords, concavity = concavity))))
            convex_poly <- st_convex_hull(st_sfc(st_multipoint(coords)))
            sqrt((st_perimeter(convex_poly) * st_area(concave_poly)) / (st_perimeter(concave_poly) * st_area(convex_poly)))
        }) %>% unlist()
    )
}

# Return cell abundance for an object.
return_cell_abd <- function(obj, cell_type_mappings = NULL) {
    if (!is.null(cell_type_mappings)) {
        obj <- RenameIdents(obj, cell_type_mappings)
        obj$cell_type <- Idents(obj)
    }

    # Relative abundance per geoloc (column-wise normalization after counting).
    cell_type_abd <- as.data.frame(table(obj$cell_type, obj$geoloc)) %>%
        group_by(Var2) %>%
        mutate(abd = Freq / sum(Freq)) %>%
        select(-Freq) %>%
        spread(Var2, abd)

    colnames(cell_type_abd)[1] <- "cell_type"
    cell_type_abd$cell_type <- as.character(cell_type_abd$cell_type)
    cell_type_abd[order(cell_type_abd$cell_type), ]
}

# Wrapper function to plot individual islets, with trimming and resizing.
make_individual_islet_plot <- function(obj,
                                       coords,
                                       palette,
                                       file_name,
                                       magick_path = magick_path_var,
                                       max_size = NULL,
                                       padding = 25,
                                       canvas_size = 8000,
                                       auto_border = TRUE,
                                       all_grey = FALSE,
                                       dark_background = FALSE,
                                       border_color = "#FFFFFF") {
    cols <- if (all_grey) {
        rep("#BEBEBE", length(table(Idents(obj))))
    } else {
        palette[names(palette) %fin% names(table(Idents(obj)))]
    }

    # Use the larger of width/height to define an islet "figure size" in pixels.
    figure_size <- max(max(coords[, 1]) - min(coords[, 1]), max(coords[, 2]) - min(coords[, 2]))

    ggsave(file_name,
        suppressWarnings(ImageDimPlot(
            obj,
            boundaries = "segmentation",
            # Scale boundary thickness inversely with size so small islets still have visible borders.
            border.size = ifelse(auto_border, 1.5 / figure_size * max_size, 1.5),
            dark.background = dark_background,
            border.color = border_color,
            col = cols
        ) +
            NoLegend() +
            theme(
                plot.margin = unit(c(0, 0, 0, 0), "cm"),
                panel.spacing = unit(0, "cm"),
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.ticks.length = unit(0, "pt")
            ) +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0))),
        width = canvas_size + padding * 2,
        height = canvas_size + padding * 2,
        units = "px",
        limitsize = FALSE,
        scale = 1,
        device = "png"
    )

    # Post-process using ImageMagick: crop padding and then resize to ensure consistent scales.
    system(paste0(
        magick_path, " mogrify -crop ",
        canvas_size, "x", canvas_size, "+", padding, "+", padding, " ",
        file_name
    ))
    system(paste0(
        magick_path, " mogrify -resize ",
        figure_size, "x", figure_size, " ",
        file_name
    ))
}

# Plot the annotated FOV. Used to select exocrine regions.
make_annotated_fov_plot <- function(obj, save_path, group_by = "menv_region") {
    coords <- as.data.frame(obj@images[[1]]$centroids@coords)
    fov_coords <- data.table(
        fov = obj$fov,
        coord_x = coords$x,
        coord_y = coords$y
    ) %>%
        group_by(fov) %>%
        summarize(
            x_mean = mean(coord_x),
            y_mean = mean(coord_y)
        ) %>%
        as.data.frame()

    # Note swapped axes to match ImageDimPlot coordinate orientation.
    ImageDimPlot(obj, group.by = group_by) +
        annotate(
            geom = "text",
            x = fov_coords$y_mean,
            y = fov_coords$x_mean,
            label = fov_coords$fov,
            color = "#FFFFFF"
        ) +
        NoLegend()

    ggsave(save_path, width = 8000, height = 6400, units = "px", limitsize = FALSE)
}

# Subset and normalize data for stellate cell signature analysis
process_stellate <- function(obj_menv, obj_exo, barcodes) {
    ScaleData(
        NormalizeData(
            subset(
                JoinLayers(
                    merge(
                        DietSeurat(
                            subset(obj_menv, cells = barcodes),
                            misc = FALSE
                        ),
                        DietSeurat(obj_exo, misc = FALSE)
                    )
                ),
                subset = cell_type %in% c("aStellate", "qStellate")
            )
        )
    )
}

# Extract barcodes from a list of islets.
extract_barcodes <- function(c_list) {
    unique(unlist(mclapply(c_list, function(islet) {
        V(islet)$barcode
    }, mc.cores = num_cores)))
}

# Rename cell_types in all islets.
rename_cell_type_graph <- function(c_list, from, to) {
    lapply(c_list, function(islet) {
        V(islet)$cell_type[V(islet)$cell_type == from] <- to
        islet
    })
}

# Annotate "QC-Dropped" cells in the raw object and graph lists based on graph comparison.
# If graph_list is NULL, the function returns the annotated Seurat object, otherwise it returns the annotated graph list.
annotate_qc_dropped <- function(obj,
                                graph_list = NULL,
                                cell_type_levels,
                                qc_thresholds = "nFeature_RNA > 50 & nCount_RNA > 100 & nCount_negprobes < 5 & qcCellsFlagged == FALSE & qcFlagsFOV == \"Pass\"") {
    barcodes_qc_dropped <- colnames(obj)[!with(obj@meta.data, eval(parse(text = qc_thresholds)))]
    obj$cell_type <- as.character(Idents(obj))
    obj$cell_type[colnames(obj) %fin% barcodes_qc_dropped & obj$cell_type == "Unclassified"] <- "QC-Dropped"
    obj$cell_type <- factor(obj$cell_type, levels = c(cell_type_levels, "QC-Dropped"))
    Idents(obj) <- obj$cell_type

    if (is.null(graph_list)) {
        obj
    } else {
        mclapply(graph_list, function(islet) {
            V(islet)$cell_type <- factor(V(islet)$cell_type, levels = c(cell_type_levels, "QC-Dropped"))
            V(islet)$cell_type[V(islet)$barcode %fin% barcodes_qc_dropped & V(islet)$cell_type == "Unclassified"] <- "QC-Dropped"
            islet
        }, mc.cores = num_cores)
    }
}

# Plot individual islets or islet clusters with a simplified wrapper function.
calculate_max_size <- function(index, obj, list_c) {
    # Sort index by the number of unique barcodes in each islet, descending.
    coords <- as.data.frame(subset(obj, cells = unique(V(list_c[[index[order(sapply(index, function(i) {
        length(unique(V(list_c[[i]])$barcode))
    }), decreasing = TRUE)][1]]])$barcode))@images[[1]]$centroids@coords)

    # Set a reference max_size from the largest islet (first after sorting) if not provided.
    max(max(coords[, 1]) - min(coords[, 1]), max(coords[, 2]) - min(coords[, 2]))
}

# Plot individual islets or islet clusters with a simplified wrapper function.
individual_islet_plot_wrapper <- function(index,
                                          obj,
                                          list_c,
                                          palette,
                                          path_prefix,
                                          max_size,
                                          magick_path = magick_path_var) {
    # Subset the object to speed up the loop.
    obj <- subset(obj, cells = extract_barcodes(list_c[index]))

    for (i in seq_along(index)) {
        obj_sub <- subset(obj, cells = unique(V(list_c[[index[i]]])$barcode))
        coords <- as.data.frame(obj_sub@images[[1]]$centroids@coords)

        make_individual_islet_plot(
            obj_sub,
            coords,
            palette,
            paste0(path_prefix, "_", i, ".png"),
            magick_path,
            max_size = max_size
        )
    }
}

# Wrapper function to add module scores for each cluster.
add_module_scores <- function(obj,
                              clusters,
                              top_features,
                              ctm_genes,
                              assay = "SCT",
                              ctrl = 50,
                              prefix = "module_score") {
    for (cluster_name in clusters) {
        # Exclude contaminations (ctm_genes) from cluster-specific signatures.
        gene_list <- setdiff(
            top_features[top_features$cluster == cluster_name, "gene", drop = TRUE],
            ctm_genes
        )

        obj <- AddModuleScore(
            object = obj,
            assay = assay,
            features = list(gene_list),
            ctrl = ctrl,
            name = paste0(prefix, "_", cluster_name)
        )
    }

    obj
}

# Create dual logFC barplot for Drop-seq and CosMx comparison.
dual_logfc_barplot <- function(sig_dropseq,
                               sig_cosmx,
                               sig_list,
                               output_file,
                               dropseq_limits = c(-3, 2),
                               cosmx_limits = c(-1.5, 1),
                               width = 5,
                               height = 3) {
    sig_logfc_df <- data.frame(
        DropSeq = sig_dropseq$avg_logFC,
        CosMx = sig_cosmx$avg_log2FC,
        row.names = sig_list
    ) %>%
        arrange(desc(DropSeq)) %>%
        mutate(
            gene = factor(rownames(.), levels = rownames(.)),
            DropSeq_orig = DropSeq,
            CosMx_orig = CosMx,
            DropSeq_clipped = pmax(dropseq_limits[1], pmin(dropseq_limits[2], DropSeq)),
            CosMx_clipped = pmax(cosmx_limits[1], pmin(cosmx_limits[2], CosMx)),
            DropSeq_label = ifelse(DropSeq < dropseq_limits[1] | DropSeq > dropseq_limits[2], sprintf("%.2f", DropSeq), NA),
            CosMx_label = ifelse(CosMx < cosmx_limits[1] | CosMx > cosmx_limits[2], sprintf("%.2f", CosMx), NA),
            color_DropSeq = ifelse(DropSeq > 0, "#FF0000", "#FFA500"),
            color_CosMx = ifelse(CosMx > 0, "#FF0000", "#FFA500"),
        )

    p1 <- ggplot(sig_logfc_df, aes(x = DropSeq_clipped, y = gene, fill = color_DropSeq)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(aes(label = DropSeq_label, x = DropSeq_clipped), hjust = ifelse(sig_logfc_df$DropSeq_orig > 0, 1.1, -0.1), size = 3, color = "white") +
        scale_fill_identity() +
        theme_classic() +
        scale_x_continuous(limits = dropseq_limits) +
        labs(x = "DropSeq", y = "")
    p2 <- ggplot(sig_logfc_df, aes(x = CosMx_clipped, y = gene, fill = color_CosMx)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(aes(label = CosMx_label, x = CosMx_clipped), hjust = ifelse(sig_logfc_df$CosMx_orig > 0, 1.1, -0.1), size = 3, color = "white") +
        scale_fill_identity() +
        theme_classic() +
        scale_x_continuous(limits = cosmx_limits) +
        labs(x = "CosMx", y = "") +
        theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank()
        )
    ggsave(output_file, p1 + p2, width = width, height = height)
}

# Plot individual islets or islet clusters while highlighting core cells.
# Save two images for one islet (and its surrounding cells), further Photoshop processing is required.
# Utility function not used in the current analysis.
multi_layer_islet_plot_wrapper <- function(index,
                                           obj,
                                           list_c,
                                           palette,
                                           col_non_islet,
                                           path_prefix,
                                           max_size = NULL,
                                           magick_path = magick_path_var,
                                           dark_background = FALSE,
                                           border_color = "#000000") {
    index <- index[order(sapply(index, function(i) {
        length(unique(V(list_c[[i]])$barcode))
    }), decreasing = TRUE)]

    for (i in seq_along(index)) {
        obj_sub <- subset(obj, cells = unique(V(list_c[[index[i]]])$barcode))
        coords <- as.data.frame(obj_sub@images[[1]]$centroids@coords)

        if (i == 1 && is.null(max_size)) max_size <- max(max(coords[, 1]) - min(coords[, 1]), max(coords[, 2]) - min(coords[, 2]))
        make_individual_islet_plot(
            obj_sub,
            coords,
            palette,
            paste0(path_prefix, "_", i, ".png"),
            magick_path,
            max_size = max_size,
            dark_background = dark_background,
            border_color = border_color
        )

        # Second pass: label cells outside c0 hull as "Non-islet" and overlay using a distinct color.
        islet_bcs <- unique(V(c0_list[[index[i]]])$barcode)
        obj_sub$cell_type <- as.character(obj_sub$cell_type)
        obj_sub$cell_type[!colnames(obj_sub) %fin% islet_bcs] <- "Non-islet"
        # Ensure factor levels preserve palette ordering and include Non-islet at the end.
        obj_sub$cell_type <- factor(obj_sub$cell_type, levels = c(names(palette)[names(palette) %in% as.character(obj_sub$cell_type)], "Non-islet"))
        Idents(obj_sub) <- obj_sub$cell_type

        make_individual_islet_plot(
            obj_sub,
            coords,
            c(palette, "Non-islet" = col_non_islet),
            paste0(path_prefix, "_", i, ".png"),
            magick_path,
            max_size = max_size,
            dark_background = dark_background,
            border_color = border_color
        )
    }
}

# Plot a grid of islet graphs to a single PNG using KK layout.
# Utility function not used in the current analysis.
plot_as_graph <- function(islet_list,
                          indices,
                          filename,
                          ncol,
                          nrow,
                          width,
                          height,
                          vertex_size) {
    png(filename, width = width, height = height)
    par(mfrow = c(nrow, ncol), mar = c(1, 1, 1, 1), oma = c(0.5, 0.5, 0.5, 0.5))
    for (i in indices) {
        # KK layout with edge weights; suppress labels for density.
        plot(islet_list[[i]],
            vertex.size = vertex_size, vertex.label = NA,
            vertex.color = V(islet_list[[i]])$color,
            vertex.label.cex = 0, layout = layout_with_kk,
            axes = FALSE, xlab = "", ylab = "",
            epsilon = 1E-10, weights = E(islet_list[[i]])$weight
        )
    }
    dev.off()
}

# Return indices of islets that contain any barcode present in the Seurat object.
# Utility function not used in the current analysis.
detect_islet <- function(list_c, obj) {
    require(purrr)

    which(unlist(lapply(list_c, function(islet) {
        purrr::reduce(
            lapply(V(islet)$barcode, function(barcode) {
                barcode %in% colnames(obj)
            }),
            `|`
        )
    })))
}

# Assign colors to vertices based on cell types.
# Utility function not used in the current analysis.
assign_node_color <- function(graph, palette) {
    V(graph)$color <- "#BEBEBE"
    for (cell_type in names(palette)) {
        V(graph)$color[V(graph)$cell_type == cell_type] <- palette[cell_type]
    }
    graph
}

# Subset and add molecule coordinates from to a Seurat object, while store a molecule-to-barcode mapping for subsetting.
# Utility function not used in the current analysis.
add_molecule_coords <- function(tx_file,
                                object,
                                num_cores = 1,
                                fov = NULL,
                                assay = "RNA",
                                pigz_path = pigz_path_var) {
    if (is.null(fov)) fov <- Images(object)[1]

    temp_dir <- file.path(tempdir(), paste0("subset_ops_", Sys.getpid()))
    if (!dir.exists(temp_dir)) dir.create(temp_dir)
    on.exit({
        unlink(temp_dir, recursive = TRUE)
    })

    writeLines(as.character(object@meta.data[, "cell"]), file.path(temp_dir, "barcodes.txt"))

    # Read first line to find column indices
    col_indices <- match(c("x_global_px", "y_global_px", "target", "cell"), strsplit(system(paste0(pigz_path, " -dc ", tx_file, " | head -n 1"), intern = TRUE), ",")[[1]])

    # Filter on disk using awk to avoid loading full file into RAM.
    # Load barcodes into array 'b', strip quotes from cell column with gsub, print if match found.
    system(paste0(
        pigz_path, " -p ", num_cores, " -dc ", tx_file, " | ",
        "awk -F, 'BEGIN{OFS=\",\"} ",
        "NR==FNR{b[$1]; next} ",
        "FNR==1 {print ", paste(paste0("$", col_indices), collapse = ","), "; next} ",
        "{ val=$", col_indices[4], "; ",
        "gsub(/\"/, \"\", val); ",
        "if (val in b) {print ", paste(paste0("$", col_indices), collapse = ","), "}}' ",
        file.path(temp_dir, "barcodes.txt"), " - | ",
        pigz_path, " -p ", num_cores, " > ", file.path(temp_dir, "filtered.csv.gz")
    ))

    mx <- fread(file.path(temp_dir, "filtered.csv.gz"), nThread = num_cores)
    setnames(mx, c("x", "y", "gene", "cell"))

    mx <- mx[gene %in% rownames(object)]

    meta_lookup <- data.table(
        barcode = rownames(object@meta.data),
        cell = as.character(object@meta.data$cell)
    )
    setkey(meta_lookup, cell)

    mx <- merge(mx, meta_lookup, by = "cell", all.x = TRUE)

    centroids <- object@images[[fov]]$centroids
    object[[fov]] <- CreateFOV(
        coords = object[[fov]]$segmentation,
        molecules = as.data.frame(mx[, .(x, y, gene)]),
        assay = assay
    )
    object[[fov]][["centroids"]] <- centroids

    object@misc$molecules <- split(mx$barcode, mx$gene)

    object
}

# Subset molecule data in a Seurat object based on a molecule-to-barcode mapping, which is not implemented in Seurat.
# Utility function not used in the current analysis.
subset_molecule <- function(object, fov = NULL) {
    if (is.null(fov)) fov <- Images(object)[1]

    # Calculate keep indices for every gene
    processed_data <- lapply(object@misc$molecules, function(barcodes_vec) {
        keep_indices <- which(barcodes_vec %in% Cells(object))

        if (length(keep_indices) == 0) {
            return(NULL)
        }

        list(
            meta = barcodes_vec[keep_indices],
            indices = keep_indices
        )
    })

    final_data <- processed_data[!vapply(processed_data, is.null, logical(1))]

    # Update molecule mappings in object@misc
    object@misc$molecules <- lapply(final_data, `[[`, "meta")

    # Extract coordinates, add gene name, and concat into a single data.table
    # Reconstruct the FOV using the concatenated dataframe
    centroids <- object@images[[fov]]$centroids
    object[[fov]] <- CreateFOV(
        coords = object[[fov]]$segmentation,
        molecules = as.data.frame(rbindlist(mapply(
            function(mol, data, gene_name) {
                coords <- mol@coords[data$indices, 1:2, drop = FALSE]
                dt <- as.data.table(coords)
                setnames(dt, c("x", "y"))
                dt[, gene := gene_name]
                dt
            },
            object[[fov]][["molecules"]][names(final_data)],
            final_data,
            names(final_data),
            SIMPLIFY = FALSE
        ))),
        assay = "RNA"
    )
    object[[fov]][["centroids"]] <- centroids

    object
}
