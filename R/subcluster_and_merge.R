#' View Subclusters for a Specified Cluster
#'
#' This function subsets a specific cluster in a Seurat object and performs subclustering.
#' It returns the subclustered Seurat object for viewing the subclustering result.
#'
#' @param seurat_obj A Seurat object.
#' @param cluster_id The ID of the cluster to be subclustered.
#' @param resolution The resolution for subclustering.
#' @param dims A numeric vector specifying the dimensions to use in FindNeighbors. Default is 1:10.
#' @return A Seurat object containing only the cells from the specified cluster, now subclustered.
#' @export
subcluster_view <- function(seurat_obj, cluster_id, resolution = 0.5, dims = 1:10) {

  # Step 1: Extract cells from the specified cluster
  cluster_cells <- WhichCells(seurat_obj, idents = cluster_id)
  cluster_subset <- subset(seurat_obj, cells = cluster_cells)

  # Step 2: Perform subclustering on the subsetted cluster
  cluster_subset <- FindVariableFeatures(cluster_subset)
  cluster_subset <- ScaleData(cluster_subset)
  cluster_subset <- RunPCA(cluster_subset)
  cluster_subset <- FindNeighbors(cluster_subset, dims = dims)  # Use user-specified dims
  cluster_subset <- FindClusters(cluster_subset, resolution = resolution)

  # Return the subclustered Seurat object for viewing
  return(cluster_subset)
}

#' Merge Subclusters into Original Seurat Object with Refined Labels
#'
#' This function takes a subclustered Seurat object, assigns new labels, and integrates the refined labels back into the original Seurat object.
#' The refined labels are stored in a new column named 'refined_clusters'.
#'
#' @param seurat_obj The original Seurat object.
#' @param subcluster_obj The subclustered Seurat object returned from `subcluster_view`.
#' @param new_labels A vector of new labels to assign to the subclusters.
#' @param new_id_column The name of the column to store combined original and refined subcluster IDs. Default is "refined_clusters".
#' @return Updated Seurat object with the refined cluster IDs stored in a new column.
#' @export
merge_subclusters <- function(seurat_obj, subcluster_obj, new_labels, new_id_column = "refined_clusters") {

  # Step 1: Apply user-defined labels to subcluster IDs
  if (length(new_labels) != length(levels(subcluster_obj@active.ident))) {
    stop("Number of new_labels must match the number of subclusters in subcluster_obj.")
  }

  # Assign new labels to subcluster_obj's active.ident
  levels(subcluster_obj@active.ident) <- new_labels

  # Step 2: Extract the refined labels and prepare for merging
  refined_labels <- as.character(subcluster_obj@active.ident)
  names(refined_labels) <- rownames(subcluster_obj@meta.data)

  # Step 3: Ensure the refined_clusters column exists in the original Seurat object
  if (!new_id_column %in% colnames(seurat_obj@meta.data)) {
    seurat_obj@meta.data[[new_id_column]] <- as.character(Idents(seurat_obj))  # Initialize with original identities
  }

  # Step 4: Map the refined labels back to the original Seurat object
  subcluster_cells <- rownames(subcluster_obj@meta.data)  # Cells from subcluster_obj
  seurat_obj@meta.data[subcluster_cells, new_id_column] <- refined_labels

  # Step 5: Ensure refined_clusters is a factor and has the correct levels
  valid_levels <- unique(seurat_obj@meta.data[[new_id_column]])
  seurat_obj@meta.data[[new_id_column]] <- factor(seurat_obj@meta.data[[new_id_column]], levels = valid_levels)

  # Step 6: Set the refined_clusters column as the active identity
  Idents(seurat_obj) <- new_id_column

  # Return the updated Seurat object
  return(seurat_obj)
}




