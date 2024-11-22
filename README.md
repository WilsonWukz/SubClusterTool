# SubClusterTool

## Overview
`SubClusterTool` is an R package designed to facilitate subclustering and integration of subclusters back into Seurat objects. It allows users to extract a specific cluster from a Seurat object, perform subclustering with custom resolutions and dimensions, and merge the refined subclusters back into the original Seurat object with user-defined labels. This tool is particularly useful for detailed single-cell RNA sequencing analysis, where deeper insights into specific clusters are needed.

## Installation

First, install the `devtools` package if you haven't already:

```r
install.packages("devtools")
```

Then, install the `SubClusterTool` package from GitHub:

```r
devtools::install_github("WilsonWukz/SubClusterTool")
```

## Usage

To use the `SubClusterTool` package, follow these steps:

### 1. Load the package and required dependencies:

```r
library(Seurat)
library(SubClusterTool)
```

### 2. Subcluster a specific cluster within a Seurat object:
Example: Subcluster cluster 0 with resolution 0.5 and using the first 10 PCA dimensions

```r
cluster_subset <- subcluster_view(seurat_obj = your_seurat_object, cluster_id = 0, resolution = 0.5, dims = 1:10)

# Visualize the subclustered result
DimPlot(cluster_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
```

### 3. Merge subclusters back into the original Seurat object:
Define new labels for the subclusters
```r
new_labels <- c("Subcluster_1", "Subcluster_2", "Subcluster_3")
```
Merge the subclusters back into the original object
```r
your_seurat_object <- merge_subclusters(
  seurat_obj = your_seurat_object,
  subcluster_obj = cluster_subset,
  new_labels = new_labels,
  new_id_column = "refined_clusters"
)
```
Visualize the updated Seurat object with refined clusters
```r
DimPlot(your_seurat_object, reduction = "umap", group.by = "refined_clusters", label = TRUE, pt.size = 0.5)
```

## Example
![image](https://github.com/user-attachments/assets/b4d441b8-addf-4dcb-b3a0-74f301e6d1ac)
Then the code should be
```r
B_subset <- subcluster_view(seurat_obj = Cell.integrated, cluster_id = "B cells", resolution = 0.6, dims = 1:10)
```
Then we will get
![image](https://github.com/user-attachments/assets/9d37f9b8-4ba1-425a-97f1-b4966febd531)
Then we try to merge them back to the Cell.intergrated, and if I want to devide the B_subsets into 2 kinds, **B0** and **B1**
```r
Cell.integrated <- merge_subclusters(
  seurat_obj = Cell.integrated,
  subcluster_obj = B_subset,
  new_labels = c("B0", "B0", "B0", "B0", "B1", "B1","B1"),
  new_id_column = "refined_clusters"
)
```
![image](https://github.com/user-attachments/assets/ace61fcf-3cea-4faa-a084-a6fcc4f7a7f2)


## Functions

### Built-in Functions
- **`subcluster_view(seurat_obj, cluster_id, resolution, dims)`**  
  Extracts a specified cluster from a Seurat object, performs subclustering, and returns the subclustered Seurat object.

- **`merge_subclusters(seurat_obj, subcluster_obj, new_labels, new_id_column)`**  
  Merges the refined subclusters back into the original Seurat object with user-defined labels. The merged labels are stored in a new column.

## Author
**Kezhao Wu**  
Maintainer: [Kezhao Wu](mailto:wilsonkwu@gmail.com)

## License
MIT

## Bug Reports
If you encounter any issues or have suggestions for improvements, please report them on the [GitHub issues page](https://github.com/WilsonWukz/SubClusterTool/issues).

## URL
[GitHub Repository](https://github.com/WilsonWukz/SubClusterTool)
