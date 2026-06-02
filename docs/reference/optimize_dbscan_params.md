# Optimize DBSCAN Parameters via Silhouette Score

Performs a grid search over eps and minPts parameter combinations for
DBSCAN clustering of cancer cell coordinates, evaluating each
combination using the average silhouette score. Returns a data frame of
results ranked by silhouette score to help select optimal parameters for
[`add_dbscan_polygons()`](add_dbscan_polygons.md).

## Usage

``` r
optimize_dbscan_params(
  seurat_obj,
  fov_name = "fov",
  cancer_col = "is_cancer",
  eps_range = seq(10, 100, by = 5),
  minPts_range = c(10, 15, 20, 30, 50),
  max_sample = 46000,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  Seurat object with spatial coordinates

- fov_name:

  Name of the field of view (default: "fov")

- cancer_col:

  Metadata column indicating cancer cells (default: "is_cancer")

- eps_range:

  Numeric vector of eps values to test (default: seq(10, 100, by = 5))

- minPts_range:

  Numeric vector of minPts values to test (default: c(10, 15, 20, 30,
  50))

- max_sample:

  Maximum number of non-noise points to sample for silhouette
  calculation to control memory usage (default: 46000)

- verbose:

  Print progress messages (default: TRUE)

## Value

A data frame with columns: eps, minPts, sil_score, n_clusters,
noise_pct, sorted by descending silhouette score. The best combination
is in the first row.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run parameter optimization
results <- optimize_dbscan_params(seurat_obj)

# View best parameters
head(results)

# Use best parameters
best <- results[1, ]
seurat_obj <- add_dbscan_polygons(seurat_obj, eps = best$eps, minPts = best$minPts)
} # }
```
