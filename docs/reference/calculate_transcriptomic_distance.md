# Calculate Transcriptomic Distance

Compares average gene expression profiles between two cell groups

## Usage

``` r
calculate_transcriptomic_distance(inside_cells, outside_cells, seurat_obj)
```

## Arguments

- inside_cells:

  Character vector of cell IDs inside region

- outside_cells:

  Character vector of cell IDs outside region

- seurat_obj:

  Seurat object

## Value

Numeric value (1 - correlation)
