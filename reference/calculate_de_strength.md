# Calculate Differential Expression Strength

Calculates mean absolute log2 fold change of top DE genes

## Usage

``` r
calculate_de_strength(inside_cells, outside_cells, seurat_obj)
```

## Arguments

- inside_cells:

  Character vector of cell IDs inside region

- outside_cells:

  Character vector of cell IDs outside region

- seurat_obj:

  Seurat object

## Value

Numeric value (mean \|log2FC\|)
