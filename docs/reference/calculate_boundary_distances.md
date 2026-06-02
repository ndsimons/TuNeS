# Calculate Boundary Distances

Calculates signed distance from each cell to polygon boundary

## Usage

``` r
calculate_boundary_distances(seurat_obj, polygons, fov_name = "fov")
```

## Arguments

- seurat_obj:

  Seurat object with spatial data

- polygons:

  sf polygon object

- fov_name:

  Name of FOV in Seurat object

## Value

Numeric vector of signed distances (negative = inside)
