# Add Polygon Membership to Seurat Object

Adds metadata column indicating if cells are inside polygons

## Usage

``` r
add_polygon_membership(seurat_obj, polygons, fov_name = "fov")
```

## Arguments

- seurat_obj:

  Seurat object

- polygons:

  sf polygon object

- fov_name:

  FOV name

## Value

Seurat object with added metadata
