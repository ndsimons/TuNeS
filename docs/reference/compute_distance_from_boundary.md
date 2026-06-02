# Compute Signed Distance From Tumor Boundary

Computes per-cell signed distance to an intratumoral boundary polygon
built from intratumoral cells using a concave hull.

## Usage

``` r
compute_distance_from_boundary(
  seurat_obj,
  polygon_col = "inside_polygon",
  intratumoral_label = TRUE,
  fov_name = "fov",
  concavity = 2
)
```

## Arguments

- seurat_obj:

  Seurat object.

- polygon_col:

  Metadata column indicating inside/outside polygon.

- intratumoral_label:

  Value in \`polygon_col\` that marks intratumoral cells.

- fov_name:

  Field of view/image name for coordinate extraction.

- concavity:

  Concavity parameter passed to \`concaveman\`.

## Value

Seurat object with \`dist_to_boundary\` and \`abs_dist_to_boundary\`
metadata.
