# Plot T Cell Hotspots on Spatial Coordinates

Plot T Cell Hotspots on Spatial Coordinates

## Usage

``` r
plot_tcell_hotspots(
  seurat_obj,
  hotspot_result,
  celltype_col = "singleR.predicted.id.brcaAtlas",
  tcell_labels = c("T Cells CD4+", "T Cells CD8+"),
  sample_label = "",
  fov_name = "fov"
)
```

## Arguments

- seurat_obj:

  Seurat object.

- hotspot_result:

  Output from \`define_tcell_hotspots\`.

- celltype_col:

  Metadata cell type column.

- tcell_labels:

  T cell labels to highlight.

- sample_label:

  Optional label for title.

- fov_name:

  Field of view/image name for coordinate extraction.

## Value

ggplot object.
