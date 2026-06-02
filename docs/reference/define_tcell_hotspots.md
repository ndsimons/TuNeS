# Define Intratumoral T Cell Hotspots

Uses DBSCAN on intratumoral CD4+/CD8+ T cells and returns hotspot
polygons, per-cell hotspot assignment, and per-cell distance to nearest
hotspot.

## Usage

``` r
define_tcell_hotspots(
  seurat_obj,
  celltype_col = "singleR.predicted.id.brcaAtlas",
  tcell_labels = c("T Cells CD4+", "T Cells CD8+"),
  polygon_col = "inside_polygon",
  intratumoral_label = TRUE,
  eps = 50,
  min_pts = 10,
  min_hotspot_size = 5,
  fov_name = "fov"
)
```

## Arguments

- seurat_obj:

  Seurat object.

- celltype_col:

  Metadata cell type column.

- tcell_labels:

  Cell type labels considered T cells.

- polygon_col:

  Metadata polygon membership column.

- intratumoral_label:

  Intratumoral value in \`polygon_col\`.

- eps:

  DBSCAN neighborhood radius.

- min_pts:

  DBSCAN minimum points.

- min_hotspot_size:

  Minimum cells to keep a hotspot cluster.

- fov_name:

  Field of view/image name for coordinate extraction.

## Value

List with hotspot polygons, per-cell metadata, and summary fields.
