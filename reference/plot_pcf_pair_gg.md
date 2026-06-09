# Plot a Single Outside-Polygon PCF Pair as a ggplot

Computes the pair correlation function (with simulation envelopes) for
one requested cell-type pair outside the tumor polygons and returns a
ggplot object suitable for embedding in R Markdown documents.

## Usage

``` r
plot_pcf_pair_gg(
  seurat_obj,
  polygons,
  type1,
  type2,
  celltype_col = "singleR.predicted.id.brcaAtlas",
  fov_name = "fov",
  min_cells = 10,
  nsim = 39,
  bw = 10,
  window_obj = NULL
)
```

## Arguments

- seurat_obj:

  Seurat object with spatial coordinates.

- polygons:

  sf polygon object representing tumor regions.

- type1:

  First cell type label (character).

- type2:

  Second cell type label (character).

- celltype_col:

  Metadata column containing cell type labels.

- fov_name:

  Field of view/image name passed to \`Seurat::GetTissueCoordinates\`.

- min_cells:

  Minimum number of cells required per cell type.

- nsim:

  Number of simulations for envelope bands.

- bw:

  Bandwidth for PCF smoothing.

- window_obj:

  Optional \`spatstat.geom::owin\` window.

## Value

A ggplot object with the observed g(r) curve and simulation envelope
bands.
