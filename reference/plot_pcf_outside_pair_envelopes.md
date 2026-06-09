# Plot Outside-Polygon PCF Envelopes for All Cell-Type Pairs

Computes and plots envelope-based pair correlation curves for all
requested cell-type pairs outside polygons. Plots can be sent to a
multi-page PDF when \`file\` is provided.

## Usage

``` r
plot_pcf_outside_pair_envelopes(
  seurat_obj,
  polygons,
  celltype_col = "singleR.predicted.id.brcaAtlas",
  fov_name = "fov",
  min_cells = 10,
  include_self_pairs = FALSE,
  ordered_pairs = FALSE,
  nsim = 10,
  bw = 10,
  window_obj = NULL,
  file = NULL,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  Seurat object with spatial coordinates.

- polygons:

  sf polygon object representing tumor regions.

- celltype_col:

  Metadata column containing cell type labels.

- fov_name:

  Field of view/image name passed to \`Seurat::GetTissueCoordinates\`.

- min_cells:

  Minimum number of cells required per cell type.

- include_self_pairs:

  Logical; include same-type pairs.

- ordered_pairs:

  Logical; if \`TRUE\`, computes ordered pairs.

- nsim:

  Number of simulations for envelopes.

- bw:

  Bandwidth passed to \`pcfcross\` in envelope evaluation.

- window_obj:

  Optional \`spatstat.geom::owin\` window.

- file:

  Optional file path for a multi-page PDF output.

- verbose:

  Logical; print progress messages.

## Value

Named list of envelope objects.
