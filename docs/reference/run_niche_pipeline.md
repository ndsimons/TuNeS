# Run Spatial Niche Pipeline for One Compartment

Computes spatial neighborhoods, clusters niche archetypes, labels
niches, and returns patient-level niche proportions for either
intratumoral or peritumoral cells.

## Usage

``` r
run_niche_pipeline(
  seurat_list,
  subset_value,
  pcr_status,
  all_celltypes,
  celltype_col = "singleR.predicted.id.brcaAtlas",
  polygon_col = "inside_polygon",
  k = 20,
  n_niches = 8,
  fov_name = "fov"
)
```

## Arguments

- seurat_list:

  Named list of Seurat objects (one per patient/sample).

- subset_value:

  Logical value in \`polygon_col\` to select compartment.

- pcr_status:

  Data frame with \`patient_id\` and \`pCR\` columns.

- all_celltypes:

  Character vector of cell type levels to use.

- celltype_col:

  Metadata column containing cell types.

- polygon_col:

  Metadata column indicating inside/outside polygon.

- k:

  Number of nearest neighbors for neighborhood composition.

- n_niches:

  Number of k-means niche archetypes.

- fov_name:

  Field of view/image name for coordinate extraction.

## Value

List with neighborhood-level data, archetypes, labels, and patient
proportions.
