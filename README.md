# stGradient: Tumor Nest Selector

<img align="right" src="https://github.com/user-attachments/assets/ab584c5a-2067-4506-971a-fc3f21d719f8" width="305" height="352">

**Interactive spatial analysis of tumor boundaries in spatial transcriptomics data**

stGradient provides tools to interactively select tumor regions and analyze genomic separation at varying distances from tumor boundaries. This package is designed for spatial transcriptomics analysis, particularly for understanding tumor microenvironment architecture.

Features

- **Interactive Polygon Selection**: Shiny app for drawing regions on spatial data
- **Distance-Based Analysis**: Calculate genomic metrics at varying distances from boundaries
- **Multiple Metrics**: Transcriptomic distance, differential expression, cell type composition
- **Visualization**: Comprehensive plotting functions for results
- **Flexible Modes**: Compare entire inside region vs boundary cells

## Installation

```r
# Install devtools if needed
install.packages("devtools")

# from GitHub
devtools::install_github("ndsimons/stGradient")
```

## Dependencies

stGradient requires the following packages:
- shiny
- plotly
- sf
- dplyr
- ggplot2
- Seurat
- tidyr
- magrittr
- dbscan
- concaveman
- cowplot

Optional helper packages used in some plotting workflows:
- patchwork
- scales
- spatstat.geom
- spatstat.explore

## Quick Start

If using shiny app to select polygons:
### Launch Interactive Polygon Selector

```r
library(stGradient)

# Launch Shiny app
launch_stGradient()

# In the app:
# 1. Enter your Seurat object name
# 2. Enter FOV name
# 3. Click "Load Data"
# 4. Draw polygons using the toolbar
# 5. Enter polygon object name
# 6. Click "Extract Polygons"
```

### Add Polygon Information to Seurat Object

```r
# Add inside/outside membership
seurat_obj <- add_polygon_membership(seurat_obj, drawn_polygons)

# Calculate boundary distances
# if using interactive polygon selection, use drawn_polygons or whatever you name it
seurat_obj$boundary_distance <- calculate_boundary_distances(seurat_obj, drawn_polygons)
```

If using DBscan to select polygons:
### Use DBscan to select and add polygons to Seurat Object

```r
auto_polygons <- plot_dbscan_polygons(seurat_obj, plot = FALSE)
seurat_obj <- add_polygon_membership(seurat_obj, auto_polygons)
seurat_obj$boundary_distance <- calculate_boundary_distances(seurat_obj, auto_polygons)

# Or add polygons directly to the Seurat object boundaries + membership in one step
# seurat_obj <- add_dbscan_polygons(seurat_obj)
# auto_polygons <- seurat_obj@images[["fov"]]@boundaries$cancer_regions
# seurat_obj$boundary_distance <- calculate_boundary_distances(seurat_obj, auto_polygons)

# customize dbscan parameters:
# seurat_obj <- add_dbscan_polygons(
#   seurat_obj,
#   eps = 35,
#   minPts = 20,
#   concavity = 2
# )
```

## Xenium Workflow Helpers

### Load Xenium Seurat Object with Stable Cell IDs

```r
library(stGradient)

seurat_obj <- load_xenium_seurat(
  rds_path = "path/to/sample_seuratObject.rds",
  celltype_col = "singleR.predicted.id.brcaAtlas",
  cancer_labels = c("Cancer LumA SC", "Cancer LumB SC", "Cancer Her2 SC")
)
```

### QC Plotting for ROI and Cell Type Overlays

```r
qc_out <- qc_plot_xenium(
  obj = seurat_obj,
  roi_col = "inside_polygon",
  celltype_col = "singleR.predicted.id.brcaAtlas",
  show = "both",
  zoom = TRUE,
  zoom_pad = 0.05,
  grey_outside_roi = TRUE
)

qc_out$p
```

### ROI Composition Plotting Helpers

```r
# Distribution inside ROI only
plot_roi_celltype_dist(
  df = qc_out$df,
  roi_col = "inside_polygon",
  celltype_col = "singleR.predicted.id.brcaAtlas",
  metric = "proportion",
  labels = "both"
)

# Stacked inside vs outside
plot_roi_inout_celltype_stack(
  df = qc_out$df,
  roi_col = "inside_polygon",
  celltype_col = "singleR.predicted.id.brcaAtlas",
  metric = "proportion"
)

# Dodged inside vs outside (with optional top_n)
plot_roi_inout_celltype_dodge(
  df = qc_out$df,
  roi_col = "inside_polygon",
  celltype_col = "singleR.predicted.id.brcaAtlas",
  metric = "count",
  top_n = 10,
  top_n_by = "inside"
)
```

## Pair Correlation Function (PCF) Outside Tumor Polygons

stGradient includes optimized helpers for `spatstat`-based PCF analysis using only cells
outside tumor polygons.

### Calculate Outside-Polygon PCF Summary at a Target Distance

```r
pcf_summary <- calculate_pcf_outside_summary(
  seurat_obj = seurat_obj,
  polygons = auto_polygons,
  distance = 200,
  min_cells = 10,
  celltype_col = "singleR.predicted.id.brcaAtlas",
  include_self_pairs = TRUE,
  ordered_pairs = FALSE
)

head(pcf_summary)
```

### Visualize Outside-Polygon PCF as a Heatmap

```r
plot_pcf_outside_heatmap(pcf_summary, midpoint = 1)
```

### Plot PCF Envelope Curves for All Cell-Type Pairs

```r
# Creates one page per pair in a PDF file
env_list <- plot_pcf_outside_pair_envelopes(
  seurat_obj = seurat_obj,
  polygons = auto_polygons,
  celltype_col = "singleR.predicted.id.brcaAtlas",
  min_cells = 10,
  include_self_pairs = FALSE,
  ordered_pairs = FALSE,
  nsim = 10,
  bw = 10,
  file = "pcf_outside_all_pairs.pdf"
)
```

Notes:
- `ordered_pairs = FALSE` computes each pair once (`A-B`), while `TRUE` computes both `A-B` and `B-A`.
- These PCF helpers analyze only cells outside polygons by design.
- The underlying PCF functions come from the `spatstat` ecosystem.

## pCR Niche Analysis Workflow (Intratumoral vs Peritumoral)

These functions package the niche workflow from the pCR analysis Rmd into reusable APIs.

```r
# Example pCR status table with required columns
# pcr_status <- data.frame(patient_id = c("P1", "P2"), pCR = c(1L, 0L))

# all_celltypes should be a stable vector of labels across samples
intra <- run_niche_pipeline(
  seurat_list = seurat_list,
  subset_value = TRUE,
  pcr_status = pcr_status,
  all_celltypes = all_celltypes,
  celltype_col = "singleR.predicted.id.brcaAtlas",
  polygon_col = "inside_polygon",
  k = 20,
  n_niches = 8
)

peri <- run_niche_pipeline(
  seurat_list = seurat_list,
  subset_value = FALSE,
  pcr_status = pcr_status,
  all_celltypes = all_celltypes
)

# Per-compartment statistical testing
intra_niches <- setdiff(colnames(intra$niche_prop_df), c("patient_id", "pCR"))
intra_X <- as.matrix(intra$niche_prop_df[, intra_niches, drop = FALSE])
intra_y <- intra$niche_prop_df$pCR

stat_intra <- run_niche_stat_tests(intra_X, intra_y, intra_niches, n_perm = 1000)
lr_intra <- run_niche_univariate_lr(intra_X, intra_y, intra_niches)

# Plot helpers
plot_niche_stacked_composition(intra$niche_prop_df)
plot_niche_grouped_bars(stat_intra, intra_X, intra_y, intra_niches, title = "Intratumoral")
plot_niche_coefficients(lr_intra, title = "Intratumoral")
```

## Spatial Distance-Gene Workflow

These functions package the distance-gene analysis from `spatialDE_v2.Rmd`.

```r
# 1) Compute signed boundary distance
seurat_obj <- compute_distance_from_boundary(
  seurat_obj,
  polygon_col = "inside_polygon",
  intratumoral_label = TRUE,
  fov_name = "fov"
)

# 2) Run gene-distance models for peritumoral region
dist_res <- find_distance_associated_genes(
  seurat_obj,
  assay = "Xenium",
  focus = "peritumoral",
  min_cells_expressing = 10,
  n_distance_bins = 10,
  max_dist_quantile = 0.95,
  p_adj_threshold = 0.05
)

# 3) Plot volcano and trend panels
dist_plots <- plot_distance_genes(dist_res, top_n = 12, sample_label = "Sample_01")
dist_plots$volcano
dist_plots$trends
```

## T-Cell Hotspot Distance-Gene Workflow

These functions package the intratumoral CD4/CD8 hotspot analysis.

```r
# 1) Identify intratumoral T-cell hotspots
hs <- define_tcell_hotspots(
  seurat_obj,
  celltype_col = "singleR.predicted.id.brcaAtlas",
  tcell_labels = c("T Cells CD4+", "T Cells CD8+"),
  polygon_col = "inside_polygon",
  intratumoral_label = TRUE,
  eps = 50,
  min_pts = 10,
  min_hotspot_size = 5
)

# 2) Plot hotspots
plot_tcell_hotspots(seurat_obj, hs, celltype_col = "singleR.predicted.id.brcaAtlas")

# 3) Model gene expression vs distance from hotspots
hs_res <- find_hotspot_distance_genes(
  seurat_obj,
  hotspot_result = hs,
  assay = "Xenium",
  focus = "perihotspot",
  min_cells_expressing = 10,
  n_distance_bins = 10,
  max_dist_quantile = 0.95,
  p_adj_threshold = 0.05,
  exclude_tcell_labels = c("T Cells CD4+", "T Cells CD8+")
)
```

### Run Distance Profile Analysis

```r
# Calculate all metrics at varying distances
results <- calculate_distance_profile(
  seurat_obj = seurat_obj,
  distance_thresholds = seq(10, 200, by = 10),
  metric = "all",
  group_by = "cellType",
  inside_mode = "distance"
)

# Extract results
distance_profile <- results$profile
celltype_contributions <- results$celltype_contributions
```

### Visualize Results

```r
# Plot transcriptomic distance
plot_distance_comparison(distance_profile)

# Plot all metrics
plot_all_metrics(distance_profile)

# Heatmap of cell type contributions
plot_celltype_heatmap(celltype_contributions)

# Cell type proportions over distance
plot_celltype_proportions(celltype_contributions, top_n = 5)
```

### Find Critical Distances

```r
# Identify distances where separation is maximized
critical_dist <- find_critical_distance(distance_profile, "transcriptomic_distance")

cat("Maximum separation at:", critical_dist$max_distance, "μm\n")
cat("Plateau begins at:", critical_dist$plateau_distance, "μm\n")
```

## Complete Workflow Example

```r
library(stGradient)
library(Seurat)

# Skip step 1 and proceed to step 2 if using DBscan for polygons
# Step 1: Select regions interactively
launch_stGradient()
seurat_obj <- add_polygon_membership(seurat_obj, drawn_polygons)

# Step 2: Generate polygons automatically (optional alternative to interactive)
# auto_polygons <- plot_dbscan_polygons(seurat_obj, plot = FALSE)
# seurat_obj <- add_polygon_membership(seurat_obj, auto_polygons)

# Step 3: Calculate signed boundary distances
seurat_obj$boundary_distance <- calculate_boundary_distances(seurat_obj, drawn_polygons)

# Step 4: Analyze both modes
results_all <- calculate_distance_profile(
  seurat_obj = seurat_obj,
  distance_thresholds = seq(10, 200, by = 10),
  metric = "all",
  group_by = "cellType",
  inside_mode = "all"
)

results_distance <- calculate_distance_profile(
  seurat_obj = seurat_obj,
  distance_thresholds = seq(10, 200, by = 10),
  metric = "all",
  group_by = "cellType",
  inside_mode = "distance"
)

# Step 5: Combine and visualize
combined_results <- rbind(
  results_all$profile,
  results_distance$profile
)

plot_distance_comparison(combined_results)
plot_celltype_heatmap(results_all$celltype_contributions)
```

## Main Functions

### Interactive Selection
- `launch_stGradient()` - Launch Shiny app for polygon drawing

### Data Processing
- `load_xenium_seurat()` - Load Xenium Seurat RDS and create stable cell mapping
- `add_dbscan_polygons()` - Generate DBscan polygons and add to seurat object
- `add_polygon_membership()` - Add inside/outside labels to Seurat object
- `calculate_boundary_distances()` - Calculate signed distances to boundary

### Metrics
- `calculate_transcriptomic_distance()` - Gene expression correlation
- `calculate_de_strength()` - Differential expression magnitude
- `calculate_composition_difference()` - Cell type distribution differences

### Analysis
- `calculate_distance_profile()` - Main analysis function
- `find_critical_distance()` - Identify key separation distances
- `calculate_pcf_outside_summary()` - PCF summary for outside-polygon cell-type pairs
- `run_niche_pipeline()` - Build spatial neighborhoods, cluster niches, and compute patient niche proportions
- `run_niche_stat_tests()` - Wilcoxon + permutation testing of niche proportions by pCR group
- `run_niche_univariate_lr()` - Univariate logistic models of pCR by niche
- `compute_distance_from_boundary()` - Signed distance from tumor boundary for each cell
- `find_distance_associated_genes()` - Distance-associated genes using linear and Spearman models
- `define_tcell_hotspots()` - DBSCAN-based intratumoral CD4/CD8 hotspot detection
- `find_hotspot_distance_genes()` - Genes associated with distance to T-cell hotspots

### Visualization
- `plot_distance_comparison()` - Compare inside modes
- `plot_all_metrics()` - Faceted view of all metrics
- `plot_celltype_heatmap()` - Cell type contribution heatmap
- `plot_celltype_proportions()` - Proportion trajectories
- `plot_dbscan_polygons()` - Plot cell coordinates with dbscan polygons
- `qc_plot_xenium()` - Side-by-side ROI and cell type spatial QC plots
- `plot_roi_celltype_dist()` - ROI-only cell type count/proportion bar plot
- `plot_roi_inout_celltype_stack()` - Stacked inside vs outside composition plot
- `plot_roi_inout_celltype_dodge()` - Dodged inside vs outside comparison plot
- `plot_pcf_outside_heatmap()` - Heatmap of outside-polygon PCF values at target distance
- `plot_pcf_outside_pair_envelopes()` - All-pairs outside-polygon PCF envelope plots
- `plot_niche_stacked_composition()` - Per-patient niche composition by pCR group
- `plot_niche_grouped_bars()` - Side-by-side pCR vs non-pCR niche means with error bars
- `plot_niche_coefficients()` - Logistic regression coefficient plot for niche predictors
- `plot_distance_genes()` - Volcano/trend plots for boundary distance-associated genes
- `plot_tcell_hotspots()` - Spatial hotspot overlay for intratumoral CD4/CD8 clusters

## Understanding Inside Modes

stGradient supports two modes for defining "inside" cells:

### Mode: "all" (default)
- **Inside**: All cells inside polygon (entire tumor core)
- **Outside**: Cells outside polygon within distance threshold
- **Use case**: "How does the entire tumor differ from the surrounding tissue?"

### Mode: "distance"
- **Inside**: Only cells inside polygon within distance threshold
- **Outside**: Cells outside polygon within distance threshold
- **Use case**: "How do cells at the tumor edge differ on both sides of the boundary?"

## Biological Interpretation

### Transcriptomic Distance
- Range: 0 (identical) to 2 (opposite)
- Higher values indicate greater gene expression differences

### Composition Difference
- Range: 0 (identical) to ~2 (completely different)
- Sum of squared differences in cell type proportions
- Independent of sample size

### DE Strength
- Mean |log2FC| of top 100 differentially expressed genes
- Higher values indicate stronger transcriptional changes

## Tips and Best Practices

1. **Drawing Polygons**: Use the freehand tool for irregular tumor shapes
2. **Distance Thresholds**: Start with `seq(10, 200, by = 10)` for most tissues
3. **Cell Types**: Ensure your `group_by` column has clean cell type labels
4. **Minimum Cells**: Analysis requires ≥10 cells per group at each distance
5. **Multiple Polygons**: You can draw multiple regions and combine them

## Troubleshooting

**"No shapes found in expected format"**
- Make sure you drew shapes before clicking "Extract Polygons"
- Try drawing a new shape

**"Insufficient cells" warnings**
- Increase distance thresholds
- Draw larger polygons
- Check that cells have the required metadata

**Plotting errors**
- Ensure you extracted both `$profile` and `$celltype_contributions` from results
- Check that inside_mode is correctly specified

## License

MIT License - see LICENSE file for details

## Contact

For questions, issues, or contributions:
- GitHub Issues: https://github.com/ndsimons/stGradient/issues
- Email: noah.simons@providence.org


stGradient was developed for analyzing tumor microenvironment architecture in spatial transcriptomics data, with a focus on understanding immune-tumor boundaries and cell type spatial distributions.
