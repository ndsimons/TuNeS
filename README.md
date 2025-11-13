# TuNeS: Tumor Nest Selector

**Interactive spatial analysis of tumor boundaries in spatial transcriptomics data**

TuNeS provides tools to interactively select tumor regions and analyze genomic separation at varying distances from tumor boundaries. This package is designed for spatial transcriptomics analysis, particularly for understanding tumor microenvironment architecture.

## Features

- **Interactive Polygon Selection**: Shiny app for drawing regions on spatial data
- **Distance-Based Analysis**: Calculate genomic metrics at varying distances from boundaries
- **Multiple Metrics**: Transcriptomic distance, differential expression, cell type composition
- **Visualization**: Comprehensive plotting functions for results
- **Flexible Modes**: Compare entire inside region vs boundary cells

## Installation

```r
# Install devtools if needed
install.packages("devtools")

# Install TuNeS from local directory
devtools::install_local("path/to/TuNeS")

# Or from GitHub (if hosted)
devtools::install_github("yourusername/TuNeS")
```

## Dependencies

TuNeS requires the following packages:
- shiny
- plotly
- sf
- dplyr
- ggplot2
- Seurat
- tidyr
- magrittr

## Quick Start

### 1. Launch Interactive Polygon Selector

```r
library(TuNeS)

# Launch Shiny app
launch_tunes()

# In the app:
# 1. Enter your Seurat object name
# 2. Enter FOV name
# 3. Click "Load Data"
# 4. Draw polygons using the toolbar
# 5. Enter polygon object name
# 6. Click "Extract Polygons"
```

### 2. Add Polygon Information to Seurat Object

```r
# Add inside/outside membership
seurat_obj <- add_polygon_membership(seurat_obj, drawn_polygons)

# Calculate boundary distances
seurat_obj$boundary_distance <- calculate_boundary_distances(seurat_obj, drawn_polygons)
```

### 3. Run Distance Profile Analysis

```r
# Calculate all metrics at varying distances
results <- calculate_distance_profile(
  seurat_obj = seurat_obj,
  distance_thresholds = seq(10, 200, by = 10),
  metric = "all",
  group_by = "cellType",
  inside_mode = "all"
)

# Extract results
distance_profile <- results$profile
celltype_contributions <- results$celltype_contributions
```

### 4. Visualize Results

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

### 5. Find Critical Distances

```r
# Identify distances where separation is maximized
critical_dist <- find_critical_distance(distance_profile, "transcriptomic_distance")

cat("Maximum separation at:", critical_dist$max_distance, "μm\n")
cat("Plateau begins at:", critical_dist$plateau_distance, "μm\n")
```

## Complete Workflow Example

```r
library(TuNeS)
library(Seurat)

# Step 1: Select regions interactively
launch_tunes()
# ... draw polygons and save as 'tumor_regions' ...

# Step 2: Process Seurat object
seurat_obj <- add_polygon_membership(seurat_obj, tumor_regions)
seurat_obj$boundary_distance <- calculate_boundary_distances(seurat_obj, tumor_regions)

# Step 3: Analyze both modes
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

# Step 4: Combine and visualize
combined_results <- rbind(
  results_all$profile,
  results_distance$profile
)

plot_distance_comparison(combined_results)
plot_celltype_heatmap(results_all$celltype_contributions)
```

## Main Functions

### Interactive Selection
- `launch_tunes()` - Launch Shiny app for polygon drawing

### Data Processing
- `add_polygon_membership()` - Add inside/outside labels to Seurat object
- `calculate_boundary_distances()` - Calculate signed distances to boundary

### Metrics
- `calculate_transcriptomic_distance()` - Gene expression correlation
- `calculate_de_strength()` - Differential expression magnitude
- `calculate_composition_difference()` - Cell type distribution differences

### Analysis
- `calculate_distance_profile()` - Main analysis function
- `find_critical_distance()` - Identify key separation distances

### Visualization
- `plot_distance_comparison()` - Compare inside modes
- `plot_all_metrics()` - Faceted view of all metrics
- `plot_celltype_heatmap()` - Cell type contribution heatmap
- `plot_celltype_proportions()` - Proportion trajectories

## Understanding Inside Modes

TuNeS supports two modes for defining "inside" cells:

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
- Typical tumor values: 0.1-0.5

### Composition Difference
- Range: 0 (identical) to ~2 (completely different)
- Sum of squared differences in cell type proportions
- Independent of sample size

### DE Strength
- Mean |log2FC| of top 100 differentially expressed genes
- Higher values indicate stronger transcriptional changes
- Typical values: 0.5-2.0

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

## Citation

If you use TuNeS in your research, please cite:

```
TuNeS: Tumor Nest Selector for Spatial Transcriptomics Analysis
[Your Name et al., Year]
```

## License

MIT License - see LICENSE file for details

## Contact

For questions, issues, or contributions:
- GitHub Issues: [repository URL]
- Email: [contact email]

## Acknowledgments

TuNeS was developed for analyzing tumor microenvironment architecture in spatial transcriptomics data, with a focus on understanding immune-tumor boundaries and cell type spatial distributions.
