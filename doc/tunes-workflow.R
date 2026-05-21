## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE,
  fig.width = 8,
  fig.height = 5
)

## ----load---------------------------------------------------------------------
library(stGradient)
library(Seurat)
library(dplyr)
library(ggplot2)

# Locate demo data both when developing and when installed.
demo_path <- system.file("demoData", "ps20_27981.rds", package = "stGradient")
if (demo_path == "") demo_path <- file.path("demoData", "ps20_27981.rds")
if (!file.exists(demo_path)) demo_path <- file.path("..", "demoData", "ps20_27981.rds")
stopifnot(file.exists(demo_path))

obj <- load_xenium_seurat(
  rds_path = demo_path,
  celltype_col = "singleR.predicted.id.brcaAtlas"
)

## ----polygons-----------------------------------------------------------------
# Add a broad cancer flag if it does not already exist.
if (!"is_cancer" %in% colnames(obj@meta.data)) {
  ct <- as.character(obj$singleR.predicted.id.brcaAtlas)
  obj$is_cancer <- grepl("Cancer|Tumor|Malignant", ct, ignore.case = TRUE)
}

# If the heuristic yields too few cancer cells, use a fallback based on top niche labels.
if (sum(obj$is_cancer, na.rm = TRUE) < 100) {
  obj$is_cancer <- obj$singleR.predicted.id.brcaAtlas %in% names(sort(table(obj$singleR.predicted.id.brcaAtlas), decreasing = TRUE))[1:2]
}

# Build polygons and add inside/outside membership.
auto_polygons <- plot_dbscan_polygons(obj, plot = FALSE)
obj <- add_polygon_membership(obj, auto_polygons)
obj$boundary_distance <- calculate_boundary_distances(obj, auto_polygons)

table(obj$inside_polygon)

## ----distance-profile---------------------------------------------------------
res <- calculate_distance_profile(
  seurat_obj = obj,
  distance_thresholds = seq(25, 125, by = 25),
  metric = "all",
  group_by = "singleR.predicted.id.brcaAtlas",
  inside_mode = "distance"
)

head(res$profile)

## ----distance-plots, fig.width=10, fig.height=8-------------------------------
plot_distance_comparison(res$profile)
plot_all_metrics(res$profile)
plot_celltype_heatmap(res$celltype_contributions)
plot_celltype_proportions(res$celltype_contributions, top_n = 8)

## ----qc-----------------------------------------------------------------------
qc_out <- qc_plot_xenium(
  obj = obj,
  roi_col = "inside_polygon",
  celltype_col = "singleR.predicted.id.brcaAtlas",
  show = "both",
  zoom = TRUE,
  grey_outside_roi = TRUE
)

qc_out$p

## ----roi-plots, fig.width=10, fig.height=8------------------------------------
plot_roi_celltype_dist(
  df = qc_out$df,
  roi_col = "inside_polygon",
  celltype_col = "singleR.predicted.id.brcaAtlas",
  metric = "proportion",
  labels = "both"
)

plot_roi_inout_celltype_stack(
  df = qc_out$df,
  roi_col = "inside_polygon",
  celltype_col = "singleR.predicted.id.brcaAtlas",
  metric = "proportion"
)

plot_roi_inout_celltype_dodge(
  df = qc_out$df,
  roi_col = "inside_polygon",
  celltype_col = "singleR.predicted.id.brcaAtlas",
  metric = "count",
  top_n = 10,
  top_n_by = "inside"
)

## ----pcf----------------------------------------------------------------------
pcf_summary <- calculate_pcf_outside_summary(
  seurat_obj = obj,
  polygons = auto_polygons,
  distance = 100,
  min_cells = 50,
  celltype_col = "singleR.predicted.id.brcaAtlas",
  include_self_pairs = TRUE,
  ordered_pairs = FALSE,
  show_progress = FALSE
)

head(pcf_summary)

## ----pcf-plot-----------------------------------------------------------------
plot_pcf_outside_heatmap(pcf_summary, midpoint = 1)

## ----pcf-cd8-lumB, fig.width=7, fig.height=5----------------------------------
plot_pcf_pair_gg(
  seurat_obj = obj,
  polygons = auto_polygons,
  type1 = "T cells CD8+",
  type2 = "Cancer LumB SC",
  celltype_col = "singleR.predicted.id.brcaAtlas",
  min_cells = 10,
  nsim = 39,
  bw = 10
)

## ----pcf-b-memory-naive, fig.width=7, fig.height=5----------------------------
plot_pcf_pair_gg(
  seurat_obj = obj,
  polygons = auto_polygons,
  type1 = "B cells Memory",
  type2 = "B cells Naive",
  celltype_col = "singleR.predicted.id.brcaAtlas",
  min_cells = 10,
  nsim = 39,
  bw = 10
)

## ----pcf-envelope, eval=FALSE-------------------------------------------------
#  # Full envelope run for all pairs (slow; saves to PDF).
#  env_list <- plot_pcf_outside_pair_envelopes(
#    seurat_obj = obj,
#    polygons = auto_polygons,
#    celltype_col = "singleR.predicted.id.brcaAtlas",
#    min_cells = 50,
#    include_self_pairs = FALSE,
#    ordered_pairs = FALSE,
#    nsim = 39,
#    bw = 10,
#    file = "pcf_outside_pairs.pdf"
#  )

## ----boundary-gene------------------------------------------------------------
obj <- compute_distance_from_boundary(
  obj,
  polygon_col = "inside_polygon",
  intratumoral_label = TRUE,
  fov_name = "fov"
)

dist_res <- find_distance_associated_genes(
  seurat_obj = obj,
  assay = "Xenium",
  focus = "peritumoral",
  min_cells_expressing = 50,
  n_distance_bins = 8,
  max_dist_quantile = 0.95,
  p_adj_threshold = 0.05
)

if (!is.null(dist_res)) {
  head(dist_res$results)
}

## ----boundary-gene-plot-------------------------------------------------------
if (!is.null(dist_res)) {
  dist_plots <- plot_distance_genes(dist_res, top_n = 12, sample_label = "Demo")
  print(dist_plots$volcano)
  if (!is.null(dist_plots$trends)) print(dist_plots$trends)
}

## ----hotspots-----------------------------------------------------------------
hs <- define_tcell_hotspots(
  seurat_obj = obj,
  celltype_col = "singleR.predicted.id.brcaAtlas",
  tcell_labels = c("T Cells CD4+", "T Cells CD8+"),
  polygon_col = "inside_polygon",
  intratumoral_label = TRUE,
  eps = 50,
  min_pts = 10,
  min_hotspot_size = 5
)

if (!is.null(hs)) {
  plot_tcell_hotspots(obj, hs, celltype_col = "singleR.predicted.id.brcaAtlas")
}

## ----hotspot-genes------------------------------------------------------------
if (!is.null(hs)) {
  hs_res <- find_hotspot_distance_genes(
    seurat_obj = obj,
    hotspot_result = hs,
    assay = "Xenium",
    focus = "perihotspot",
    min_cells_expressing = 50,
    n_distance_bins = 8,
    max_dist_quantile = 0.95,
    p_adj_threshold = 0.05,
    exclude_tcell_labels = c("T Cells CD4+", "T Cells CD8+"),
    celltype_col = "singleR.predicted.id.brcaAtlas"
  )

  if (!is.null(hs_res)) head(hs_res$results)
}

## ----niche-setup--------------------------------------------------------------
set.seed(42)
all_cells <- colnames(obj)
patient_ids <- c("demo_A", "demo_B", "demo_C", "demo_D")

split_idx <- split(sample(all_cells), rep(patient_ids, each = floor(length(all_cells) / length(patient_ids))))

seurat_list <- lapply(patient_ids, function(pid) subset(obj, cells = split_idx[[pid]]))
names(seurat_list) <- patient_ids

pcr_status <- data.frame(
  patient_id = patient_ids,
  pCR = c(1L, 0L, 1L, 0L),
  stringsAsFactors = FALSE
)

all_celltypes <- sort(unique(as.character(obj$singleR.predicted.id.brcaAtlas)))
all_celltypes <- all_celltypes[!grepl("^\\d+$", all_celltypes)]

## ----niche-run----------------------------------------------------------------
intra <- run_niche_pipeline(
  seurat_list = seurat_list,
  subset_value = TRUE,
  pcr_status = pcr_status,
  all_celltypes = all_celltypes,
  celltype_col = "singleR.predicted.id.brcaAtlas",
  polygon_col = "inside_polygon",
  k = 10,
  n_niches = 4
)

intra_niches <- setdiff(colnames(intra$niche_prop_df), c("patient_id", "pCR"))
intra_X <- as.matrix(intra$niche_prop_df[, intra_niches, drop = FALSE])
intra_y <- intra$niche_prop_df$pCR

stat_intra <- run_niche_stat_tests(intra_X, intra_y, intra_niches, n_perm = 100)
lr_intra <- run_niche_univariate_lr(intra_X, intra_y, intra_niches)

head(stat_intra)
head(lr_intra)

## ----niche-plots--------------------------------------------------------------
plot_niche_stacked_composition(intra$niche_prop_df)
plot_niche_grouped_bars(stat_intra, intra_X, intra_y, intra_niches, title = "Intratumoral pseudo-cohort")
plot_niche_coefficients(lr_intra, title = "Intratumoral pseudo-cohort")

## ----shiny-app, eval=FALSE----------------------------------------------------
#  # Interactive polygon drawing app (manual usage)
#  launch_stGradient()

