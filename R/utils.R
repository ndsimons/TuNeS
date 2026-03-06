#' Parse SVG Path to Coordinates
#'
#' Internal function to parse plotly SVG paths
#'
#' @param path_string SVG path string from plotly
#' @return Matrix of x,y coordinates
#' @keywords internal
parse_svg_path <- function(path_string) {
  path_clean <- gsub("M|Z", "", path_string)
  
  coords_str <- strsplit(path_clean, "L")[[1]]
  coords_str <- coords_str[coords_str != ""]
  
  coords_list <- lapply(coords_str, function(coord) {
    vals <- as.numeric(strsplit(trimws(coord), ",")[[1]])
    if (length(vals) == 2) return(vals)
    return(NULL)
  })
  
  coords_list <- coords_list[!sapply(coords_list, is.null)]
  
  if (length(coords_list) > 0) {
    coords_matrix <- do.call(rbind, coords_list)
    colnames(coords_matrix) <- c("x", "y")
    
    if (!all(coords_matrix[1,] == coords_matrix[nrow(coords_matrix),])) {
      coords_matrix <- rbind(coords_matrix, coords_matrix[1,])
    }
    
    return(coords_matrix)
  }
  
  return(NULL)
}

#' Add Polygon Membership to Seurat Object
#'
#' Adds metadata column indicating if cells are inside polygons
#'
#' @param seurat_obj Seurat object
#' @param polygons sf polygon object
#' @param fov_name FOV name
#' @return Seurat object with added metadata
#' @export
add_polygon_membership <- function(seurat_obj, polygons, fov_name = "fov") {
  coords <- Seurat::GetTissueCoordinates(seurat_obj, image = fov_name)
  coords_sf <- sf::st_as_sf(coords, coords = c("x", "y"))
  
  cells_inside <- sf::st_intersects(coords_sf, polygons, sparse = FALSE)
  is_inside <- apply(cells_inside, 1, any)
  
  seurat_obj$inside_polygon <- is_inside
  
  return(seurat_obj)
}

#' Automatically detect cancer cell regions and create polygons
#'
#' Uses DBSCAN clustering on cancer cell coordinates to identify spatial regions,
#' then creates concave hull polygons around each cluster. Adds polygons to the
#' Seurat object's FOV.
#'
#' @param seurat_obj Seurat object with spatial coordinates
#' @param fov_name Name of the field of view (default: "fov")
#' @param cancer_col Metadata column indicating cancer cells (default: "is_cancer")
#' @param eps DBSCAN epsilon parameter - max distance between points (default: 35)
#' @param minPts DBSCAN minimum points per cluster (default: 20)
#' @param concavity Concave hull concavity parameter, lower = tighter fit (default: 2)
#' @param length_threshold Edge length threshold for concave hull (default: 0)
#'
#' @return Seurat object with cancer region polygons added to FOV
#' @export
#' 
#' @examples
#' seurat_obj <- detect_cancer_regions(seurat_obj, eps = 50, minPts = 30)
#'
add_dbscan_polygons <- function(seurat_obj,
                                   fov_name = "fov",
                                   cancer_col = "is_cancer",
                                   eps = 35,
                                   minPts = 20,
                                   concavity = 2,
                                   length_threshold = 0) {
  
  # Check required packages
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required. Install with: install.packages('dbscan')")
  }
  if (!requireNamespace("concaveman", quietly = TRUE)) {
    stop("Package 'concaveman' is required. Install with: install.packages('concaveman')")
  }
  
  # Get tissue coordinates
  coords <- Seurat::GetTissueCoordinates(seurat_obj, image = fov_name)
  
  # Check if cancer column exists
  if (!cancer_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste0("Column '", cancer_col, "' not found in metadata"))
  }
  
  coords$is_cancer <- seurat_obj@meta.data[[cancer_col]]
  
  # Filter to cancer cells only
  cancer_coords <- coords %>% dplyr::filter(is_cancer == TRUE)
  
  if (nrow(cancer_coords) == 0) {
    stop("No cancer cells found in metadata")
  }
  
  # DBSCAN clustering
  message("Clustering cancer cells with DBSCAN...")
  clusters <- dbscan::dbscan(cancer_coords[, c("x", "y")], 
                             eps = eps, 
                             minPts = minPts)
  cancer_coords$cluster <- clusters$cluster
  
  # Remove noise points (cluster 0)
  cancer_coords <- cancer_coords %>% dplyr::filter(cluster > 0)
  
  if (nrow(cancer_coords) == 0) {
    stop("No clusters found. Try adjusting eps or minPts parameters")
  }
  
  n_clusters <- length(unique(cancer_coords$cluster))
  message(paste0("Found ", n_clusters, " cancer regions"))
  
  # Create polygons for each cluster
  polygon_list <- list()
  for (cluster_id in unique(cancer_coords$cluster)) {
    cluster_coords <- cancer_coords %>% dplyr::filter(cluster == cluster_id)
    cluster_matrix <- as.matrix(cluster_coords[, c("x", "y")])
    
    # Create concave hull
    concave <- concaveman::concaveman(cluster_matrix,
                                      concavity = concavity,
                                      length_threshold = length_threshold)
    
    # Convert to sf polygon
    polygon_list[[cluster_id]] <- sf::st_polygon(list(concave)) %>%
      sf::st_sfc() %>%
      sf::st_sf(region_id = cluster_id)
  }
  
  # Combine all polygons
  auto_polygons <- do.call(rbind, polygon_list)
  
  # Add polygons to Seurat object FOV
  message("Adding polygons to Seurat object...")
  #seurat_obj@images[[fov_name]]@boundaries$cancer_regions <- auto_polygons
  seurat_obj <- add_polygon_membership(seurat_obj,auto_polygons)
  message(paste0("Successfully added ", n_clusters, " cancer region polygons"))
  
  return(list(
    seurat_obj = seurat_obj,
    dbscan_polygons = auto_polygons
  ))
}

#' Plot DBSCAN-based cancer region polygons
#'
#' Uses DBSCAN clustering on cancer cell coordinates to identify spatial regions,
#' creates concave hull polygons around each cluster, and optionally visualizes them.
#'
#' @param seurat_obj Seurat object with spatial coordinates
#' @param fov Name of the field of view (default: "fov")
#' @param eps DBSCAN epsilon parameter - max distance between points (default: 35)
#' @param minPts DBSCAN minimum points per cluster (default: 20)
#' @param concavity Concave hull concavity parameter, lower = tighter fit (default: 2)
#' @param length_threshold Edge length threshold for concave hull (default: 0)
#' @param sample_name Sample name for plot title
#' @param plot Whether to generate visualization plots (default: TRUE)
#'
#' @return sf object containing polygon geometries for each cancer region
#' @export
#' 
#' @importFrom Seurat GetTissueCoordinates
#' @importFrom sf st_as_sf st_polygon st_sfc st_sf
#' @importFrom dbscan dbscan
#' @importFrom concaveman concaveman
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot geom_sf aes scale_color_manual scale_x_reverse scale_y_reverse theme_minimal theme element_rect element_blank element_text ggtitle
#' @importFrom cowplot plot_grid
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' polygons <- plotDBSCANpolygons(seurat_obj, sample_name = "Patient1")
#' 
#' # Custom parameters without plotting
#' polygons <- plotDBSCANpolygons(
#'   seurat_obj, 
#'   eps = 50, 
#'   minPts = 30,
#'   sample_name = "Patient1",
#'   plot = FALSE
#' )
#' }
#'
plot_dbscan_polygons <- function(seurat_obj, 
                                fov = "fov", 
                                eps = 35, 
                                minPts = 20, 
                                concavity = 2,
                                length_threshold = 0,
                                sample_name = NULL,
                                plot = TRUE) {
  
  # Check required packages
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required. Install with: install.packages('dbscan')")
  }
  if (!requireNamespace("concaveman", quietly = TRUE)) {
    stop("Package 'concaveman' is required. Install with: install.packages('concaveman')")
  }
  
  # Get coordinates and cancer annotation
  coords <- Seurat::GetTissueCoordinates(seurat_obj, image = fov)
  coords$is_cancer <- seurat_obj$is_cancer
  
  # Convert to sf object
  coords_sf <- sf::st_as_sf(coords, coords = c("x", "y"))
  
  # Filter cancer cells
  cancer_coords <- coords %>% dplyr::filter(is_cancer == TRUE)
  
  if (nrow(cancer_coords) == 0) {
    stop("No cancer cells found in seurat_obj$is_cancer")
  }
  
  # DBSCAN clustering
  clusters <- dbscan::dbscan(cancer_coords[, c("x", "y")], eps = eps, minPts = minPts)
  cancer_coords$cluster <- clusters$cluster
  
  # Remove noise points (cluster 0)
  cancer_coords <- cancer_coords %>% dplyr::filter(cluster > 0)
  
  if (nrow(cancer_coords) == 0) {
    warning("No clusters found. Try adjusting eps or minPts parameters")
    return(NULL)
  }
  
  # Create polygons for each cluster
  polygon_list <- list()
  for (cluster_id in unique(cancer_coords$cluster)) {
    cluster_coords <- cancer_coords %>% dplyr::filter(cluster == cluster_id)
    cluster_matrix <- as.matrix(cluster_coords[, c("x", "y")])
    
    # Create concave hull
    concave <- concaveman::concaveman(cluster_matrix, 
                                      concavity = concavity,
                                      length_threshold = length_threshold)
    
    # Convert to sf polygon
    polygon_list[[cluster_id]] <- sf::st_polygon(list(concave)) %>%
      sf::st_sfc() %>%
      sf::st_sf(region_id = cluster_id)
  }
  
  # Combine all polygons
  auto_polygons <- do.call(rbind, polygon_list)
  
  # Visualize if requested
  if (plot) {
    if (is.null(sample_name)) {
      sample_name <- "Sample"
    }
    
    p1 <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = coords_sf, ggplot2::aes(color = is_cancer), size = 0.005, alpha = 0.5) +
      ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "lightgray")) +
      ggplot2::scale_x_reverse() +
      ggplot2::scale_y_reverse() +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "black"),
                     plot.background = ggplot2::element_rect(fill = "black"),
                     panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(color = "white")) +
      ggplot2::ggtitle(paste0(sample_name, ": Xenium cancer cell annotations"))
    
    p2 <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = coords_sf, ggplot2::aes(color = is_cancer), size = 0.005, alpha = 0.5) +
      ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "lightgray")) +
      ggplot2::geom_sf(data = auto_polygons, 
                       fill = NA,
                       color = "blue",
                       linewidth = 1) +
      ggplot2::scale_x_reverse() +
      ggplot2::scale_y_reverse() +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "black"),
                     plot.background = ggplot2::element_rect(fill = "black"),
                     panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(color = "white")) +
      ggplot2::ggtitle(paste0(sample_name, ": cluster polygons"))
    
    print(cowplot::plot_grid(p1, p2))
  }
  
  # Return the polygons
  return(auto_polygons)
}
