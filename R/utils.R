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
