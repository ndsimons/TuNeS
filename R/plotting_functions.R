#' Plot Distance Profile Comparison
#'
#' Plots transcriptomic distance for different inside modes
#'
#' @param distance_profile_combined Combined data frame from both modes
#' @return ggplot object
#' @importFrom magrittr %>%
#' @export
plot_distance_comparison <- function(distance_profile_combined) {
  ggplot2::ggplot(distance_profile_combined, 
                  ggplot2::aes(x = distance, y = transcriptomic_distance, 
                               color = inside_mode, linetype = inside_mode)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = c("all" = "steelblue", "distance" = "darkorange"),
                                labels = c("all" = "All inside vs boundary outside",
                                           "distance" = "Boundary inside vs boundary outside")) +
    ggplot2::scale_linetype_manual(values = c("all" = "solid", "distance" = "dashed"),
                                   labels = c("all" = "All inside vs boundary outside",
                                              "distance" = "Boundary inside vs boundary outside")) +
    ggplot2::labs(
      title = "Transcriptomic Distance: Comparison of Inside Modes",
      subtitle = "How does the definition of 'inside' affect separation metrics?",
      x = "Distance from boundary (μm)",
      y = "Transcriptomic distance (1 - correlation)",
      color = "Inside definition",
      linetype = "Inside definition"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

#' Plot All Metrics Faceted
#'
#' Creates faceted plot of all separation metrics
#'
#' @param distance_profile_combined Combined data frame from both modes
#' @return ggplot object
#' @importFrom magrittr %>%
#' @export
plot_all_metrics <- function(distance_profile_combined) {
  distance_profile_long <- distance_profile_combined %>%
    tidyr::pivot_longer(cols = c(transcriptomic_distance, composition_difference, de_strength),
                        names_to = "metric",
                        values_to = "value")
  
  ggplot2::ggplot(distance_profile_long, 
                  ggplot2::aes(x = distance, y = value, 
                               color = inside_mode, linetype = inside_mode)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 1) +
    ggplot2::scale_color_manual(values = c("all" = "steelblue", "distance" = "darkorange"),
                                labels = c("all" = "All inside",
                                           "distance" = "Boundary inside")) +
    ggplot2::scale_linetype_manual(values = c("all" = "solid", "distance" = "dashed"),
                                   labels = c("all" = "All inside",
                                              "distance" = "Boundary inside")) +
    ggplot2::labs(
      title = "Genomic Separation: All Metrics Comparison",
      subtitle = "Effect of inside definition on separation metrics",
      x = "Distance from boundary (μm)",
      y = "Separation metric value",
      color = "Inside definition",
      linetype = "Inside definition"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

#' Plot Cell Type Contributions Heatmap
#'
#' Shows which cell types drive differences at each distance
#'
#' @param celltype_contributions Data frame from calculate_distance_profile
#' @return ggplot object
#' @export
plot_celltype_heatmap <- function(celltype_contributions) {
  ggplot2::ggplot(celltype_contributions, 
                  ggplot2::aes(x = distance, y = celltype, fill = contribution)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "white", high = "red", midpoint = 0) +
    ggplot2::labs(
      title = "Cell Type Contribution to Composition Difference",
      subtitle = "Which cell types drive inside/outside differences at each distance?",
      x = "Distance from boundary (μm)",
      y = "Cell Type",
      fill = "Squared\ndifference"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
}

#' Plot Cell Type Proportions
#'
#' Shows proportion trajectories for top cell types
#'
#' @param celltype_contributions Data frame from calculate_distance_profile
#' @param top_n Number of top cell types to show
#' @return ggplot object
#' @importFrom magrittr %>%
#' @export
plot_celltype_proportions <- function(celltype_contributions, top_n = 5) {
  # Find peak distance
  total_by_distance <- celltype_contributions %>%
    dplyr::group_by(distance) %>%
    dplyr::summarise(total_contrib = sum(contribution), .groups = "drop")
  
  peak_distance <- total_by_distance$distance[which.max(total_by_distance$total_contrib)]
  
  # Get top cell types at peak
  top_celltypes <- celltype_contributions %>%
    dplyr::filter(distance == peak_distance) %>%
    dplyr::arrange(dplyr::desc(contribution)) %>%
    dplyr::slice(1:top_n) %>%
    dplyr::pull(celltype)
  
  # Prepare data
  celltype_props_long <- celltype_contributions %>%
    dplyr::filter(celltype %in% top_celltypes) %>%
    tidyr::pivot_longer(cols = c(inside_prop, outside_prop),
                        names_to = "region",
                        values_to = "proportion") %>%
    dplyr::mutate(region = ifelse(region == "inside_prop", "Inside", "Outside"))
  
  ggplot2::ggplot(celltype_props_long, 
                  ggplot2::aes(x = distance, y = proportion, color = region, linetype = region)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::facet_wrap(~celltype, scales = "free_y") +
    ggplot2::scale_color_manual(values = c("Inside" = "red", "Outside" = "blue")) +
    ggplot2::labs(
      title = "Cell Type Proportions: Inside vs Outside",
      subtitle = paste("Top", top_n, "most different cell types"),
      x = "Distance from boundary (μm)",
      y = "Proportion",
      color = "Region",
      linetype = "Region"
    ) +
    ggplot2::theme_minimal()
}
