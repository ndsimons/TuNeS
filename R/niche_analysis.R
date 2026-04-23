#' Run Spatial Niche Pipeline for One Compartment
#'
#' Computes spatial neighborhoods, clusters niche archetypes, labels niches, and
#' returns patient-level niche proportions for either intratumoral or peritumoral
#' cells.
#'
#' @param seurat_list Named list of Seurat objects (one per patient/sample).
#' @param subset_value Logical value in `polygon_col` to select compartment.
#' @param pcr_status Data frame with `patient_id` and `pCR` columns.
#' @param all_celltypes Character vector of cell type levels to use.
#' @param celltype_col Metadata column containing cell types.
#' @param polygon_col Metadata column indicating inside/outside polygon.
#' @param k Number of nearest neighbors for neighborhood composition.
#' @param n_niches Number of k-means niche archetypes.
#' @param fov_name Field of view/image name for coordinate extraction.
#' @return List with neighborhood-level data, archetypes, labels, and patient proportions.
#' @export
run_niche_pipeline <- function(seurat_list,
                               subset_value,
                               pcr_status,
                               all_celltypes,
                               celltype_col = "singleR.predicted.id.brcaAtlas",
                               polygon_col = "inside_polygon",
                               k = 20,
                               n_niches = 8,
                               fov_name = "fov") {
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required. Install with: install.packages('FNN')")
  }

  neighborhood_list <- lapply(names(seurat_list), function(pid) {
    obj <- seurat_list[[pid]]
    meta <- obj@meta.data

    if (!polygon_col %in% colnames(meta)) {
      return(NULL)
    }

    keep <- which(meta[[polygon_col]] == subset_value)
    if (length(keep) < k + 1) {
      return(NULL)
    }

    meta_sub <- meta[keep, , drop = FALSE]
    tissue_coords <- Seurat::GetTissueCoordinates(obj, image = fov_name)

    if ("cell" %in% colnames(tissue_coords)) {
      rownames(tissue_coords) <- tissue_coords$cell
      tissue_coords$cell <- NULL
    }

    tissue_coords <- tissue_coords[rownames(meta_sub), , drop = FALSE]
    coords <- as.matrix(tissue_coords[, c("x", "y"), drop = FALSE])

    ct <- factor(meta_sub[[celltype_col]], levels = all_celltypes)
    valid_cells <- !is.na(ct)

    knn_out <- FNN::get.knn(coords, k = k)
    neighbor_ct <- matrix(
      0,
      nrow = nrow(coords),
      ncol = length(all_celltypes),
      dimnames = list(NULL, all_celltypes)
    )

    for (i in seq_len(nrow(coords))) {
      nn_idx <- knn_out$nn.index[i, ]
      nn_types <- ct[nn_idx]
      nn_types <- nn_types[!is.na(nn_types)]
      if (length(nn_types) > 0) {
        tab <- table(nn_types) / length(nn_types)
        neighbor_ct[i, names(tab)] <- as.numeric(tab)
      }
    }

    data.frame(
      patient_id = pid,
      cell_barcode = rownames(meta_sub)[valid_cells],
      cell_type = as.character(ct[valid_cells]),
      neighbor_ct[valid_cells, , drop = FALSE],
      check.names = FALSE
    )
  })

  all_neighborhoods <- dplyr::bind_rows(neighborhood_list)
  if (nrow(all_neighborhoods) == 0) {
    stop("No neighborhoods available. Check polygon labels and minimum cell counts.")
  }

  comp_mat <- as.matrix(all_neighborhoods[, all_celltypes, drop = FALSE])
  km <- stats::kmeans(comp_mat, centers = n_niches, nstart = 25, iter.max = 100)
  all_neighborhoods$niche <- km$cluster

  archetype_profiles <- all_neighborhoods %>%
    dplyr::group_by(niche) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(all_celltypes), mean), n_cells = dplyr::n(), .groups = "drop")

  arch_mat <- as.matrix(archetype_profiles[, all_celltypes, drop = FALSE])

  niche_labels <- sapply(seq_len(n_niches), function(i) {
    sorted <- sort(arch_mat[i, ], decreasing = TRUE)
    top_types <- names(sorted)[1:min(3, length(sorted))]
    top_props <- sorted[1:min(3, length(sorted))]

    if (top_props[1] > 0.40) {
      label <- top_types[1]
      if (length(top_props) > 1 && top_props[2] > 0.10) {
        label <- paste0(label, " + ", top_types[2])
      }
    } else if (length(top_props) > 1 && top_props[1] > 0.15 && top_props[2] > 0.15) {
      label <- paste0(top_types[1], " + ", top_types[2])
    } else {
      label <- paste0("Mixed: ", paste(top_types[1:min(2, length(top_types))], collapse = " / "))
    }
    paste0(label, " niche")
  })
  names(niche_labels) <- as.character(seq_len(n_niches))

  if (any(duplicated(niche_labels))) {
    for (d in unique(niche_labels[duplicated(niche_labels)])) {
      idx <- which(niche_labels == d)
      niche_labels[idx] <- paste0(niche_labels[idx], " (", idx, ")")
    }
  }

  all_neighborhoods$niche_label <- niche_labels[as.character(all_neighborhoods$niche)]

  niche_prop_df <- all_neighborhoods %>%
    dplyr::group_by(patient_id, niche_label) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
    dplyr::mutate(proportion = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-n) %>%
    tidyr::pivot_wider(names_from = niche_label, values_from = proportion, values_fill = 0)

  niche_prop_df <- dplyr::inner_join(
    niche_prop_df,
    pcr_status %>% dplyr::select(patient_id, pCR),
    by = "patient_id"
  )

  list(
    neighborhoods = all_neighborhoods,
    archetype_profiles = archetype_profiles,
    arch_mat = arch_mat,
    niche_labels = niche_labels,
    niche_prop_df = niche_prop_df
  )
}

#' Run Niche Differential Tests Between pCR Groups
#'
#' @param X Numeric matrix/data.frame of patient niche proportions.
#' @param y Integer/binary response vector (1 = pCR, 0 = non-pCR).
#' @param niche_names Character vector of column names in `X` to test.
#' @param n_perm Number of permutations for permutation p-values.
#' @return Data frame of per-niche test results.
#' @export
run_niche_stat_tests <- function(X, y, niche_names, n_perm = 1000) {
  dplyr::bind_rows(lapply(niche_names, function(niche) {
    vals_pcr <- X[y == 1, niche]
    vals_nonpcr <- X[y == 0, niche]

    if (stats::sd(c(vals_pcr, vals_nonpcr)) == 0) {
      return(data.frame(
        niche = niche,
        mean_pCR = mean(vals_pcr),
        mean_non_pCR = mean(vals_nonpcr),
        diff = 0,
        wilcox_p = NA_real_,
        perm_p = NA_real_,
        direction = "no difference",
        stringsAsFactors = FALSE
      ))
    }

    wt <- stats::wilcox.test(vals_pcr, vals_nonpcr, exact = FALSE)
    obs_diff <- mean(vals_pcr) - mean(vals_nonpcr)
    all_vals <- c(vals_pcr, vals_nonpcr)
    n1 <- length(vals_pcr)

    perm_diffs <- replicate(n_perm, {
      shuf <- sample(all_vals)
      mean(shuf[1:n1]) - mean(shuf[(n1 + 1):length(shuf)])
    })
    perm_p <- (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_perm + 1)

    data.frame(
      niche = niche,
      mean_pCR = mean(vals_pcr),
      mean_non_pCR = mean(vals_nonpcr),
      diff = obs_diff,
      wilcox_p = wt$p.value,
      perm_p = perm_p,
      direction = ifelse(obs_diff > 0, "enriched in pCR", "enriched in non-pCR"),
      stringsAsFactors = FALSE
    )
  })) %>%
    dplyr::mutate(
      wilcox_padj = stats::p.adjust(wilcox_p, method = "BH"),
      perm_padj = stats::p.adjust(perm_p, method = "BH")
    ) %>%
    dplyr::arrange(perm_p)
}

#' Run Univariate Logistic Regression for Each Niche
#'
#' @param X Numeric matrix/data.frame of patient niche proportions.
#' @param y Integer/binary response vector (1 = pCR, 0 = non-pCR).
#' @param niche_names Character vector of column names in `X` to test.
#' @return Data frame with per-niche logistic regression statistics.
#' @export
run_niche_univariate_lr <- function(X, y, niche_names) {
  dplyr::bind_rows(lapply(niche_names, function(niche) {
    fit <- stats::glm(y ~ X[, niche], family = stats::binomial)
    s <- summary(fit)
    coef_row <- s$coefficients[2, , drop = FALSE]

    data.frame(
      niche = niche,
      estimate = coef_row[1, "Estimate"],
      std_error = coef_row[1, "Std. Error"],
      z_value = coef_row[1, "z value"],
      p_value = coef_row[1, "Pr(>|z|)"],
      odds_ratio = exp(coef_row[1, "Estimate"]),
      direction = ifelse(coef_row[1, "Estimate"] > 0, "higher in pCR", "higher in non-pCR"),
      stringsAsFactors = FALSE
    )
  })) %>%
    dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(p_value)
}

#' Plot Patient Niche Stacked Composition
#'
#' @param niche_prop_df Output `niche_prop_df` from `run_niche_pipeline`.
#' @return ggplot object.
#' @export
plot_niche_stacked_composition <- function(niche_prop_df) {
  niche_cols <- setdiff(colnames(niche_prop_df), c("patient_id", "pCR"))

  niche_prop_df %>%
    dplyr::mutate(pCR_group = ifelse(pCR == 1, "pCR", "non-pCR")) %>%
    dplyr::arrange(dplyr::desc(pCR), patient_id) %>%
    dplyr::mutate(patient_id = factor(patient_id, levels = patient_id)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(niche_cols), names_to = "niche", values_to = "proportion") %>%
    ggplot2::ggplot(ggplot2::aes(x = patient_id, y = proportion, fill = niche)) +
    ggplot2::geom_bar(stat = "identity", width = 0.85) +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    ggplot2::facet_grid(~pCR_group, scales = "free_x", space = "free_x") +
    ggplot2::labs(x = NULL, y = "Proportion", fill = "Niche", title = "Niche composition by patient") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "bottom",
      panel.grid.major.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

#' Plot Grouped Niche Bars for pCR vs non-pCR
#'
#' @param stat_res Output from `run_niche_stat_tests`.
#' @param X Numeric matrix/data.frame of patient niche proportions.
#' @param y Integer/binary response vector (1 = pCR, 0 = non-pCR).
#' @param niche_names Character vector of niche columns.
#' @param title Plot title.
#' @return ggplot object.
#' @export
plot_niche_grouped_bars <- function(stat_res, X, y, niche_names, title = "Niche group comparison") {
  sem_df <- dplyr::bind_rows(lapply(niche_names, function(niche) {
    data.frame(
      niche = niche,
      group = c("pCR", "non-pCR"),
      sem = c(
        stats::sd(X[y == 1, niche]) / sqrt(sum(y == 1)),
        stats::sd(X[y == 0, niche]) / sqrt(sum(y == 0))
      ),
      stringsAsFactors = FALSE
    )
  }))

  grouped_df <- stat_res %>%
    dplyr::select(niche, mean_pCR, mean_non_pCR, perm_padj) %>%
    tidyr::pivot_longer(cols = c(mean_pCR, mean_non_pCR), names_to = "group", values_to = "mean_prop") %>%
    dplyr::mutate(
      group = dplyr::recode(group, mean_pCR = "pCR", mean_non_pCR = "non-pCR"),
      sig_label = dplyr::case_when(
        perm_padj < 0.01 ~ "**",
        perm_padj < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::left_join(sem_df, by = c("niche", "group"))

  ggplot2::ggplot(grouped_df, ggplot2::aes(x = stats::reorder(niche, -mean_prop), y = mean_prop, fill = group)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(0.8), width = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_prop - sem, ymax = mean_prop + sem),
      position = ggplot2::position_dodge(0.8), width = 0.2, linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = grouped_df %>% dplyr::filter(group == "pCR", sig_label != ""),
      ggplot2::aes(x = niche, y = mean_prop + sem + 0.02, label = sig_label),
      position = ggplot2::position_dodge(0.8),
      color = "#D4537E", size = 5, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = c(pCR = "#1D9E75", `non-pCR` = "#D85A30")) +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(x = NULL, y = "Mean proportion", fill = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, size = 8), legend.position = "top")
}

#' Plot Niche Logistic Coefficients
#'
#' @param univar_res Output from `run_niche_univariate_lr`.
#' @param title Plot title.
#' @return ggplot object.
#' @export
plot_niche_coefficients <- function(univar_res, title = "Niche logistic coefficients") {
  ggplot2::ggplot(
    univar_res,
    ggplot2::aes(x = stats::reorder(niche, estimate), y = estimate, fill = ifelse(estimate > 0, "pCR", "non-pCR"))
  ) +
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = estimate - 1.96 * std_error, ymax = estimate + 1.96 * std_error),
      width = 0.2, linewidth = 0.4
    ) +
    ggplot2::scale_fill_manual(values = c(pCR = "#1D9E75", `non-pCR` = "#D85A30")) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Log odds ratio (95% CI)", fill = "Associated with", title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "top")
}
