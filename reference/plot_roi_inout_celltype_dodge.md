# Plot Dodged Cell Type Comparison Inside vs Outside ROI

Creates a horizontal dodged bar plot comparing inside and outside ROI
counts or proportions by cell type.

## Usage

``` r
plot_roi_inout_celltype_dodge(
  df,
  roi_col,
  celltype_col,
  metric = c("count", "proportion"),
  roi_colors = c(`FALSE` = "grey80", `TRUE` = "red"),
  inside_label = "Inside ROI",
  outside_label = "Outside ROI",
  order_by = c("total", "inside", "none"),
  top_n = NULL,
  top_n_by = c("total", "inside"),
  other_label = "Other",
  keep_other = TRUE,
  legend = TRUE,
  legend_title = "ROI status",
  legend_position = "right",
  title = NULL,
  xlab = "Cell type",
  ylab = NULL,
  labels = c("none", "count", "percent", "both"),
  label_size = 3,
  label_hjust = -0.1,
  label_pad_mult = 0.15,
  label_accuracy = 0.1
)
```

## Arguments

- df:

  Data frame with ROI and cell type columns

- roi_col:

  ROI indicator column name

- celltype_col:

  Cell type column name

- metric:

  One of "count" or "proportion"

- roi_colors:

  Named colors for \`FALSE\` and \`TRUE\`

- inside_label:

  Label for inside ROI group

- outside_label:

  Label for outside ROI group

- order_by:

  One of "total", "inside", "none"

- top_n:

  Optional number of top cell types to keep

- top_n_by:

  One of "total" or "inside"

- other_label:

  Label for collapsed non-top groups

- keep_other:

  Logical; keep collapsed other group

- legend:

  Logical; show legend

- legend_title:

  Legend title

- legend_position:

  Legend position

- title:

  Optional plot title

- xlab:

  X-axis label

- ylab:

  Optional y-axis label

- labels:

  One of "none", "count", "percent", "both"

- label_size:

  Label text size

- label_hjust:

  Label horizontal adjustment

- label_pad_mult:

  Extra axis expansion for labels

- label_accuracy:

  Percent label rounding accuracy

## Value

ggplot object
