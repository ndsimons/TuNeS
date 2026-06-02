# Plot Cell Type Distribution Inside ROI

Creates a horizontal bar plot of cell type counts or proportions among
ROI cells.

## Usage

``` r
plot_roi_celltype_dist(
  df,
  roi_col,
  celltype_col,
  celltype_colors = NULL,
  metric = c("count", "proportion"),
  order_by = c("value", "palette", "none"),
  legend = FALSE,
  title = "Cell type distribution inside ROI",
  xlab = "Cell type",
  ylab = NULL,
  labels = c("none", "count", "percent", "both"),
  label_size = 3,
  label_hjust = -0.1,
  label_accuracy = 0.1,
  label_pad_mult = 0.15
)
```

## Arguments

- df:

  Data frame with ROI and cell type columns

- roi_col:

  ROI indicator column name

- celltype_col:

  Cell type column name

- celltype_colors:

  Optional named cell type color vector

- metric:

  One of "count" or "proportion"

- order_by:

  One of "value", "palette", "none"

- legend:

  Logical; show legend

- title:

  Plot title

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

- label_accuracy:

  Percent label rounding accuracy

- label_pad_mult:

  Extra axis expansion for labels

## Value

ggplot object
