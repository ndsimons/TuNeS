# Plot Grouped Niche Bars for pCR vs non-pCR

Plot Grouped Niche Bars for pCR vs non-pCR

## Usage

``` r
plot_niche_grouped_bars(
  stat_res,
  X,
  y,
  niche_names,
  title = "Niche group comparison"
)
```

## Arguments

- stat_res:

  Output from \`run_niche_stat_tests\`.

- X:

  Numeric matrix/data.frame of patient niche proportions.

- y:

  Integer/binary response vector (1 = pCR, 0 = non-pCR).

- niche_names:

  Character vector of niche columns.

- title:

  Plot title.

## Value

ggplot object.
