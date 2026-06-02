# Plot Outside-Polygon PCF Heatmap

Converts PCF summary output into a symmetric matrix and visualizes it as
a heatmap.

## Usage

``` r
plot_pcf_outside_heatmap(pcf_summary, midpoint = 1)
```

## Arguments

- pcf_summary:

  Output from \`calculate_pcf_outside_summary\`.

- midpoint:

  Midpoint value for diverging color scale.

## Value

ggplot object.
