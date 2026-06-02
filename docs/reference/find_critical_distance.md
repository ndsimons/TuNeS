# Find Critical Distance

Identifies distance at which separation is maximized

## Usage

``` r
find_critical_distance(
  distance_profile,
  metric_column = "transcriptomic_distance"
)
```

## Arguments

- distance_profile:

  Data frame from calculate_distance_profile

- metric_column:

  Column name of metric to analyze

## Value

List with max_distance, max_value, and plateau_distance
