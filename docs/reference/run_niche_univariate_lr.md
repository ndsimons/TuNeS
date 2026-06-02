# Run Univariate Logistic Regression for Each Niche

Run Univariate Logistic Regression for Each Niche

## Usage

``` r
run_niche_univariate_lr(X, y, niche_names)
```

## Arguments

- X:

  Numeric matrix/data.frame of patient niche proportions.

- y:

  Integer/binary response vector (1 = pCR, 0 = non-pCR).

- niche_names:

  Character vector of column names in \`X\` to test.

## Value

Data frame with per-niche logistic regression statistics.
