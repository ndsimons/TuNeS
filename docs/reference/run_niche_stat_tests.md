# Run Niche Differential Tests Between pCR Groups

Run Niche Differential Tests Between pCR Groups

## Usage

``` r
run_niche_stat_tests(X, y, niche_names, n_perm = 1000)
```

## Arguments

- X:

  Numeric matrix/data.frame of patient niche proportions.

- y:

  Integer/binary response vector (1 = pCR, 0 = non-pCR).

- niche_names:

  Character vector of column names in \`X\` to test.

- n_perm:

  Number of permutations for permutation p-values.

## Value

Data frame of per-niche test results.
