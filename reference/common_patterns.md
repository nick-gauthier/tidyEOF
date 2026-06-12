# Compute Common EOF Patterns from Multiple Sources

Performs joint PCA on multiple datasets that share a spatial domain.
Each dataset is anomalized with its own climatology, then concatenated
along time for joint PCA. The resulting shared spatial patterns get
source-specific amplitudes and climatologies, enabling CCA coupling and
cross-source prediction.

## Usage

``` r
common_patterns(
  datasets,
  k = 4,
  scale = TRUE,
  rotate = FALSE,
  monthly = FALSE,
  weight = TRUE,
  irlba_threshold = 5e+05
)
```

## Arguments

- datasets:

  Named list of stars objects sharing the same spatial grid. Names
  become the source identifiers used for extraction.

- k:

  Number of EOF modes to retain

- scale:

  Logical, whether to scale anomalies by standard deviation (default
  TRUE)

- rotate:

  Logical, whether to apply varimax rotation (default FALSE)

- monthly:

  Logical, whether to use monthly climatology (default FALSE)

- weight:

  Logical, whether to apply area weighting (default TRUE)

- irlba_threshold:

  Minimum data elements to trigger IRLBA (default 500000)

## Value

A \`common_patterns\` S3 object. Source-specific patterns are extracted
with \`\$\` or \`\[\[\` using source names (e.g., \`cpat\$era\`). Each
extracted element is a standard \`patterns\` object with shared EOFs but
source-specific climatology and amplitudes.

## Details

Each source must currently contain a single variable; combining
\`common_patterns()\` with multivariate (multi-attribute) input is not
yet supported.

## Examples

``` r
if (FALSE) { # \dontrun{
cpat <- common_patterns(
  list(era = era_coarse, phyda = phyda_coarse),
  k = 11, scale = TRUE
)

# Extract source-specific patterns
cpat$era     # patterns object with ERA climatology + amplitudes
cpat$phyda   # patterns object with PHYDA climatology + amplitudes

# Use with existing couple/predict workflow
fine_pat <- patterns(era_fine, k = 13)
coupled <- couple(cpat$era, fine_pat, k = 9)
predict(coupled, era_new)

# Cross-source prediction
predict(coupled, phyda_new, predictor_patterns = cpat$phyda)
} # }
```
