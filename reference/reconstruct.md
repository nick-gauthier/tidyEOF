# Reconstruct spatial field from EOF amplitudes

Converts PC amplitudes back into a full spatial-temporal field by
multiplying by the EOF patterns and adding back the climatology.

## Usage

``` r
reconstruct(target_patterns, amplitudes = NULL)
```

## Arguments

- target_patterns:

  A patterns object containing EOFs and climatology

- amplitudes:

  Amplitudes to use for reconstruction. Can be: - NULL (default): uses
  original amplitudes from target_patterns - tibble: with time column
  and PC columns - stars object: will be projected onto patterns first

## Value

A stars object with reconstructed spatial-temporal data

## Details

For bounded variables like precipitation, the reconstructed field may
contain small negative values due to EOF truncation. To clamp these, use
\`mutate(result, across(everything(), ~pmax(.x, 0 \* .x)))\` (the \`0 \*
.x\` trick preserves units).

For multivariate patterns the return is a multi-attribute stars object,
with each variable's climatology and units restored.
