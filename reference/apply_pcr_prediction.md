# Apply PCR Prediction Transform

Internal function that maps predictor PC amplitudes to predicted
predictand PC amplitudes via the fitted OLS coefficient matrix.

## Usage

``` r
apply_pcr_prediction(new_amplitudes, pcr)
```

## Arguments

- new_amplitudes:

  Tibble with \`time\` and predictor PC amplitudes

- pcr:

  The \`pcr\` slot of a coupled object: \`coefficients\`, \`xcenter\`,
  \`ycenter\`

## Value

Tibble with \`time\` and predicted response PC amplitudes
