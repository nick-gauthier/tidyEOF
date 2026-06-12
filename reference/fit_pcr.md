# Fit the OLS coefficient matrix for PCR coupling

Ordinary least squares of the predictand PC amplitudes on the predictor
PC amplitudes. The EOF truncation already performed the dimension
reduction, so this is plain (optionally centered) multivariate
regression. The fit is rank-safe: a rank-deficient predictor matrix
(e.g. \`k_pred \>= n_times\` or collinear amplitudes) aborts rather than
returning \`NA\` coefficients.

## Usage

``` r
fit_pcr(pred_amps, resp_amps, center = TRUE)
```

## Arguments

- pred_amps:

  Predictor amplitude matrix (time x k_pred)

- resp_amps:

  Predictand amplitude matrix (time x k_resp)

- center:

  Logical; center both sides before fitting (default TRUE)

## Value

List with \`coefficients\` (k_pred x k_resp), \`xcenter\`, \`ycenter\`
(each a named numeric vector when centered, or \`FALSE\`)
