# Multivariate EOF/CCA Support — Design

**Date:** 2026-06-11
**Status:** Approved

## Motivation

Downscaling temperature and precipitation with two independent EOF/CCA stacks
can produce physically inconsistent predictions (e.g., warm-wet where the
joint climate is warm-dry). Combined (multivariate) EOF analysis — joint PCA
on multiple variables concatenated along space — extracts modes of covariability
and gives each mode a single shared amplitude time series, so CCA-based
downscaling predicts all variables coherently.

Mechanically this is the dual of `common_patterns()`: that function
concatenates anomalies along **time** (same variable, multiple sources →
shared EOFs, per-source amplitudes); multivariate EOF concatenates along
**space** (same times, multiple variables → per-variable spatial loadings,
shared amplitudes).

## Decisions

| Decision | Choice |
|---|---|
| Input shape | Multi-attribute stars object only (shared grid + times). Variables on different grids are out of scope. |
| Normalization | Multivariate input **requires `scale = TRUE`** (per-pixel standardization via climatological sd). Block normalization is not implemented in v1. |
| CV metrics | Per-variable columns plus a pooled score. The plain metric name (`rmse`, `cor_spatial`, `cor_temporal`) is the pooled score; multivariate runs add `rmse_<var>` etc. |
| Architecture | Attribute-aware core (Approach A). Univariate is the single-block special case; one code path throughout. |
| EOF attribute naming | `patterns$eofs` attributes are named by variable name in **all** cases, replacing the current `"weight"` name for univariate (small breaking change, approved). |
| `common_patterns()` | Keeps its single-attribute gate in v1. Multivariate × multi-source composition is future work. |

## 1. API and constraints

- `patterns(dat, ...)` accepts a multi-attribute stars object. Shared grid
  and time steps are guaranteed by the stars data model.
- If `length(dat) > 1` and `scale = FALSE`, abort with class
  `tidyeof_multivariate_scale` and a message explaining that PCA is
  variance-driven and mixed units require standardization.
- `rotate` and `monthly` work unchanged. Varimax on concatenated loadings is
  mathematically valid; the monthly climatology machinery gains per-attribute
  support (see §3).
- The full stack works with multivariate patterns on either side of the
  coupling: `couple()`, `predict.coupled_patterns()`, `reconstruct()`,
  `project_patterns()`, `tune_eof()`, `prep_cv_folds()`/`tune_cca()`, plot
  methods.
- `couple()` and `apply_cca_prediction()` need **no changes**: they operate
  purely on amplitude tibbles.
- `prep_cv_folds(common_with = ...)` with a multivariate predictor remains
  blocked by `common_patterns()`'s gate; document this limitation.

### Documented caveats (not code)

- Precipitation is skewed; users may want to sqrt- or log-transform before
  analysis. Document in `patterns()` and the vignette.
- Per-pixel standardization is undefined where climatological sd ≈ 0; such
  cells are dropped automatically via the non-finite valid-pixel check (§2).

## 2. Data model

**Concatenated space.** The space-time matrix becomes `time × (V·p)` where V
is the number of variables and p the number of grid cells. Each variable is
flattened with the existing spatial ordering; blocks are cbind'd
variable-major (all of variable 1's cells, then variable 2's, in
`names(dat)` order).

**`patterns` object changes:**

- `block_map` (new): named list mapping variable name → integer column range
  in the concatenated space (e.g., `list(tmean = 1:2601, ppt = 2602:5202)`).
- `eofs`: stars object with one attribute per variable (each `x, y, PC` or
  `geometry, PC`), attributes named by variable. Univariate output is also
  named by its variable (was `"weight"`).
- `valid_pixels`: indexes the concatenated space. The validity check is
  extended from `anyNA` to **any non-finite value**, so Inf/NaN cells produced
  by sd ≈ 0 standardization are dropped like NA cells. Per-variable NA masks
  need no special handling.
- `proj_matrix`: unchanged shape `(n_valid × k)` over concatenated valid
  columns.
- `climatology`: `mean`/`sd` are multi-attribute stars objects.
- `units`, `names`: already per-variable (no change).

## 3. Core flow

**Single-convention helpers.** The flatten/unflatten column-ordering
convention lives in exactly two shared helpers:

- `flatten_time_space()` flattens every attribute and cbinds; returns the
  matrix plus block map and spatial metadata.
- `eof_loading_matrix()` (generalized) stacks the `eofs` attributes back into
  the concatenated `(V·p × k)` matrix.

Every current consumer of `patterns$eofs[[1]]` — `reconstruct()`,
`flip_patterns()`/`compute_eof_signs()`, `get_canonical_patterns()`,
`evaluate_eof_fold()` — switches to the helper so no other code embeds the
ordering.

**Per-function changes:**

- `check_single_attribute()`: removed from `patterns()`/`project_patterns()`
  entry points (retained in `common_patterns()`). Replaced by the
  scale-requirement check and the attribute-match check below.
- `get_climatology()` / `get_anomalies()`: the annual path is expected to
  work via stars attribute-wise arithmetic (verify with tests); the monthly
  path (`apply_monthly_climatology`, monthly branch of `get_climatology`)
  currently flattens `dat[1]` only and gains a per-attribute loop.
- `get_eofs()`: drops the `[1]` subsetting; computes `area_weights()` once
  per grid and replicates the vector once per block; PCA, rotation,
  eigenvalue/North machinery unchanged; loadings are split by `block_map`
  into per-variable arrays and rebuilt as the multi-attribute `eofs` stars
  object.
- `project_patterns()`: validates that `names(newdata)` is a permutation of
  `patterns$names` (reorders internally; aborts with
  `tidyeof_attribute_mismatch` otherwise). Flattening, `valid_pixels`
  subsetting, and `%*% proj_matrix` are unchanged. Area weighting of
  anomalies applies per attribute.
- `reconstruct()`: computes anomalies in concatenated space, splits columns
  by `block_map`, rebuilds each variable via `matrix_to_spacetime()`, and
  combines into a multi-attribute stars object; climatology restoration and
  the existing per-variable units loop apply attribute-wise.
- `[.patterns` (k-truncation used by CV): slices the PC dimension across all
  `eofs` attributes; `proj_matrix` column slicing unchanged.
- Plot methods: EOF maps render one panel row per variable with separate
  fill scales (per-variable ggplot objects assembled with the existing
  patchwork layout machinery). Amplitude plots are unchanged — modes share
  one time series.

## 4. CV metrics

`compute_spatial_metrics()` computes each requested metric per variable block
and a pooled score:

- **Pooled RMSE**: per-variable RMSE normalized by that variable's observed
  standard deviation over the comparison set, then RMS-combined across
  variables (a normalized RMSE; raw-unit pooling across °C and mm is
  meaningless).
- **Pooled correlations**: mean of per-variable correlations (already
  unitless).

Column naming: the plain name (`rmse`, …) is always the pooled score, so
`summarize_cv()`, `tune_cca()`, and `tune_eof()` work unchanged and optimize
the pooled score by default. Multivariate runs add suffixed columns
(`rmse_tmean`, `rmse_ppt`, …) for per-variable inspection; univariate output
is byte-identical to today. `evaluate_eof_fold()`'s speckled holdout samples
hidden cells across the concatenated space and reports per-variable +
pooled metrics using the fold patterns' `block_map`.

## 5. Error handling

| Condition | Behavior |
|---|---|
| Multivariate input with `scale = FALSE` | `tidyeof_multivariate_scale` abort with hint to set `scale = TRUE` |
| `newdata` attributes ≠ training variables | `tidyeof_attribute_mismatch` abort (same names, any order, is accepted and reordered) |
| Multivariate input to `common_patterns()` | existing `tidyeof_multiple_attributes` abort; documented as future work |
| Non-finite cells after standardization | silently excluded via `valid_pixels` (same as NA today) |

## 6. Testing

- **Univariate invariance** (the regression guard): V=1 through the new code
  path reproduces current results numerically for `patterns()`,
  `project_patterns()`, `reconstruct()`, and CV metric values; only the EOF
  attribute name changes (`"weight"` → variable name).
- **Analytic sanity**: joint EOF of a duplicated variable
  `c(a = x, b = x)` yields identical loading blocks; full-rank round-trip
  `reconstruct(patterns(dat, k = max))` recovers the input within tolerance.
- **End-to-end**: univariate coarse predictor → multivariate fine predictand
  through `couple()`/`predict()`/`tune_cca()`; multivariate predictor side;
  monthly mode; differing NA masks per variable; rotation.
- **Metrics**: multivariate runs emit pooled + suffixed columns; univariate
  emits exactly today's columns.
- **Test data**: a helper synthesizes a second attribute from the existing
  PRISM temperature test field (nonlinear transform + noise), keeping
  `inst/testdata` unchanged. A real PRISM `ppt` extract for the vignette is
  optional follow-up at the author's discretion.

## Out of scope (future work)

- Variables on different grids (list-of-stars interface).
- Block normalization (`scale = FALSE` multivariate via domain-mean sd).
- Multivariate sources in `common_patterns()` (space × time concatenation).
