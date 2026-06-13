# Re-align a training centering vector to projected amplitude columns by name

Centering vectors from cancor()/fit_pcr() are named by the PC labels
seen during training. When the projected amplitudes carry the same
labels this is an order-safe no-op; when a label cannot be matched, name
indexing returns NA and would silently corrupt every prediction, so we
abort with a clear error instead. \`cols = NULL\` (an unnamed matrix
product) falls through to positional centering, preserving prior
behavior.

## Usage

``` r
reindex_centering(center, cols, what)
```
