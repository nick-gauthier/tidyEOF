# Evaluate one tuning argument with the global RNG temporarily seeded

Runs \`code\` after \`set.seed(seed)\` and restores the previous RNG
state on exit, so cross-validation masks are reproducible without
disturbing the caller's random stream.

## Usage

``` r
with_seed(seed, code)
```
