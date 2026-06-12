# Flatten a time-indexed stars object into a matrix

Internal helper to turn a stars object (one or more attributes) with a
\`time\` dimension into a matrix of shape time x (V \* space) along with
the metadata needed to reconstruct the original spatial layout. When the
object has multiple attributes they are concatenated variable-major: all
columns for the first attribute, then all columns for the second, etc.
The returned \`block_map\` is a named list mapping each attribute name
to its column indices, and \`n_space\` is the number of spatial cells
per attribute.

## Usage

``` r
flatten_time_space(dat)
```

## Arguments

- dat:

  A stars object with a \`time\` dimension

## Value

A list containing the flattened matrix, block_map (named list of column
index ranges per attribute), n_space (spatial cells per attribute),
spatial dimension names, their shape, and coordinate values
