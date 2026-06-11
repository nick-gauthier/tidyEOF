# Check that a stars object has exactly one attribute

EOF analysis operates on a single variable; silently using the first
attribute of a multi-attribute object would hide the others.

## Usage

``` r
check_single_attribute(
  x,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
)
```

## Arguments

- x:

  A stars object

- arg:

  Argument name for error messages

- call:

  Calling environment for error messages
