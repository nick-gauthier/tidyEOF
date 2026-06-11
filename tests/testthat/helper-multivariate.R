# Build a two-variable stars object from the single-variable PRISM test data.
# The second attribute is a deterministic nonlinear transform of temperature
# (no RNG, so helpers stay reproducible without touching the seed).
make_multivar <- function(x) {
  vals <- units::drop_units(x[[1]])
  idx <- array(seq_along(vals), dim = dim(vals))
  ppt <- exp(0.15 * vals) + 3 * sin(idx / 7)
  out <- x
  out$ppt <- units::set_units(ppt, "mm")
  out
}
