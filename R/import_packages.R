#' @import stars dplyr tidyr cubelyr
#' @importFrom stats screeplot cancor cor cor.test p.adjust pgamma prcomp predict sd setNames time varimax
#' @importFrom utils head data
#' @importFrom sf st_crs
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_vline geom_col facet_wrap scale_fill_distiller coord_quickmap coord_sf theme_void theme_bw labs theme element_text unit margin
#' @importFrom scico scale_fill_scico
#' @importFrom glue glue glue_collapse
#' @importFrom cli cli_h1 cli_h2 cli_text cli_progress_step cli_progress_done cli_abort cli_inform cli_warn
#' @importFrom rlang abort caller_arg caller_env inherits_any .data
#' @importFrom tibble as_tibble tibble
NULL

# Suppress R CMD check notes for tidyverse NSE variables
utils::globalVariables(c(
  ".", "PC", "amplitude", "correlation", "cumulative", "cumvar_line",
  "data", "error", "fold", "hi", "k", "k_cca", "k_pred", "k_resp",
  "low", "multiplet", "percent", "scale_factor", "separated", "std.dev",
  "value"
))
