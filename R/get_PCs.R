#' Get PCs
#'
#' @param pc_object
#' @param k
#'
#' @return
#' @export
#'
#' @examples
get_pcs <- function(pc_object, k){
  tidy(pc_object, 'samples') %>%
    mutate(year = as.numeric(as.character(row))) %>%
    select(-row) %>%
    filter(PC <= k) %>%
    mutate(PC = as.character(PC)) %>%
    rename(amplitude = value)
}
