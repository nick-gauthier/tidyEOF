#' Preprocessing functions
#'
#' @param x
#' @param var
#' @param flip
#' @param regrid
#' @param monthly
#'
#' @return
#' @export
#'
#' @examples
preprocess <- function(x, var, bbox, flip = FALSE, regrid = FALSE, monthly = FALSE) {
  maps <- brick(x, varname = var)

  indices <- getZ(maps) %>% # find the index for April 1st
    {if(monthly) str_detect(., '-03-') else str_detect(., '-04-01')} %>%
    which()

  subset(maps, indices) %>%
    {if(flip) rotate(.) else .} %>%
    crop(bbox) %>%
    {if(regrid) aggregate(., fact = 2, na.rm = TRUE) else .} # resample to lower res to speed up analysis
}

preprocess_prism <- function(x, bbox, var = 'SWE'){
  maps <- brick(x, varname = var)

  indices <- getZ(maps) %>% # find the indices for March
    str_detect('-03-') %>%
    which()

  raster::subset(maps, indices) %>%
    crop(bbox) %>%
    mean() %>%
    aggregate(fact = 2, fun = "mean") # resample to lower res to speed up analysis
}

snow_only <- function(dat){
  dat %>%
    group_by(x, y) %>%
    mutate(test = all(near(SWE, 0) | is.na(SWE))) %>%
    filter(test == FALSE) %>%
    ungroup() %>%
    dplyr::select(-test)
}
