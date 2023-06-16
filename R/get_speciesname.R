#' Function to get species name from GBIF file name
#'
#' @param x a GBIF file name
#'
#' @return the name of the species for the specified GBIF file name
#' @export
#'
#' @examples
get_speciesname = function(x) {
  x <- strsplit(x, 'occurencesGBIF_')
  x <- x[[1]][2]
  x <- strsplit(x, '.csv')
  x <- unlist(x)
  return(x)
}