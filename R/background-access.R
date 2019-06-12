#' Get the background paths
#'
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
get_background_paths <- function(n = 1){
  flatten(spagi_universe[1:n])
}

