
#' @param x a network representation as a an adyacency matrix, edgelist,
#' or a network object
#' @importFrom network as.network
#' @importFrom methods is
#' @noRd

.DataInit <- function(x){
  if(is(x) == "network"){
    x
  }
  if(is(x) == "matrix"){
    x <- as.network(x, loops = TRUE)
  }
  return(x)
}
