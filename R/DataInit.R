
#' @param x a network representation as a an adyacency matrix, edgelist,
#' or a network object
#' @importFrom network as.network
#' @noRd

.DataInit <- function(x){
  if(class(x) == "network"){
    x
  }
  if(class(x) == "matrix"){
    x <- as.network(x, loops = TRUE)
  }
  return(x)
}
