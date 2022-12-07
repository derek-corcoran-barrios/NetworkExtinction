
#' @param x a network representation as a an adyacency matrix, edgelist,
#' or a network object
#' @importFrom network as.network
#' @importFrom methods is
#' @noRd

.DataInit <- function(x){
  if(is(x)[1] == "network"){
    x
  }
  if(is(x)[1] == "matrix"){
    x <- network::as.network(x,
                             matrix.type='adjacency',
                             # loops = TRUE,
                             directed = TRUE,
                             ignore.eval = FALSE,
                             names.eval = 'weight')
  }
  return(x)
}
