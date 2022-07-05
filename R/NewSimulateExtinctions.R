SimulateExtinctions <- function(Network, Method, Order = NULL,
                                clust.method = "cluster_infomap",
                                IS = 0,
                                Rewiring = FALSE, RewiringDist = NULL,
                                verbose = TRUE){
  Network <- .DataInit(x = Network)

  if(!is.null(Order)){Method <- "Ordered"}

  '%ni%'<- Negate('%in%')
  if(Method %ni% c("Mostconnected", "Ordered")) stop('Choose the right method. See ?SimulateExtinction.')

  if(Method == "Mostconnected"){
    edgelist <- as.matrix.network.edgelist(Network,matrix.type="edgelist") #Prey - Predator
    Conected <- data.frame(ID = 1:network.size(Network), Grado = degree(edgelist, c("total")))
    Conected <- arrange(Conected, desc(Grado))
    DF <- ExtinctionOrder(Network = Network, Order = Conected$ID, clust.method = clust.method,
                          IS = IS, Rewiring = Rewiring, RewiringDist = RewiringDist, verbose = verbose)
  }
  if(Method == "Ordered"){
    DF <- ExtinctionOrder(Network = Network, Order = Order, clust.method = clust.method,
                          IS = IS, Rewiring = Rewiring, RewiringDist = RewiringDist, verbose = verbose)
  }

  return(DF)
}

