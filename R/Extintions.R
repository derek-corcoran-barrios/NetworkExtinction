#' Extinctions analysis for trophic networks
#'
#' The SimulateExtinctions function, can be used to test how the order of species
#' extinctions might affect the stability of the network by comparing  The extintion history
#' and checking for secondary extinctions.
#'
#' @param Method a character with the options Mostconnected, Oredered, or Random
#' @param Network a network representation as a an adyacency matrix, edgelist,
#' or a network object
#' @param Order this should be NULL, unless using the Ordered method, in that case
#' it should be a vector with the order of extinctions by ID
#' @param clust.method a character with the options cluster_edge_betweenness, cluster_spinglass,
#' cluster_label_prop or cluster_infomap, defaults to cluster_infomap
#' @return exports data frame with the characteristics of the network after every
#' extintion. The resulting data frame contains 11 columns that incorporate the
#' topological index, the secondary extinctions, predation release, and total extinctions of the network
#' in each primary extinction.
#'
#' @details When method is Mostconnected, it takes a network and it calculates wich node is the most connected
#' of the network, using total degree. Then remove the most connected node,
#' and calculates the the topological indexes of the network and the number of
#' secundary extintions (how many species have indegree 0, without considered
#' primary producers). After that, remove the nodes that were secondarily extinct
#' in the previous step and recalculate which is the new most connected
#' node and so on, until the number of links in the network is zero.
#'
#' When method is Ordered, it takes a network, and extinguishes nodes using a custom order,
#' then it calculates the secondary extinctions and plots the accumulated
#' secondary extinctions.
#'
#' When clust.method = cluster_edge_betweenness computes the network modularity using cluster_edge_betweenness methods from igraph to detect communities
#' When clust.method = cluster_spinglass computes the network modularity using cluster_spinglass methods from igraph to detect communities, here the number of spins are equal to the nerwork size
#' When clust.method = cluster_label_prop computes the network modularity using cluster_label_prop methods from igraph to detect communities
#' When clust.method = cluster_infomap computes the network modularity using cluster_infomap methods from igraph to detect communities, here the number of nb.trials are equal to the nerwork size
#'
#' @examples
#' # Mostconnected example
#' data("net")
#' SimulateExtinctions(Network = net, Method = "Mostconnected", clust.method = "cluster_infomap")
#' #first Ordered example
#' data("net")
#' SimulateExtinctions(Network = net, Order = c(1,2,3,4,5,6,7,8,9,10), Method = "Ordered" , clust.method = "cluster_infomap")
#' #Second Ordered example
#' data("net")
#' SimulateExtinctions(Network = net, Order = c(2,8,9), Method = "Ordered", clust.method = "cluster_infomap")

#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author M. Isidora Ávila-Thieme <msavila@uc.cl>
#' @export


SimulateExtinctions <- function(Network, Method, Order = NULL, clust.method = "cluster_infomap"){
  Network <- .DataInit(x = Network)

  '%ni%'<- Negate('%in%')
  if(Method %ni% c("Mostconnected", "Ordered")) stop('Choose the right method. See ?SimulateExtinction.')

  if(Method == "Mostconnected"){
    DF <- .Mostconnected(Network = Network, clust.method = clust.method)
  }
  if(Method == "Ordered"){
    DF <- .ExtinctionOrder(Network = Network, Order = Order, clust.method = clust.method)
  }
  return(DF)
}

#' Extinctions analysis from most connected to less connected nodes in the network
#'
#' It takes a network and it calculates wich node is the most connected
#' of the network, using total degree. Then remove the most connected node,
#' and calculates the the topological indexes of the network and the number of
#' secundary extintions (how many species have indegree 0, without considered
#' primary producers). After that, remove the nodes that were secondarily extinct
#' in the previous step and recalculate which is the new most connected
#' node and so on, until the number of links in the network is zero.
#'
#
#' @param Network a trophic network of class network
#' @return exports data frame with the characteristics of the network after every
#' extintion. The resulting data frame contains 11 columns that incorporate the
#' topological index, the secondary extinctions, predation release, and total extinctions of the network
#' in each primary extinction.
#' @importFrom network as.matrix.network.edgelist
#' @importFrom network delete.vertices
#' @importFrom network network.density
#' @importFrom network network.edgecount
#' @importFrom network network.size
#' @importFrom sna degree
#' @importFrom stats complete.cases
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author M. Isidora Ávila-Thieme <msavila@uc.cl>
#' @seealso [NetworkExtinction::ExtinctionOrder()]
#' @export


Mostconnected <- function(Network){
  Grado <- NULL
  Network <- .DataInit(x = Network)
  edgelist <- as.matrix.network.edgelist(Network,matrix.type="edgelist") #Prey - Predator
  Conected <- data.frame(ID = 1:network.size(Network), Grado = degree(edgelist, c("total")))
  Conected <- arrange(Conected, desc(Grado))
  Conected1<- c(Conected$ID)
  indegreebasenet <- degree(Network, cmode = "indegree")
  indegreebasenetzeros <- sum(degree(Network, cmode = "indegree") == 0)
  indegreetopnetzeros <- sum(degree(Network, cmode = "outdegree") == 0)
  Producers <- (1:length(degree(Network, cmode = "indegree")))[degree(Network, cmode = "indegree") == 0]
  TopPredators <- (1:length(degree(Network, cmode = "outdegree")))[degree(Network, cmode = "outdegree") == 0]
  DF <- data.frame(Spp = rep(NA, network.size(Network)), S = rep(NA, network.size(Network)), L = rep(NA, network.size(Network)), C = rep(NA, network.size(Network)), Link_density = rep(NA, network.size(Network)),SecExt = rep(NA,network.size(Network)), Pred_release = rep(NA,network.size(Network)), Iso_nodes =rep (NA,network.size(Network)))

  Secundaryext <- c()
  Predationrel <- c()
  accExt <- c()
  totalExt <- c()
  FinalExt <- list()
  Conected3 <- c()

  ####LOOP####
  for (i in 1:network.size(Network)){
    #esta lista tiene el mismo orden que conected 1, hay que
    #volver a hacer la red y calcular el grado
    if (length(accExt)==0){
      Temp <- Network
      DF$Spp[i] <- Conected1[i]
      delete.vertices(Temp, c(DF$Spp[1:i]))
    }
    if (length(accExt)>0){
      Temp <- Network
      Temp <- delete.vertices(Temp, c(accExt))
      edgelist <- as.matrix.network.edgelist(Temp,matrix.type="edgelist")
      Conected2 <- data.frame(ID = 1:network.size(Temp), Grado = degree(edgelist, c("total")))
      Conected2 <- arrange(Conected2, desc(Grado))
      for(j in sort(accExt)){
        Conected2$ID <- ifelse(Conected2$ID < j, Conected2$ID, Conected2$ID + 1)
      }

      DF$Spp[i] <- Conected2$ID[1]
      Temp <- Network

      delete.vertices(Temp, unique(c(c(DF$Spp[1:i]),accExt)))
    }

    DF$S[i] <- network.size(Temp)
    DF$L[i] <- network.edgecount(Temp)
    DF$C[i] <- network.density(Temp)
    DF$Link_density [i] <- DF$L[i]/DF$S[i]
    SecundaryextTemp <- (1:length(degree(Temp, cmode = "indegree")))[degree(Temp, cmode = "indegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      SecundaryextTemp <- ifelse(SecundaryextTemp < j, SecundaryextTemp, SecundaryextTemp + 1)
    }
    Secundaryext <- SecundaryextTemp
    Secundaryext <- Secundaryext[!(Secundaryext %in% Producers)]
    DF$SecExt[i]<- length(Secundaryext)

    PredationrelTemp <- (1:length(degree(Temp, cmode = "outdegree")))[degree(Temp, cmode = "outdegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      PredationrelTemp <- ifelse(PredationrelTemp < j, PredationrelTemp, PredationrelTemp + 1)
    }
    Predationrel <- PredationrelTemp
    Predationrel <- Predationrel[!(Predationrel %in% TopPredators)]
    DF$Pred_release[i]<- length(Predationrel)

    DF$Iso_nodes[i] <- sum(degree(Temp) == 0)
    print(i)
    FinalExt[[i]] <-(Secundaryext)
    accExt <- append(accExt, DF$Spp[1:i])
    accExt <- unique(append(accExt,Secundaryext))

    if (DF$L[i] == 0) break
  }
  DF <- DF[complete.cases(DF),]
    DF$AccSecExt<- cumsum(DF$SecExt)
  DF$NumExt <- 1:nrow(DF)
  DF$TotalExt <- DF$AccSecExt + DF$NumExt
  class(DF) <- c("data.frame", "Mostconnected")
  .Deprecated("SimulateExtinctions")
  return(DF)
}

#' Extinctions analysis from most connected to less connected nodes in the network
#'
#' It takes a network and it calculates wich node is the most connected
#' of the network, using total degree. Then remove the most connected node,
#' and calculates the the topological indexes of the network and the number of
#' secundary extintions (how many species have indegree 0, without considered
#' primary producers). After that, remove the nodes that were secondarily extinct
#' in the previous step and recalculate which is the new most connected
#' node and so on, until the number of links in the network is zero.
#'
#
#' @param Network a trophic network of class network
#' @param clust.method a character with the options cluster_edge_betweenness, cluster_spinglass, cluster_label_prop or cluster_infomap
#' @return exports data frame with the characteristics of the network after every
#' extintion. The resulting data frame contains 11 columns that incorporate the
#' topological index, the secondary extinctions, predation release, and total extinctions of the network
#' in each primary extinction.
#' @importFrom network as.matrix.network.edgelist
#' @importFrom network delete.vertices
#' @importFrom network network.density
#' @importFrom network network.edgecount
#' @importFrom network network.size
#' @importFrom sna degree
#' @importFrom stats complete.cases
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr relocate
#' @importFrom network as.matrix.network.adjacency
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph cluster_edge_betweenness
#' @importFrom igraph cluster_spinglass
#' @importFrom igraph cluster_label_prop
#' @importFrom igraph cluster_infomap
#' @importFrom igraph modularity
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author M. Isidora Ávila-Thieme <msavila@uc.cl>
#' @seealso [NetworkExtinction::ExtinctionOrder()]

.Mostconnected <- function(Network, clust.method = "cluster_infomap"){
  Grado <- NULL
  Network <- Network
  edgelist <- as.matrix.network.edgelist(Network,matrix.type="edgelist") #Prey - Predator
  Conected <- data.frame(ID = 1:network.size(Network), Grado = degree(edgelist, c("total")))
  Conected <- arrange(Conected, desc(Grado))
  Conected1<- c(Conected$ID)
  indegreebasenet <- degree(Network, cmode = "indegree")
  indegreebasenetzeros <- sum(degree(Network, cmode = "indegree") == 0)
  indegreetopnetzeros <- sum(degree(Network, cmode = "outdegree") == 0)
  Producers <- (1:length(degree(Network, cmode = "indegree")))[degree(Network, cmode = "indegree") == 0]
  TopPredators <- (1:length(degree(Network, cmode = "outdegree")))[degree(Network, cmode = "outdegree") == 0]
  DF <- data.frame(Spp = rep(NA, network.size(Network)), S = rep(NA, network.size(Network)), L = rep(NA, network.size(Network)), C = rep(NA, network.size(Network)), Link_density = rep(NA, network.size(Network)),SecExt = rep(NA,network.size(Network)), Pred_release = rep(NA,network.size(Network)), Iso_nodes =rep (NA,network.size(Network)))

  Secundaryext <- c()
  Predationrel <- c()
  accExt <- c()
  totalExt <- c()
  FinalExt <- list()
  Conected3 <- c()

  ####LOOP####
  for (i in 1:network.size(Network)){
    #esta lista tiene el mismo orden que conected 1, hay que
    #volver a hacer la red y calcular el grado
    if (length(accExt)==0){
      Temp <- Network
      DF$Spp[i] <- Conected1[i]
      delete.vertices(Temp, c(DF$Spp[1:i]))
    }
    if (length(accExt)>0){
      Temp <- Network
      Temp <- delete.vertices(Temp, c(accExt))
      edgelist <- as.matrix.network.edgelist(Temp,matrix.type="edgelist")
      Conected2 <- data.frame(ID = 1:network.size(Temp), Grado = degree(edgelist, c("total")))
      Conected2 <- arrange(Conected2, desc(Grado))
      for(j in sort(accExt)){
        Conected2$ID <- ifelse(Conected2$ID < j, Conected2$ID, Conected2$ID + 1)
      }

      DF$Spp[i] <- Conected2$ID[1]
      Temp <- Network

      delete.vertices(Temp, unique(c(c(DF$Spp[1:i]),accExt)))
    }

    DF$S[i] <- network.size(Temp)
    DF$L[i] <- network.edgecount(Temp)
    DF$C[i] <- network.density(Temp)
    DF$Link_density [i] <- DF$L[i]/DF$S[i]

    Networkclass = class(Temp)

    if (Networkclass[1] == "matrix"){
      netgraph = graph_from_adjacency_matrix(Temp, mode = "directed", weighted = NULL)
    }

    if (Networkclass[1] == "network"){
      net = as.matrix.network.adjacency(Temp)
      netgraph = graph_from_adjacency_matrix(net, mode = "directed", weighted = NULL)
    }

    if (clust.method == "cluster_edge_betweenness"){
      Membership = suppressWarnings(cluster_edge_betweenness(netgraph, weights=NULL, directed = FALSE, edge.betweenness = TRUE,
                                            merges = TRUE, bridges = TRUE, modularity = TRUE, membership = TRUE))
    } else if (clust.method == "cluster_spinglass"){
      spins = network.size(Temp)
      Membership = suppressWarnings(cluster_spinglass(netgraph, spins=spins)) #spins could be the Richness
    }else if (clust.method == "cluster_label_prop"){
      Membership = suppressWarnings(cluster_label_prop(netgraph, weights = NULL, initial = NULL,
                                      fixed = NULL))
    }else if (clust.method == "cluster_infomap"){
      nb.trials = network.size(Temp)
      Membership = suppressWarnings(cluster_infomap(netgraph, e.weights = NULL, v.weights = NULL,
                                   nb.trials = nb.trials, modularity = TRUE))

    } else stop('Select a valid method for clustering. ?SimulateExtinction')

    DF$Modularity[i] <- suppressMessages(modularity(Membership))

    SecundaryextTemp <- (1:length(degree(Temp, cmode = "indegree")))[degree(Temp, cmode = "indegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      SecundaryextTemp <- ifelse(SecundaryextTemp < j, SecundaryextTemp, SecundaryextTemp + 1)
    }
    Secundaryext <- SecundaryextTemp
    Secundaryext <- Secundaryext[!(Secundaryext %in% Producers)]
    DF$SecExt[i]<- length(Secundaryext)

    PredationrelTemp <- (1:length(degree(Temp, cmode = "outdegree")))[degree(Temp, cmode = "outdegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      PredationrelTemp <- ifelse(PredationrelTemp < j, PredationrelTemp, PredationrelTemp + 1)
    }
    Predationrel <- PredationrelTemp
    Predationrel <- Predationrel[!(Predationrel %in% TopPredators)]
    DF$Pred_release[i]<- length(Predationrel)

    DF$Iso_nodes[i] <- sum(degree(Temp) == 0)
    print(i)
    FinalExt[[i]] <-(Secundaryext)
    accExt <- append(accExt, DF$Spp[1:i])
    accExt <- unique(append(accExt,Secundaryext))

    if (DF$L[i] == 0) break
  }
  DF <- DF[complete.cases(DF),]
  DF$AccSecExt<- cumsum(DF$SecExt)
  DF$NumExt <- 1:nrow(DF)
  DF$TotalExt <- DF$AccSecExt + DF$NumExt
  DF <- relocate(DF, Modularity, .after = Link_density)
  class(DF) <- c("data.frame", "SimulateExt")
  return(DF)
}

#' Extinctions analysis from custom order
#'
#' It takes a network, and extinguishes nodes using a custom order,
#' then it calculates the secondary extinctions and plots the accumulated
#' secondary extinctions.
#'
#' @param Network a network of class network
#' @param Order Vector with the order of extinctions by ID
#' @return exports data frame with the characteristics of the network after every
#' extintion, and a graph with the mean and 95% interval
#'
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom network as.matrix.network.edgelist
#' @importFrom network delete.vertices
#' @importFrom network network.edgecount
#' @importFrom network network.size
#' @importFrom sna degree
#' @importFrom stats complete.cases
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author M. Isidora Ávila-Thieme <msavila@uc.cl>
#' @export

ExtinctionOrder <- function(Network, Order){
  Grado <- NULL
  Network <- .DataInit(x = Network)
  edgelist <- as.matrix.network.edgelist(Network,matrix.type="edgelist") #Prey - Predator
  Conected <- data.frame(ID = 1:network.size(Network), Grado = degree(edgelist, c("total")))

  Conected1<-  Order

  indegreebasenet <- degree(Network, cmode = "indegree")

  indegreebasenetzeros <- sum(degree(Network, cmode = "indegree") == 0)
  indegreetopnetzeros <- sum(degree(Network, cmode = "outdegree") == 0)

  Producers <- (1:length(degree(Network, cmode = "indegree")))[degree(Network, cmode = "indegree") == 0]
  TopPredators <- (1:length(degree(Network, cmode = "outdegree")))[degree(Network, cmode = "outdegree") == 0]

  DF <- data.frame(Spp = rep(NA, network.size(Network)), S = rep(NA, network.size(Network)), L = rep(NA, network.size(Network)), C = rep(NA, network.size(Network)),  SecExt = rep(NA,network.size(Network)), Pred_release = rep(NA,network.size(Network)))

  Secundaryext <- c()
  Predationrel <- c()
  accExt <- c()
  totalExt <- c()
  FinalExt <- list()
  Conected3 <- c()

  for (i in 1:length(Order)){

    if (length(accExt)==0){
      Temp <- Network
      DF$Spp[i] <- Conected1[i]
      delete.vertices(Temp, c(DF$Spp[1:i]))
    }
    if (length(accExt)>0){
      Temp <- Network
      Temp <- delete.vertices(Temp, c(accExt))
      edgelist <- as.matrix.network.edgelist(Temp,matrix.type="edgelist")
      Conected2 <- data.frame(ID =1:network.size(Temp), Grado = degree(edgelist, c("total")))
      for(j in sort(accExt)){
        Conected2$ID <- ifelse(Conected2$ID < j, Conected2$ID, Conected2$ID + 1)
      }

      DF$Spp[i] <- Conected1[i]
      Temp <- Network

      delete.vertices(Temp, unique(c(c(DF$Spp[1:i]),accExt)))
    }

    DF$S[i] <- network.size(Temp)
    DF$L[i] <- network.edgecount(Temp)
    DF$C[i] <- network.density(Temp)

    SecundaryextTemp <- (1:length(degree(Temp, cmode = "indegree")))[degree(Temp, cmode = "indegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      SecundaryextTemp <- ifelse(SecundaryextTemp < j, SecundaryextTemp, SecundaryextTemp + 1)
    }
    Secundaryext <- SecundaryextTemp
    Secundaryext <- Secundaryext[!(Secundaryext %in% Producers)]
    DF$SecExt[i]<- length(Secundaryext)

    PredationrelTemp <- (1:length(degree(Temp, cmode = "outdegree")))[degree(Temp, cmode = "outdegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      PredationrelTemp <- ifelse(PredationrelTemp < j, PredationrelTemp, PredationrelTemp + 1)
    }
    Predationrel <- PredationrelTemp
    Predationrel <- Predationrel[!(Predationrel %in% TopPredators)]
    DF$Pred_release[i]<- length(Predationrel)

    message(i)
    FinalExt[[i]] <-(Secundaryext)
    accExt <- append(accExt, DF$Spp[1:i])
    accExt <- unique(append(accExt,Secundaryext))

    if (DF$L[i] == 0) break
  }
  DF <- DF[complete.cases(DF),]
  DF$AccSecExt <- cumsum(DF$SecExt)
  DF$NumExt <- 1:nrow(DF)
  DF$TotalExt <- DF$AccSecExt + DF$NumExt
  class(DF) <- c("data.frame", "ExtinctionOrder")
  .Deprecated("SimulateExtinctions")
  return(DF)
}

#' Extinctions analysis from custom order
#'
#' It takes a network, and extinguishes nodes using a custom order,
#' then it calculates the secondary extinctions and plots the accumulated
#' secondary extinctions.
#'
#' @param Network a network of class network
#' @param Order Vector with the order of extinctions by ID
#' @return exports data frame with the characteristics of the network after every
#' extintion, and a graph with the mean and 95% interval
#'
#' @importFrom dplyr relocate
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom network as.matrix.network.edgelist
#' @importFrom network delete.vertices
#' @importFrom network network.edgecount
#' @importFrom network network.size
#' @importFrom sna degree
#' @importFrom stats complete.cases
#' @importFrom network as.matrix.network.adjacency
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph cluster_edge_betweenness
#' @importFrom igraph cluster_spinglass
#' @importFrom igraph cluster_label_prop
#' @importFrom igraph cluster_infomap
#' @importFrom igraph modularity
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author M. Isidora Ávila-Thieme <msavila@uc.cl>

.ExtinctionOrder <- function(Network, Order, clust.method = "cluster_infomap"){
  Grado <- NULL
  edgelist <- as.matrix.network.edgelist(Network,matrix.type="edgelist") #Prey - Predator
  Conected <- data.frame(ID = 1:network.size(Network), Grado = degree(edgelist, c("total")))

  Conected1<-  Order

  indegreebasenet <- degree(Network, cmode = "indegree")

  indegreebasenetzeros <- sum(degree(Network, cmode = "indegree") == 0)
  indegreetopnetzeros <- sum(degree(Network, cmode = "outdegree") == 0)

  Producers <- (1:length(degree(Network, cmode = "indegree")))[degree(Network, cmode = "indegree") == 0]
  TopPredators <- (1:length(degree(Network, cmode = "outdegree")))[degree(Network, cmode = "outdegree") == 0]

  DF <- data.frame(Spp = rep(NA, network.size(Network)), S = rep(NA, network.size(Network)), L = rep(NA, network.size(Network)), C = rep(NA, network.size(Network)), Link_density = rep(NA, network.size(Network)),SecExt = rep(NA,network.size(Network)), Pred_release = rep(NA,network.size(Network)), Iso_nodes =rep (NA,network.size(Network)))

  Secundaryext <- c()
  Predationrel <- c()
  accExt <- c()
  totalExt <- c()
  FinalExt <- list()
  Conected3 <- c()

  for (i in 1:length(Order)){

    if (length(accExt)==0){
      Temp <- Network
      DF$Spp[i] <- Conected1[i]
      delete.vertices(Temp, c(DF$Spp[1:i]))
    }
    if (length(accExt)>0){
      Temp <- Network
      Temp <- delete.vertices(Temp, c(accExt))
      edgelist <- as.matrix.network.edgelist(Temp,matrix.type="edgelist")
      Conected2 <- data.frame(ID =1:network.size(Temp), Grado = degree(edgelist, c("total")))
      for(j in sort(accExt)){
        Conected2$ID <- ifelse(Conected2$ID < j, Conected2$ID, Conected2$ID + 1)
      }

      DF$Spp[i] <- Conected1[i]
      Temp <- Network

      delete.vertices(Temp, unique(c(c(DF$Spp[1:i]),accExt)))
    }

    DF$S[i] <- network.size(Temp)
    DF$L[i] <- network.edgecount(Temp)
    DF$C[i] <- network.density(Temp)
    DF$Link_density [i] <- DF$L[i]/DF$S[i]

    Networkclass = class(Temp)

    if (Networkclass[1] == "matrix"){
      netgraph = graph_from_adjacency_matrix(Temp, mode = "directed", weighted = NULL)
    }

    if (Networkclass[1] == "network"){
      net = as.matrix.network.adjacency(Temp)
      netgraph = suppressMessages(graph_from_adjacency_matrix(net, mode = "directed", weighted = NULL))
    }

    if (clust.method == "cluster_edge_betweenness"){
      Membership = suppressWarnings(cluster_edge_betweenness(netgraph, weights=NULL, directed = TRUE, edge.betweenness = TRUE,
                                            merges = TRUE, bridges = TRUE, modularity = TRUE, membership = TRUE))
    } else if (clust.method == "cluster_spinglass"){
      spins = 107#network.size(Temp)
      Membership = suppressWarnings(cluster_spinglass(netgraph, spins=spins)) #spins could be the Richness
    }else if (clust.method == "cluster_label_prop"){
      Membership = suppressWarnings(cluster_label_prop(netgraph, weights = NULL, initial = NULL,
                                      fixed = NULL))
    }else if (clust.method == "cluster_infomap"){
      nb.trials = 107#network.size(Temp)
      Membership = suppressWarnings(cluster_infomap(netgraph, e.weights = NULL, v.weights = NULL,
                                   nb.trials = nb.trials, modularity = TRUE))

    } else if (clust.method == "none"){
      Membership = NA
    }else stop('Select a valid method for clustering. ?SimulateExtinction')
    if(is.na(Membership)){
      DF$Modularity[i] <- NA
    }else{
      DF$Modularity[i] <- suppressWarnings(modularity(Membership))
    }

    SecundaryextTemp <- (1:length(degree(Temp, cmode = "indegree")))[degree(Temp, cmode = "indegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      SecundaryextTemp <- ifelse(SecundaryextTemp < j, SecundaryextTemp, SecundaryextTemp + 1)
    }
    Secundaryext <- SecundaryextTemp
    Secundaryext <- Secundaryext[!(Secundaryext %in% Producers)]
    DF$SecExt[i]<- length(Secundaryext)

    PredationrelTemp <- (1:length(degree(Temp, cmode = "outdegree")))[degree(Temp, cmode = "outdegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      PredationrelTemp <- ifelse(PredationrelTemp < j, PredationrelTemp, PredationrelTemp + 1)
    }
    Predationrel <- PredationrelTemp
    Predationrel <- Predationrel[!(Predationrel %in% TopPredators)]
    DF$Pred_release[i]<- length(Predationrel)
    DF$Iso_nodes[i] <- sum(degree(Temp) == 0)

    message(i)
    FinalExt[[i]] <-(Secundaryext)
    accExt <- append(accExt, DF$Spp[1:i])
    accExt <- unique(append(accExt,Secundaryext))

    if (DF$L[i] == 0) break
  }
  DF <- DF[complete.cases(DF),]
  DF$AccSecExt <- cumsum(DF$SecExt)
  DF$NumExt <- 1:nrow(DF)
  DF$TotalExt <- DF$AccSecExt + DF$NumExt
  DF <- relocate(DF, Modularity, .after = Link_density)
  class(DF) <- c("data.frame", "ExtinctionOrder")
  return(DF)
}

#' Random extinction
#'
#' Generates a null model by generating random extinction histories and calculating
#' the mean and standard deviation of the accumulated secondary extinctions developed
#' by making n random extinction histories
#'
#' @param Network a trophic network of class network
#' @param nsim number of simulations
#' @param parallel if TRUE, it will use parallel procesing, if FALSE (default) it will run
#' sequentially
#' @param ncores number of cores to use if using parallel procesing
#' @param Record logical, if TRUE, records every simulation and you can read the
#' raw results in the object FullSims
#' @param plot logical if true, will add a graph to the results
#' @return exports data frame with the characteristics of the network after every
#' extintion, and a graph with the mean and 95% interval
#' @examples
#' #first example
#' data("More_Connected")
#' RandomExtinctions(Network = More_Connected, nsim = 20)
#'
#' # Using parallel procesing
#' ## Detect your number of cores divide by 2
#' \dontrun{
#' cores <- ceiling(parallel::detectCores()/2)
#'
#' RandomExtinctions(Network = More_Connected, nsim = 20, parallel = TRUE, ncores = cores)
#' }
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr summarise
#' @importFrom foreach `%dopar%`
#' @importFrom foreach foreach
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom magrittr "%>%"
#' @importFrom network network.size
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom scales muted
#' @importFrom stats sd
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author M. Isidora Ávila-Thieme <msavila@uc.cl>
#' @export

RandomExtinctions <- function(Network, nsim = 10, parallel = FALSE, ncores, Record = F, plot = F){
  NumExt <- sd <- AccSecExt <- AccSecExt_95CI <- AccSecExt_mean <- Lower <- NULL
  network <- .DataInit(x = Network)

  if(parallel){
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    sims <- foreach(i=1:nsim, .packages = "NetworkExtinction")%dopar%{
      sims <- try(.ExtinctionOrder(Network = network, Order = sample(1:network.size(network)), clust.method = "none"), silent = T)
      try({sims$simulation <- i}, silent = T)
      sims
    }
    stopCluster(cl)
  }else{
    sims <- list()
    for(i in 1:nsim){
      sims[[i]] <- try(.ExtinctionOrder(Network = network, Order = sample(1:network.size(network))), silent = T)
      try({sims[[i]]$simulation <- i}, silent = T)
      message(paste("Simulation", i, "of", nsim, "ready"))
    }
  }

  cond <- sapply(sims, function(x) "data.frame" %in% class(x))
  cond <- c(1:length(cond))[cond]
  sims <- sims[cond]
  sims <- do.call(rbind, sims)
  if(Record == TRUE){
    FullSims <- sims
  }

  sims <- sims %>% group_by(NumExt) %>% summarise(AccSecExt_95CI = 1.96*sd(AccSecExt), AccSecExt_mean = mean(AccSecExt)) %>% mutate(Upper = AccSecExt_mean + AccSecExt_95CI, Lower = AccSecExt_mean - AccSecExt_95CI, Lower = ifelse(Lower < 0, 0, Lower))
  if(plot == T){
    g <- ggplot(sims, aes_string(x = "NumExt", y = "AccSecExt_mean")) + geom_ribbon(aes_string(ymin = "Lower", ymax = "Upper"), fill = muted("red")) + geom_line() + ylab("Acc. Secondary extinctions") + xlab("Primary extinctions") + theme_bw()
    g
  }

  if(Record == T & plot == T){
    return(list(sims = sims, graph = g, FullSims = FullSims))
  }else if(Record == F & plot == T){
    return(list(sims = sims, graph = g))
  }else if(Record == F & plot == F){
    return(sims)
  }else if(Record == T & plot == F){
    return(list(sims = sims, FullSims = FullSims))
  }
}

#' Comparison of Null hypothesis with other extinction histories
#'
#' It compares an object genrated either by the Mostconnected or ExtinctionOrder functions
#' with a null hypothesis generated by the RandomExtinctions function it is important that
#' RandomExtinctions is in plot = T.
#'
#' @param Nullmodel an object generated by the RandomExtinctions
#' @param Hypothesis Extinction history generated by the Mostconnected or ExtinctionOrder
#' fuction
#' @return a plot comparing the expected value of secondary extinctions originated at random
#' with the observed extinction history.
#'
#' @examples
#' data("net")
#' History <- SimulateExtinctions(Network = net, Method = "Mostconnected")
#'
#' NullHyp <- RandomExtinctions(Network = net, nsim = 100, plot = TRUE)
#'
#' CompareExtinctions(Nullmodel = NullHyp, Hypothesis = History)
#' @importFrom broom tidy
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 theme_bw
#' @importFrom scales muted
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author M. Isidora Ávila-Thieme <msavila@uc.cl>
#' @export

CompareExtinctions <- function(Nullmodel, Hypothesis){
  if(class(Hypothesis)[2] == "SimulateExt"){
    NumExt <- sd <- AccSecExt <- AccSecExt_mean <-NULL
    if(class(Nullmodel)[1] == "list"){
      g <- Nullmodel$graph + geom_line(aes(color = "blue"))
      g <- g + geom_point(data = Hypothesis, aes(y = AccSecExt), color = "black") + geom_line(data = Hypothesis, aes(y = AccSecExt, color = "black")) + scale_color_manual(name = "Comparison",values =c("black", "blue"), label = c("Observed","Null hypothesis"))
    } else {
      g <- ggplot(Nullmodel, aes(x = NumExt, y = AccSecExt_mean)) + geom_ribbon(aes_string(ymin = "Lower", ymax = "Upper"), fill = muted("red")) + geom_line() + ylab("Acc. Secondary extinctions") + xlab("Primary extinctions") + theme_bw()
      g <- g + geom_point(data = Hypothesis, aes(y = AccSecExt), color = "black") + geom_line(data = Hypothesis, aes(y = AccSecExt, color = "black")) + scale_color_manual(name = "Comparison",values =c("black", "blue"), label = c("Observed","Null hypothesis"))
      g
    }

    g

    return(g)
  }
  if(class(Hypothesis)[2] %in% c("Mostconnected", "ExtinctionOrder")){
  NumExt <- sd <- AccSecExt <- AccSecExt_mean <-NULL
  g <- Nullmodel$graph + geom_line(aes(color = "blue"))
  g <- g + geom_point(data = Hypothesis, aes(y = AccSecExt), color = "black") + geom_line(data = Hypothesis, aes(y = AccSecExt, color = "black")) + scale_color_manual(name = "Comparison", values =c("black", "blue"), label = c("Observed","Null hypothesis"))
  g
  return(g)
  }
  else{
    message("Hipothesis not of class Mostconnected or ExtinctionOrder")
  }
}
