#' Extinctions analysis for ecological networks
#'
#' The SimulateExtinctions function, can be used to test how the order of species
#' extinctions, species-dependency on existing interaction strength, and rewiring potential might affect the stability of the network by comparing  The extinction history
#' and checking for secondary extinctions.
#'
#' @param Network a network representation as a an adjacency matrix, edgelist,
#' or a network object
#' @param Method a character with the options Mostconnected, Leastconnected and Ordered
#' @param Order a numeric vector indexing order of primary extinctions. For Method = Mostconnected or Leastconnected Order must be NULL. If Order is not NULL, Method is internally forced to be Ordered.
#' @param NetworkType a character with the options Trophic and Mutualistic - is used to calculate secondary extinctions.
#' @param clust.method a character with the options cluster_edge_betweenness,
#' cluster_label_prop or cluster_infomap, defaults to cluster_infomap
#' @param IS either numeric or a named vector of numerics. Identifies the threshold of relative interaction strength which species require to not be considered secondarily extinct (i.e. IS = 0.3 leads to removal of all nodes which lose 70 percent of their interaction strength in the Network argument). If a named vector, names must correspond to vertex names in Network argument.
#' @param Rewiring either a function or a named vector of functions. Signifies how rewiring probabilities are calculated from the RewiringDist argument. If FALSE, no rewiring is carried out.
#' @param RewiringDist a numeric matrix of NxN dimension (N... number of nodes in Network). Contains, for example, phylogenetic or functional trait distances between nodes in Network which are used by the Rewiring argument to calculate rewiring probabilities. If Rewiring == function(x){x}, this matrix is expected to contain probabilities of a connection being present between species-pairs.
#' @param RewiringProb a numeric which identifies the threshold at which to assume rewiring potential is met.
#' @param verbose Logical. Whether to report on function progress or not.
#' @param forceFULL Logical. Whether to continue removal of nodes after initial order has been depleted. This will force the simulations to execute extinctions and check for secondary extinctions/new links until the network does not change anylonger.
#' @return exports list containing a data frame with the characteristics of the network after every extinction and a network object containing the final network. The resulting data frame contains 11 columns that incorporate the topological index, the secondary extinctions, predation release, and total extinctions of the network in each primary extinction.
#' @details When method is Mostconnected, the function takes the network and calculates which node is the most connected of the network, using total degree. Then remove the most connected node, and calculates the the topological indexes of the network and the number of secondary extinctions. This process is repeated until the entire network has gone extinct. When method is Leastconnected, this process prioritises nodes of lowest degree.
#'
#' When method is Ordered, it takes a network, and extinguishes nodes using a custom order, then it calculates the secondary extinctions and plots the accumulated secondary extinctions.
#'
#' When NetworkType = Trophic, secondary extinctions only occur for any predator, but not producers. If NetworkType = Mutualistic, secondary extinctions occur for all species in the network.
#'
#' When clust.method = cluster_edge_betweenness computes the network modularity using cluster_edge_betweenness methods from igraph to detect communities
#' When clust.method = cluster_label_prop computes the network modularity using cluster_label_prop methods from igraph to detect communities
#' When clust.method = cluster_infomap computes the network modularity using cluster_infomap methods from igraph to detect communities, here the number of nb.trials are equal to the network size
#' @examples
#' # Mostconnected example
#' data("net")
#' SimulateExtinctions(Network = net, Method = "Mostconnected",
#' clust.method = "cluster_infomap")
#'
#' #first Ordered example
#' data("net")
#' SimulateExtinctions(Network = net, Order = c(1,2,3,4,5,6,7,8,9,10),
#' Method = "Ordered" , clust.method = "cluster_infomap")
#'
#'  #Second Ordered example
#' data("net")
#' SimulateExtinctions(Network = net, Order = c(2,8,9),
#' Method = "Ordered", clust.method = "cluster_infomap")
#'
#' #Network-Dependency Example
#' data("net")
#' SimulateExtinctions(Network = net, Order = c(2,8), IS = 0.3,
#' Method = "Ordered", clust.method = "cluster_infomap")
#'
#'  #Rewiring
#' data("net")
#' data(dist)
#' SimulateExtinctions(Network = net, Order = c(2,8), IS = 0.3,
#' # assuming an exponential decline in rewiring potential
#' # as values in RewiringDist increase
#' Rewiring = function(x){1-pexp(x, rate = 1/0.5)},
#' RewiringDist = dist, # distance matrix
#' RewiringProb = 0.2, # low threshold for rewiring potential
#' Method = "Ordered", clust.method = "cluster_infomap")
#'
#' #Rewiring, assuming dist contains probabilities
#' #' data("net")
#' data(dist)
#' SimulateExtinctions(Network = net, Order = c(2,8), IS = 0.3,
#' Rewiring = function(x){x}, # no changes to the RewiringDist object means
#' RewiringDist = dist, RewiringProb = 0.2,
#' Method = "Ordered", clust.method = "cluster_infomap")
#'
#' ## mutualistic network example
#' data(mutual)
#' # tallying of first-order secondary extinctions only
#' SimulateExtinctions(Network = mutual, Order = 3, NetworkType = "Mutualistic",
#' IS = 1, forceFULL = FALSE)
#' # tallying of all secondary extinctions until network contains no
#' #more potential secondary extinctions
#' SimulateExtinctions(Network = mutual, Order = 3, NetworkType = "Mutualistic",
#' IS = 1, forceFULL = TRUE)
#' @importFrom dplyr desc
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author M. Isidora Ávila-Thieme <msavila@uc.cl>
#' @author Erik Kusch <erik.kusch@bio.au.dk>
#' @export
SimulateExtinctions <- function(Network, Method, Order = NULL,
                                NetworkType = "Trophic", clust.method = "cluster_infomap",
                                IS = 0,
                                Rewiring = FALSE, RewiringDist, RewiringProb = 0.5,
                                verbose = TRUE, forceFULL = FALSE){
  Network <- .DataInit(x = Network)
  if(!NetworkType %in% c("Trophic", "Mutualistic")){stop("Please specify NetworkType as either 'Trophic' or 'Mutualistic'")}

  if(!is.null(Order)){Method <- "Ordered"}

  '%ni%'<- Negate('%in%')
  if(Method %ni% c("Mostconnected", "Leastconnected", "Ordered")) stop('Choose the right method. See ?SimulateExtinction.')

  edgelist <- network::as.matrix.network.edgelist(Network,matrix.type="edgelist") #Prey - Predator
  if(Method == "Mostconnected" | Method == "Leastconnected"){
    # if(NetworkType == "Trophic"){
    #   Conected <- as.numeric(names(sort(table(edgelist[,1]), decreasing = TRUE)))
    # }else{
    Grado <- NULL
    Conected <- data.frame(ID = 1:network::network.size(Network), Grado = sna::degree(edgelist, c("total")))
    if(Method == "Mostconnected"){
      Conected <- dplyr::arrange(Conected, desc(Grado))$ID
      RecalcConnect = 1
    }
    if(Method == "Leastconnected"){
      Conected <- dplyr::arrange(Conected[Conected$Grado != 0,], Grado)$ID
      RecalcConnect = -1
    }

    # }
    DF <- ExtinctionOrder(Network = Network, Order = Conected, clust.method = clust.method,
                          IS = IS, Rewiring = Rewiring, RewiringDist = RewiringDist,
                          verbose = verbose, RewiringProb = RewiringProb, NetworkType = NetworkType,
                          RecalcConnect = RecalcConnect, forceFULL = forceFULL)
  }
  if(Method == "Ordered"){
    DF <- ExtinctionOrder(Network = Network, Order = Order, clust.method = clust.method,
                          IS = IS, Rewiring = Rewiring, RewiringDist = RewiringDist,
                          verbose = verbose, RewiringProb = RewiringProb, NetworkType = NetworkType, forceFULL = forceFULL)
  }

  return(DF)
}

#' Extinctions analysis from custom order
#'
#' This function takes a network and eliminates nodes using a custom order. Subsequently, secondary extinctions are tallied up. Secondary extinction severity can be targeted by manipulating the node-dependency on network edges (IS) and node-rewiring potential upon loss of links (Rewiring).
#'
#' @param Network a network representation as a an adjacency matrix, edgelist, or a network object
#' @param Order a numeric vector indexing order of primary extinctions. For Method = Mostconnected Order must be NULL. If Order is not NULL, Method is internally forced to be Ordered.
#' @param NetworkType a character with the options Trophic and Mutualistic - is used to calculate secondary extinctions.
#' @param clust.method a character with the options cluster_edge_betweenness,
#' cluster_label_prop or cluster_infomap, defaults to cluster_infomap
#' @param IS either numeric or a named vector of numerics. Identifies the threshold of relative interaction strength which species require to not be considered secondarily extinct (i.e. IS = 0.3 leads to removal of all nodes which lose 70percent of their interaction strength in the Network argument). If a named vector, names must correspond to vertex names in Network argument.
#' @param Rewiring either a function or a named vector of functions. Signifies how rewiring probabilities are calculated from the RewiringDist argument. If FALSE, no rewiring is carried out.
#' @param RewiringDist a numeric matrix of NxN dimension (N... number of nodes in Network). Contains, for example, phylogenetic or functional trait distances between nodes in Network which are used by the Rewiring argument to calculate rewiring probabilities. If Rewiring == function(x){x}, this matrix is expected to contain probabilities of a connection being present between species-pairs.
#' @param RewiringProb a numeric which identifies the threshold at which to assume rewiring potential is met.
#' @param verbose Logical. Whether to report on function progress or not.
#' @param RecalcConnect Logical or Numeric. Whether to recalculate connectedness of each node following each round of extinction simulation and subsequently update extinction order with newly mostconnected nodes.
#' @param forceFULL Logical. Whether to continue removal of nodes after initial order has been depleted. This will force the simulations to execute extinctions and check for secondary extinctions/new links until the network does not change anylonger.
#' @return exports list containing a data frame with the characteristics of the network after every extinction and a network object containing the final network. The resulting data frame contains 11 columns that incorporate the topological index, the secondary extinctions, predation release, and total extinctions of the network in each primary extinction.
#' @details When NetworkType = Trophic, secondary extinctions only occur for any predator, but not producers. If NetworkType = Mutualistic, secondary extinctions occur for all species in the network.
#'
#' When clust.method = cluster_edge_betweenness computes the network modularity using cluster_edge_betweenness methods from igraph to detect communities
#' When clust.method = cluster_label_prop computes the network modularity using cluster_label_prop methods from igraph to detect communities
#' When clust.method = cluster_infomap computes the network modularity using cluster_infomap methods from igraph to detect communities, here the number of nb.trials are equal to the network size
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
#' @importFrom network network.density
#' @importFrom network get.vertex.attribute
#' @importFrom network get.edge.attribute
#' @importFrom igraph as.undirected
#' @importFrom igraph E
#' @importFrom sna degree
#' @importFrom stats complete.cases
#' @importFrom network as.matrix.network.adjacency
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph cluster_edge_betweenness
#' @importFrom igraph cluster_label_prop
#' @importFrom igraph cluster_infomap
#' @importFrom igraph modularity
#' @importFrom stats na.omit
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom dplyr desc
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author M. Isidora Ávila-Thieme <msavila@uc.cl>
#' @author Erik Kusch <erik.kusch@bio.au.dk>
#' @export
ExtinctionOrder <- function(Network, Order, NetworkType = "Trophic", clust.method = "cluster_infomap",
                            IS = 0,
                            Rewiring = FALSE, RewiringDist, RewiringProb = 0.5,
                            verbose = TRUE,
                            RecalcConnect = FALSE, forceFULL = FALSE
){
  if(!NetworkType %in% c("Trophic", "Mutualistic")){stop("Please specify NetworkType as either 'Trophic' or 'Mutualistic'")}
  # Setting up Objects for function run ++++++++++ ++++++++++ ++++++++++ ++++++++++
  Link_density <- Modularity <- Grado <- NULL
  Network <- .DataInit(x = Network)
  edgelist <- network::as.matrix.network.edgelist(Network,matrix.type="edgelist") #Prey - Predator
  Conected <- data.frame(ID = 1:network::network.size(Network), Grado = sna::degree(edgelist, c("total")))
  Conected <- dplyr::arrange(Conected, desc(Grado))
  Conected1 <- Order

  ## Interaction Strength Loss Preparations ++++++++++ ++++++++++
  if(length(IS )== 1){ # when the same dependency is to be applied for all species
    IS <- rep(IS, network::network.size(Network)) # repeat IS argument for all species
    names(IS) <- network::get.vertex.attribute(Network, "vertex.names") # assign species names to IS arguments
  }else{
    if(sum(!(network::get.vertex.attribute(Network, "vertex.names") %in% names(IS))) != 0){stop("Please ensure that the names of the nodes in your Network are matched by names of the elements in your IS argument vector.")}
  }
  ## Rewiring Preparations ++++++++++ ++++++++++
  if(!isFALSE(Rewiring)){ # if rewiring has been specified
    if(!exists("RewiringDist")){stop("To execute rewiring simulations, you need to specify the RewiringDist argument as well as the Rewiring argument.")}
    diag(RewiringDist)<- NA # set diagonal to NA as one can't rewire to oneself
    if(length(Rewiring) == 1){ # when the same Rewiring happen for all species
      fun <- deparse1(Rewiring) # turn function into string
      Rewiring <- rep(fun, network::network.size(Network)) # repeat function for each species
      names(Rewiring) <- network::get.vertex.attribute(Network, "vertex.names") # assign names of species in network to rewiring functions
    }
    if(sum(!(network::get.vertex.attribute(Network, "vertex.names") %in% names(Rewiring))) != 0){stop("Please ensure that the names of the nodes in your Network are matched by names of the elements in your Rewiring argument vector.")}
  }

  # Base net calculations ++++++++++ ++++++++++ ++++++++++ ++++++++++
  ## base interaction strengths per node ++++++++++ ++++++++++
  options(warn=-1) # turn warning off
  Weight_mat <- net <- network::as.matrix.network.adjacency(Network, attrname = "weight")
  options(warn=0) # turn warnings on
  if(sum(IS) != 0){
    # if(sum(network::get.edge.attribute(Network, "weight"), na.rm = TRUE) == 0){
    #   stop("Either your network does not contain any edges with weights or your network does not have the edge attribute `weight` required for calculation of extinctions based on relative interaction strength loss.")
    # }
    netgraph <- suppressMessages(igraph::graph_from_adjacency_matrix(net, weighted = TRUE))
    if(NetworkType == "Trophic"){
      strengthbaseout <- igraph::strength(netgraph, mode = "out")
      strengthbasein <- igraph::strength(netgraph, mode = "in")
    }
    if(NetworkType == "Mutualistic"){
      strengthbasenet <- igraph::strength(netgraph)
    }
  }

  ## identification of producers and top predators ++++++++++ ++++++++++
  indegreebasenet <- sna::degree(Network, cmode = "indegree")
  indegreebasenetzeros <- sum(sna::degree(Network, cmode = "indegree") == 0)
  indegreetopnetzeros <- sum(sna::degree(Network, cmode = "outdegree") == 0)
  Producers <- network::get.vertex.attribute(Network, "vertex.names")[sna::degree(Network, cmode = "indegree") == 0]
  TopPredators <- network::get.vertex.attribute(Network, "vertex.names")[sna::degree(Network, cmode = "outdegree") == 0]

  ## output object ++++++++++ ++++++++++
  DF <- data.frame(Spp = rep(NA, length(Order)),
                   S = rep(NA, length(Order)),
                   L = rep(NA, length(Order)),
                   C = rep(NA, length(Order)),
                   Link_density = rep(NA, length(Order)),
                   SecExt = rep(NA,length(Order)),
                   Pred_release = rep(NA,length(Order)),
                   Iso_nodes = rep (NA,length(Order)),
                   Modularity = rep (NA,length(Order)))
  Secundaryext <- c()
  Predationrel <- c()
  accExt <- c()
  totalExt <- c()
  FinalExt <- list()
  Conected3 <- c()

  # Sequential extinction simulation ++++++++++ ++++++++++ ++++++++++ ++++++++++
  if(verbose){ProgBar <- txtProgressBar(max = length(Order), style = 3)}
  primskip <- c()
  i <- 1
  while(i <= length(Order)){
    # for (i in 1:length(Order)){
    if(is.na(Conected1[i])){
      if(verbose){setTxtProgressBar(ProgBar, i)}
      i <- i+1
      next()
    }
    # print(i)

    if(i > nrow(DF)){
      DF <- rbind(DF, data.frame(
        Spp = NA,
        S = NA,
        L = NA,
        C = NA,
        Link_density = NA,
        SecExt = NA,
        Pred_release = NA,
        Iso_nodes = NA,
        Modularity = NA)
      )
    }

    ### creating temporary network + deleting vertices if they have been set to go extinct ++++++++++ ++++++++++
    if(length(accExt)==0){ # on first iteration
      Temp <- Network
      DF$Spp[i] <- Conected1[i]
      network::delete.vertices(Temp, c(DF$Spp[1:i]))

      if(network::network.size(Temp) < 1){
        if(verbose){setTxtProgressBar(ProgBar, length(Order))}
        warning(paste("All nodes in your network went extinct before all primary extinctions were simulated. This happened at extinction step", i-1, "out of", length(Order)))
        break
      }

      if(dim(edgelist)[1] == 0){
        if(verbose){setTxtProgressBar(ProgBar, length(Order))}
        warning(paste("Your network became completely unconnected before all primary extinctions were simulated. This happened at extinction step", i-1, "out of", length(Order)))
        break
      }

    }
    if(length(accExt)>0){ # on any subsequent iteration
      Temp <- Network
      Temp <- network::delete.vertices(Temp, c(accExt))
      edgelist <- network::as.matrix.network.edgelist(Temp,matrix.type="edgelist")

      if(network::network.size(Temp) < 1){
        if(verbose){setTxtProgressBar(ProgBar, length(Order))}
        warning(paste("All nodes in your network went extinct before all primary extinctions were simulated. This happened at extinction step", i-1, "out of", length(Order)))
        break
      }

      if(dim(edgelist)[1] == 0){
        if(verbose){setTxtProgressBar(ProgBar, length(Order))}
        warning(paste("Your network became completely unconnected before all primary extinctions were simulated. This happened at extinction step", i-1, "out of", length(Order)))
        break
      }

      if(RecalcConnect != FALSE){
        Conected2 <- data.frame(ID =
                                  # get.vertex.attribute(Temp, "vertex.names"),
                                  1:network::network.size(Temp),
                                Grado = sna::degree(edgelist, c("total")))
        if(RecalcConnect == 1){
          Conected2 <- arrange(Conected2, desc(Grado))
        }
        if(RecalcConnect == -1){
          Conected2 <- arrange(Conected2[Conected2$Grado != 0,], Grado)
        }
        for(j in sort(accExt)){
          Conected2$ID <- ifelse(Conected2$ID < j, Conected2$ID, Conected2$ID + 1)
        }
        DF$Spp[i] <- Conected2$ID[1]
      }else{

        if(Conected1[i] %in% accExt){
          primskip <- c(primskip, Conected1[i])
          Conected1 <- Conected1[-i]
          if(verbose){setTxtProgressBar(ProgBar, length(Order))}
          if(is.na(Conected1[i])){break()}
        }

        DF$Spp[i] <- Conected1[i]
      }

      Temp <- Network
      network::delete.vertices(Temp, unique(c(c(DF$Spp[1:i]),accExt)))
    }

    ## network metrics to output object  ++++++++++ ++++++++++
    DF$S[i] <- network::network.size(Temp)
    DF$L[i] <- network::network.edgecount(Temp)
    DF$C[i] <- network::network.density(Temp)
    DF$Link_density[i] <- DF$L[i]/DF$S[i]

    ## premature complete annihilation message ++++++++++ ++++++++++
    if(i > 1){
      if(DF$L[i-1] == 0){
        if(verbose){setTxtProgressBar(ProgBar, length(Order))}
        warning(paste("Your network become completely unconnected before all primary extinctions were simulated. This happened at extinction step", i-1, "out of", length(Order)))
        break
      }
    }

    ## calculating modularity ++++++++++ ++++++++++
    Networkclass = class(Temp)
    if (Networkclass[1] == "matrix"){
      netgraph = igraph::graph_from_adjacency_matrix(Temp, mode = "directed", weighted = TRUE)
    }
    if (Networkclass[1] == "network"){
      net = network::as.matrix.network.adjacency(Temp)
      netgraph = suppressMessages(igraph::graph_from_adjacency_matrix(net, mode = "directed", weighted = TRUE))
    }
    if (clust.method == "cluster_edge_betweenness"){
      Membership = suppressWarnings(cluster_edge_betweenness(netgraph, weights = igraph::E(
        igraph::as.undirected(netgraph)
      )$weight, directed = TRUE, edge.betweenness = TRUE,
      merges = TRUE, bridges = TRUE, modularity = TRUE, membership = TRUE))
    }else if (clust.method == "cluster_label_prop"){
      Membership = suppressWarnings(cluster_label_prop(netgraph, weights = igraph::E(igraph::as.undirected(netgraph))$weight, initial = NULL,
                                                       fixed = NULL))
    }else if (clust.method == "cluster_infomap"){
      nb.trials = network::network.size(Temp)
      Membership = suppressWarnings(igraph::cluster_infomap(igraph::as.undirected(netgraph),
                                                            e.weights = igraph::E(
                                                              igraph::as.undirected(netgraph)
                                                            )$weight,
                                                            v.weights = NULL,
                                                            nb.trials = nb.trials,
                                                            modularity = TRUE))

    } else if (clust.method == "none"){
      Membership = NA
    }else stop('Select a valid method for clustering. ?SimulateExtinction')
    DF$Modularity[i] <- Membership$modularity

    ## rewiring ++++++++++ ++++++++++
    accExt <- unique(append(accExt, DF$Spp[1:i]))
    if(!isFALSE(Rewiring)){
      ### identify rewiring potential ++++++++++
      Rewiring_df <- data.frame(Direction = NA,
                                Species = NA,
                                NewPartner = NA,
                                LostPartner = NA,
                                IS = NA)
      Rewiring_df <- na.omit(Rewiring_df)
      #### loop over all deleted vertices and the connections lost because of their exclusion
      for(Iter_PrimaryExt in 1:length(accExt)){
        # Iter_PrimaryExt = 1
        LostPartner <- network::get.vertex.attribute(Network, "vertex.names")[accExt[Iter_PrimaryExt]] # name of primary extinction species
        LostISCol <- Weight_mat[, LostPartner] # lost interaction strength with nodes now slated for secondary extinction
        LostISRow <- Weight_mat[LostPartner, ]
        Lost_df <- data.frame(LostIS = c(LostISCol, LostISRow),
                              Direction = rep(c(1,2), c(length(LostISCol), length(LostISRow))),
                              names = c(names(LostISCol), names(LostISRow))
        )
        Lost_df <- Lost_df[Lost_df$LostIS != 0, ]
        if(nrow(Lost_df)!=0){
          for(Iter_LostIS in 1:nrow(Lost_df)){ ## looping over all species that were linked to the current primary extinction
            # Iter_LostIS = 1
            if(Rewiring[which(names(Rewiring) == Lost_df$names[Iter_LostIS])] == "function (x)  {     x }"){
              LostPartnerSim <- eval(str2lang(Rewiring[which(names(Rewiring) == Lost_df$names[Iter_LostIS])]))(RewiringDist[,Lost_df$names[Iter_LostIS]]) # probability of rewiring to each node
            }else{
              LostPartnerSim <- eval(str2lang(Rewiring[which(names(Rewiring) == Lost_df$names[Iter_LostIS])]))(RewiringDist[,LostPartner]) # probability of rewiring to each node in network given rewiring function and species similraity
            }
            names(LostPartnerSim) <- colnames(RewiringDist)
            RewiringCandidates <- LostPartnerSim[LostPartnerSim > RewiringProb & names(LostPartnerSim) %in% network::get.vertex.attribute(Temp, "vertex.names")] # rewiring probability for nodes still in temporary network and having a higher rewiring probability than 0.5
            RewiredPartner <- names(which.max(RewiringCandidates)) # most likely rewiring partner
            if(!is.null(RewiredPartner)){ # if a rewired partner has been found
              Rewiring_df <- rbind(Rewiring_df,
                                   data.frame(Direction = Lost_df$Direction[Iter_LostIS],
                                              Species = Lost_df$names[Iter_LostIS],
                                              NewPartner = RewiredPartner,
                                              LostPartner = LostPartner,
                                              IS = Lost_df$LostIS[Iter_LostIS])
              )
            }
          }
        }
      }

      ### shift rewired interaction strengths ++++++++++
      if(nrow(Rewiring_df) != 0){
        #### shift interaction weights in Weight_mat
        for(Iter_Rewiring in 1:nrow(Rewiring_df)){
          # Iter_Rewiring = 1
          ## assigning shifted interaction strength
          ColSpec <- Rewiring_df[Iter_Rewiring,4-Rewiring_df[Iter_Rewiring,"Direction"]]
          RowSpec <- Rewiring_df[Iter_Rewiring,1+Rewiring_df[Iter_Rewiring,"Direction"]]
          Weight_mat[RowSpec, ColSpec] <- Weight_mat[RowSpec, ColSpec] + Rewiring_df[Iter_Rewiring,"IS"]
          ## deleting shiften interaction strength

          ColLost <- ifelse(Rewiring_df[Iter_Rewiring, "Direction"] == 1,
                            Rewiring_df[Iter_Rewiring, "LostPartner"],
                            Rewiring_df[Iter_Rewiring, "Species"])
          RowLost <- ifelse(Rewiring_df[Iter_Rewiring, "Direction"] == 1,
                            Rewiring_df[Iter_Rewiring, "Species"],
                            Rewiring_df[Iter_Rewiring, "LostPartner"])
          Weight_mat[RowLost, ColLost] <- 0
        }
        #### establishing rewired network and deleting primary extinction nodes
        Network <- as.network(Weight_mat, matrix.type = "adjacency", ignore.eval=FALSE, names.eval='weight')
        Temp <- Network
        network::delete.vertices(Temp, unique(c(c(DF$Spp[1:i]),accExt)))
      }
    }

    ## identify secondary extinctions ++++++++++ ++++++++++
    ### Relative Interaction Strength loss ++++++++++
    if(sum(IS) == 0){
      SecundaryextNames <- network::get.vertex.attribute(Temp, "vertex.names")[which(sna::degree(Temp) == 0)]
      Secundaryext <- match(SecundaryextNames, network::get.vertex.attribute(Network, "vertex.names"))
    }else{
      if(NetworkType == "Trophic"){
        AbsISin <- igraph::strength(suppressMessages(igraph::graph_from_adjacency_matrix(
          network::as.matrix.network.adjacency(Temp, attrname = "weight"),
          weighted = TRUE)
        ), mode = "in")
        RelISRemainIn <-  AbsISin / strengthbasein[names(strengthbasein) %in% network::get.vertex.attribute(Temp, "vertex.names")]
        SecundaryextNames <- names(which(AbsISin == 0 | RelISRemainIn < IS[match(names(RelISRemainIn), names(IS))]))
        SecundaryextNames <- SecundaryextNames[!(SecundaryextNames %in% Producers)]
        Secundaryext <- match(SecundaryextNames, network::get.vertex.attribute(Network, "vertex.names"))
      }

      if(NetworkType == "Mutualistic"){
        AbsIS <- igraph::strength(suppressMessages(igraph::graph_from_adjacency_matrix(
          network::as.matrix.network.adjacency(Temp, attrname = "weight"),
          weighted = TRUE)
        ))
        RelISloss <-  AbsIS / strengthbasenet[names(strengthbasenet) %in% network::get.vertex.attribute(Temp, "vertex.names")]
        SecundaryextNames <- names(which(AbsIS == 0 | RelISloss < IS[match(names(RelISloss), names(IS))]))
        Secundaryext <- match(SecundaryextNames, network::get.vertex.attribute(Network, "vertex.names"))
      }
    }

    ### for trophic networks ++++++++++
    if(NetworkType == "Trophic"){
      MidPredExt <- network::get.vertex.attribute(Temp, "vertex.names")[sna::degree(Temp, cmode = "indegree") == 0]
      MidPredExt <- match(MidPredExt[!(MidPredExt %in% Producers)], network::get.vertex.attribute(Network, "vertex.names"))
      SecundaryextTrue <- unique(c(SecundaryextNames[!(SecundaryextNames %in% as.character(Producers))],
                                   MidPredExt))
      Secundaryext <- match(SecundaryextTrue, network::get.vertex.attribute(Network, "vertex.names"))
      DF$SecExt[i] <- length(Secundaryext)
      DF$Pred_release[i] <- length(SecundaryextNames[!(SecundaryextNames %in% as.character(TopPredators))])
    }
    ### for mutualistic networks ++++++++++
    if(NetworkType == "Mutualistic"){
      DF$SecExt[i] <- length(Secundaryext)
    }
    DF$Iso_nodes[i] <- sum(sna::degree(Temp) == 0)

    ## Return of objects ++++++++++ ++++++++++
    FinalExt[[i]] <- Secundaryext
    accExt <- append(accExt, DF$Spp[1:i])
    accExt <- unique(append(accExt,Secundaryext))
    if(verbose){setTxtProgressBar(ProgBar, i)}

    ## final Temp deletion ++++++++++ ++++++++++
    if(length(accExt)>0){
      Temp <- Network
      Temp <- network::delete.vertices(Temp, c(accExt))
      edgelist <- network::as.matrix.network.edgelist(Temp,matrix.type="edgelist")
      # DF$Spp[i] <- Conected1[i]
      Temp <- Network
      network::delete.vertices(Temp, unique(c(c(DF$Spp[1:i]),accExt)))
    }

    if(forceFULL){
      Conected1 <- Order <- unique(c(Order, Secundaryext))
    }
    i <- i+1
  }


  # return of final data objects ++++++++++ ++++++++++
  DF <- DF[!is.na(DF$Spp),]
  DF$AccSecExt <- cumsum(DF$SecExt)
  DF$NumExt <- 1:nrow(DF)
  DF$TotalExt <- DF$AccSecExt + DF$NumExt
  DF <- relocate(DF, Modularity, .after = Link_density)
  class(DF) <- c("data.frame", "SimulateExt")

  # return of robustness metrics



  if(length(primskip)!= 0){warning(paste("Primary extinctions of", paste(primskip, collapse = ", "), "skipped due to their prior extinction as secondary extinctions."))}

  NetworkSize <- network::network.size(Network)
  FracTotalExt <- DF$TotalExt/NetworkSize
  positionR50 <- which(round((FracTotalExt),1) >= 0.5)
  positionR100 <- which(round((FracTotalExt),1) >= 1)

  if (length(positionR50) > 0){
    position <- positionR50[1]
    R50 <- round((position/NetworkSize), 2)
  } else{
    nodeContribution <- 1/NetworkSize
    HighestValue <- max(FracTotalExt, na.rm = TRUE)
    NremovalslackingtoR50<- (0.5-HighestValue)/nodeContribution
    ObservedNremovals <- which(FracTotalExt == HighestValue)
    R50 <- round(((ObservedNremovals+NremovalslackingtoR50)/NetworkSize),2)
  }

  if (length(positionR100) > 0){
    position <- positionR100[1]
    R100 <- round((position/NetworkSize), 2)
  } else{
    nodeContribution <- 1/NetworkSize
    HighestValue <- max(FracTotalExt, na.rm = TRUE)
    NremovalslackingtoR50<- (1-HighestValue)/nodeContribution
    ObservedNremovals <- which(FracTotalExt == HighestValue)
    R100 <- round(((ObservedNremovals+NremovalslackingtoR50)/NetworkSize),2)
  }

  # Your function code here

  return(list(sims = DF,
              R50 = R50,
              R100 = R100,
              Network = Temp))
}

#' Random extinction
#'
#' Generates a null model by generating random extinction histories and calculating the mean and standard deviation of the accumulated secondary extinctions developed by making n random extinction histories.
#'
#' @param Network a network representation as a an adjacency matrix, edgelist,
#' or a network object
#' @param nsim numeric, number of simulations
#' @param Record logical, if TRUE, records every simulation and you can read the
#' raw results in the object FullSims
#' @param plot logical if TRUE, will add a graph to the results
#' @param SimNum numeric, how many nodes to register for primary extinction. By default sets all of them.
#' @param NetworkType a character with the options Trophic and Mutualistic - is used to calculate secondary extinctions.
#' @param clust.method a character with the options cluster_edge_betweenness,
#' cluster_label_prop or cluster_infomap, defaults to cluster_infomap
#' @param parallel if TRUE, it will use parallel procesing, if FALSE (default) it will run
#' sequentially
#' @param ncores numeric, number of cores to use if using parallel procesing
#' @param IS either numeric or a named vector of numerics. Identifies the threshold of relative interaction strength which species require to not be considered secondarily extinct (i.e. IS = 0.3 leads to removal of all nodes which lose 70 precent of their interaction strength in the Network argument). If a named vector, names must correspond to vertex names in Network argument.
#' @param Rewiring either a function or a named vector of functions. Signifies how rewiring probabilities are calculated from the RewiringDist argument. If FALSE, no rewiring is carried out.
#' @param RewiringDist a numeric matrix of NxN dimension (N... number of nodes in Network). Contains, for example, phylogenetic or functional trait distances between nodes in Network which are used by the Rewiring argument to calculate rewiring probabilities. If Rewiring == function(x){x}, this matrix is expected to contain probabilities of a connection being present between species-pairs.
#' @param RewiringProb a numeric which identifies the threshold at which to assume rewiring potential is met.
#' @param verbose Logical. Whether to report on function progress or not.
#' @param forceFULL Logical. Whether to continue removal of nodes after initial order has been depleted. This will force the simulations to execute extinctions and check for secondary extinctions/new links until the network does not change anylonger.
#' @return exports list containing a data frame with the characteristics of the network after every extinction, a network object containing the final network, and a graph with the mean and 95percent interval. The resulting data frame contains 11 columns that incorporate the topological index, the secondary extinctions, predation release, and total extinctions of the network in each primary extinction.
#' @details
#'
#' "Note: When using the pre-defined order of nodes for primary
#' removals option in the random extinction scenario, it is possible
#' that some of the species in the predefined order may be lost as
#' secondary extinctions. As such, they should not be counted as
#' primary removals. For example, if a network has five species {A,B,
#' ,D,E} and a pre-defined "random" order of removal {C,A,B,E,D} with
#' removal of C causing the additional loss of A and removal of B
#' causing the additional loss of E and D, only two primary removals
#' (C and B) would be required for total network collapse, even though
#' the algorithm would terminate at the third element of the removal
#' vector, i.e., {C,A,B}."
#'
#' When NetworkType = Trophic, secondary extinctions only occur for any predator, but not producers. If NetworkType = Mutualistic, secondary extinctions occur for all species in the network.
#'
#' When clust.method = cluster_edge_betweenness computes the network modularity using cluster_edge_betweenness methods from igraph to detect communities
#' When clust.method = cluster_label_prop computes the network modularity using cluster_label_prop methods from igraph to detect communities
#' When clust.method = cluster_infomap computes the network modularity using cluster_infomap methods from igraph to detect communities, here the number of nb.trials are equal to the network size
#'
#' @examples
#' #first example
#' \dontrun{
#' data("More_Connected")
#' RandomExtinctions(Network = More_Connected, nsim = 20)
#'
#' # Using parallel procesing
#' ## Detect your number of cores divide by 2
#'
#' cores <- ceiling(parallel::detectCores()/2)
#'
#' RandomExtinctions(Network = More_Connected, nsim = 20, parallel = TRUE, ncores = cores)
#' }
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr n
#' @importFrom dplyr summarise
#' @importFrom dplyr relocate
#' @importFrom dplyr everything
#' @importFrom foreach `%dopar%`
#' @importFrom foreach foreach
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 element_blank
#' @importFrom magrittr "%>%"
#' @importFrom network network.size
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel clusterExport
#' @importFrom scales muted
#' @importFrom stats sd
#' @importFrom stats na.omit
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom patchwork wrap_plots
#' @importFrom patchwork plot_annotation
#' @importFrom doSNOW registerDoSNOW


#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author M. Isidora Ávila-Thieme <msavila@uc.cl>
#' @author Erik Kusch <erik.kusch@bio.au.dk>
#' @export
RandomExtinctions <- function(Network, nsim = 10,
                              Record = FALSE, plot = FALSE,
                              SimNum = NULL,
                              NetworkType = "Trophic", clust.method = "cluster_infomap",
                              parallel = FALSE, ncores,
                              IS = 0,
                              Rewiring = FALSE, RewiringDist = NULL, RewiringProb = 0.5,
                              verbose = TRUE, forceFULL = FALSE){
  if(!NetworkType %in% c("Trophic", "Mutualistic")){stop("Please specify NetworkType as either 'Trophic' or 'Mutualistic'")}
  ## setting up objects
  NumExt <- sd <- AccSecExt <- AccSecExt_95CI <- AccSecExt_mean <- Lower <- Upper <-NULL
  network <- .DataInit(x = Network)
  if(is.null(SimNum)){
    SimNum <- network::network.size(network)
  }

  ## simulations
  if(verbose & !parallel){ProgBar <- txtProgressBar(max = nsim, style = 3)}
  if(parallel){
    cl <- makeCluster(ncores)
    # registerDoParallel(cl)
    registerDoSNOW(cl)
    parallel::clusterExport(cl,
                            varlist = c("network", "SimNum", "IS", "Rewiring", "RewiringDist", "RewiringProb", "SimNum", "forceFULL"),
                            envir = environment()
    )
    pb <- txtProgressBar(max = nsim, style = 3)
    progress <- function(i) setTxtProgressBar(pb, i)
    opts <- list(progress = progress)

    sims <- foreach(i=1:nsim, .options.snow = opts, .packages = c("NetworkExtinction", "dplyr"))%dopar%{
      Order <- sample(1:network::network.size(network), size = SimNum)
      sims1 <- try(ExtinctionOrder(Network = network,
                                   Order = Order,
                                   IS = IS, NetworkType = NetworkType,
                                   Rewiring = Rewiring, RewiringDist = RewiringDist,
                                   verbose = FALSE, RewiringProb = RewiringProb), silent = TRUE)
      try({sims1$sims$simulation <- i}, silent = TRUE)
      return(sims1)
    }
    stopCluster(cl)
  }else{
    sims <- list()
    for(i in 1:nsim){
      Order <- sample(1:network::network.size(network), size = SimNum)
      sims[[i]] <- try(ExtinctionOrder(Network = network, Order = Order,
                                       IS = IS, NetworkType = NetworkType,
                                       Rewiring = Rewiring, RewiringDist = RewiringDist,
                                       verbose = FALSE, RewiringProb = RewiringProb), silent = TRUE)
      try({sims[[i]]$sims$simulation <- i}, silent = TRUE)
      if(verbose){setTxtProgressBar(ProgBar, i)}
    }
  }

  r50values <- c()
  r100values <- c()
  for(i in 1:nsim){
    r50values[i] <- sims[[i]]$R50
    r100values[i] <-sims[[i]]$R100
  }

  R50mean <- mean(r50values)
  R50CI <- 1.96*sd(r50values)
  R50result <- c(R50mean,R50CI)
  names(R50result) <- c("Mean", "CI")

  R100mean <- mean(r100values)
  R100CI <- 1.96*sd(r100values)
  R100result <- c(R100mean,R100CI)
  names(R50result) <- c("Mean", "CI")

  ## extract objects
  temps <- lapply(sims, "[[", 4)
  sims <- lapply(sims, "[[", 1)
  cond <- sapply(sims, function(x) "data.frame" %in% class(x))
  cond <- c(1:length(cond))[cond]
  sims <- sims[cond]
  sims <- do.call(rbind, sims)
  if(Record == TRUE){
    FullSims <- sims
  }


  sims <- sims[!is.na(sims$SecExt), ] %>% dplyr::group_by(NumExt) %>% summarise(AccSecExt_95CI = 1.96*sd(AccSecExt), AccSecExt_mean = mean(AccSecExt), nsim = dplyr::n()) %>% mutate(Upper = AccSecExt_mean + AccSecExt_95CI, Lower = AccSecExt_mean - AccSecExt_95CI, Lower = ifelse(Lower < 0, 0, Lower)) %>% dplyr::relocate(nsim, .after = dplyr::everything())

  ## plot output
  if(plot == TRUE){
    g <- ggplot(sims, aes(x = NumExt, y = AccSecExt_mean)) + geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = scales::muted("red")) + geom_line() + ylab("Acc. Secondary extinctions") + xlab("Primary extinctions") + theme_bw() + ggplot2::theme(axis.title.x = ggplot2::element_blank())
    h <- ggplot(sims, aes(x = NumExt, y = nsim/max(nsim))) + geom_line() + theme_bw() + labs(y = "Prop", x = "Primary extinctions")
    I <- patchwork::wrap_plots(g, h, ncol = 1, guides = "keep", heights = c(3,1)) + patchwork::plot_annotation(tag_levels = 'A')


    print(I)
  }
  message("Note: If using the pre-defined order of nodes for primary removals option, please note that some species in the predefined order may be lost as secondary extinctions and should not be counted as primary removals.")

  ## object output
  if(Record == T & plot == T){
    return(list(sims = sims, graph = g, FullSims = FullSims, nets = temps, R50result = R50result, R100result = R100result))
  }else if(Record == F & plot == T){
    return(list(sims = sims, graph = g, nets = temps, R50result = R50result, R100result = R100result))
  }else if(Record == F & plot == F){
    return(list(sims = sims, nets = temps, R50result = R50result, R100result = R100result))
  }else if(Record == T & plot == F){
    return(list(sims = sims, FullSims = FullSims, nets= temps, R50result = R50result, R100result = R100result))
  }
}

#' Comparison of Null hypothesis with other extinction histories
#'
#' It compares an object generated either by the Mostconnected or ExtinctionOrder functions
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
#' \dontrun{
#' data("Less_Connected")
#' History <- SimulateExtinctions(Network = Less_Connected, Method = "Mostconnected")
#' NullHyp <- RandomExtinctions(Network = Less_Connected, nsim = 100)
#' CompareExtinctions(Nullmodel = NullHyp, Hypothesis = History)
#' }
#' @importFrom broom tidy
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 theme
#' @importFrom patchwork wrap_plots
#' @importFrom patchwork plot_annotation


#' @importFrom scales muted
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author M. Isidora Ávila-Thieme <msavila@uc.cl>
#' @export

CompareExtinctions <- function(Nullmodel, Hypothesis){
  if(class(Hypothesis$sims)[2] == "SimulateExt"){
    NumExt <- sd <- AccSecExt <- AccSecExt_mean <- Lower <- Upper <- nsim<-NULL
    g <- ggplot(Nullmodel$sims, aes(x = NumExt, y = AccSecExt_mean)) + geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = scales::muted("red")) + geom_line(aes(color = "blue")) + ylab("Acc. Secondary extinctions") + xlab("Primary extinctions") + theme_bw() + ggplot2::theme(axis.title.x = ggplot2::element_blank())
    g <- g + geom_point(data = Hypothesis$sims, aes(y = AccSecExt), color = "black") + geom_line(data = Hypothesis$sims, aes(y = AccSecExt, color = "black")) + scale_color_manual(name = "Comparison",values =c("black", "blue"), label = c("Observed","Null hypothesis"))
    h <- ggplot(Nullmodel$sims, aes(x = NumExt, y = nsim/max(nsim))) + geom_line() + theme_bw() + labs(y = "Prop", x = "Primary extinctions")
    I <- patchwork::wrap_plots(g, h, ncol = 1, guides = "keep", heights = c(3,1)) + patchwork::plot_annotation(tag_levels = 'A')
    I
    return(I)
  }
  if(class(Hypothesis$sims)[2] %in% c("Mostconnected", "ExtinctionOrder")){
    NumExt <- sd <- AccSecExt <- AccSecExt_mean <-NULL
    g <- ggplot(Nullmodel$sims, aes(x = NumExt, y = AccSecExt_mean)) + geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = scales::muted("red")) + geom_line(aes(color = "blue")) + ylab("Acc. Secondary extinctions") + xlab("Primary extinctions") + theme_bw() + ggplot2::theme(axis.title.x = ggplot2::element_blank())
    g <- g + geom_point(data = Hypothesis, aes(y = AccSecExt), color = "black") + geom_line(data = Hypothesis, aes(y = AccSecExt, color = "black")) + scale_color_manual(name = "Comparison", values =c("black", "blue"), label = c("Observed","Null hypothesis"))
    I <- patchwork::wrap_plots(g, h, ncol = 1, guides = "keep", heights = c(3,1))
    I
    return(I)
  }
  else{
    message("Hypothesis not of class Mostconnected or ExtinctionOrder")
  }
}
