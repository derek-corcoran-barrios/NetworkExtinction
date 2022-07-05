ExtinctionOrder <- function(Network, Order, IS = 0,
                            Rewiring = NULL, RewiringDist = NULL,
                            verbose = TRUE, clust.method = "cluster_infomap"){
  ## Setting up Objects for function run
  Link_density <- Modularity <- Grado <- NULL
  Network <- .DataInit(x = Network)
  edgelist <- as.matrix.network.edgelist(Network,matrix.type="edgelist") #Prey - Predator
  Conected <- data.frame(ID = 1:network.size(Network), Grado = degree(edgelist, c("total")))
  Conected1 <- Order
  if(length(IS )== 1){
    IS <- rep(IS, network.size(Network))
    names(IS) <- get.vertex.attribute(Network, "vertex.names")
  }

  ## Base net calculations
  ### identify base interaction strengths per node
  if(sum(IS) != 0){
    if(sum(get.edge.attribute(Network, "weight"), na.rm = TRUE) == 0){
      stop("Either your network does not contain any edges with weights or your network does not have the edge attribute `weight` required for calculation of extinctions based on relative interaction strength loss.")
    }
    net <- as.matrix.network.adjacency(Network, attrname = "weight")
    netgraph <- suppressMessages(graph_from_adjacency_matrix(net, weighted = TRUE))
    strengthbasenet <- igraph::strength(netgraph)
  }
  Weight_mat <- as.matrix.network.adjacency(Network, attrname = "weight")

  ### identification of producers and top predators
  indegreebasenet <- degree(Network, cmode = "indegree")
  indegreebasenetzeros <- sum(degree(Network, cmode = "indegree") == 0)
  indegreetopnetzeros <- sum(degree(Network, cmode = "outdegree") == 0)
  Producers <- (1:length(degree(Network, cmode = "indegree")))[degree(Network, cmode = "indegree") == 0]
  TopPredators <- (1:length(degree(Network, cmode = "outdegree")))[degree(Network, cmode = "outdegree") == 0]

  ### output object
  DF <- data.frame(Spp = rep(NA, length(Order)),
                   S = rep(NA, length(Order)),
                   L = rep(NA, length(Order)),
                   C = rep(NA, length(Order)),
                   Link_density = rep(NA, length(Order)),
                   SecExt = rep(NA,length(Order)),
                   Pred_release = rep(NA,length(Order)),
                   Iso_nodes =rep (NA,length(Order)))
  Secundaryext <- c()
  Predationrel <- c()
  accExt <- c()
  totalExt <- c()
  FinalExt <- list()
  Conected3 <- c()

  ## Rewiring
  if(!isFALSE(Rewiring)){
    ### default decay parameter
    diag(RewiringDist)<- NA
    if(is.null(Rewiring)){
      dist_10 <- quantile(RewiringDist, na.rm = TRUE, 0.1) # 10% distance quantile
      decay <- 1/as.numeric(dist_10) # assuming mean rewiring capability lies at dist_10
      Rewiring <- function(x){1-pexp(x, rate = decay)}
    }

    if(length(Rewiring )== 1){
      fun <- deparse1(Rewiring)
      Rewiring <- rep(fun, network.size(Network))
      names(Rewiring) <- get.vertex.attribute(Network, "vertex.names")
    }
    # plot_seq <- seq(from = 0, to = max(RewiringDist[RewiringDist != 0], na.rm = TRUE), length = 1e3)
    # plot(plot_seq,
    #      eval(str2lang(Rewiring[1]))(x)
    #      )
  }


  ## Sequential extinction simulation
  if(verbose){ProgBar <- txtProgressBar(max = length(Order), style = 3)}
  for (i in 1:length(Order)){
    # print(i)

    ### creating temporary network representations and deleting vertices if they have been set to go extinct
    if (length(accExt)==0){
      Temp <- Network
      DF$Spp[i] <- Conected1[i]
      delete.vertices(Temp, c(DF$Spp[1:i]))
    }
    if (length(accExt)>0){
      Temp <- Network
      Temp <- delete.vertices(Temp, c(accExt))
      edgelist <- as.matrix.network.edgelist(Temp,matrix.type="edgelist")
      # Conected2 <- data.frame(ID =1:network.size(Temp), Grado = degree(edgelist, c("total")))
      # for(j in sort(accExt)){
      #   Conected2$ID <- ifelse(Conected2$ID < j, Conected2$ID, Conected2$ID + 1)
      # }

      DF$Spp[i] <- Conected1[i]
      Temp <- Network

      delete.vertices(Temp, unique(c(c(DF$Spp[1:i]),accExt)))

    }
    # print(Temp)

    ### network metrics to output object
    DF$S[i] <- network.size(Temp)
    DF$L[i] <- network.edgecount(Temp)
    DF$C[i] <- network.density(Temp)
    DF$Link_density[i] <- DF$L[i]/DF$S[i]

    if(i > 1 ){
      if(DF$L[i-1] == 0){
        if(verbose){setTxtProgressBar(ProgBar, length(Order))}
        warning(paste("All species in network went extinct through secondary extinction before all primary extinctions were simulated. This happened at extinction step", i-1, "out of", length(Order)))
        break
      }
    }


    ### calculating modularity
    Networkclass = class(Temp)
    if (Networkclass[1] == "matrix"){
      netgraph = graph_from_adjacency_matrix(Temp, mode = "directed", weighted = TRUE)
    }

    if (Networkclass[1] == "network"){
      net = as.matrix.network.adjacency(Temp)
      netgraph = suppressMessages(graph_from_adjacency_matrix(net, mode = "directed", weighted = TRUE))
    }

    if (clust.method == "cluster_edge_betweenness"){
      Membership = suppressWarnings(cluster_edge_betweenness(netgraph, weights = TRUE, directed = TRUE, edge.betweenness = TRUE,
                                                             merges = TRUE, bridges = TRUE, modularity = TRUE, membership = TRUE))
    } else if (clust.method == "cluster_spinglass"){
      spins = 107#network.size(Temp)
      Membership = suppressWarnings(cluster_spinglass(netgraph, spins=spins)) #spins could be the Richness
    }else if (clust.method == "cluster_label_prop"){
      Membership = suppressWarnings(cluster_label_prop(netgraph, weights = TRUE, initial = NULL,
                                                       fixed = NULL))
    }else if (clust.method == "cluster_infomap"){
      nb.trials = 107#network.size(Temp)
      Membership = suppressWarnings(cluster_infomap(as.undirected(netgraph),
                                                    e.weights = E(netgraph)$weight,
                                                    v.weights = NULL,
                                                    nb.trials = nb.trials,
                                                    modularity = TRUE))

    } else if (clust.method == "none"){
      Membership = NA
    }else stop('Select a valid method for clustering. ?SimulateExtinction')
    # if(is.na(Membership)){
    DF$Modularity[i] <- Membership$modularity
    # }else{
    #   DF$Modularity[i] <- suppressWarnings(modularity(Membership))
    # }

    ### identifying secondary extinctions
    #### Producers
    SecundaryextTemp <- (1:length(degree(Temp, cmode = "indegree")))[degree(Temp, cmode = "indegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      SecundaryextTemp <- ifelse(SecundaryextTemp < j, SecundaryextTemp, SecundaryextTemp + 1)
    }
    Secundaryext <- SecundaryextTemp
    Secundaryext <- Secundaryext[!(Secundaryext %in% Producers)]
    DF$SecExt[i] <- length(Secundaryext)

    #### Predators
    PredationrelTemp <- (1:length(degree(Temp, cmode = "outdegree")))[degree(Temp, cmode = "outdegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      PredationrelTemp <- ifelse(PredationrelTemp < j, PredationrelTemp, PredationrelTemp + 1)
    }
    Predationrel <- PredationrelTemp
    Predationrel <- Predationrel[!(Predationrel %in% TopPredators)]
    DF$Pred_release[i]<- length(Predationrel)
    DF$Iso_nodes[i] <- sum(degree(Temp) == 0)

    ### rewiring
    accExt <- append(accExt, DF$Spp[1:i])
    if(!isFALSE(Rewiring)){
      Rewiring_df <- data.frame(Direction = NA,
                                Species = NA,
                                NewPartner = NA,
                                LostPartner = NA,
                                IS = NA)
      Rewiring_df <- na.omit(Rewiring_df)
      #### loop over all deleted vertices and the connections lost because of their exclusion
      for(Iter_PrimaryExt in 1:length(accExt)){
        # Iter_PrimaryExt = 1
        LostPartner <- get.vertex.attribute(Network, "vertex.names")[accExt[Iter_PrimaryExt]] # name of primary extinction species
        LostIS <- Weight_mat[, LostPartner] # lost interaction strength with nodes now slated for secondary extinction
        Direction <- 1 # identify column-driven loss
        if(sum(abs(LostIS))==0){# if the primary species is not an animal, LostIS will be filled with 0s, so we need to look for LosTIS in other orientiation in Weight_mat
          LostIS <- Weight_mat[LostPartner, ]
          Direction <- 2 # identify row-driven loss
        }


        for(Iter_LostIS in 1:length(LostIS)){ ## looping over all species that were linked to the current primary extinction
          # Iter_LostIS = 1
          LostPartnerSim <- eval(str2lang(Rewiring[Iter_LostIS]))(dist_mat[,LostPartner]) # probability of rewiring too each node in network given rewiring function and species similraity
          RewiringCandidates <- LostPartnerSim[LostPartnerSim > 0.3 & names(LostPartnerSim) %in% get.vertex.attribute(Temp, "vertex.names")] # rewiring probability for nodes still in temporary network and having a higher rewiring probability than 0.3
          RewiredPartner <- names(which.max(RewiringCandidates)) # most likely rewiring partner
          if(!is.null(RewiredPartner)){ # if a rewired partner has been found
            Rewiring_df <- rbind(Rewiring_df,
                                 data.frame(Direction = Direction,
                                            Species = names(LostIS[Iter_LostIS]),
                                            NewPartner = RewiredPartner,
                                            LostPartner = LostPartner,
                                            IS = LostIS[Iter_LostIS])
            )
          }
        }
      }


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
      delete.vertices(Temp, unique(c(c(DF$Spp[1:i]),accExt)))
    }

    #### Relative Interaction Strength loss
    if(sum(IS) == 0){
      Secundaryext <- get.vertex.attribute(Temp, "vertex.names")[which(degree(Temp) == 0)]
      Secundaryext <- match(Secundaryext, get.vertex.attribute(Network, "vertex.names"))
    }else{
      AbsIS <- igraph::strength(suppressMessages(graph_from_adjacency_matrix(
        as.matrix.network.adjacency(Temp, attrname = "weight"),
        weighted = TRUE)
      ))
      RelISloss <-  AbsIS / strengthbasenet[names(strengthbasenet) %in% get.vertex.attribute(Temp, "vertex.names")]
      Secundaryext <- which(AbsIS == 0 | RelISloss < IS[match(names(RelISloss), names(IS))])
      Secundaryext <- match(names(Secundaryext), get.vertex.attribute(Network, "vertex.names"))
    }
    DF$SecExt[i] <- length(Secundaryext)

    ### Return of objects
    FinalExt[[i]] <- (Secundaryext)
    accExt <- append(accExt, DF$Spp[1:i])
    accExt <- unique(append(accExt,Secundaryext))
    if(verbose){setTxtProgressBar(ProgBar, i)}

  }
  DF <- DF[complete.cases(DF),]
  DF$AccSecExt <- cumsum(DF$SecExt)
  DF$NumExt <- 1:nrow(DF)
  DF$TotalExt <- DF$AccSecExt + DF$NumExt
  DF <- relocate(DF, Modularity, .after = Link_density)
  class(DF) <- c("data.frame", "SimulateExt")
  return(list(sims = DF,
              Network = Temp))
}
