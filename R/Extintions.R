#' Extinction analyses from more to less conected specie in a network
#'
#' It takes a network and it calculates wich species is the most conected
#' of the network, then it extinguishes that species, and calculates the
#' secundary extintions. After that it recalculates the conections and
#' repeats the process until every species becomes extinct.
#'
#' @param Network a trophic network of class network
#' @return exports data frame with the characteristics of the network after every
#' extintion. The resulting data frame contains seven columns
#' @examples
#' data("net")
#' Mostconnected(Network = net)
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
#' @author Isidora Avila <msavila@uc.cl>
#' @seealso [NetworkExtinction::ExtinctionOrder()]
#' @export


Mostconnected <- function(Network){
  Grado <- NULL
  Network <- Network
  edgelist <- as.matrix.network.edgelist(Network,matrix.type="edgelist") #Prey - Predator
  Conected <- data.frame(ID = 1:network.size(Network), Grado = degree(edgelist, c("total")))
  Conected <- arrange(Conected, desc(Grado))
  Conected1<- c(Conected$ID)
  indegreebasenet <- degree(Network, cmode = "indegree")
  indegreebasenetzeros <- sum(degree(Network, cmode = "indegree") == 0)
  Producers <- (1:length(degree(Network, cmode = "indegree")))[degree(Network, cmode = "indegree") == 0]
  DF <- data.frame(Spp = rep(NA, network.size(Network)), nodesS = rep(NA, network.size(Network)), linksS = rep(NA, network.size(Network)), Conectance = rep(NA, network.size(Network)), LinksPerSpecies = rep(NA, network.size(Network)),Secondary_extinctions = rep(NA,network.size(Network)), aislate_nodes =rep (NA,network.size(Network)))

  Secundaryext <- c()
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

    DF$nodesS[i] <- network.size(Temp)
    DF$linksS[i] <- network.edgecount(Temp)
    DF$Conectance[i] <- network.density(Temp)
    DF$LinksPerSpecies [i] <- DF$linksS[i]/DF$nodesS[i]
    SecundaryextTemp <- (1:length(degree(Temp, cmode = "indegree")))[degree(Temp, cmode = "indegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      SecundaryextTemp <- ifelse(SecundaryextTemp < j, SecundaryextTemp, SecundaryextTemp + 1)
    }
    Secundaryext <- SecundaryextTemp
    Secundaryext <- Secundaryext[!(Secundaryext %in% Producers)]
    DF$Secondary_extinctions[i]<- length(Secundaryext)
    DF$aislate_nodes[i] <- sum(degree(Temp) == 0)
    print(i)
    FinalExt[[i]] <-(Secundaryext)
    accExt <- append(accExt, DF$Spp[1:i])
    accExt <- unique(append(accExt,Secundaryext))

    if (DF$linksS[i] == 0) break
  }
  DF <- DF[complete.cases(DF),]
    DF$AccSecondaryExtinction<- cumsum(DF$Secondary_extinctions)
  DF$NumExt <- 1:nrow(DF)
  DF$TotalExt <- DF$AccSecondaryExtinction + DF$NumExt
  class(DF) <- c("data.frame", "Mostconnected")
  return(DF)
}


#' Extinction analysis by custom order
#'
#' It takes a network extinguishes species using a given order, then
#' it calculates the secundary extintions.
#'
#' @param Network a trophic network of class network
#' @param Order Vector with the order of extinctions by ID
#' @return exports data frame with the characteristics of the network after every
#' extintion, and a graph with the mean and 95% interval
#' @examples
#' #first example
#' data("net")
#' ExtinctionOrder(Network = net, Order = c(1,2,3,4,5,6,7,8,9,10))
#' #Second example
#' data("net")
#' ExtinctionOrder(Network = net, Order = c(2,8,9))
#'
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom network as.matrix.network.edgelist
#' @importFrom network delete.vertices
#' @importFrom network network.edgecount
#' @importFrom network network.size
#' @importFrom sna degree
#' @importFrom stats complete.cases
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Isidora Avila <msavila@uc.cl>
#' @export

ExtinctionOrder <- function(Network, Order){
  Grado <- NULL
  edgelist <- as.matrix.network.edgelist(Network,matrix.type="edgelist") #Prey - Predator
  Conected <- data.frame(ID = 1:network.size(Network), Grado = degree(edgelist, c("total")))

  Conected1<-  Order

  indegreebasenet <- degree(Network, cmode = "indegree")
  indegreebasenetzeros <- sum(degree(Network, cmode = "indegree") == 0)
  Producers <- (1:length(degree(Network, cmode = "indegree")))[degree(Network, cmode = "indegree") == 0]
  DF <- data.frame(Spp = rep(NA, network.size(Network)), nodesS = rep(NA, network.size(Network)), linksS = rep(NA, network.size(Network)), Conectance = rep(NA, network.size(Network)),  Secondary_extinctions = rep(NA,network.size(Network)))

  Secundaryext <- c()
  accExt <- c()
  totalExt <- c()
  FinalExt <- list()
  Conected3 <- c()

  ####LOOP####
  for (i in 1:length(Order)){
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
      Conected2 <- data.frame(ID =1:network.size(Temp), Grado = degree(edgelist, c("total")))
      for(j in sort(accExt)){
        Conected2$ID <- ifelse(Conected2$ID < j, Conected2$ID, Conected2$ID + 1)
      }

      DF$Spp[i] <- Conected1[i]
      Temp <- Network

      delete.vertices(Temp, unique(c(c(DF$Spp[1:i]),accExt)))
    }

    DF$nodesS[i] <- network.size(Temp)
    DF$linksS[i] <- network.edgecount(Temp)
    DF$Conectance[i] <- network.density(Temp)
    SecundaryextTemp <- (1:length(degree(Temp, cmode = "indegree")))[degree(Temp, cmode = "indegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      SecundaryextTemp <- ifelse(SecundaryextTemp < j, SecundaryextTemp, SecundaryextTemp + 1)
    }
    Secundaryext <- SecundaryextTemp
    Secundaryext <- Secundaryext[!(Secundaryext %in% Producers)]
    DF$Secondary_extinctions[i]<- length(Secundaryext)
    message(i)
    FinalExt[[i]] <-(Secundaryext)
    accExt <- append(accExt, DF$Spp[1:i])
    accExt <- unique(append(accExt,Secundaryext))

    if (DF$linksS[i] == 0) break
  }
  DF <- DF[complete.cases(DF),]
  DF$AccSecondaryExtinction <- cumsum(DF$Secondary_extinctions)
  DF$NumExt <- 1:nrow(DF)
  DF$TotalExt <- DF$AccSecondaryExtinction + DF$NumExt
  G <- ggplot(DF, aes_string(x = "NumExt", y = "AccSecondaryExtinction")) + geom_line() + ylab("Secondary extinctions") + xlab("number of exctinctions") + theme_classic()
  Results <- list(DF= DF, Graph = G)
  class(Results) <- c("ExtinctionOrder")
  return(Results)
}

#' Neutral model
#'
#' It takes a network extinguishes species using a given order, then
#' it calculates the secundary extintions.
#'
#' @param Network a trophic network of class network
#' @param nsim number of simulations
#' @return exports data frame with the characteristics of the network after every
#' extintion, and a graph with the mean and 95% interval
#' @examples
#' #first example
#' data("net")
#' ExtinctionOrder(Network = net, Order = c(1,2,3,4,5,6,7,8,9,10))
#' #Second example
#' data("net")
#' ExtinctionOrder(Network = net, Order = c(2,8,9))
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom magrittr "%>%"
#' @importFrom network network.size
#' @importFrom scales muted
#' @importFrom stats sd
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Isidora Avila <msavila@uc.cl>
#' @export

RandomExtinctions <- function(Network, nsim = 10){
  NumExt <- sd <- AccSecondaryExtinction <- NULL
  network <- Network
  sims <- list()
  for(i in 1:nsim){
    sims[[i]] <- try(ExtinctionOrder(Network = network, Order = sample(1:network.size(network)))$DF)
    sims[[i]]$simulation <- i
    message(paste("Simulation", i, "of", nsim, "ready"))
  }
  cond <- sapply(sims, function(x) class(x) != "try-error")
  sims <- sims[cond]
  sims <- do.call(rbind, sims)
  sims <- sims %>% group_by(NumExt) %>% summarise(SdAccSecondaryExtinction = sd(AccSecondaryExtinction), AccSecondaryExtinction = mean(AccSecondaryExtinction))
  g <- ggplot(sims, aes_string(x = "NumExt", y = "AccSecondaryExtinction")) + geom_ribbon(aes_string(ymin = "AccSecondaryExtinction - SdAccSecondaryExtinction", ymax = "AccSecondaryExtinction + SdAccSecondaryExtinction"), fill = muted("red")) + geom_line() + theme_classic()
  g
  return(list(sims = sims, graph = g))
}

#' Compare extinction hypothesis
#'
#' It takes a network extinguishes species using a given order, then
#' it calculates the secundary extintions.
#'
#' @param Nullmodel an object generated by the RandomExtinctions
#' @param Hypothesis Extinction history generated by the Mostconected or ExtinctionOrder
#' fuction
#' @return a plot comparing the Nullmodel and the observed model, and a goodness of
#' fit test showing if there are significant differences between the null model and
#' the observed extinction history.
#' @examples
#' data("net")
#' History <- Mostconnected(Network = net)
#'
#' NullHyp <- RandomExtinctions(Network = net, nsim = 100)
#'
#' CompareExtinctions(Nullmodel = NullHyp, Hypothesis = History)
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom stats chisq.test
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Isidora Avila <msavila@uc.cl>
#' @export

CompareExtinctions <- function(Nullmodel, Hypothesis){
  if(class(Hypothesis)[1] == "ExtinctionOrder"){
    NumExt <- sd <- AccSecondaryExtinction <- NULL
    g <- Nullmodel$graph
    g <- g + geom_point(data = Hypothesis$DF) + geom_line(data = Hypothesis$DF, lty = 2)
    g
    Test <- chisq.test(x = Hypothesis$DF$AccSecondaryExtinction, y = Nullmodel$sims$AccSecondaryExtinction[1:length(Hypothesis$DF$AccSecondaryExtinction)])
    return(list(Test = Test, graph = g))
  }
  if(class(Hypothesis)[2] == "Mostconnected"){
  NumExt <- sd <- AccSecondaryExtinction <- NULL
  g <- Nullmodel$graph
  g <- g + geom_point(data = Hypothesis) + geom_line(data = Hypothesis, lty = 2)
  g
  Test <- chisq.test(x = Hypothesis$AccSecondaryExtinction, y = Nullmodel$sims$AccSecondaryExtinction[1:length(Hypothesis$AccSecondaryExtinction)])
  return(list(Test = Test, graph = g))
  }
  else{
    message("Hipothesis not of class Mostconnected or ExtinctionOrder")
  }
}
