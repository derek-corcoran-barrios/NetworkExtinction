#' Extnction analyses from more to less conected
#'
#' It takes a network and it calculates wich species is the most conected
#' of the network, then it extinguishes that species, and calculates the
#' secundary extintions. After that it recalculates the conections and
#' repeats the process until every species becomes extinct.
#'
#' @param Network a trophic network of class network
#' @return exports data frame with the characteristics of the network after every
#' extintion
#' @examples
#' data("net")
#' Mostconnected(Network = net)
#' @importFrom network as.matrix.network.edgelist
#' @importFrom network delete.vertices
#' @importFrom network network.edgecount
#' @importFrom network network.size
#' @importFrom sna degree
#' @importFrom stats complete.cases
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Isidora Avila <msavila@uc.cl>
#' @export

Mostconnected <- function(Network){
  Grado <- NULL
  edgelist <- as.matrix.network.edgelist(Network,matrix.type="edgelist") #Prey - Predator
  Conected <- data.frame(ID = 1:network.size(Network), Grado = degree(edgelist, c("total")))
  Conected <- arrange(Conected, desc(Grado))
  Conected1<- c(Conected$ID)
  indegreebasenet <- degree(Network, cmode = "indegree")
  indegreebasenetzeros <- sum(degree(Network, cmode = "indegree") == 0)
  Producers <- (1:length(degree(Network, cmode = "indegree")))[degree(Network, cmode = "indegree") == 0]
  DF <- data.frame(Spp = rep(NA, network.size(Network)), nodesS = rep(NA, network.size(Network)), linksS = rep(NA, network.size(Network)),  indegreecero = rep(NA,network.size(Network)))

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
    SecundaryextTemp <- (1:length(degree(Temp, cmode = "indegree")))[degree(Temp, cmode = "indegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      SecundaryextTemp <- ifelse(SecundaryextTemp < j, SecundaryextTemp, SecundaryextTemp + 1)
    }
    Secundaryext <- SecundaryextTemp
    Secundaryext <- Secundaryext[!(Secundaryext %in% Producers)]
    DF$indegreecero[i]<- length(Secundaryext)
    print(i)
    FinalExt[[i]] <-(Secundaryext)
    accExt <- append(accExt, DF$Spp[1:i])
    accExt <- unique(append(accExt,Secundaryext))

    if (DF$linksS[i] == 0) break
  }
  DF <- DF[complete.cases(DF),]
  class(DF) <- c("NetworkTopology", "data.frame")
  return(DF)
}

#' Extnction analyses from more to less conected in a network
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
#' MostconnectedExp(Network = net)
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
#' @seealso [NetworkExtintion::ExtinctionOrder()]
#' @export


MostconnectedExp <- function(Network){
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
#' data("net")
#' ExtinctionOrder(Network = net, Order = c(1,2,3,4,5,6,7,8,9,10))
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 ggplot
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
  DF <- data.frame(Spp = rep(NA, network.size(Network)), nodesS = rep(NA, network.size(Network)), linksS = rep(NA, network.size(Network)),  indegreecero = rep(NA,network.size(Network)))

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
    SecundaryextTemp <- (1:length(degree(Temp, cmode = "indegree")))[degree(Temp, cmode = "indegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      SecundaryextTemp <- ifelse(SecundaryextTemp < j, SecundaryextTemp, SecundaryextTemp + 1)
    }
    Secundaryext <- SecundaryextTemp
    Secundaryext <- Secundaryext[!(Secundaryext %in% Producers)]
    DF$indegreecero[i]<- length(Secundaryext)
    print(i)
    FinalExt[[i]] <-(Secundaryext)
    accExt <- append(accExt, DF$Spp[1:i])
    accExt <- unique(append(accExt,Secundaryext))

    if (DF$linksS[i] == 0) break
  }
  DF <- DF[complete.cases(DF),]
  DF$Cumsec <- cumsum(DF$indegreecero)
  DF$NumExt <- 1:nrow(DF)
  DF$TotalExt <- DF$Cumsec + DF$NumExt
  G <- ggplot(DF, aes_string(x = "NumExt", y = "Cumsec")) + geom_line() + ylab("Secondary extinctions") + xlab("number of exctinctions")
  return(list(DF= DF, Graph = G))
}
