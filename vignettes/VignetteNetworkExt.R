## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE---------------------------------------------------------
#  install.packages(NetworkExtinction)
#  library(NetworkExtinction)

## ------------------------------------------------------------------------
#create a matrix
a<- matrix(c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), nrow=10, ncol=10)

#the matrix of the data(net) used below is the transpose of a
a<- t(a) #transpose matrix a
a
#create a network
library(network)
net <- as.network(a, loops = TRUE)
net

## ------------------------------------------------------------------------
library(NetworkExtinction)
data("net")
MostconnectedExp(Network = net)

## ------------------------------------------------------------------------
data("net")
history <- MostconnectedExp(Network = net)
ExtinctionPlot(History = history, Variable = "AccSecondaryExtinction")

#"Fig 1. The graph shows the number of accumulated secondary extinctions that occur when removing species from the most connected to the least connected"

## ------------------------------------------------------------------------
data("net")
ExtinctionOrder(Network = net, Order = c(2,4,7))

## ------------------------------------------------------------------------
data(net)
RandomExtinctions(Network= net, nsim= 5)

## ------------------------------------------------------------------------
data("net")
degree_distribution(net, name = "Test")

## ---- fig.cap="Fig 1. This graph shows something"------------------------
plot(1:10)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

