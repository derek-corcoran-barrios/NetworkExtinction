## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.width=6, fig.height=4) 

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
Mostconnected(Network = net)

## ---- fig.cap="Fig 1. The graph shows the number of accumulated secondary extinctions that occur when removing species from the most connected to the least connected"----
data("net")
history <- Mostconnected(Network = net)
ExtinctionPlot(History = history, Variable = "AccSecondaryExtinction")


## ------------------------------------------------------------------------
data("net")
ExtinctionOrder(Network = net, Order = c(2,4,7))

## ---- message=FALSE------------------------------------------------------
data(net)
RandomExtinctions(Network= net, nsim= 10)

## ----message=FALSE, warning=FALSE----------------------------------------
data("net")
History <- ExtinctionOrder(Network = net, Order = c(1,2,3,4,5,6,7,8,9,10))

set.seed(2)
NullHyp <- RandomExtinctions(Network = net, nsim = 100)

Comparison <- CompareExtinctions(Nullmodel = NullHyp, Hypothesis = History)

## ------------------------------------------------------------------------
Comparison$graph

## ------------------------------------------------------------------------
Comparison$Test

## ------------------------------------------------------------------------
data("net")
degree_distribution(net, name = "Test")

