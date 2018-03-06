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
##To import the matrix from excel:

#library(readxl)
#MatrixB <- read_excel(HERE YOU WRITE THE DIRECTORY WHERE YOU EXCEL MATRIX ARE)
#View(MatrixB)

##Once the matrix was already imported, you need to build a network object, to do this:

#library(network)
#net <- as.network(MatrixB, loops = TRUE)
#net 

## ------------------------------------------------------------------------
#To create the matrix a
a<- matrix(c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), nrow=10, ncol=10)
a

#Note that, the previous matrix a, the consumers are in the rows and the resources in the columns. For this reason, you need to  transpose the matrix a. 

#To transpose, this is the form:
a<- t(a) #transpose matrix a
a

#Once the matrix is ready, you need to build a network object, to do this:
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
RandomExtinctions(Network= net, nsim= 50)

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
data(net)
history <- Mostconnected(Network = net)
ExtinctionPlot(History = history)

# To specify the variable to be ploted in the y axis, you need to
# add the name of the variable that you want plot
ExtinctionPlot(History = history, Variable = "LinksPerSpecies")


## ------------------------------------------------------------------------
data("chilean_intertidal")
degree_distribution(chilean_intertidal, name = "Test")

