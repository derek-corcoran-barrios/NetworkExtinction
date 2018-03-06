## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.width=6, fig.height=4, message = FALSE) 

## ---- eval=FALSE---------------------------------------------------------
#  install.packages(NetworkExtinction)
#  library(NetworkExtinction)

## ------------------------------------------------------------------------
a<- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0),nrow=10, ncol=10)

a

## ------------------------------------------------------------------------
library(network)
net <- as.network(a, loops = TRUE)
net

## ---- eval=FALSE---------------------------------------------------------
#  library(NetworkExtinction)
#  data("net")
#  Mostconnected(Network = net)

## ---- echo=FALSE, message=FALSE------------------------------------------
library(NetworkExtinction)
data("net")
knitr::kable(Mostconnected(Network = net), caption = "Table 1: The resulting dataframe of the Mostconnected function")

## ---- fig.cap="Figure 3. The graph shows the number of accumulated secondary extinctions that occur when removing species from the most to the least connected species"----
data("net")
history <- Mostconnected(Network = net)
ExtinctionPlot(History = history, Variable = "AccSecondaryExtinction")


## ---- eval=FALSE---------------------------------------------------------
#  data("net")
#  ExtinctionOrder(Network = net, Order = c(2,4,7))

## ---- echo=FALSE---------------------------------------------------------
data("net")
knitr::kable(ExtinctionOrder(Network = net, Order = c(2,4,7))$DF, caption = "Table 2: The resulting dataframe of the ExtinctionOrder function")

## ---- echo=FALSE, fig.cap= "Figure 4. The graph shows the number of accumulated secondary extinctions that occur when removing species in a custom order. In this example species 2 is removed followed by 4 and lastly species 7 is removed"----
data("net")
ExtinctionOrder(Network = net, Order = c(2,4,7))$Graph

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

