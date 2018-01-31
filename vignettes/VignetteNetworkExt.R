## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(NetworkExtintion)
data("net")
MostconnectedExp(Network = net)

## ------------------------------------------------------------------------
data("net")
ExtinctionOrder(Network = net, Order = c(2,4,7))

## ------------------------------------------------------------------------
data("net")
degree_distribution(net, name = "Test")

## ---- fig.cap="Fig 1. This graph shows something"------------------------
plot(1:10)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

