#' Degree distribution of the network
#'
#' This function calculates the degree distribution of the network. First it
#' fits exponential, power law and truncated power law distribution models,
#' and calculates the AIC values to select the best fit, and finally it plots
#' the degree distribution in a log log scale showing the three fitted models
#' mentioned above against the observed distribution.
#'
#'
#' @param Network a trophic network of class network
#' @param name a categorical variable that represent
#' the distribution model
#' @return exports three principal results:
#' 1. A list with network degree distribution values and with the value of each fit model
#' 2. A list with each model results and AIC of the distribution models
#' 3. A Ghraph of the degree distribution with the models adjust
#' In DDvalues, k represent the degree of the network and cumulative
#' the probability that each specie could be have this degree (pk).
#' Observation: In the graph, the zero values are not represented but this result are incorporate in the DF result
#'
#'@examples
#'library(NetworkExtinction)
#'data("net")
#'degree_distribution(net, name = "Test")
#'
#'@importFrom sna degree
#'@importFrom stats nls
#'@importFrom broom glance
#'@importFrom dplyr filter
#'@importFrom dplyr arrange
#'@importFrom dplyr full_join
#'@importFrom tidyr gather
#'@importFrom magrittr "%>%"
#'@importFrom ggplot2 aes_string
#'@importFrom ggplot2 geom_line
#'@importFrom ggplot2 geom_point
#'@importFrom ggplot2 scale_x_log10
#'@importFrom ggplot2 scale_y_log10
#'@importFrom ggplot2 theme_classic
#'@importFrom stats predict
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author M.Isidora Avila Thieme <msavila@uc.cl>
#' @export


degree_distribution <- function(Network, name){
  AIC <- Cumulative <- Exp <- fit <- model <- truncated <- NULL
  totaldegree<- degree(Network)
  K <- 0:max(totaldegree)
  For.Graph<- data.frame(K = K, Cumulative = NA, Scenario = name)
  for(i in 1:length(K)){
    For.Graph$Cumulative[i] <- sum(totaldegree>K[i])/length(totaldegree)
  }

  exp.model <- nls(Cumulative~exp(-K/y),start= list(y=0.1), data = For.Graph)

  For.Graph$Exp <- predict(exp.model)
  Summs.exp <- glance(exp.model)
  Summs.exp$model <- "Exponential"

  power <- filter(For.Graph, K != 0 & Cumulative != 0)
  powerlaw.model <- nls(Cumulative~K^-y, start= list(y=0), data = power)
  power$power <- predict(powerlaw.model)
  For.Graph <- full_join(For.Graph, power)
  Summs.power <- glance(powerlaw.model)
  Summs.power$model <- "Power"

  truncated.powerlaw.model <- nls(Cumulative~(K^-y)*(exp(-K/y)), start = list(y=1), data = power)
  power$truncated <- predict(truncated.powerlaw.model)
  For.Graph <- full_join(For.Graph, power)
  Summs.truncated <- glance(truncated.powerlaw.model)
  Summs.truncated$model <- "truncated"

  Summs <- full_join(Summs.exp, Summs.power)
  Summs <- full_join(Summs, Summs.truncated)
  Summs <- arrange(Summs, AIC)

  DF2 <- For.Graph %>% filter(K != 0 & Cumulative != 0) %>% gather(key = model, value = fit, Exp, power, truncated)

  g <- ggplot(DF2, aes_string(x = "K", y = "Cumulative")) + geom_line() + geom_point()+ scale_x_log10() + scale_y_log10() + theme_classic() + geom_line(aes_string(y ="fit", color = "model"))

  g

  return(list(DDvalues = For.Graph, models = Summs, graph = g))
}


