#' Plots the extinctions history of a network
#'
#' It takes a NetworkTopology class object and plots the network index
#' after every extinction
#'
#' @param History a NetworkTopology object obtained from the Mostconnected function
#' or the ExtinctionOrder function
#' @param Variable the variable of the NetworkTopology object that you want as a y variable
#' @return A plot of number of extinctions in the x axis vs the choosen variable in the Y axis
#' @examples
#' # If you don't specify the y variable it will plot the secondary extinctions
#' # by default
#' data("net")
#' history <- SimulateExtinctions(Network = net, Method = "Mostconnected")
#' ExtinctionPlot(History = history$sims)
#' # You can also specify the variable to be ploted in the y axis
#' ExtinctionPlot(History = history$sims, Variable = "Link_density")
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author M. Isidora Ávila-Thieme <msavila@uc.cl>
#' @seealso [NetworkExtintion::ExtinctionOrder()]
#' @export

ExtinctionPlot <- function(History, Variable = "AccSecExt"){
  History$X <- 1:nrow(History)
  ggplot(History, aes_string(x = "X", y = Variable)) + geom_line() + theme_bw() + ylab(Variable) + xlab("Primary extinctions")
}
