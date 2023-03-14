#' A toymodel trophic network
#'
#' A trophic network with 10 species where the first four species are primery producters
#'
#' @format a network
#'
"net"

#' A toymodel distance matrix
#'
#' A distance matrix used for demonstration of rewiring capabilities
#'
#' @format a distance matrix
#'
"dist"

#' The binaryfoodweb of the intertidal zone in central chile
#'
#' A trophic network with 107 species present in the  intertidal zone of central Chile.
#' The food web was reconstructed from the Kefi et al. 2015
#
#' @format a network
#' @references Kefi, Sonia, Eric L. Berlow, Evie A. Wieters, Lucas N. Joppa, Spencer A. Wood, Ulrich Brose, and Sergio A. Navarrete. "Network structure beyond food webs: mapping non trophic and trophic interactions on Chilean rocky shores." Ecology 96, no. 1 (2015.
"chilean_intertidal"

#' The weighted foodweb of the intertidal zone in central chile
#'
#' A trophic network with 107 species present in the  intertidal zone of central Chile.
#' The food web was reconstructed from the Kefi et al. 2015
#
#' @format a network
#' @references Kefi, Sonia, Eric L. Berlow, Evie A. Wieters, Lucas N. Joppa, Spencer A. Wood, Ulrich Brose, and Sergio A. Navarrete. "Network structure beyond food webs: mapping non trophic and trophic interactions on Chilean rocky shores." Ecology 96, no. 1 (2015.
"chilean_weighted"

#' The potential foodweb of the intertidal zone in central chile
#'
#' A trophic network with 107 species present in the  intertidal zone of central Chile.
#' The food web was reconstructed from the Kefi et al. 2015
#
#' @format a network
"chilean_potential"

#' A densely connected foodweb
#'
#' A trophic network with 30 species and 222 trophic interactions.
#' This foodweb has a connectance of 0.3
#
#' @format a network
#' @seealso \code{\link{Less_Connected}}
"More_Connected"

#' A sparsely connected foodweb
#'
#' A network with 30 species and 47 interactions.
#' This network has a connectance of 0.03
#' @format a network
#' @seealso \code{\link{More_Connected}}
"Less_Connected"

#' A mutualistic web
#'
#' A network with 10 species (5 basal and 5 of higher order)
#'
#' @format a network
"mutual"
