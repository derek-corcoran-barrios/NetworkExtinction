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
#' @param scale a character stating if the graph is on a log-log scale
#' ("LogLog") or arithmetic scale ("arithmetic"), defaults to arithmetic
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
#'data("chilean_intertidal")
#'DegreeDistribution(chilean_intertidal)
#'
#'@importFrom sna degree
#'@importFrom stats nls
#'@importFrom broom augment
#'@importFrom broom glance
#'@importFrom dplyr arrange
#'@importFrom dplyr bind_rows
#'@importFrom dplyr case_when
#'@importFrom dplyr filter
#'@importFrom dplyr full_join
#'@importFrom dplyr group_split
#'@importFrom dplyr mutate
#'@importFrom dplyr select
#'@importFrom tidyr pivot_longer
#'@importFrom magrittr "%>%"
#'@importFrom ggplot2 aes
#'@importFrom ggplot2 geom_line
#'@importFrom ggplot2 geom_point
#'@importFrom ggplot2 scale_x_log10
#'@importFrom ggplot2 scale_y_log10
#'@importFrom ggplot2 theme_bw
#'@importFrom ggplot2 ylim
#'@importFrom MASS fitdistr
#'@importFrom purrr map
#'@importFrom purrr reduce
#'@importFrom purrr discard
#'@importFrom stats ks.test
#'@importFrom stats glm
#'@importFrom stats logLik
#'@importFrom stats predict
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author M.Isidora Avila Thieme <msavila@uc.cl>
#' @export


DegreeDistribution <- function(Network, scale = "arithmetic"){
  AIC <- Cumulative <- Exp <- fit <- model <- LogPower <- logLik <- BIC <- Power <- Normal.Resid <- LogExp <- family <- AICcNorm <- NULL
  Network <- .DataInit(Network)
  totaldegree<- degree(Network)
  K <- 0:max(totaldegree)
  For.Graph<- data.frame(K = K, Cumulative = NA)
  for(i in 1:length(K)){
    For.Graph$Cumulative[i] <- sum(totaldegree>K[i])/length(totaldegree)
  }
  For.Graph <- For.Graph %>% mutate(LogK = log(K), LogCum = log(Cumulative))


  #exponential model nls
  tryCatch({
    exp.model <- nls(Cumulative~exp(K*lambda+ c),start= list(lambda=0.1, c = 0), data = For.Graph)
    For.Graph$Exp <- predict(exp.model)
    Summs.exp <- glance(exp.model)
    Summs.exp$model <- "Exp"
    Summs.exp$Normal.Resid <- ifelse(tidy(ks.test(augment(exp.model)$.resid,y='pnorm',alternative='two.sided'))$p.value < 0.05, "No", "Yes")
    Summs.exp$family <- "Exponential"
    Summs.exp$AICcNorm <- logLik(MASS::fitdistr(augment(exp.model)$.resid, "normal"))[1]
    Summs.exp$AICcNorm <- (4 - 2*Summs.exp$AICcNorm) + (12/(nrow(augment(exp.model)) - 1))
    Params.exp <- tidy(exp.model)
    Params.exp$model <- "Exp"

  }, error = function(e) {
    # display error message when there is an error
    message("Could not fit Exponential model: ", e$message)
  })
  if(!exists("Summs.exp")) {
    Summs.exp <- NULL
  }
  if(!exists("Params.exp")) {
    Params.exp <- NULL
  }

  #exponential model Log
  tryCatch({
    power <- filter(For.Graph, K != 0 & Cumulative != 0)
    logexp.model <- glm(LogCum ~ K, data = power)
    power$LogExp <- exp(predict(logexp.model))
    Summs.logexp <- glance(logexp.model)
    Summs.logexp$model <- "LogExp"
    Summs.logexp$Normal.Resid <- ifelse(tidy(ks.test(augment(logexp.model)$.resid,y='pnorm',alternative='two.sided'))$p.value < 0.05, "No", "Yes")
    Summs.logexp$family <- "Exponential"
    Summs.logexp$AICcNorm <- logLik(MASS::fitdistr(augment(logexp.model)$.resid, "normal"))[1]
    Summs.logexp$AICcNorm <- (4 - 2*Summs.logexp$AICcNorm) + (12/(nrow(augment(logexp.model)) - 1))
    Params.logexp <- tidy(logexp.model)
    Params.logexp$model <- "LogExp"}, error = function(e) {
      # display error message when there is an error
      message("Could not fit Logexponential model: ", e$message)
    })
  if(!exists("Summs.logexp")) {
    Summs.logexp <- NULL
  }
  if(!exists("Params.logexp")) {
    Params.logexp <- NULL
  }

  #logpowerlaw
  tryCatch({
    logpower.model <- glm(LogCum ~ I(log(K)), data = power)
    power$LogPower <- exp(predict(logpower.model))
    For.Graph <- full_join(For.Graph, power)
    Summs.logpower <- glance(logpower.model)
    Summs.logpower$model <- "LogPower"
    Summs.logpower$Normal.Resid <- ifelse(tidy(ks.test(augment(logpower.model)$.resid,y='pnorm',alternative='two.sided'))$p.value < 0.05, "No", "Yes")
    Summs.logpower$family <- "PowerLaw"
    Summs.logpower$AICcNorm <- logLik(MASS::fitdistr(augment(logpower.model)$.resid, "normal"))[1]
    Summs.logpower$AICcNorm <- (4 - 2*Summs.logpower$AICcNorm) + (12/(nrow(augment(logpower.model)) - 1))
    Params.logpower <- tidy(logpower.model)
    Params.logpower$model <- "LogPower"}, error = function(e) {
      # display error message when there is an error
      message("Could not fit Logpower model: ", e$message)
    })
  if(!exists("Summs.logpower")) {
    Summs.logpower <- NULL
  }
  if(!exists("Params.logpower")) {
    Params.logpower <- NULL
  }

  #powerlaw
  tryCatch({
    power.model <- nls(Cumulative~a*K^y, start= list(y=0, a = 1), data = power)
    power$Power <- predict(power.model)
    For.Graph <- full_join(For.Graph, power)
    Summs.power <- glance(power.model)
    Summs.power$model <- "Power"
    Summs.power$Normal.Resid <- ifelse(tidy(ks.test(augment(power.model)$.resid,y='pnorm',alternative='two.sided'))$p.value < 0.05, "No", "Yes")
    Summs.power$family <- "PowerLaw"
    Summs.power$AICcNorm <- logLik(MASS::fitdistr(augment(power.model)$.resid, "normal"))[1]
    Summs.power$AICcNorm <- (4 - 2*Summs.power$AICcNorm) + (12/(nrow(augment(power.model)) - 1))
    Params.power <- tidy(power.model)
    Params.power$model <- "Power"

  }, error = function(e) {
    # display error message when there is an error
    message("Could not fit Power model: ", e$message)
  })
  if(!exists("Summs.power")) {
    Summs.power <- NULL
  }
  if(!exists("Params.power")) {
    Params.power <- NULL
  }


  #all together
  Summs <- list(Summs.exp, Summs.power, Summs.logexp, Summs.logpower) %>%
    purrr::discard(is.null) %>%
    purrr::reduce(full_join) %>%
    select(logLik, AIC, BIC, model, Normal.Resid, family, AICcNorm) %>%      arrange(Normal.Resid, AIC) %>%
    dplyr::select(logLik, AIC, BIC, model, Normal.Resid, family) %>%
    dplyr::filter(AIC > -Inf, AIC < Inf)


  params <- list(Params.logpower, Params.power, Params.logexp, Params.exp) %>% purrr::discard(is.null) %>%
    purrr::reduce(bind_rows) %>%
    dplyr::filter(model %in% Summs$model) %>%
    mutate(term = case_when(term == "y" ~ "Beta",
                            term == "a" ~ "c",
                            term == "(Intercept)" ~ "c",
                            term == "I(log(K))" ~ "Beta",
                            term == "lambda" ~ "Lambda",
                            term == "K" ~ "Lambda",
                            TRUE ~ term))
  DF2 <- For.Graph %>%
    filter(K != 0 & Cumulative != 0) %>%
    tidyr::pivot_longer(
      cols = c(Exp, Power, LogExp, LogPower),
      names_to = "model",
      values_to = "fit"
    ) %>% dplyr::filter(model %in% Summs$model)

  g <- ggplot(DF2, aes(x = K, y = Cumulative)) + geom_line() + geom_point()+ theme_bw() + geom_line(aes(y =fit, color = model)) + ylim(c(0,1))

  if(scale == "LogLog"){
    g <- g  + scale_x_log10() + scale_y_log10(breaks=c(0, .001,.01,1))
  }

  g

  return(list(DDvalues = For.Graph, models = Summs, graph = g, params = params))
}


