

test_that("ExtinctionPlot works", {
  data("net")
  expect_warning(history <- SimulateExtinctions(Network = net, Method = "Mostconnected"))
  expect_s3_class(ExtinctionPlot(History = history$sims), "ggplot")
})
