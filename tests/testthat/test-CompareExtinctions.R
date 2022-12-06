test_that("CompareExtinctions works", {
  skip_on_cran()
  data("Less_Connected")
  expect_warning(History <- SimulateExtinctions(Network = Less_Connected, Method = "Mostconnected"))
  capture_warnings(NullHyp <- RandomExtinctions(Network = Less_Connected, nsim = 100))
  capture_warnings(Compare <- CompareExtinctions(Nullmodel = NullHyp, Hypothesis = History))
  expect_s3_class(Compare, "ggplot")
})
