test_that("RandomExtinctions works", {
  skip_on_cran()
  capture_warnings(DF <- RandomExtinctions(Network = More_Connected, nsim = 20)$sims)
  expect_s3_class(DF, "data.frame")
})
