test_that("ExtinctionOrder works", {
  capture_warnings(DF <- ExtinctionOrder(net, Order = 1:8, NetworkType = "Trophic"))
  expect_s3_class(DF$sims, "data.frame")
})
