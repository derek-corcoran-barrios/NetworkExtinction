test_that("DegreeDistribution works", {
  data("chilean_intertidal")
  capture_messages(DD <- DegreeDistribution(chilean_intertidal))
  expect_s3_class(DD$DDvalues, "data.frame")
})
