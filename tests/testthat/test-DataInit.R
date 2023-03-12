
test_that("DataInit() returns a network object", {
  # Create a matrix object
  m <- matrix(c(0,1,1,0), nrow=2, ncol=2)

  # Call DataInit()
  n <- .DataInit(m)

  # Test that the output is a network object
  expect_true(is(n)[1] == "network")
})
