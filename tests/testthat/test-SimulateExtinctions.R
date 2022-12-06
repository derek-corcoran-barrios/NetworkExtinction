test_that("SimulateExtinction works", {
  skip_on_cran()
  capture_warnings(DF <- SimulateExtinctions(Network = NetworkExtinction::net, Method = "Mostconnected",
                                             clust.method = "cluster_infomap")$sims)
  expect_s3_class(DF, "data.frame")
  capture_warnings(DF1 <- SimulateExtinctions(Network = NetworkExtinction::net, Method = "Ordered",
                                             Order = 1:8,
                                             clust.method = "cluster_infomap")$sims)
  expect_s3_class(DF1, "data.frame")
  capture_warnings(DF2 <- SimulateExtinctions(Network = NetworkExtinction::net, Method = "Ordered",
                                              Order = 1:8,
                                              clust.method = "cluster_edge_betweenness")$sims)
  expect_s3_class(DF2, "data.frame")
  capture_warnings(DF3 <- SimulateExtinctions(Network = NetworkExtinction::net, Method = "Ordered",
                                              Order = 1:8,
                                              clust.method = "cluster_label_prop")$sims)
  expect_s3_class(DF3, "data.frame")
  capture_warnings(Mostconnected_Thresh <- SimulateExtinctions(Network = chilean_intertidal, Method = "Mostconnected",IS = 0.5))
  expect_s3_class(Mostconnected_Thresh$sims, "data.frame")

  data(chilean_intertidal)
  data(chilean_potential)
  capture_warnings(Mostconnected_Rewiring <- SimulateExtinctions(Network = chilean_intertidal, Method = "Mostconnected", Rewiring = function(x){x}, RewiringDist = chilean_potential, RewiringProb = 0.5))
  expect_s3_class(Mostconnected_Rewiring$sims, "data.frame")

})
