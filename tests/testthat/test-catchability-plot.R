test_that("catchability plots connect inputs without smoothing overshoot", {
  q <- c(0, 0.001, 0.01, 0.5, 0.99, 1)
  m <- list(model_type = "ACL", len_mid = seq_along(q),
            obj = list(env = list(data = list(log_q = log(q)))))
  p <- plot_compare_selectivity(m, m, model1_name = "one", model2_name = "two")
  expect_true(inherits(p$layers[[1]]$geom, "GeomLine"))
  plotted <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(sort(plotted$y), sort(rep(q, 2)))
  expect_true(all(plotted$y >= 0 & plotted$y <= 1))
})
