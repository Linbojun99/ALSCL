test_that("ggsci settings propagate and explicit overrides survive", {
  old <- options(acl.theme = NULL); on.exit(options(old))
  expect_identical(acl_theme("palette"), "npg")
  expect_equal(acl_theme("line_color"), ggsci::pal_npg()(10)[4])
  f <- list(report=list(Rec=1:4), year=2000:2003)
  expect_equal(ggplot2::ggplot_build(plot_recruitment(f))$data[[1]]$colour,
               rep(acl_theme("line_color"),4))
  acl_theme_set(palette="jama", line_color="purple")
  expect_identical(acl_theme("line_color"), "purple")
  expect_equal(acl_theme("compare_colors"), ggsci::pal_jama()(7)[1:2])
  expect_true(all(ggplot2::ggplot_build(plot_recruitment(f,line_color="orange"))$data[[1]]$colour=="orange"))
  expect_identical(acl_theme_set()$palette, "jama")
  expect_error(acl_theme_set(palette="typo"), "arg")
  acl_theme_reset()
  expect_identical(acl_theme("palette"), "npg")
})

test_that("retrospective colors follow terminal periods and allow overrides", {
  old <- options(acl.theme = NULL); on.exit(options(old))
  d <- expand.grid(Year=2000:2002, RetrospectiveYear=c("2019","2018","2017"),Variable="SSB")
  d$Value <- seq_len(nrow(d))
  r <- list(results=d,last_points=d[d$Year==2002,],rho_text=data.frame(Variable="SSB",Rho=.01))
  p <- plot_retro(r, palette="npg", ylab="SSB units")
  expect_equal(ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")$map(c("2019","2018","2017")), ggsci::pal_npg()(3))
  expect_equal(p$labels$y,"SSB units")
  r$results <- d[nrow(d):1,]
  expect_equal(ggplot2::ggplot_build(plot_retro(r))$plot$scales$get_scales("colour")$map(c("2017","2019")), ggsci::pal_npg()(3)[c(3,1)])
  p <- plot_retro(r,colors=c("2017"="black","2019"="red","2018"="blue"))
  expect_equal(ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")$map(c("2019","2018","2017")),c("red","blue","black"))
  expect_error(plot_retro(r,colors="red"),"every retrospective period")
  expect_error(plot_retro(r,colors=c(a="red",b="blue",c="green")),"Named colors")
  expect_length(ALSCL:::.acl_palette("npg",15),15)
  expect_false(anyNA(ALSCL:::.acl_palette("npg",15)))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("continuous and ridge palettes retain legacy options", {
  old <- options(acl.theme = NULL); on.exit(options(old))
  for (pal in c("npg","aaas","nejm","lancet","jama")) {
    acl_theme_set(palette=pal)
    expect_no_error(ALSCL:::.acl_ridges_palette(NULL,20))
    expect_no_error(ALSCL:::.acl_fill_continuous(pal))
  }
  expect_no_error(ALSCL:::.acl_fill_continuous("inferno"))
  expect_no_error(ALSCL:::.acl_ridges_palette("viridis",20))
  expect_no_error(ALSCL:::.acl_ridges_palette(c("navy","white"),20))
})


test_that("ridges keep ordered gradient colors independently of categorical themes", {
  old <- options(acl.theme = NULL); on.exit(options(old))
  years <- c(2, 10, 11)
  survey <- matrix(log(seq_len(12)), nrow=4)
  f <- list(report=list(Elog_index=survey),
            obj=list(env=list(.data=list(logN_at_len=survey))),
            len_mid=seq(10,40,10), year=years)
  expected <- ggplot2::scale_fill_viridis_d()$palette(3)
  # The palette must remain sequential for both implicit and explicit NULL calls.
  for (pal in c("npg", "jama")) {
    acl_theme_set(palette=pal)
    expect_equal(ALSCL:::.acl_ridges_palette(NULL,3)$palette(3), expected)
    expect_equal(ALSCL:::.acl_ridges_palette("viridis",3)$palette(3), expected)
    expect_no_warning(ggplot2::ggplot_build(plot_ridges(f)))
  }
  expect_identical(formals(plot_ridges)$palette, "viridis")
  expect_equal(ALSCL:::.acl_ridges_palette("cividis",3)$palette(3),
               ggplot2::scale_fill_viridis_d(option="E")$palette(3))
  expect_equal(ALSCL:::.acl_ridges_palette(c("#EEEEEE","#003366"),3)$palette(3),
               c("#EEEEEE","#7790AA","#003366"))
})
