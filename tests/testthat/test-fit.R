test_that("both public fitting APIs return finite usable results and consistent diagnostics", {
  for(model in c("ACL","ALSCL")) {
    x<-fit_fixture(model);m<-x$result
    expect_true(is.finite(m$opt$objective));expect_equal(m$convergence_code,0L)
    expect_lt(m$max_abs_gradient,1e-3)
    expect_true(isTRUE(m$pdHess));expect_true(all(is.finite(m$est_std)))
    expect_equal(m$obj$fn(m$opt$par),m$opt$objective,tolerance=1e-8)
    expect_identical(rownames(m$par_low_up),names(m$opt$par))
    expect_equal(nrow(diagnose_model(x$input$data.CatL,m)),16L)
    for(fun in list(plot_CatL,plot_recruitment,plot_SSB,plot_biomass,plot_abundance,
                    plot_catch,plot_fishing_mortality,plot_residuals,plot_VB,plot_pla)) {
      expect_s3_class(ggplot2::ggplot_build(fun(m)),"ggplot_built")
    }
  }
})
test_that("parallel starts agree with serial fits without jittering fixed effects", {
  x<-fit_fixture("ACL",ncores=2);y<-fit_fixture("ACL",ncores=1)
  expect_equal(x$result$opt$objective,y$result$opt$objective,tolerance=1e-6)
  expect_equal(x$result$report$Linf,60,tolerance=1e-12)
  expect_length(x$result$starts,2)
  x<-fit_fixture("ALSCL",ncores=2);y<-fit_fixture("ALSCL",ncores=1)
  expect_equal(x$result$opt$objective,y$result$opt$objective,tolerance=1e-6)
})
test_that("retrospective wrappers preserve positional parameters and socket results", {
  x <- fit_fixture("ACL")
  args <- c(list(nyear=2),x$input,list(parameters=x$parameters,map=x$map,silent=TRUE))
  serial <- do.call(retro_acl,args)
  parallel <- do.call(retro_acl,c(args,list(ncores=2)))
  expect_equal(parallel$rho_text,serial$rho_text,tolerance=1e-6)
  positional <- do.call(retro_acl,c(list(2),unname(x$input),list(x$parameters),
                                  list(map=x$map,silent=TRUE)))
  expect_equal(positional$rho_text,serial$rho_text,tolerance=1e-6)
})
test_that("sim_acl reads saved simulations and honors bounds-only calls", {
  x <- fit_fixture("ACL");f <- x$input;m <- x$result
  path <- tempfile("simulation-");dir.create(path)
  sim.data <- list(SN_at_len=t(as.matrix(f$data.CatL[-1])),nyear=length(m$year),
    nage=length(m$age),ages=m$age,growth_step=1,len_mid=m$len_mid,len_border=m$len_border,
    weight=as.matrix(f$data.wgt[-1]),mat=as.matrix(f$data.mat[-1]),
    q_surv_L50=f$sel_L50,q_surv_L95=f$sel_L95)
  save(sim.data,file=file.path(path,"sim_rep4"))
  fitted <- sim_acl(4,path,path,parameters=x$parameters,map=x$map,M=f$M)[[1]]
  expect_null(fitted$error)
  expect_equal(fitted$opt$objective,m$opt$objective,tolerance=1e-6)
  fixed <- x$parameters;fixed$mean_log_R <- NULL
  testthat::local_mocked_bindings(create_parameters=function(model_type, species=NULL,
      parameters=NULL,parameters.L=NULL,parameters.U=NULL) {
    expect_null(parameters);expect_equal(parameters.L,list(mean_log_R=0))
    stop("bounds reached by name")
  },.package="ALSCL")
  fitted <- sim_acl(4,path,path,parameters.L=list(mean_log_R=0),M=f$M)[[1]]
  expect_match(fitted$error,"bounds reached by name")
})
