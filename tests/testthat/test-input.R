test_that("midpoints preserve unequal widths and decimals", {
  f<-fixture_frames();x<-do.call(ALSCL:::.acl_prepare_data,f)
  expect_equal(x$bins$mid,c(10,seq(22.5,57.5,5)))
  b<-ALSCL:::.acl_length_bins(c("<20","20-25",">25"))
  expect_equal(b$mid,c(17.5,22.5,27.5));expect_equal(b$lower,c(-Inf,20,25));expect_equal(b$upper,c(20,25,Inf))
  b<-ALSCL:::.acl_length_bins(as.character(c(10,20,35)),len_mid=c(11,21,36),len_border=c(15,27.5))
  expect_equal(b$mid,c(11,21,36))
  expect_error(ALSCL:::.acl_length_bins(c("0-10","11-20")),"contiguous")
  expect_error(ALSCL:::.acl_length_bins(c("20","10")),"increasing")
})
test_that("quarterly ages agree with simulation and missing cells are safe", {
  f<-fixture_frames();f$rec.age<-.25;f$nage<-20L;f$growth_step<-.25
  f$data.CatL[1,2]<-NA_real_;f$data.CatL[2,2]<-0
  x<-do.call(ALSCL:::.acl_prepare_data,f)
  expect_equal(x$data$age,initialize_params(species="tuna")$ages)
  expect_equal(x$data$na_matrix[1:2,1],c(0,0));expect_true(all(is.finite(x$data$logN_at_len)))
  expect_error(do.call(ALSCL:::.acl_prepare_data,c(f,list(zero_action="error"))),"Zero")
  f$data.CatL[1,2]<- -1;expect_error(do.call(ALSCL:::.acl_prepare_data,f),"nonnegative")
  f<-fixture_frames();f$data.wgt<-f$data.wgt[9:1,];expect_error(do.call(ALSCL:::.acl_prepare_data,f),"match")
})
test_that("legacy parameter calls, named bounds and maps are aligned", {
  p<-create_parameters(NULL,list(log_vbk=log(.15)),list(log_vbk=log(.3)))
  expect_equal(p$parameters.L$log_vbk,log(.15));expect_equal(p$parameters.U$log_vbk,log(.3))
  expect_equal(create_parameters(list(mean_log_R=6))$parameters$mean_log_R,6)
  expect_equal(create_parameters(parameters.L=list(log_vbk=log(.15)))$parameters$log_vbk,log(.2))
  mapping<-generate_map(list(log_Linf=factor(NA)))
  expect_setequal(names(mapping),c("log_std_log_F","logit_log_F_y","logit_log_F_a","t0","log_Linf"))
  expect_false("t0" %in% names(generate_map(list(t0=NULL))))
  free<-setdiff(names(p$parameters),names(mapping));v<-unlist(p$parameters)[free]
  b<-ALSCL:::.acl_bounds(v,p)
  expect_identical(names(b$lower),names(v));expect_equal(unname(b$lower["log_vbk"]),log(.15))
  expect_error(create_parameters(parameters=list(unknown=1)),"known")
})
test_that("one requested optimization pass performs one call", {
  count<-0L
  local_mocked_bindings(nlminb=function(...) {count<<-count+1L;list(par=1,objective=0)},.package="stats")
  ALSCL:::.acl_optimize(list(fn=function(x)0,gr=function(x)0),1,-Inf,Inf,1L,list())
  expect_equal(count,1L)
  ALSCL:::.acl_optimize(list(fn=function(x)0,gr=function(x)0),1,-Inf,Inf,3L,list())
  expect_equal(count,4L)
})
test_that("empty overrides preserve defaults and unnamed ALSCL maps are rejected", {
  expect_equal(generate_map(list()),generate_map())
  expect_equal(create_parameters(parameters=list()),create_parameters())
  expect_error(ALSCL:::.acl_map(create_parameters_alscl()$parameters,list(factor(NA)),"ALSCL"),"named list")
})
