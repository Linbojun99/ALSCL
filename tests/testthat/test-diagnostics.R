test_that("Mohn rho averages terminal-year differences over all peels", {
  f<-fixture_frames()
  local_mocked_bindings(run_acl=function(data.CatL,...) {
    n<-ncol(data.CatL)-1L;v<-rep(100,n)
    if(n<5) v[n]<-if(n==4)120 else 140
    list(year=2000+seq_len(n)-1L,report=list(B=v,Rec=v,SSB=v,N=v))
  },.package="ALSCL")
  x<-do.call(retro_model,c(list(nyear=2),f))
  expect_equal(x$rho_text$Rho,rep(.3,4))
  expect_error(do.call(retro_model,c(list(nyear=4),f)),"at least two")
})
test_that("diagnostics use maximum absolute gradients and valid observation pairs", {
  f<-fixture_frames();d<-f$data.CatL
  d[-1]<-matrix(rep(10:14,each=9),9,5)
  m<-list(report=list(Elog_index=log(as.matrix(d[-1])-2)),opt=list(par=1:2,objective=10,convergence=0L),
          bound_hit=FALSE,final_outer_mgc=c(-.4,.001,.2),pdHess=TRUE)
  x<-diagnostic_metrics(d,m)
  expect_equal(x$Value[x$Metric=="MASE"],2)
  a<-diagnose_model(d,m);expect_equal(a$Value[a$Metric=="Final Outer mgc"],.4)
  m$final_outer_mgc<-.5
  expect_equal(diagnose_model(d,m)$Value[15],.5)
  d[1,2]<-NA;d[2,2]<-0
  expect_true(all(is.finite(diagnostic_metrics(d,m)$Value)))
  d[-1]<-10;m$report$Elog_index[,]<-log(10)
  x<-diagnostic_metrics(d,m)
  expect_true(is.na(x$Value[x$Metric=="MASE"]))
  expect_true(is.na(x$Value[x$Metric=="Rsquared"]))
})
test_that("growth confidence intervals include Linf and parameter covariances", {
  v<-matrix(c(.04,.006,.006,.01),2,dimnames=list(c("log_Linf","log_vbk"),c("log_Linf","log_vbk")))
  m<-list(report=list(Linf=60,vbk=.2,t0=-.1),vcov=v,pdHess=TRUE)
  age<-c(1,5,20);x<-ALSCL:::.acl_growth_interval(m,age)
  numerical<-sapply(seq_len(2),function(j) {
    a<-b<-c(log(60),log(.2));a[j]<-a[j]+1e-5;b[j]<-b[j]-1e-5
    (exp(a[1])*(1-exp(-exp(a[2])*(age+.1)))-exp(b[1])*(1-exp(-exp(b[2])*(age+.1))))/2e-5
  })
  expect_equal(x$se,sqrt(rowSums((numerical%*%v)*numerical)),tolerance=1e-7)
  m$vcov<-NULL;expect_null(ALSCL:::.acl_growth_interval(m,age))
})
