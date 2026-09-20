if (nzchar(Sys.getenv("ALSCL_TEST_CACHE"))) options(ALSCL.tmb.cache = Sys.getenv("ALSCL_TEST_CACHE"))
fixture_frames <- function() {
  labels <- c("0-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60")
  data <- data.frame(LengthBin=labels, matrix(seq_len(45)+10,9,5), check.names=FALSE)
  names(data) <- c("LengthBin",as.character(2000:2004))
  w <- m <- data; w[-1] <- 1; m[-1] <- .5
  list(data.CatL=data,data.wgt=w,data.mat=m,rec.age=1,nage=4,M=.2,sel_L50=20,sel_L95=30)
}
cpp_fixture <- function(model="ALSCL", cv=.2, quarterly=TRUE) {
  sp <- initialize_params(species="tuna")
  L<-length(sp$len_mid); A<-6L; Y<-4L
  p<-create_parameters(model_type=tolower(model),species="tuna")$parameters
  p$log_cv_len<-log(cv);p$mean_log_F<-log(.2)
  p$dev_log_R<-rep(0,Y);p$dev_log_F<-array(0,c(if(model=="ACL")A else L,Y));p$dev_log_N0<-rep(0,A-1L)
  d<-list(logN_at_len=matrix(log(100),L,Y),na_matrix=matrix(1,L,Y),log_q=rep(0,L),
    len_border=seq(15,120,5),age=if(quarterly)seq(.25,1.5,.25) else 1:A,Y=Y,A=A,L=L,
    weight=matrix(seq_len(L)/10,L,Y),mat=matrix(.5,L,Y),M=.2,
    len_mid=sp$len_mid,len_lower=c(-Inf,seq(15,120,5)),len_upper=c(seq(15,120,5),Inf),growth_step=if(quarterly).25 else 1)
  info<-ALSCL:::.acl_compile(model)
  obj<-TMB::MakeADFun(d,p,DLL=info$dll_name,silent=TRUE)
  list(obj=obj,data=d,parameters=p,info=info)
}
fit_fixture <- function(model="ACL", ncores=1) {
  f<-fixture_frames(); f$data.CatL[[1]]<-as.character(seq(15,55,5));f$data.wgt[[1]]<-f$data.mat[[1]]<-f$data.CatL[[1]]
  p<-create_parameters(model_type=tolower(model))$parameters
  p$log_Linf<-log(60);p$log_vbk<-log(.2);p$mean_log_R<-log(100);p$log_init_Z<-log(.5)
  if(model=="ACL") {p$log_std_index<-log(.2);p$log_std_log_N0<-log(.1);p$log_std_log_R<-log(.1)} else {
    p$log_sigma_index<-log(.2);p$log_sigma_log_N0<-log(.1);p$log_sigma_log_R<-log(.1)
  }
  prep<-do.call(ALSCL:::.acl_prepare_data,f); d<-prep$data
  if(model=="ALSCL") {d$len_mid<-prep$bins$mid;d$len_lower<-prep$bins$lower;d$len_upper<-prep$bins$upper;d$growth_step<-1}
  pars<-p;pars$dev_log_R<-rep(0,d$Y);pars$dev_log_F<-array(0,c(if(model=="ACL")d$A else d$L,d$Y));pars$dev_log_N0<-rep(0,d$A-1L)
  info<-ALSCL:::.acl_compile(model)
  truth<-TMB::MakeADFun(d,pars,DLL=info$dll_name,silent=TRUE)
  f$data.CatL[-1]<-exp(truth$report(truth$par)$Elog_index)
  mapping<-lapply(p,function(x)factor(NA));mapping$mean_log_R<-NULL
  args<-c(f,list(parameters=p,map=mapping,ncores=ncores,silent=TRUE))
  result<-do.call(if(model=="ACL")run_acl else run_alscl,args)
  list(result=result,input=f,parameters=p,map=mapping)
}
