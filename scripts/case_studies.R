# 两个完整模拟案例，重新拟合并绘图 / Refit and plot both simulation cases.
# 从仓库根目录运行 / Run from the repository root: Rscript scripts/case_studies.R
library(ALSCL)
library(ggplot2)
library(patchwork)
source('scripts/simulation_workflow.R', encoding='UTF-8')
run_case_studies <- function(out='docs') {
  cache <- Sys.getenv('ALSCL_MANUAL_CACHE')
  if(nzchar(cache)) options(ALSCL.tmb.cache=cache)
  old <- acl_theme(); on.exit(do.call(acl_theme_set,old),add=TRUE)
  acl_theme_set(palette='npg',base_theme='theme_bw')
  cols <- setNames(acl_theme('compare_colors'),c('ACL','ALSCL'))
  fig <- file.path(out,'figures','cases'); dir.create(fig,recursive=TRUE,showWarnings=FALSE)
  save_plot <- function(id,p,height=5) ggsave(file.path(fig,paste0(id,'.png')),p,width=10,height=height,dpi=160,bg='white')
  statuses <- list()
  for(sp in c('flatfish','tuna')) {
    message('CASE ',sp)
    folder <- file.path(out,'data','cases',sp)
    z <- readRDS(file.path(folder,'truth.rds')); pa<-z$parameters; sm<-z$truth
    inp <- simulation_to_tables(sm,pa); bio <- sim_cal(pa)
    # 固定已知生物参数，保留条件估计 / Conditional estimation with known biology.
    pA <- list(log_Linf=log(pa$Linf),log_vbk=log(pa$vbk),t0=pa$t0,
      log_cv_len=log(pa$cv_L),log_std_log_F=log(pa$F_sd),
      log_std_log_N0=log(pa$std_logN0),logit_log_R=qlogis(pa$R_ar))
    pB <- list(log_Linf=log(pa$Linf),log_vbk=log(pa$vbk),log_t0=log(pa$t0),
      log_cv_len=log(pa$cv_L),log_cv_grow=log(pa$cv_inc),log_sigma_log_F=log(pa$F_sd),
      log_sigma_log_N0=log(pa$std_logN0),logit_log_R=qlogis(pa$R_ar))
    mA <- lapply(pA,function(v)factor(NA)); mB <- lapply(pB,function(v)factor(NA))
    mB$logit_log_F_l <- mB$logit_log_F_y <- factor(NA)
    common <- c(inp,list(rec.age=pa$rec.age,nage=pa$nage,M=pa$M,
      sel_L50=pa$q_surv_L50,sel_L95=pa$q_surv_L95,len_mid=pa$len_mid,
      len_border=pa$len_border[-c(1,length(pa$len_border))],growth_step=pa$growth_step,
      train_times=2,ncores=1,silent=TRUE))
    a <- do.call(run_acl,c(common,list(parameters=pA,map=mA)))
    b <- do.call(run_alscl,c(common,list(parameters=pB,map=mB)))
    status <- function(f,model) data.frame(Case=sp,Model=model,Code=f$convergence_code,
      MaxGradient=f$max_abs_gradient,pdHess=f$pdHess,Boundary=f$bound_hit)
    statuses[[sp]] <- rbind(status(a,'ACL'),status(b,'ALSCL'))
    years <- a$year
    obs <- data.frame(Year=rep(years,each=length(pa$len_mid)),
      Length=rep(pa$len_mid,length(years)),Value=as.vector(t(sm$SN_at_len)))
    save_plot(paste0(sp,'_data'),ggplot(obs,aes(Year,Length,fill=log10(Value)))+
      geom_tile()+scale_fill_gradient(low=acl_theme('low_col'),high=acl_theme('high_col'))+
      labs(title=paste(sp,'synthetic survey'),fill='log10 index')+theme_bw(),4)
    truth <- data.frame(Year=rep(years,3),Quantity=rep(c('B','SSB','Rec'),each=length(years)),
      Truth=c(sm$TB,sm$SSB,sm$Rec))
    estimates <- rbind(data.frame(Year=rep(years,3),Quantity=truth$Quantity,
      Value=c(a$report$B,a$report$SSB,a$report$Rec),Model='ACL'),
      data.frame(Year=rep(years,3),Quantity=truth$Quantity,
      Value=c(b$report$B,b$report$SSB,b$report$Rec),Model='ALSCL'))
    write.csv(truth,file.path(folder,'truth_timeseries.csv'),row.names=FALSE)
    write.csv(estimates,file.path(folder,'estimated_timeseries.csv'),row.names=FALSE)
    save_plot(paste0(sp,'_truth'),ggplot(estimates,aes(Year,Value,color=Model))+
      geom_line()+scale_color_manual(values=cols)+
      geom_line(data=truth,aes(Year,Truth),inherit.aes=FALSE,linetype=2)+
      facet_wrap(~Quantity,scales='free_y',ncol=1)+labs(y=NULL,subtitle='Dashed: generating truth; solid: conditional estimates')+theme_bw(),6)
    growth <- data.frame(Age=pa$ages,Length=VB_func(pa$Linf,pa$vbk,pa$t0,pa$ages))
    prob <- rbind(data.frame(Length=pa$len_mid,Value=mat_func(pa$mat_L50,pa$mat_L95,pa$len_mid),Type='Maturity'),
      data.frame(Length=pa$len_mid,Value=bio$q_surv,Type='Survey q'))
    save_plot(paste0(sp,'_biology'),(ggplot(growth,aes(Age,Length))+geom_line(color=cols[1])+theme_bw())/
      (ggplot(prob,aes(Length,Value,color=Type))+geom_line()+
       scale_color_manual(values=setNames(unname(cols),c('Maturity','Survey q')))+
       labs(y='Probability')+theme_bw()),6)
    if(sp=='tuna') {
      G <- b$report$G
      stopifnot(max(abs(colSums(G)-1))<1e-10,all(is.finite(G)),all(G>=0))
      g <- expand.grid(Destination=pa$len_mid,Source=pa$len_mid);g$Probability<-as.vector(G)
      save_plot('tuna_growth_transition',ggplot(g,aes(Source,Destination,fill=Probability))+
        geom_tile()+coord_equal()+scale_fill_gradient(low=acl_theme('low_col'),high=acl_theme('high_col'))+
        labs(title='ALSCL growth transition G',subtitle='Columns sum to one; no shrinkage')+theme_bw(),6)
    }
  }
  checks <- do.call(rbind,statuses);print(checks)
  write.csv(checks,file.path(out,'results','case_convergence.csv'),row.names=FALSE)
  stopifnot(all(checks$Code==0),all(checks$MaxGradient<.001),all(checks$pdHess),!any(checks$Boundary))
  invisible(checks)
}
if(sys.nframe()==0L) run_case_studies()
