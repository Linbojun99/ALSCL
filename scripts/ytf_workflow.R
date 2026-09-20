# GitHub 图谱：所有拟合图使用内置 YTF_example / All fitted plots use YTF_example.
# 从仓库根目录运行 / Run from the repository root:
# Rscript scripts/ytf_workflow.R
library(ALSCL)
library(ggplot2)
library(patchwork)

run_ytf_guide <- function(out = "docs", npeels = 3L) {
  cache <- Sys.getenv("ALSCL_MANUAL_CACHE")
  if (nzchar(cache)) options(ALSCL.tmb.cache = cache)
  for (d in c("figures/ytf", "results", "data"))
    dir.create(file.path(out, d), recursive=TRUE, showWarnings=FALSE)
  data("YTF_example", package="ALSCL")
  x <- YTF_example
  inputs <- x[c("data.CatL", "data.wgt", "data.mat")]
  common <- c(inputs, x$fit_args, list(train_times=2, ncores=1, silent=TRUE))
  # 固定生物参数的教学拟合 / Conditional teaching fits with known biology.
  message("Fitting YTF ACL")
  a <- do.call(run_acl, c(common, x$fit_config$acl))
  message("Fitting YTF ALSCL")
  b <- do.call(run_alscl, c(common, x$fit_config$alscl))
  dat <- inputs$data.CatL
  status <- function(f, name) data.frame(Model=name, Code=f$convergence_code,
    MaxGradient=f$max_abs_gradient, pdHess=f$pdHess, Boundary=f$bound_hit)
  convergence <- rbind(status(a,"ACL"),status(b,"ALSCL"))
  print(convergence)
  stopifnot(all(convergence$Code==0), all(convergence$MaxGradient<.001),
            all(convergence$pdHess), !any(convergence$Boundary))
  write.csv(convergence,file.path(out,"results","convergence.csv"),row.names=FALSE)
  for (name in names(inputs)) write.csv(inputs[[name]],
    file.path(out,"data",paste0(name,".csv")),row.names=FALSE)
  write.csv(diagnose_model(dat,a),file.path(out,"results","diagnostics_ACL.csv"),row.names=FALSE)
  write.csv(diagnose_model(dat,b),file.path(out,"results","diagnostics_ALSCL.csv"),row.names=FALSE)
  cmp <- compare_models(a,b,dat)
  for (nm in c("summary","growth","fit_metrics","correlation"))
    write.csv(cmp[[nm]],file.path(out,"results",paste0("comparison_",nm,".csv")),row.names=FALSE)
  # 返回图和数据的示例 / Access plotting data explicitly.
  tab <- plot_abundance(b,type="NL",return_data=TRUE)
  stopifnot(is.list(tab),!is.null(tab$plot))
  old_theme <- acl_theme()
  on.exit(do.call(acl_theme_set, old_theme),add=TRUE)
  acl_theme_set(base_theme="bw",font_family="sans",title_size=12,
    axis_title_size=10,axis_text_size=9,strip_text_size=9,
    palette="npg")
  registry <- list()
  render <- function(id,zh,en,expr,height=4.6,note_zh="",note_en="") {
    call <- substitute(expr); warnings <- character()
    message("FIGURE ",id)
    p <- withCallingHandlers(eval(call,parent.frame()), warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
    })
    withCallingHandlers(ggsave(file.path(out,"figures","ytf",paste0(id,".png")),
      p,width=11,height=height,dpi=160,bg="white"),warning=function(w){
        warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
      })
    registry[[id]] <<- list(id=id,zh=zh,en=en,
      code=paste(deparse(call,width.cutoff=90),collapse="\n"),
      note_zh=note_zh,note_en=note_en,warnings=unique(warnings),status="PASS")
    jsonlite::write_json(registry,file.path(out,"results","figure_manifest.json"),
      auto_unbox=TRUE,pretty=TRUE,null="null")
  }
for(typ in c('length','year')) for(ex in c(FALSE,TRUE)) {
 render(paste0('CatL_',typ,'_',ex),paste('调查拟合',typ,if(ex)'原尺度' else '对数尺度'),paste('Survey fit',typ,if(ex)'original scale' else 'log scale'),
  plot_CatL(b,type=typ,exp_transform=ex,facet_ncol=4,x_breaks=if(typ=='length')c(2000,2010,2015)else c(10,30,50),point_size=1,line_size=.6,ylab=if(ex)'Survey index' else 'Log survey index'),height=8.5,
  note_zh=if(ex)'原尺度曲线为 exp(Elog_index)，即对数正态中位数。'else if(typ=='length')'点为观测对数，曲线为拟合的 Elog_index。'else '朱红色为观测对数，深蓝色为拟合的 Elog_index。',note_en=if(ex)'Original-scale curves are lognormal medians, exp(Elog_index).'else if(typ=='length')'Points are log observations; curves show fitted Elog_index.'else 'Vermilion curves are log observations; navy curves show fitted Elog_index.')
}
for(fun in c('plot_abundance','plot_biomass','plot_SSB','plot_catch')) {
 types<-switch(fun,plot_abundance=c('N','NA','NL'),plot_biomass=c('B','BA','BL'),plot_SSB=c('SSB','SBA','SBL'),plot_catch=c('CN','CNA','CNL'))
 for(typ in types) render(paste0(fun,'_',typ),paste(switch(fun,plot_abundance='丰度',plot_biomass='生物量',plot_SSB='产卵生物量',plot_catch='捕获尾数'),typ),paste(switch(fun,plot_abundance='Abundance',plot_biomass='Biomass',plot_SSB='Spawning biomass',plot_catch='Catch numbers'),typ),
  do.call(fun,list(b,type=typ,se=TRUE,facet_ncol=3,line_size=.6)),height=if(nchar(typ)==1 || typ %in% c('CN','SSB'))4 else 8.5,
  note_zh=paste('ALSCL 示例；A 表示年龄组，L 表示长度组。有可用标准误时显示近似 95% 区间。',if(fun=='plot_catch')'捕获尾数是模型推算的渔业捕获，不是调查指数。'else '总量图显示各组之和；分面图的纵轴按组调整。'),note_en=paste('ALSCL example: A denotes age, L length. Approximate 95% intervals appear where standard errors are available.',if(fun=='plot_catch')'These are model-derived fishery catches, not survey indices.'else 'Totals sum across groups; faceted y axes adjust to each group.'))
}
render('recruitment','补充量','Recruitment',plot_recruitment(b,se=TRUE),note_zh='补充量对应最小年龄组；受设定的补充年龄和时间步长影响。',note_en='Recruitment enters the youngest age class; its interpretation depends on recruitment age and time step.')
render('SSB_Rec','亲体与补充关系','Spawning biomass and recruitment',plot_SSB_Rec(b,age_at_recruitment=1),note_zh='横轴与补充量错开 1 个观测步长。这里只画散点，不拟合资源补充函数。',note_en='Recruitment is shifted by one observation step. This is a scatter plot, not a fitted stock recruitment curve.')
render('VB','生长曲线与区间','Growth curve and interval',plot_VB(b,age_range=c(1,15),se=TRUE),note_zh='生长参数在本例中固定，区间退化为曲线；释放参数后区间才反映估计不确定性，不是个体长度分布。',note_en='Growth is fixed in this example, so its interval collapses. With estimated growth, intervals reflect parameter uncertainty, not individual length variation.')
render('pla','年龄长度转换','Age length conversion',plot_pla(b)+scale_x_discrete(labels=as.character(1:15)),note_zh='色值是长度组给定年龄的概率；不是增长转移矩阵 G。年龄标签为组序号。',note_en='Cells are length-bin probabilities conditional on age, not the growth transition matrix G. Age labels are class indices.')
for(model in c('ACL','ALSCL'))for(typ in c('year','age','length')) {
 fit<-if(model=='ACL')a else b
 render(paste('F',model,typ,sep='_'),paste(model,'捕捞死亡率',typ),paste(model,'fishing mortality',typ),plot_fishing_mortality(fit,type=typ,se=TRUE,facet_ncol=4,line_size=.6),height=8.5,
 note_zh=if(model=='ACL' && typ=='length')'当前 ACL 的 length 分支与 age 分支相同，仍是年龄曲线，不能作为长度别 F。' else '年度示例的 F 为每年瞬时率；季度数据得到每季度率。不同模型的原生 F 维度不同。',
 note_en=if(model=='ACL' && typ=='length')'In this version, ACL length is an alias of age; the output is still age-specific F.' else 'F is an instantaneous rate per model step. Native F dimensions differ between the models.')
}
for(typ in c('length','year')) render(paste0('residuals_',typ),paste('对数残差',typ),paste('Log residuals',typ),plot_residuals(b,type=typ,facet_ncol=4,x_breaks=if(typ=='length')c(2000,2010,2015)else c(10,30,50),line_size=.6),height=8.5,note_zh='观测对数减预测对数；未经标准差标准化。平滑线用于发现结构，不是显著性检验。',note_en='Observed minus fitted log index, without standardization. Smooth curves reveal structure, not significance.')
for(typ in c('R','F'))for(lg in c(TRUE,FALSE))render(paste0('deviation_',typ,'_',lg),paste('过程偏差',typ,lg),paste('Process deviations',typ,lg),plot_deviance(b,type=typ,log=lg,se=TRUE,facet_ncol=3,point_size=1.2,line_size=.5)+scale_x_continuous(breaks=c(2000,2010,2019),expand=expansion(mult=.04)),height=if(typ=='F')8.5 else 4,
note_zh=if(lg)'对数过程偏差以 0 为参照。这不是似然偏差统计量。' else '指数转换后以 1 为参照，区间不对称。这不是原尺度残差。',note_en=if(lg)'Log process deviations use zero as reference; these are not likelihood deviance statistics.' else 'Exponentiated deviations use one as reference and asymmetric intervals; these are not original-scale residuals.')
render('ridges','长度组成山脊图','Length composition ridges',plot_ridges(b),height=8.5,note_zh='每期分别归一化；左为观测，右为拟合。不能由此读出总丰度趋势。需当前会话中的 obj。',note_en='Each period is normalized separately: observations left, fit right. This does not show total abundance trends. A live obj is required.')
render('compare_ts','种群时间序列比较','Population time series comparison',plot_compare_ts(a,b,se=TRUE,ncol=2),height=5.7,note_zh='只比较共同输出与相同时间窗；自由纵轴不能直接比较不同量的振幅。',note_en='Compare common quantities over the same time window. Free y axes preclude amplitude comparisons across quantities.')
render('compare_F','死亡率矩阵比较','Fishing mortality matrices',plot_compare_F(a,b)+scale_x_continuous(breaks=c(2000,2010,2019)),height=5,note_zh='ACL 纵轴是年龄，ALSCL 纵轴是长度；颜色可辅助观察，但两行不逐格对应。',note_en='ACL uses age and ALSCL length on the vertical axis; cells do not correspond one to one.')
render('compare_residuals','残差联合诊断','Combined residual diagnostics',plot_compare_residuals(a,b,dat),height=5.5,note_zh='包括直方图、QQ 图、年度残差和长度组箱线图。年度误差棒是标准差，不是置信区间。',note_en='Includes histograms, QQ plots, annual residuals and length-bin boxplots. Annual error bars represent SD, not confidence intervals.')
# 调整组合图排版，保留全部参数 / Wrap subtitles without changing estimates.
growth_display <- function(p) {
 p[[1]] <- p[[1]] + labs(subtitle=paste(strsplit(p[[1]]$labels$subtitle,' | ',fixed=TRUE)[[1]],collapse='\n')) + theme(legend.text=element_text(size=7))
 p[[2]] <- p[[2]] + scale_x_continuous(breaks=c(1,5,10,15))
 p
}
# 参考生长点也是模拟数据 / Reference growth points are synthetic too.
set.seed(42); ref<-data.frame(Age=rep(1:15,each=5)); ref$Length<-VB_func(60,.2,1/60,ref$Age)+rnorm(nrow(ref),0,2)
write.csv(ref,file.path(out,'data','synthetic_growth_reference.csv'),row.names=FALSE)
render('compare_growth','生长与外部参照接口','Growth and reference-data interface',growth_display(plot_compare_growth(a,b,age_range=c(1,15),ref_data=ref,ref_name='Synthetic reference',nls_start=list(Linf=60,k=.2,t0=1/60))),height=8.5,note_zh='灰色参照曲线来自模拟 Age 与 Length；不代表论文实测生长。',note_en='The reference curve uses synthetic Age and Length data, not empirical measurements from the paper.')
render('compare_CatL','两模型调查拟合','Survey fits from both models',plot_compare_CatL(a,b,dat,years=c(2000,2005,2010,2015),ncol=2),height=5,note_zh='显示 4 个指定年份；省略 years 将显示全部年份。点为调查数据，线为拟合中位数。',note_en='Four selected years are shown. Omit years to show all periods. Points are observations and curves fitted medians.')
render('compare_metrics','拟合指标比较','Fit metric comparison',plot_compare_metrics(a,b,dat,ncol=2)&scale_y_continuous(expand=expansion(mult=c(0,.3))),height=5.5,note_zh='当前数据由年龄模型生成。图示差异不能证明任一模型普遍较优；IC 还需要相同数据与似然口径。',note_en='The generating model is age based. This comparison cannot establish universal superiority; information criteria also require comparable data and likelihoods.')
render('compare_selectivity','固定调查可捕性','Fixed survey catchability',plot_compare_selectivity(a,b),note_zh='当前函数展示输入 q 的连接线，不是估计的渔业选择性。两模型使用同一输入，曲线重合。',note_en='This function connects the input survey q; it does not estimate fishery selectivity. Identical inputs produce overlapping curves.')
for(method in c('apical','mean'))render(paste0('compare_F_',method),paste('总体 F 比较',method),paste('Summary F comparison',method),plot_compare_annual_F(a,b,method=method)+labs(title=if(method=='apical')'Apical fishing mortality'else 'Mean fishing mortality'),note_zh='apical 取最大值，mean 取组间算术平均，均非丰度加权。函数名称含 annual，但不会自动把季度率年化。',note_en='Apical is the maximum; mean is an unweighted group average. Despite its name, the function does not annualize quarterly rates.')
# 回溯逐次删除末端时间步 / Peel successive terminal time steps.
retro_a <- do.call(retro_acl, c(inputs, x$fit_args, x$fit_config$acl,
  list(nyear=npeels, train_times=2, ncores=1, silent=TRUE)))
retro_b <- do.call(retro_alscl, c(inputs, x$fit_args, x$fit_config$alscl,
  list(nyear=npeels, train_times=2, ncores=1, silent=TRUE)))
render('retro_ACL','ACL 回溯分析','ACL retrospective analysis',
  plot_retro(retro_a,facet_col=2,rho_digits=3,rho_size=3,point_size=1,line_size=.6),height=6)
render('retro_ALSCL','ALSCL 回溯分析','ALSCL retrospective analysis',
  plot_retro(retro_b,facet_col=2,rho_digits=3,rho_size=3,point_size=1,line_size=.6),height=6)
write.csv(retro_a$rho_text,file.path(out,'results','rho_ACL.csv'),row.names=FALSE)
write.csv(retro_b$rho_text,file.path(out,'results','rho_ALSCL.csv'),row.names=FALSE)

  # 以模拟真值检查估计的实际含义 / Compare with the generating truth.
  years <- a$year
  truth <- data.frame(Year=rep(years,3),Quantity=rep(c("B","SSB","Rec"),each=length(years)),
                      Truth=c(x$truth$TB,x$truth$SSB,x$truth$Rec))
  estimates <- rbind(data.frame(Year=rep(years,3),Quantity=truth$Quantity,
    Value=c(a$report$B,a$report$SSB,a$report$Rec),Model="ACL"),
    data.frame(Year=rep(years,3),Quantity=truth$Quantity,
    Value=c(b$report$B,b$report$SSB,b$report$Rec),Model="ALSCL"))
  render("truth","模拟真值与条件估计","Simulation truth and conditional estimates",
    ggplot(estimates,aes(Year,Value,color=Model))+geom_line()+
      scale_color_manual(values=setNames(acl_theme("compare_colors"),c("ACL","ALSCL")))+
      geom_line(data=truth,aes(Year,Truth),inherit.aes=FALSE,linetype=2)+
      facet_wrap(~Quantity,scales="free_y",ncol=1)+theme_bw()+labs(y=NULL),height=6,
    note_zh="黑虚线为生成真值；固定生物参数，且模拟与拟合的补充过程不同。",
    note_en="Black dashed curves are generating truth. Biology is fixed; recruitment processes differ.")
  write.csv(truth,file.path(out,"results","truth.csv"),row.names=FALSE)
  write.csv(estimates,file.path(out,"results","estimates.csv"),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(out,"results","sessionInfo.txt"))
  message("COMPLETE: ",length(registry)," YTF figures")
  invisible(list(acl=a,alscl=b,retro_acl=retro_a,retro_alscl=retro_b,inputs=inputs))
}
if (sys.nframe() == 0L) run_ytf_guide()
