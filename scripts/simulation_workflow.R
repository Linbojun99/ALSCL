# 两条模拟路径与批量拟合 / Two simulation pathways and batch fitting.
# source("scripts/simulation_workflow.R", encoding="UTF-8")
library(ALSCL)

simulation_to_tables <- function(sm, pa, first_year=2000) {
  years <- first_year + (seq_len(sm$nyear)-1)*pa$growth_step
  frame <- function(m) {
    d <- data.frame(LengthBin=as.character(pa$len_mid),m,check.names=FALSE)
    names(d) <- c("LengthBin",as.character(years)); d
  }
  list(data.CatL=frame(t(sm$SN_at_len)),data.wgt=frame(sm$weight),data.mat=frame(sm$mat))
}

run_simulation_examples <- function(out="simulation_examples",fit_batch=FALSE) {
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  # 简单年度模型：直接返回三张表 / Simple annual model returns fit-ready tables.
  quick <- simulate_example_data(years=2000:2019, seed=42,
    Linf=60, vbk=.2, t0=1/60, M=.2, nage=15,
    bin_breaks=seq(5,51,2), L50_sel=15, L95_sel=20,
    L50_mat=35, L95_mat=40, cv_catch=.2,
    save_csv=TRUE,output_dir=file.path(out,"quick"))
  # 完整物种模拟：保留真值、预热期与重复文件 / Full simulation retains truth.
  cases <- list()
  for (sp in c("flatfish","tuna")) {
    folder <- file.path(out,sp)
    pa <- initialize_params(species=sp,observation_error="independent")
    bio <- sim_cal(pa)
    sm <- sim_data(bio,pa,iter_range=4:5,return_iter=4,output_dir=folder)
    inputs <- simulation_to_tables(sm,pa)
    for (n in names(inputs)) utils::write.csv(inputs[[n]],
      file.path(folder,paste0(n,".csv")),row.names=FALSE)
    saveRDS(list(parameters=pa,truth=sm),file.path(folder,"truth.rds"))
    cases[[sp]] <- list(parameters=pa,truth=sm,inputs=inputs)
  }
  # 文件由 save 写入，使用 load / sim_rep files use save/load, not RDS.
  e <- new.env(parent=emptyenv())
  load(file.path(out,"flatfish","sim_rep4"),envir=e)
  stopifnot(is.list(e$sim.data))
  if (fit_batch) {
    data(YTF_example,package="ALSCL")
    batch <- do.call(sim_acl,c(list(iter_range=4:5,
      sim_data_path=file.path(out,"flatfish"),output_dir=file.path(out,"batch_acl"),
      M=cases$flatfish$parameters$M,ncores=1,train_times=2),YTF_example$fit_config$acl))
    # 保存失败原因；不能只计算成功重复的均值 / Retain failures in summaries.
    checks <- lapply(seq_along(batch),function(i) {
      f<-batch[[i]]
      if (!is.null(f$error)) return(data.frame(Seed=i+3,Code=NA,Gradient=NA,
        pdHess=FALSE,Boundary=NA,Error=f$error))
      data.frame(Seed=i+3,Code=f$convergence_code,Gradient=f$max_abs_gradient,
        pdHess=f$pdHess,Boundary=f$bound_hit,Error="")
    })
    utils::write.csv(do.call(rbind,checks),file.path(out,"batch_checks.csv"),row.names=FALSE)
  }
  invisible(list(quick=quick,cases=cases))
}
if (sys.nframe()==0L) run_simulation_examples()
