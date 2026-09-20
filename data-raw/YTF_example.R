# 从仓库根目录运行 / Run from the repository root after installing ALSCL.
# YTF 是参数列表；此脚本创建可直接拟合的教学数据 / YTF is a parameter list.
library(ALSCL)
data("YTF", package = "ALSCL")

build_ytf_example <- function() {
  # 固定算法与种子，退出时恢复随机状态 / Pin and restore the RNG state.
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  folder <- tempfile("ytf-simulation-")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  # 全部共有生物参数取自 YTF，特别是 std_SN=0.1 / Reuse YTF, including std_SN=0.1.
  shared <- intersect(names(YTF), names(formals(initialize_params)))
  overrides <- YTF[setdiff(shared, "nyear")]
  pa <- do.call(initialize_params, c(list(
    species = "flatfish", nyear = 80 + YTF$nyear, burn_in = 80,
    first.year = 2000, growth_step = 1, observation_error = "independent"
  ), overrides))
  bio <- sim_cal(pa)
  truth <- sim_data(bio, pa, iter_range = 4L, return_iter = 4L,
                    output_dir = folder)
  years <- 2000 + seq_len(truth$nyear) - 1
  frame <- function(x) {
    ans <- data.frame(LengthBin = as.character(pa$len_mid), x,
                      check.names = FALSE)
    names(ans) <- c("LengthBin", as.character(years))
    ans
  }
  # 条件拟合固定已知生物参数；不用于证明可识别性 / Conditional teaching fits.
  p_acl <- list(log_Linf = log(pa$Linf), log_vbk = log(pa$vbk),
    t0 = pa$t0, log_cv_len = log(pa$cv_L), log_std_log_F = log(pa$F_sd),
    log_std_log_N0 = log(pa$std_logN0), logit_log_R = qlogis(pa$R_ar))
  p_alscl <- list(log_Linf = log(pa$Linf), log_vbk = log(pa$vbk),
    log_t0 = log(pa$t0), log_cv_len = log(pa$cv_L),
    log_cv_grow = log(pa$cv_inc), log_sigma_log_F = log(pa$F_sd),
    log_sigma_log_N0 = log(pa$std_logN0), logit_log_R = qlogis(pa$R_ar))
  fixed <- function(x) setNames(rep(list(factor(NA)), length(x)), names(x))
  m_acl <- fixed(p_acl)
  m_alscl <- fixed(p_alscl)
  m_alscl$logit_log_F_l <- factor(NA)
  m_alscl$logit_log_F_y <- factor(NA)

  list(
    data.CatL = frame(t(truth$SN_at_len)),
    data.wgt = frame(truth$weight), data.mat = frame(truth$mat),
    parameters = pa, truth = truth,
    fit_args = list(rec.age = pa$rec.age, nage = pa$nage, M = pa$M,
      sel_L50 = pa$q_surv_L50, sel_L95 = pa$q_surv_L95,
      len_mid = pa$len_mid,
      len_border = pa$len_border[-c(1L, length(pa$len_border))],
      growth_step = 1),
    fit_config = list(acl = list(parameters = p_acl, map = m_acl),
                     alscl = list(parameters = p_alscl, map = m_alscl)),
    provenance = list(
      kind = "Synthetic teaching data; not observations from the paper",
      species = "Yellowtail flounder (Limanda ferruginea)",
      parameter_source = "ALSCL::YTF with flatfish operating-model defaults",
      code_version = "2.0.0", seed = 4L,
      rng = c("Mersenne-Twister", "Inversion", "Rejection"),
      generating_model = "age_based", burn_in_years = 80,
      observation_years = years, observation_log_sd = pa$std_SN,
      observation_error = "independent",
      time_labels = "2000:2019 are illustrative labels, not sampling dates",
      doi = "10.1111/faf.12673"))
}

YTF_example <- build_ytf_example()
dir.create("data", showWarnings = FALSE)
save(YTF_example, file = "data/YTF_example.rda", compress = "xz", version = 2)
message("Saved data/YTF_example.rda: 23 length bins x 20 annual observations")
