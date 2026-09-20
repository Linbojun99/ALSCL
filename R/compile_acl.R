#' Compile and load the ACL model in a writable session cache
#' @param openmp Logical. Request OpenMP, with a serial fallback.
#' @return Paths and DLL name for the compiled model.
#' @keywords internal
compile_and_load_acl <- function(openmp = FALSE) .acl_compile("ACL", openmp)
.acl_compile <- function(model, openmp = FALSE) {
  # Binary TMB installations do not necessarily install their LinkingTo headers.
  if (!requireNamespace("RcppEigen", quietly = TRUE)) stop("Install RcppEigen to compile model templates.")
  src <- system.file("extdata", paste0(model, ".cpp"), package = "ALSCL")
  header <- system.file("extdata", "model_math.hpp", package = "ALSCL")
  if (!nzchar(src) || !nzchar(header)) stop("Model sources are missing; reinstall the ALSCL package.")
  key_file <- tempfile(); on.exit(unlink(key_file), add = TRUE)
  writeLines(c(unname(tools::md5sum(c(src, header))), R.version.string,
               as.character(utils::packageVersion("TMB")),
               as.character(utils::packageVersion("RcppEigen")), as.character(openmp)), key_file)
  key <- substr(unname(tools::md5sum(key_file)), 1L, 12L); dll <- paste0(model, "_", key)
  cache <- file.path(getOption("ALSCL.tmb.cache", file.path(tempdir(), "ALSCL-tmb")), dll); dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  cache <- normalizePath(cache, winslash = "/", mustWork = TRUE)
  cpp <- file.path(cache, paste0(dll, ".cpp")); lib <- TMB::dynlib(file.path(cache, dll))
  if (!file.exists(lib)) {
    file.copy(src, cpp, overwrite = TRUE); file.copy(header, file.path(cache, "model_math.hpp"), overwrite = TRUE)
    # Let make see a plain filename, avoiding Windows backslashes and spaces in paths.
    previous_dir <- setwd(cache); on.exit(setwd(previous_dir), add = TRUE)
    if (openmp) {
      ok <- tryCatch({TMB::compile(basename(cpp), flags = "-O2", openmp = TRUE); TRUE}, error = function(e) FALSE)
      if (!ok) {unlink(c(lib, sub("\\.cpp$", ".o", cpp))); TMB::compile(basename(cpp), flags = "-O2", openmp = FALSE)}
    } else TMB::compile(basename(cpp), flags = "-O2", openmp = FALSE)
  }
  loaded <- vapply(getLoadedDLLs(), function(x) x[["path"]], character(1))
  loaded <- loaded[file.exists(loaded)] # R's built-in "base" entry is not a file.
  if (!normalizePath(lib, winslash = "/") %in% normalizePath(loaded, winslash = "/")) dyn.load(lib)
  invisible(list(cpp_path = cpp, dll_path = lib, dll_name = dll))
}
#' Unload a model dynamic library
#' @param dll_path Path to the library. Existing TMB objects must no longer be used.
#' @keywords internal
unload_acl <- function(dll_path) {
  paths <- vapply(getLoadedDLLs(), function(x) x[["path"]], character(1))
  path <- normalizePath(dll_path, mustWork = FALSE)
  if (path %in% paths) dyn.unload(path)
  invisible(NULL)
}
