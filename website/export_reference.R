# Extract signatures from source, without installing ALSCL or running any models.
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else "website/.build/api"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
exports <- sub("^export\\(([^)]+)\\)$", "\\1",
  grep("^export\\(", readLines("NAMESPACE"), value = TRUE))
env <- new.env(parent = baseenv())
sources <- list()
for (path in list.files("R", "\\.R$", full.names = TRUE)) {
  for (expr in parse(path, keep.source = FALSE)) {
    if (is.call(expr) && is.symbol(expr[[1]]) && as.character(expr[[1]]) %in% c("<-", "=") &&
        is.symbol(expr[[2]]) && is.call(expr[[3]]) &&
        identical(expr[[3]][[1]], as.name("function"))) {
      name <- as.character(expr[[2]])
      eval(expr, env)
      sources[[name]] <- path
    }
  }
}
records <- lapply(exports, function(name) {
  stopifnot(exists(name, env, inherits = FALSE))
  f <- get(name, env)
  defaults <- lapply(formals(f), function(x) paste(deparse(x, width.cutoff = 500), collapse = " "))
  # An empty deparse is the required-argument marker; NULL is a real default.
  signature <- paste(deparse(args(f), width.cutoff = 80), collapse = "\n")
  signature <- sub("^function", name, sub("[[:space:]]*NULL$", "", signature))
  list(name = name, defaults = defaults, signature = signature, source = sources[[name]])
})
jsonlite::write_json(records, file.path(out, "functions.json"), auto_unbox = TRUE, pretty = TRUE)
for (path in list.files("man", "\\.Rd$", full.names = TRUE)) {
  rd <- tools::parse_Rd(path)
  aliases <- unlist(lapply(rd, function(x) if (identical(attr(x, "Rd_tag"), "\\alias")) paste(x, collapse = "")))
  if (any(aliases %in% c(exports, "YTF", "YTF_example", "example_data"))) {
    tools::Rd2HTML(rd, out = file.path(out, paste0(basename(path), ".html")),
      no_links = TRUE)
  }
}
cat("Extracted", length(records), "public function signatures and R help pages.\n")
