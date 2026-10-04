# =============================================================================
# Repository path helper — source at the top of every analysis script.
#
#   source("code/_paths.R")        # run from the repository root
#
# No absolute paths are used anywhere. The repository root is resolved as:
#   1. $TBCIRCRNA_ROOT, if set and valid (optional override);
#   2. the working directory, when it looks like the repository root;
#   3. otherwise the directory containing this script, walking upwards.
#
# Layout
#   config/     analysis parameters, locked panel, locked thresholds
#   data/       INPUT placeholders — download the deposited matrices here
#   metadata/   de-identified sample dictionary and cohort summary
#   results/    lightweight derived results (fixed, match the manuscript)
#   env/        R and package versions captured from the analysis environment
#   docs/       workflow, data availability, reproducibility
#   output/     default write directory for re-runs
# =============================================================================

.paths_find_root <- function() {
  env <- Sys.getenv("TBCIRCRNA_ROOT", unset = "")
  if (nzchar(env) && dir.exists(env)) return(normalizePath(env, winslash = "/"))

  looks_like_root <- function(p) {
    dir.exists(file.path(p, "code")) && dir.exists(file.path(p, "config"))
  }

  wd <- normalizePath(".", winslash = "/", mustWork = TRUE)
  if (looks_like_root(wd)) return(wd)

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  start <- if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
  } else {
    wd
  }

  cur <- start
  for (i in seq_len(8)) {
    if (looks_like_root(cur)) return(normalizePath(cur, winslash = "/"))
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }

  stop("Repository root not found. Run from the repository root ",
       "or set TBCIRCRNA_ROOT to the cloned repository path.", call. = FALSE)
}

ROOT <- .paths_find_root()

CODE   <- file.path(ROOT, "code")
CONFIG <- file.path(ROOT, "config")
DATA   <- file.path(ROOT, "data")

DATA_META    <- file.path(DATA, "metadata")
DATA_QN      <- file.path(DATA, "qn")
DATA_QN_RAW  <- file.path(DATA_QN, "raw")
DATA_BSJ     <- file.path(DATA, "bsj")
DATA_BENCH   <- file.path(DATA, "benchmark")

RESULTS <- file.path(ROOT, "results")
DOCS    <- file.path(ROOT, "docs")
ENV     <- file.path(ROOT, "env")
OUT     <- file.path(ROOT, "output")

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Locate a required analysis input and fail with an actionable message.
# `what` is a short description; `where` is the GEO / public source of the file.
paths_require <- function(path, what, where) {
  if (file.exists(path)) return(invisible(normalizePath(path, winslash = "/")))
  stop(sprintf(
    "Required input file not found.\n  File : %s\n  What : %s\n  Where: %s\n  Hint : place the downloaded file at the path above (relative to %s).",
    sub(paste0("^", ROOT, "/?"), "", path), what, where, ROOT), call. = FALSE)
}

# Read a TSV that may start with "#" comment lines (data.table::fread has no
# comment.char argument).
paths_read_tsv <- function(path, ...) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]
  data.table::fread(text = paste(lines, collapse = "\n"), ...)
}

message("[paths] ROOT = ", ROOT)
