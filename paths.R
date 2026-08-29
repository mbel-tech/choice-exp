# paths.R -- project root discovery for the choice-exp repository.
#
# Every script here resolves file locations against PROJECT_ROOT rather than a
# hardcoded drive letter, so the repository runs from wherever it is cloned.
#
# Usage, from anywhere inside the repository:
#
#     source("paths.R")                       # if the working directory is the root
#     source(file.path(PROJECT_ROOT, "..."))  # thereafter
#
# Resolution order:
#   1. the CHOICE_EXP_ROOT environment variable, if set
#   2. walking up from the working directory for the .choice-exp-root sentinel
#
# Failure is loud and immediate. A silently wrong root would write output to the
# wrong place, which is worse than not running at all.

.choice_exp_find_root <- function(start = getwd()) {
  env <- Sys.getenv("CHOICE_EXP_ROOT", unset = "")
  if (nzchar(env)) {
    if (!dir.exists(env)) {
      stop("CHOICE_EXP_ROOT is set to '", env, "' but no such directory exists.")
    }
    return(normalizePath(env, winslash = "/", mustWork = TRUE))
  }

  d <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(d, ".choice-exp-root"))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) {
      stop(
        "Could not locate the choice-exp repository root.\n",
        "  Searched upward from : ", normalizePath(start, winslash = "/"), "\n",
        "  Looking for sentinel : .choice-exp-root\n",
        "  Fix: setwd() to a directory inside the repository, or set the\n",
        "       CHOICE_EXP_ROOT environment variable to the repository root."
      )
    }
    d <- parent
  }
}

PROJECT_ROOT <- .choice_exp_find_root()
