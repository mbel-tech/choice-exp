# config.R -- machine-specific locations that are NOT part of this repository.
#
# These are deliberately separate from paths.R. PROJECT_ROOT describes where the
# code is; these describe what else the machine happens to have. Each falls back
# to the value used on the original analysis machine, so nothing changes for the
# author, and each can be overridden by an environment variable.

if (!exists("PROJECT_ROOT")) {
  stop("config.R needs PROJECT_ROOT.\n",
       "  Run source(\"paths.R\") from the repository root first.",
       call. = FALSE)
}

# Pandoc, used to render .docx tables and reports.
# Original machine: D:/tools/pandoc-3.10.1/pandoc.exe
PANDOC <- Sys.getenv("CHOICE_EXP_PANDOC",
                     unset = "D:/tools/pandoc-3.10.1/pandoc.exe")

# Raw idtracker.ai session data (~14 GB). Not in this repository and not in the
# Zenodo deposit; it is primary tracking output held in backup. Only the STEP1/
# STEP2 stages need it -- the statistical pipeline reads derived indicators.
PIPELINE_DATA_DIR <- Sys.getenv("CHOICE_EXP_DATA_DIR", unset = "")

if (!nzchar(PIPELINE_DATA_DIR)) {
  message(
    "PIPELINE_DATA_DIR is unset. Stages that read raw idtracker.ai sessions ",
    "will not run.\n  Set CHOICE_EXP_DATA_DIR to the checked_sessions folder ",
    "if you need them."
  )
}

# The Zenodo deposit: four CSVs, a data dictionary and the source workbook.
# Not in this repository -- download the DOI named in README.md and unpack it
# here. The easy_scripts mini-scripts read their data through it.
#
# NOT the same thing as PIPELINE_DATA_DIR above, despite the similar names.
# That one is ~14 GB of raw idtracker.ai sessions and only STEP1/STEP2 read it.
# This one is ~200 kB of derived indicators, and is what every *reported*
# analysis actually needs. If you are trying to reproduce a number from the
# paper, this is the one you want.
DATA_ROOT <- Sys.getenv("CHOICE_EXP_DATA_ROOT",
                        unset = file.path(PROJECT_ROOT, "zenodo_dataset"))

# Which source the easy_scripts loader prefers when both are present:
#   "auto"     -- deposit first, then the pipeline-side CSVs (default)
#   "deposit"  -- the Zenodo CSVs only; fail if absent
#   "internal" -- the pipeline-side easy_scripts_*.csv only; fail if absent
DATA_SOURCE <- Sys.getenv("CHOICE_EXP_DATA_SOURCE", unset = "auto")
