# config.R -- machine-specific locations that are NOT part of this repository.
#
# These are deliberately separate from paths.R. PROJECT_ROOT describes where the
# code is; these describe what else the machine happens to have. Each falls back
# to the value used on the original analysis machine, so nothing changes for the
# author, and each can be overridden by an environment variable.

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
