# =============================================================================
# build_zenodo_deposit.R
# -----------------------------------------------------------------------------
# Builds: D:/CHOICE R SCRIPTS/zenodo_dataset/
#
#   README.md
#   CHECKLIST.md
#   data_dictionary.csv
#   pref_fish.csv                 80 rows   subject table + assay crosswalk
#   pref_trials.csv               48 rows   trial x interval behaviour
#   pref_cortisol.csv             80 rows   plasma cortisol
#   pref_monoamines_long.csv     888 rows   37 fish x 4 regions x 6 analytes
#   source_workbook.xlsx                    provenance copy
#
# SCOPE. The preference experiment only (control vs exercise choice). The
# observational study is not part of this deposit.
#
# The behaviour and endocrine tables are read back out of the workbook that
# build_datasets_workbook.R just wrote, so the deposit and the workbook cannot
# disagree. Fish-level metadata and the missing-value reasons come from the two
# upstream files the workbook does not carry in full.
#
# CONVENTIONS (Zenodo deposit handout, section 4)
#   plain CSV, UTF-8, comma separated, "." decimal, one header row
#   ISO 8601 dates throughout
#   lowercase_underscore column names; no spaces, parentheses or slashes
#   missing values are empty cells, with the reason in a missing_reason column
#   units live in a units column or in the column name, never in the cell
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(readxl); library(tibble)
})

.ROOT <- PROJECT_ROOT
.PIPE <- file.path(.ROOT, "choice R pipeline")
.OUT  <- file.path(.ROOT, "zenodo_dataset")
.WB   <- file.path(.PIPE, "Behavioural and neuroendocrine correlates datasets.xlsx")
.ENDO <- file.path(.ROOT, "choice exp for claude mono and cortisol/claude output REVISED")
.CORT_CLEAN <- file.path(.ENDO, "b3_cortisol_clean.csv")
.MONO_LONG  <- file.path(.ENDO, "Data/monoamine_long.csv")
.CV_TABLES  <- file.path(.ENDO, "assay_cv_tables.xlsx")

.log <- function(...) cat("[deposit] ", ..., "\n", sep = "")
dir.create(.OUT, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(.WB), file.exists(.CORT_CLEAN), file.exists(.MONO_LONG))

# CSVs are written with write_csv: it emits "" for NA, which is what the handout
# asks for (empty cell + a reason in the paired missing_reason column).
.write <- function(d, name) {
  p <- file.path(.OUT, name)
  readr::write_csv(d, p, na = "")
  .log(sprintf("%-28s %4d rows x %2d cols", name, nrow(d), ncol(d)))
  invisible(p)
}

# -----------------------------------------------------------------------------
# 1. pref_trials.csv -- trial x interval behaviour (48 rows)
# -----------------------------------------------------------------------------
beh <- as.data.frame(readxl::read_excel(.WB, sheet = "Preference_exp_behavior"))
stopifnot(nrow(beh) == 48L)

pref_trials <- beh %>%
  mutate(trial_date = format(as.Date(trial_date, format = "%d.%m.%Y"), "%Y-%m-%d"),
         interval   = as.integer(interval),
         phys_trial = as.integer(phys_trial)) %>%
  arrange(phys_trial, interval)
stopifnot(!anyNA(pref_trials$trial_date),
          length(unique(pref_trials$trial_date)) == 8L)
.write(pref_trials, "pref_trials.csv")

# The design, one row per trial. Used to give every fish its behavioural trial.
design <- pref_trials %>%
  distinct(phys_trial, trial_seq, tank, fish_density, treatment, motor_side, trial_date) %>%
  mutate(trial_key = paste0(tank, "_", fish_density))
stopifnot(nrow(design) == 16L, !anyDuplicated(design$trial_key))

# -----------------------------------------------------------------------------
# 2. pref_fish.csv -- the subject table and the assay crosswalk (80 rows)
# -----------------------------------------------------------------------------
# THE CROSSWALK. The two assays were long thought to be unlinkable, because the
# distributed workbook identified cortisol samples as "B3_01" and monoamine
# samples as the integer 1. They are the same fish: monoamine sample_id N is
# cortisol B3_<N padded to two digits>. Verified against (tank, trial,
# treatment, sex) -- 36 of the 37 monoamine fish match a cortisol row exactly,
# with zero mismatches, and the 37th (B3_19) is a fish that was assayed for
# monoamines but returned no cortisol value. Both identifiers are carried here
# so no reader has to rediscover this.
cort_raw <- readr::read_csv(.CORT_CLEAN, show_col_types = FALSE)
stopifnot(nrow(cort_raw) == 80L)

mono_raw <- readr::read_csv(.MONO_LONG, show_col_types = FALSE)
mono_ids <- mono_raw %>%
  distinct(sample_id, monoamine_sample_id = sample_n) %>%
  mutate(sample_id = as.character(sample_id),
         monoamine_sample_id = as.integer(monoamine_sample_id))
stopifnot(nrow(mono_ids) == 37L)
# the crosswalk claim, asserted rather than assumed
stopifnot(all(mono_ids$sample_id == sprintf("B3_%02d", mono_ids$monoamine_sample_id)))

fish <- cort_raw %>%
  transmute(
    cortisol_sample_id = as.character(sample_id),
    fish_id   = sprintf("PREF_%03d", as.integer(sub("^B3_", "", sample_id))),
    trial_key = as.character(trial),
    sex       = ifelse(is.na(sex) | sex == "NA", NA_character_, as.character(sex)),
    weight_g  = suppressWarnings(as.numeric(wt)),
    length_cm = suppressWarnings(as.numeric(lenght)),
    video_id  = as.character(video_id),
    batch     = suppressWarnings(as.integer(batch)),
    sampling_date = format(as.Date(date), "%Y-%m-%d")
  ) %>%
  left_join(mono_ids, by = c("cortisol_sample_id" = "sample_id")) %>%
  left_join(design,   by = "trial_key") %>%
  transmute(fish_id, cortisol_sample_id, monoamine_sample_id,
            phys_trial, trial_seq, trial_key, tank, fish_density,
            treatment, motor_side, sex, weight_g, length_cm,
            video_id, batch, sampling_date) %>%
  arrange(fish_id)

stopifnot(nrow(fish) == 80L, !anyNA(fish$phys_trial),
          sum(!is.na(fish$monoamine_sample_id)) == 37L,
          all(table(fish$phys_trial) == 5L))
.write(fish, "pref_fish.csv")

.key <- fish %>% select(fish_id, cortisol_sample_id, monoamine_sample_id,
                        phys_trial, treatment)

# -----------------------------------------------------------------------------
# 3. pref_cortisol.csv -- one row per sampled fish (80 rows, 78 assayed)
# -----------------------------------------------------------------------------
# plate is TEXT: it holds "Evg" alongside the numeric plate IDs 5-8.
pref_cortisol <- cort_raw %>%
  transmute(
    cortisol_sample_id = as.character(sample_id),
    plate       = as.character(plate),
    plate_date  = format(as.Date(plate_date), "%Y-%m-%d"),
    cortisol_plasma = suppressWarnings(as.numeric(plasma_cortisol)),
    units       = "ng_per_ml",
    cv_intra_pct = suppressWarnings(as.numeric(cv_intra_pct)),
    cv_source   = ifelse(is.na(cv_source) | cv_source == "NA", NA_character_,
                         as.character(cv_source))
  ) %>%
  left_join(.key, by = "cortisol_sample_id") %>%
  mutate(missing_reason = ifelse(is.na(cortisol_plasma), "no_value_returned",
                                 NA_character_)) %>%
  transmute(fish_id, cortisol_sample_id, phys_trial, treatment,
            plate, plate_date, cortisol_plasma, units,
            cv_intra_pct, cv_source, missing_reason) %>%
  arrange(fish_id)

stopifnot(nrow(pref_cortisol) == 80L,
          sum(!is.na(pref_cortisol$cortisol_plasma)) == 78L,
          !any(is.na(pref_cortisol$cortisol_plasma) & is.na(pref_cortisol$missing_reason)))
.write(pref_cortisol, "pref_cortisol.csv")

# -----------------------------------------------------------------------------
# 4. pref_monoamines_long.csv -- fish x region x analyte
# -----------------------------------------------------------------------------
# MISSING VALUES. The source carries a free-text note per affected sample; those
# notes are recoded here into a controlled vocabulary. Nothing is guessed: where
# the source records no reason the code is "unspecified", not an invented one.
#
# DEVIATION FROM HANDOUT SECTION 3.3, deliberate. The handout says the turnover
# ratios must not be rows because they are derived. They are derived -- verified
# here to 3e-14 -- but they are also ANALYSED cells: the 24 analyte x region
# cells the paper reports are {5-HT, 5-HIAA, 5-HIAA/5-HT, DA, DOPAC, DOPAC/DA}
# x {DM, POA, VV, VD}. Dropping them would leave 8 of the 24 reported cells
# absent from the deposit. They are kept as rows and marked role = "derived" in
# the data dictionary, with the formula stated in the README.
.reason <- function(note) {
  n <- tolower(ifelse(is.na(note), "", as.character(note)))
  dplyr::case_when(
    grepl("protein concentration", n)          ~ "no_protein_quantification",
    grepl("amplifier", n) & grepl("peak", n)   ~ "assay_failure_amplifier",
    TRUE                                       ~ "unspecified")
}

.MONO_UNITS <- c("5-HT" = "ng_per_mg_tissue", "5-HIAA" = "ng_per_mg_tissue",
                 "DA" = "ng_per_mg_tissue",  "DOPAC"  = "ng_per_mg_tissue",
                 "5-HIAA/5-HT" = "ratio", "DOPAC/DA" = "ratio")
.IS_RATIO   <- c("5-HIAA/5-HT", "DOPAC/DA")
.RATIO_PARTS <- list("5-HIAA/5-HT" = c("5-HIAA", "5-HT"),
                     "DOPAC/DA"    = c("DOPAC", "DA"))

mono <- mono_raw %>%
  transmute(cortisol_sample_id = as.character(sample_id),
            monoamine_sample_id = as.integer(sample_n),
            region  = as.character(area),
            analyte = as.character(analyte),
            concentration = suppressWarnings(as.numeric(value)),
            note = as.character(notes))
stopifnot(setequal(unique(mono$region), c("DM", "POA", "VV", "VD")),
          setequal(unique(mono$analyte), names(.MONO_UNITS)),
          !"NE" %in% unique(mono$analyte))

# ratio == quotient, asserted before we call the ratios derived
.chk <- mono %>%
  select(cortisol_sample_id, region, analyte, concentration) %>%
  tidyr::pivot_wider(names_from = analyte, values_from = concentration)
.dev <- c(abs(.chk[["5-HIAA/5-HT"]] - .chk[["5-HIAA"]] / .chk[["5-HT"]]),
          abs(.chk[["DOPAC/DA"]]    - .chk[["DOPAC"]]  / .chk[["DA"]]))
stopifnot(max(.dev, na.rm = TRUE) < 1e-9)
.log("turnover ratios verified as exact quotients (max deviation ",
     format(max(.dev, na.rm = TRUE), digits = 3), ")")

# a ratio is empty whenever either of its parts is empty; say so rather than
# repeating the part's reason
.present <- mono %>% filter(!is.na(concentration)) %>%
  distinct(cortisol_sample_id, region, analyte) %>% mutate(have = TRUE)
.part_missing <- function(sid, reg, an) {
  pr <- .RATIO_PARTS[[an]]
  !all(paste(sid, reg, pr) %in% paste(.present$cortisol_sample_id,
                                      .present$region, .present$analyte))
}

pref_mono <- mono %>%
  rowwise() %>%
  mutate(missing_reason = if (!is.na(concentration)) NA_character_
                          else if (analyte %in% .IS_RATIO &&
                                   .part_missing(cortisol_sample_id, region, analyte))
                                 "derived_from_missing_component"
                          else .reason(note)) %>%
  ungroup() %>%
  mutate(units = unname(.MONO_UNITS[analyte])) %>%
  left_join(.key %>% select(fish_id, cortisol_sample_id, phys_trial, treatment),
            by = "cortisol_sample_id") %>%
  transmute(fish_id, cortisol_sample_id, monoamine_sample_id,
            phys_trial, treatment, region, analyte,
            concentration, units, missing_reason) %>%
  arrange(monoamine_sample_id,
          match(region, c("DM", "POA", "VV", "VD")),
          match(analyte, names(.MONO_UNITS)))

stopifnot(nrow(pref_mono) == 37L * 4L * 6L,
          !anyNA(pref_mono$fish_id),
          !any(is.na(pref_mono$concentration) & is.na(pref_mono$missing_reason)),
          !any(!is.na(pref_mono$concentration) & !is.na(pref_mono$missing_reason)))
.write(pref_mono, "pref_monoamines_long.csv")

.log("missing measurement cells by reason: ",
     paste(sprintf("%s=%d", names(table(pref_mono$missing_reason)),
                   as.integer(table(pref_mono$missing_reason))), collapse = ", "))

# -----------------------------------------------------------------------------
# 5. source_workbook.xlsx
# -----------------------------------------------------------------------------
file.copy(.WB, file.path(.OUT, "source_workbook.xlsx"), overwrite = TRUE)
.log("source_workbook.xlsx copied")

# -----------------------------------------------------------------------------
# 6. Cross-file integrity gate
# -----------------------------------------------------------------------------
errs <- character(0)
E <- function(...) errs <<- c(errs, sprintf(...))

if (!all(pref_cortisol$fish_id %in% fish$fish_id)) E("pref_cortisol has a fish_id not in pref_fish")
if (!all(pref_mono$fish_id     %in% fish$fish_id)) E("pref_monoamines has a fish_id not in pref_fish")
if (!all(fish$phys_trial       %in% pref_trials$phys_trial)) E("pref_fish has a phys_trial not in pref_trials")
if (length(unique(pref_trials$phys_trial)) != 16L) E("pref_trials does not hold 16 trials")
if (!all(table(pref_trials$phys_trial) == 3L))     E("pref_trials: not every trial has 3 intervals")

for (f in list.files(.OUT, pattern = "[.]csv$", full.names = TRUE)) {
  nm <- names(readr::read_csv(f, n_max = 0, show_col_types = FALSE))
  bad <- grep("[ ()/]|[A-Z]", nm, value = TRUE)
  if (length(bad)) E("%s: column name(s) with a space, bracket, slash or capital: %s",
                     basename(f), paste(bad, collapse = ", "))
  if (length(grep("_agg$", nm))) E("%s: trial-aggregate column present", basename(f))
  if (length(grep("^ne__|^(BS|HY|TEL|OT)_|FISH_ID|trial_2nd_exp", nm)))
    E("%s: NE or observational-study column present", basename(f))
}
if (length(errs)) stop("Deposit integrity FAILED:\n  - ",
                       paste(errs, collapse = "\n  - "), call. = FALSE)
.log("Integrity: OK")

# objects the documentation builder needs
DEPOSIT <- list(trials = pref_trials, fish = fish,
                cortisol = pref_cortisol, monoamines = pref_mono,
                out = .OUT, wb = .WB)
