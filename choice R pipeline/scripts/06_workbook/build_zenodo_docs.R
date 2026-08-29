# =============================================================================
# build_zenodo_docs.R
# -----------------------------------------------------------------------------
# Writes data_dictionary.csv into D:/CHOICE R SCRIPTS/zenodo_dataset/, reading
# the CSVs that build_zenodo_deposit.R has already written plus the workbook
# legend sheets. Run it AFTER build_zenodo_deposit.R.
#
# It also prints the assay-QC figures that README.md quotes, so those numbers
# can be re-checked against their source rather than trusted. README.md and
# CHECKLIST.md are prose about a frozen deposit and are maintained by hand.
#
# Everything the repository actually knows is filled in. Everything it does not
# is written as a literal [TO FILL] marker, never guessed -- a deposit that
# states an ethics licence number or an LOQ nobody checked is worse than one
# that admits the gap.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(readxl); library(tibble)
})

.ROOT <- PROJECT_ROOT
.PIPE <- file.path(.ROOT, "choice R pipeline")
.OUT  <- file.path(.ROOT, "zenodo_dataset")
.WB   <- file.path(.OUT, "source_workbook.xlsx")
.CV   <- file.path(.ROOT, "choice exp for claude mono and cortisol/claude output REVISED/assay_cv_tables.xlsx")
.log  <- function(...) cat("[docs] ", ..., "\n", sep = "")
TOFILL <- "[TO FILL]"

trials <- readr::read_csv(file.path(.OUT, "pref_trials.csv"),          show_col_types = FALSE)
fish   <- readr::read_csv(file.path(.OUT, "pref_fish.csv"),            show_col_types = FALSE)
cort   <- readr::read_csv(file.path(.OUT, "pref_cortisol.csv"),        show_col_types = FALSE)
mono   <- readr::read_csv(file.path(.OUT, "pref_monoamines_long.csv"), show_col_types = FALSE)

# -----------------------------------------------------------------------------
# 1. data_dictionary.csv
# -----------------------------------------------------------------------------
beh_leg <- as.data.frame(readxl::read_excel(.WB, sheet = "Preference_exp_behavior_legend"))

# role is stricter than the workbook's `type`: every transformed index is
# "derived", which is what tells a reanalyst which columns are raw measurement.
.DERIVED <- c("alr_flow", "alr_high", "alr_medium", "alr_low",
              "clr_high", "clr_medium", "clr_low", "clr_calm",
              "jacobs_d_flow", "jacobs_d_high", "jacobs_d_medium", "jacobs_d_low")
.QC      <- c("obs_seconds", "n_frames", "cv_intra_pct", "cv_source")
.IDENT   <- c("fish_id", "cortisol_sample_id", "monoamine_sample_id", "phys_trial")

.role_of <- function(col, fallback = "metadata") {
  if (col %in% .DERIVED) "derived"
  else if (col %in% .QC) "qc"
  else if (col %in% .IDENT) "identifier"
  else fallback
}

.LEVELS <- c(
  treatment    = "control; exercise choice",
  tank         = "27; 28; 29; 31",
  fish_density = "4; 8; 12; 16",
  motor_side   = "FT; FD",
  interval     = "1; 2; 3",
  phys_trial   = "1-16",
  trial_seq    = "1-16",
  sex          = "F; M",
  region       = "DM; POA; VV; VD",
  analyte      = "5-HT; 5-HIAA; 5-HIAA/5-HT; DA; DOPAC; DOPAC/DA",
  plate        = "5; 6; 7; 8; Evg",
  units        = "ng_per_ml",
  cv_source    = "myassays; recomputed",
  batch        = "3",
  missing_reason = "assay_failure_amplifier; no_protein_quantification; derived_from_missing_component; unspecified; no_value_returned"
)

# --- pref_trials.csv: reuse the workbook legend, which is generated from the
#     same indicator map the analysis uses, so the two cannot drift.
dd_trials <- tibble(file = "pref_trials.csv", column = names(trials)) %>%
  left_join(beh_leg %>% transmute(column, wb_type = type, indicator_name, unit, description),
            by = "column") %>%
  mutate(
    role = vapply(seq_along(column), function(i)
      .role_of(column[i], if (is.na(wb_type[i]) || wb_type[i] == "unknown")
                            "metadata" else wb_type[i]), ""),
    levels = unname(.LEVELS[column]),
    description = ifelse(is.na(description) | description == "Needs review",
                         TOFILL, description),
    unit = ifelse(is.na(unit), "", unit),
    indicator_name = ifelse(is.na(indicator_name), "", indicator_name)) %>%
  select(file, column, role, indicator_name, unit, levels, description)

# --- the three deposit-only tables --------------------------------------------
.row <- function(file, column, role, unit, description, indicator_name = "")
  tibble(file = file, column = column, role = role,
         indicator_name = indicator_name, unit = unit,
         levels = unname(.LEVELS[column]), description = description)

dd_fish <- bind_rows(
  .row("pref_fish.csv", "fish_id", "identifier", "text",
       "Unique fish identifier, PREF_001 to PREF_080. Primary key for the preference experiment. The numeric part matches the cortisol sample number."),
  .row("pref_fish.csv", "cortisol_sample_id", "identifier", "text",
       "Cortisol sample code, B3_01 to B3_80. Joins to pref_cortisol.csv."),
  .row("pref_fish.csv", "monoamine_sample_id", "identifier", "integer",
       "Monoamine sample number as recorded by the chromatography run, 1 to 74. Empty for the 43 fish assayed for cortisol only. It equals the numeric part of cortisol_sample_id: monoamine sample N is cortisol B3_N. Joins to pref_monoamines_long.csv."),
  .row("pref_fish.csv", "phys_trial", "identifier", "integer",
       "Behavioural trial this fish came from, 1 to 16. Joins to pref_trials.csv. This is the manuscript's trial numbering."),
  .row("pref_fish.csv", "trial_seq", "metadata", "integer",
       "The same trial in chronological running order. Differs from phys_trial: trials 1 and 2, 13 and 14, and 15 and 16 are swapped between the two numberings."),
  .row("pref_fish.csv", "trial_key", "metadata", "text",
       "Experimental unit coded '<tank>_<number of fish>', e.g. '27_16'. The one key shared by every table in the study; maps one-to-one onto phys_trial."),
  .row("pref_fish.csv", "tank", "metadata", "factor", "Holding tank the school came from."),
  .row("pref_fish.csv", "fish_density", "metadata", "count (fish)",
       "Number of fish in the school for that trial."),
  .row("pref_fish.csv", "treatment", "metadata", "factor", "Experimental condition."),
  .row("pref_fish.csv", "motor_side", "metadata", "factor",
       "Side of the arena on which the flow motor was mounted. Held constant for the four trials of a holding tank and reversed between tanks, so it is perfectly confounded with tank."),
  .row("pref_fish.csv", "sex", "metadata", "factor",
       "Fish sex determined at dissection. Empty for one fish (PREF_080) where it was not recorded."),
  .row("pref_fish.csv", "weight_g", "metadata", "g", "Body mass at sampling."),
  .row("pref_fish.csv", "length_cm", "metadata", "cm", "Fork length at sampling."),
  .row("pref_fish.csv", "video_id", "metadata", "text",
       "Sampling video the fish was taken from, B3_N1 or B3_N2. Two fish were sampled per video."),
  .row("pref_fish.csv", "batch", "metadata", "integer",
       "Experimental batch. All fish in this deposit are batch 3."),
  .row("pref_fish.csv", "sampling_date", "metadata", "date (ISO 8601)",
       "Date the fish was sampled. Equal to the trial date for every fish in this deposit.")
)

dd_cort <- bind_rows(
  .row("pref_cortisol.csv", "fish_id", "identifier", "text", "Joins to pref_fish.csv."),
  .row("pref_cortisol.csv", "cortisol_sample_id", "identifier", "text", "Cortisol sample code."),
  .row("pref_cortisol.csv", "phys_trial", "identifier", "integer",
       "Behavioural trial, 1 to 16. Joins to pref_trials.csv at trial level."),
  .row("pref_cortisol.csv", "treatment", "metadata", "factor", "Experimental condition."),
  .row("pref_cortisol.csv", "plate", "metadata", "text",
       "Assay plate. TEXT, not numeric: it holds 'Evg' alongside the numeric plate IDs 5 to 8. 'Evg' is a fifth plate run separately; see the README."),
  .row("pref_cortisol.csv", "plate_date", "metadata", "date (ISO 8601)", "Date the plate was run."),
  .row("pref_cortisol.csv", "cortisol_plasma", "indicator", "ng/mL",
       "Plasma cortisol concentration. Empty for the two fish that returned no value.",
       indicator_name = "Plasma cortisol"),
  .row("pref_cortisol.csv", "units", "metadata", "text",
       "Unit of cortisol_plasma, constant across the file."),
  .row("pref_cortisol.csv", "cv_intra_pct", "qc", "percent",
       "Intra-assay coefficient of variation across the sample's replicate wells. Empty where fewer than two wells were in range."),
  .row("pref_cortisol.csv", "cv_source", "qc", "factor",
       "Whether cv_intra_pct came from the assay software ('myassays') or was recomputed from the raw well readings."),
  .row("pref_cortisol.csv", "missing_reason", "metadata", "factor",
       "Reason cortisol_plasma is empty. Empty when a value is present.")
)

dd_mono <- bind_rows(
  .row("pref_monoamines_long.csv", "fish_id", "identifier", "text", "Joins to pref_fish.csv."),
  .row("pref_monoamines_long.csv", "cortisol_sample_id", "identifier", "text",
       "Carried so a reader can pair a monoamine row with its cortisol row without going through pref_fish.csv."),
  .row("pref_monoamines_long.csv", "monoamine_sample_id", "identifier", "integer",
       "Monoamine sample number from the chromatography run."),
  .row("pref_monoamines_long.csv", "phys_trial", "identifier", "integer",
       "Behavioural trial, 1 to 16. Joins to pref_trials.csv at trial level."),
  .row("pref_monoamines_long.csv", "treatment", "metadata", "factor", "Experimental condition."),
  .row("pref_monoamines_long.csv", "region", "metadata", "factor",
       "Brain region. DM dorsomedial pallium; POA preoptic area; VV ventral part of the ventral telencephalon; VD dorsal part of the ventral telencephalon."),
  .row("pref_monoamines_long.csv", "analyte", "metadata", "factor",
       "Measured compound or turnover ratio. 5-HT serotonin; 5-HIAA its catabolite; DA dopamine; DOPAC its catabolite. The two ratios are derived, not measured."),
  .row("pref_monoamines_long.csv", "concentration", "indicator", "see units column",
       "Tissue concentration for the four measured analytes; a unitless quotient for the two ratios. Empty where the value is missing.",
       indicator_name = "Tissue monoamine"),
  .row("pref_monoamines_long.csv", "units", "metadata", "text",
       "ng_per_mg_tissue for the four measured analytes, 'ratio' for the two turnover ratios."),
  .row("pref_monoamines_long.csv", "missing_reason", "metadata", "factor",
       "Reason concentration is empty. Empty when a value is present.")
) %>% mutate(levels = ifelse(column == "units", "ng_per_mg_tissue; ratio", levels))

data_dictionary <- bind_rows(dd_trials, dd_fish, dd_cort, dd_mono) %>%
  mutate(across(everything(), ~ ifelse(is.na(.x), "", .x)))

# every column of every deposited CSV must appear exactly once
.expected <- c(paste("pref_trials.csv", names(trials)),
               paste("pref_fish.csv", names(fish)),
               paste("pref_cortisol.csv", names(cort)),
               paste("pref_monoamines_long.csv", names(mono)))
.got <- paste(data_dictionary$file, data_dictionary$column)
if (!setequal(.expected, .got))
  stop("data_dictionary does not cover every column exactly once.\n  missing: ",
       paste(setdiff(.expected, .got), collapse = ", "),
       "\n  extra:   ", paste(setdiff(.got, .expected), collapse = ", "), call. = FALSE)
if (anyDuplicated(.got)) stop("data_dictionary has duplicate rows", call. = FALSE)

readr::write_csv(data_dictionary, file.path(.OUT, "data_dictionary.csv"), na = "")
.log("data_dictionary.csv  ", nrow(data_dictionary), " rows covering 4 files")
.n_tofill <- sum(data_dictionary$description == TOFILL)
if (.n_tofill > 0) .log("  ", .n_tofill, " description(s) still marked ", TOFILL)

# -----------------------------------------------------------------------------
# 2. Assay QC figures, read from the source rather than typed in
# -----------------------------------------------------------------------------
.cvs <- as.data.frame(readxl::read_excel(.CV, sheet = "intra_assay_summary"))
.cv_line <- paste(sprintf("plate %s: %.3g%%", .cvs$plate, as.numeric(.cvs$median_cv)),
                  collapse = "; ")
.std <- as.data.frame(readxl::read_excel(.CV, sheet = "inter_assay_cv_standards"))
.std <- .std[!is.na(suppressWarnings(as.numeric(.std$cv_inter_pct))), ]
.std_line <- sprintf("%.3g%% to %.3g%% across the %d back-fitted standards",
                     min(as.numeric(.std$cv_inter_pct)), max(as.numeric(.std$cv_inter_pct)),
                     nrow(.std))

.miss_tab <- mono %>% filter(is.na(concentration)) %>% count(missing_reason, name = "n")
.miss_line <- paste(sprintf("`%s` (%d cells)", .miss_tab$missing_reason, .miss_tab$n),
                    collapse = ", ")
.miss_reg <- mono %>% filter(is.na(concentration)) %>% count(region, name = "n")
.miss_reg_line <- paste(sprintf("%s %d", .miss_reg$region, .miss_reg$n), collapse = ", ")

# Printed, not written: these numbers belong in README.md prose, and a stray
# dotfile in the deposit folder would be uploaded to Zenodo along with the data.
.log("QC figures for README.md:")
for (l in c(.cv_line, .std_line, .miss_line, .miss_reg_line)) .log("  ", l)
