# scripts/07_cross_validate/compare_deposit_to_internal.R
# -----------------------------------------------------------------------------
# Tier 1 of the deposit-repointing verification. Loads each dataset from BOTH
# sources through easy_scripts/_data_access.R and proves they are the same data.
#
# Run before touching any mini-script: everything downstream assumes the mapping
# in _data_access.R is correct, and this is what establishes that.
#
# On tolerances. The residuals below are the 15-significant-digit round trip of
# readr::write_csv, not error -- e.g. alr_flow is 1.39041969184572 in the
# deposit and 1.3904196918457186 internally. They are pinned at the measured
# magnitude rather than set to a comfortable round number, so that a genuine
# regression cannot hide underneath a generous threshold.
#
# switches_per_session is deliberately NOT given a tolerance. It genuinely
# differs by up to 1.4e-2 because build_datasets_workbook.R restores it to the
# integer count it always was. A tolerance loose enough to swallow that would
# also swallow a real regression in every other column, so it gets its own
# status instead. It backs no reported outcome.
#
# Usage:
#   source("paths.R")
#   source(file.path(PROJECT_ROOT,
#          "choice R pipeline/scripts/07_cross_validate/compare_deposit_to_internal.R"))
# -----------------------------------------------------------------------------

if (!exists("PROJECT_ROOT")) stop("Run source(\"paths.R\") first.", call. = FALSE)
source(file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts/_data_access.R"))

.TOL <- c(
  logit_flow            = 1e-12,
  lr_high               = 1e-12,
  lr_medium             = 1e-12,
  lr_low                = 1e-12,
  mean_nnd_cm           = 1e-12,
  mean_iid_cm           = 1e-12,
  mean_centroid_spd_cm  = 1e-12,
  mean_hull_area_cm2    = 1e-10,
  zone_flux_per_session = 1e-11,
  prop_flow             = 1e-14,
  prop_calm             = 1e-14,
  obs_seconds           = 0,
  n_frames              = 0,
  bouts_per_min         = 1e-11,
  longest_flow_bout_s   = 1e-11,
  longest_calm_bout_s   = 1e-11,
  mean_flow_bout_s      = 1e-11,
  mean_calm_bout_s      = 1e-11
)
.EXPECTED_DIFF <- c(switches_per_session = 0.05)

.rows <- list()
.add <- function(check, status, detail = "") {
  .rows[[length(.rows) + 1L]] <<- data.frame(check = check, status = status,
                                             detail = detail,
                                             stringsAsFactors = FALSE)
  cat(sprintf("  %-42s %-20s %s\n", check, status, detail))
}

cat("\n=== Tier 1: deposit vs pipeline CSV ===\n\n")

# ---- behaviour --------------------------------------------------------------
cat("BEHAVIOUR\n")
dep <- load_behaviour_dataset(source = "deposit",  quiet = TRUE)
int <- load_behaviour_dataset(source = "internal", quiet = TRUE)

.add("behaviour: deposit rows", if (nrow(dep) == 48L) "PASS" else "FAIL", nrow(dep))
.add("behaviour: internal rows", if (nrow(int) == 48L) "PASS" else "FAIL", nrow(int))

key <- function(d) paste(d$phys_trial_id, d$timepoint, sep = "|")
.add("behaviour: key sets identical",
     if (setequal(key(dep), key(int))) "PASS" else "FAIL")

dep <- dep[order(key(dep)), ]; int <- int[order(key(int)), ]

for (col in names(.TOL)) {
  in_dep <- col %in% names(dep); in_int <- col %in% names(int)
  if (in_dep && !in_int) {
    # The five BOUT-engine outcomes live in the deposit but not in
    # easy_scripts_dataset.csv -- the workbook joins them from a separate
    # engine. Not a defect: the deposit is strictly richer here.
    .add(paste0("behaviour: ", col), "DEPOSIT-ONLY",
         "absent from easy_scripts_dataset.csv (BOUT engine)")
    next
  }
  if (!in_dep || !in_int) {
    .add(paste0("behaviour: ", col), "FAIL",
         paste0("deposit=", in_dep, " internal=", in_int))
    next
  }
  a <- as.numeric(dep[[col]]); b <- as.numeric(int[[col]])
  if (!identical(is.na(a), is.na(b))) {
    .add(paste0("behaviour: ", col), "FAIL", "NA pattern differs"); next
  }
  d <- max(abs(a - b), na.rm = TRUE)
  if (!is.finite(d)) d <- 0
  .add(paste0("behaviour: ", col),
       if (d <= .TOL[[col]]) "PASS" else "FAIL",
       sprintf("max|delta| = %.3e  (tol %.0e)", d, .TOL[[col]]))
}

for (col in names(.EXPECTED_DIFF)) {
  a <- as.numeric(dep[[col]]); b <- as.numeric(int[[col]])
  d <- max(abs(a - b), na.rm = TRUE)
  .add(paste0("behaviour: ", col),
       if (d <= .EXPECTED_DIFF[[col]]) "EXPECTED-DIFFERENCE" else "FAIL",
       sprintf("max|delta| = %.3e  (deposit rounds to integer count; backs no reported outcome)", d))
}

for (col in c("treatment", "tank", "motor_side", "trial_date")) {
  .add(paste0("behaviour: ", col, " (character)"),
       if (identical(as.character(dep[[col]]), as.character(int[[col]])))
         "PASS" else "FAIL")
}
.add("behaviour: fish_density_f numeric",
     if (is.numeric(dep$fish_density_f)) "PASS" else "FAIL",
     paste(class(dep$fish_density_f), collapse = "/"))

# ---- endocrine --------------------------------------------------------------
cat("\nENDOCRINE\n")
edep <- load_endocrine_dataset(source = "deposit",  quiet = TRUE)
eint <- load_endocrine_dataset(source = "internal", quiet = TRUE)

.add("endocrine: deposit rows",  "INFO", nrow(edep))
.add("endocrine: internal rows", "INFO", nrow(eint))
.add("endocrine: row counts match",
     if (nrow(edep) == nrow(eint)) "PASS" else "FAIL",
     paste0(nrow(edep), " vs ", nrow(eint)))

ekey <- function(d) paste(d$sample_id, ifelse(is.na(d$area), "", d$area),
                          d$analyte, sep = "|")
.add("endocrine: key sets identical",
     if (setequal(ekey(edep), ekey(eint))) "PASS" else "FAIL")

# Join on the FULL key. Monoamines carry exactly four rows per
# (sample_id, analyte) -- one per brain region -- so joining on
# (sample_id, analyte) alone produces a 4x cartesian product and compares
# unrelated regions against each other.
.kcol <- function(d) paste(d$sample_id, ifelse(is.na(d$area), "", d$area),
                           d$analyte, sep = "|")
edep$.k <- .kcol(edep); eint$.k <- .kcol(eint)
.cols <- c(".k","area","value","tank","trial","treatment","sex")
m <- merge(edep[, .cols], eint[, .cols], by = ".k",
           suffixes = c(".dep", ".int"), all = FALSE)
m$analyte <- sub("^[^|]*\\|[^|]*\\|", "", m$.k)
.add("endocrine: join is 1:1",
     if (nrow(m) == nrow(edep)) "PASS" else "FAIL",
     paste0(nrow(m), " joined rows vs ", nrow(edep), " deposit rows"))
mono <- m[m$analyte != "cort", ]
cort <- m[m$analyte == "cort", ]

dc <- max(abs(cort$value.dep - cort$value.int), na.rm = TRUE)
.add("endocrine: cortisol value", if (dc == 0) "PASS" else "FAIL",
     sprintf("max|delta| = %.3e (exact expected)", dc))
dm <- max(abs(mono$value.dep - mono$value.int), na.rm = TRUE)
.add("endocrine: monoamine value", if (dm <= 1e-10) "PASS" else "FAIL",
     sprintf("max|delta| = %.3e  (tol 1e-10)", dm))

for (col in c("tank", "trial", "treatment", "sex")) {
  ok <- identical(as.character(m[[paste0(col, ".dep")]]),
                  as.character(m[[paste0(col, ".int")]]))
  .add(paste0("endocrine: ", col), if (ok) "PASS" else "FAIL")
}
.add("endocrine: treatment levels",
     if (setequal(unique(edep$treatment), c("control", "treat"))) "PASS" else "FAIL",
     paste(sort(unique(edep$treatment)), collapse = ", "))

# ---- verdict ----------------------------------------------------------------
res <- do.call(rbind, .rows)
n_fail <- sum(res$status == "FAIL")
cat("\n---\n")
cat(sprintf("checks: %d   PASS: %d   FAIL: %d   EXPECTED-DIFFERENCE: %d\n",
            nrow(res), sum(res$status == "PASS"), n_fail,
            sum(res$status == "EXPECTED-DIFFERENCE")))
if (n_fail) {
  cat("\nFAILURES:\n"); print(res[res$status == "FAIL", ], row.names = FALSE)
  stop("Tier 1 failed. Do not proceed to the mini-scripts.", call. = FALSE)
}
cat("Tier 1 GREEN - the mapping in _data_access.R is sound.\n\n")
invisible(res)
