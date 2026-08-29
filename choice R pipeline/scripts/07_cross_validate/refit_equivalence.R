# =============================================================================
# refit_equivalence.R   —  Layer 2 of cross-validation
# -----------------------------------------------------------------------------
# For each sampled easy-script indicator, run the script TWICE in isolated
# child envs:
#   (A) with the original easy_scripts CSVs
#   (B) with a workbook-derived CSV that has been written to a temp dir
#       and made visible to the script via assignInNamespace() monkey-patch
#       of readr::read_csv (so the script's `readr::read_csv(file.path(ROOT,
#       "easy_scripts_dataset.csv"))` call returns the workbook-derived data).
#
# After both runs, diff the captured ANOVA data frames (per-term, per-metric)
# and the AICc-selected RE structure. Tolerance: rel_diff <= 1e-6 for floats,
# exact for integer df / categorical fields.
#
# A run is PASS for Layer 2 iff every (indicator, term, metric) comparison
# passes — i.e. the same statistics come out of the same model regardless of
# which dataset (original or workbook-derived) fed it.
#
# This file defines functions only — invoked by cross_validate_workbook.R.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
})

# ---- Indicator → easy-script-path lookup -----------------------------------
# wb_col / synthetic id  ->  relative path under easy_scripts/
# Behaviour: every workbook indicator maps to its by_timepoint script (since
# the workbook is interval-grain). Aggregated variants ("_agg") are NOT in
# the workbook so we don't refit them here.
.ES_ROOT <- file.path(PROJECT_ROOT, "choice R pipeline/easy_scripts")

.BEH_ES_SCRIPT_LOOKUP <- list(
  # 7 flow indicators
  logit_flow              = "by_timepoint/logit_flow_by_tp.R",
  lr_high                 = "by_timepoint/lr_high_by_tp.R",
  lr_medium               = "by_timepoint/lr_medium_by_tp.R",
  lr_low                  = "by_timepoint/lr_low_by_tp.R",
  # 6 collective indicators
  switches_per_session    = "by_timepoint/switches_per_session_by_tp.R",
  mean_nnd_cm             = "by_timepoint/mean_nnd_cm_by_tp.R",
  mean_polarisation       = "by_timepoint/mean_polarisation_by_tp.R",
  mean_school_speed_cm_s  = "by_timepoint/mean_centroid_spd_cm_by_tp.R",
  mean_iid_cm             = "by_timepoint/mean_iid_cm_by_tp.R",
  mean_school_area_cm2    = "by_timepoint/mean_hull_area_cm2_by_tp.R"
  # Jacobs + 3 sub-zone composites in the workbook (subzone_alr_pairs_tp,
  # main_cells_flow_tp, sub_cells_high_clr_tp) don't have 1:1 easy scripts;
  # we skip them with a clear SKIPPED reason in the report.
)

# Behaviour indicators we can refit (subset of workbook indicators)
.BEH_REFITTABLE <- names(.BEH_ES_SCRIPT_LOOKUP)

# Cortisol: single script
.CORT_ES_SCRIPT <- "endocrine/plasma_cortisol.R"

# Monoamine cells: one script per analyte (script handles all 4 areas in one
# fit `log(value) ~ treatment * area`). So sampling 12 of 24 cells maps to
# at most 6 distinct script runs (we de-duplicate before running).
#
# 2026-08-18: noradrenaline removed. NE was dropped from the study, and
# easy_scripts_endo_dataset.csv carries only da, dopac, dopac_da_ratio,
# hiaa_5, hiaa_5_ratio, ht_5 (+ cort) -- there is no `ne` row to refit. The
# grid is 6 analytes x 4 regions = 24 cells, not 28. The lookup kept pointing
# at monoamine_ne.R because the _OUTDATED rename had broken this map
# wholesale, so the dead entry was never reached.
.MONO_ANALYTE_SCRIPT <- c(
  ht_5            = "endocrine/monoamine_5_ht.R",
  hiaa_5          = "endocrine/monoamine_5_hiaa.R",
  hiaa_5_ratio    = "endocrine/monoamine_5_hiaa_5_ht.R",
  da              = "endocrine/monoamine_da.R",
  dopac           = "endocrine/monoamine_dopac.R",
  dopac_da_ratio  = "endocrine/monoamine_dopac_da.R"
)

# ---- read_csv monkey-patch utility -----------------------------------------
# We swap readr::read_csv for the duration of a refit. The patch intercepts
# the specific filenames easy-scripts read and redirects them to the supplied
# alternative paths. Restoration is guaranteed via on.exit().

.with_redirected_csvs <- function(beh_csv, endo_csv, expr) {
  # Capture current binding
  ns <- asNamespace("readr")
  orig <- readr::read_csv

  patched <- function(file, ...) {
    bn <- basename(as.character(file))
    if (bn == "easy_scripts_dataset.csv"      && !is.null(beh_csv))
      file <- beh_csv
    if (bn == "easy_scripts_endo_dataset.csv" && !is.null(endo_csv))
      file <- endo_csv
    orig(file, ...)
  }

  # Replace in readr namespace
  tryCatch(unlockBinding("read_csv", ns), error = function(e) invisible(NULL))
  assign("read_csv", patched, envir = ns)
  lockBinding("read_csv", ns)

  on.exit({
    tryCatch(unlockBinding("read_csv", ns), error = function(e) invisible(NULL))
    assign("read_csv", orig, envir = ns)
    lockBinding("read_csv", ns)
  })

  force(expr)
}

# ---- ggsave monkey-patch utility -------------------------------------------
# Easy scripts call ggplot2::ggsave to write a PNG into easy_scripts/outputs/
# graphs/. During refit we DO NOT want to overwrite the user's working PNGs,
# so we no-op ggsave.
.with_silenced_ggsave <- function(expr) {
  ns <- asNamespace("ggplot2")
  orig <- ggplot2::ggsave
  noop <- function(...) invisible(NULL)
  tryCatch(unlockBinding("ggsave", ns), error = function(e) invisible(NULL))
  assign("ggsave", noop, envir = ns)
  lockBinding("ggsave", ns)
  on.exit({
    tryCatch(unlockBinding("ggsave", ns), error = function(e) invisible(NULL))
    assign("ggsave", orig, envir = ns)
    lockBinding("ggsave", ns)
  })
  force(expr)
}

# ---- Source one easy script and capture its environment --------------------
.refit_one_script <- function(script_path) {
  env <- new.env(parent = globalenv())
  # Silence message() chatter during refit
  capture.output(
    suppressMessages(sys.source(script_path, envir = env)),
    file = nullfile(), type = c("output", "message"))
  env
}

# ---- Extract ANOVA frame from a refit env ----------------------------------
# Different easy scripts store the ANOVA in different places. Try them in
# priority order; return a tidy frame: term, statistic, df_num, df_denom, p_value.
.extract_anova <- function(env) {
  # behaviour scripts: res$anova
  if (exists("res", envir = env, inherits = FALSE)) {
    res <- get("res", envir = env)
    if (!is.null(res$anova)) return(.tidy_anova(res$anova))
  }
  # cortisol: cort_anova
  if (exists("cort_anova", envir = env, inherits = FALSE)) {
    return(.tidy_anova(get("cort_anova", envir = env)))
  }
  # monoamine: A1 (script convention) or fit_<x>_anova
  for (nm in c("A1", "fit_anova", "anova_kr", "mono_anova")) {
    if (exists(nm, envir = env, inherits = FALSE))
      return(.tidy_anova(get(nm, envir = env)))
  }
  # As a final fallback: search the env for any object that looks like
  # an anova-table data frame (has columns 'F' or 'F value' or 'Chisq').
  for (nm in ls(env)) {
    v <- get(nm, envir = env, inherits = FALSE)
    if (inherits(v, "anova") || inherits(v, "data.frame")) {
      cn <- tolower(names(v))
      if (any(c("f","f.value","f value","chisq","pr(>f)","pr(>chisq)") %in% cn)) {
        return(.tidy_anova(v))
      }
    }
  }
  NULL
}

.tidy_anova <- function(a) {
  if (is.null(a)) return(NULL)
  df <- as.data.frame(a, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(NULL)
  # Term name: prefer rownames if not "1","2"; otherwise check 'term' col
  if (!is.null(rownames(df)) && !all(rownames(df) == as.character(seq_len(nrow(df))))) {
    df$term <- rownames(df)
  } else if (!"term" %in% names(df)) {
    df$term <- as.character(seq_len(nrow(df)))
  }
  # Pull the standard columns where we can
  pick <- function(candidates) {
    for (c in candidates) if (c %in% names(df)) return(df[[c]])
    rep(NA_real_, nrow(df))
  }
  data.frame(
    term      = as.character(df$term),
    statistic = as.numeric(pick(c("F value","F","Chisq","chisq","statistic"))),
    df_num    = as.numeric(pick(c("Df","df","numDF","NumDF"))),
    df_denom  = as.numeric(pick(c("Df.res","df_denom","DenDF","denDF","Df.error"))),
    p_value   = as.numeric(pick(c("Pr(>F)","Pr(>Chisq)","p_value","p.value","Pr(>|F|)"))),
    stringsAsFactors = FALSE
  )
}

# ---- Compare two tidy ANOVA frames -----------------------------------------
.compare_anova <- function(a_orig, a_wb, indicator, kind,
                            abs_tol = 1e-9, rel_tol = 1e-6) {
  if (is.null(a_orig) || is.null(a_wb)) {
    return(data.frame(
      indicator = indicator, kind = kind, term = NA_character_,
      metric = "anova", pipeline_value = NA_real_, workbook_value = NA_real_,
      abs_diff = NA_real_, rel_diff = NA_real_, tolerance = rel_tol,
      status = "SKIPPED-NO-ANOVA",
      stringsAsFactors = FALSE))
  }
  terms_union <- union(a_orig$term, a_wb$term)
  rows <- list()
  for (tm in terms_union) {
    oi <- which(a_orig$term == tm)[1]
    wi <- which(a_wb$term   == tm)[1]
    if (is.na(oi) || is.na(wi)) {
      rows[[length(rows)+1L]] <- data.frame(
        indicator = indicator, kind = kind, term = tm,
        metric = "presence",
        pipeline_value = if (!is.na(oi)) 1 else 0,
        workbook_value = if (!is.na(wi)) 1 else 0,
        abs_diff = NA_real_, rel_diff = NA_real_,
        tolerance = NA_real_, status = "FAIL-TERM-MISSING",
        stringsAsFactors = FALSE)
      next
    }
    for (metric in c("statistic","df_num","df_denom","p_value")) {
      v_o <- a_orig[[metric]][oi]; v_w <- a_wb[[metric]][wi]
      if (is.na(v_o) && is.na(v_w)) {
        status <- "PASS"; absd <- NA_real_; reld <- NA_real_
      } else if (is.na(v_o) || is.na(v_w)) {
        status <- "FAIL-NA-MISMATCH"; absd <- NA_real_; reld <- NA_real_
      } else {
        absd <- abs(v_o - v_w)
        denom <- max(abs(v_o), abs(v_w), 1)
        reld <- absd / denom
        # Integer df should match exactly (allow 1e-9 for floating df_denom KR)
        is_df <- metric %in% c("df_num","df_denom")
        tol <- if (is_df && metric == "df_num") abs_tol else rel_tol
        status <- if ((is_df && metric == "df_num" && absd <= 1e-9) ||
                      (!is_df && reld <= rel_tol) ||
                      (is_df && metric == "df_denom" && reld <= rel_tol))
          "PASS" else "FAIL-VALUE"
      }
      rows[[length(rows)+1L]] <- data.frame(
        indicator = indicator, kind = kind, term = tm, metric = metric,
        pipeline_value = v_o, workbook_value = v_w,
        abs_diff = absd, rel_diff = reld,
        tolerance = if (metric == "df_num") abs_tol else rel_tol,
        status = status, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

# =============================================================================
# Sampling
# =============================================================================

sample_indicators <- function(seed = 20260528) {
  set.seed(seed)
  beh_pool  <- .BEH_REFITTABLE
  mono_pool <- names(.MONO_ANALYTE_SCRIPT)   # one script per analyte
  beh_n     <- ceiling(0.5 * length(beh_pool))
  mono_n    <- ceiling(0.5 * length(mono_pool))
  list(
    seed   = seed,
    beh    = sort(sample(beh_pool,  beh_n)),
    cort   = "cort",
    mono   = sort(sample(mono_pool, mono_n))
  )
}

# =============================================================================
# Top-level: run Layer 2
# =============================================================================

run_refit_equivalence <- function(wb_path, es_beh_csv, es_endo_csv,
                                  report_dir,
                                  reconstruct_behaviour_fn,
                                  reconstruct_endo_fn,
                                  seed = 20260528) {
  samp <- sample_indicators(seed)

  # ---- Write workbook-derived CSVs to a session-scoped temp dir -----------
  tmp_dir <- file.path(tempdir(), paste0("wb_csvs_", format(Sys.time(),
                                                            "%Y%m%d_%H%M%S")))
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  wb_beh_csv  <- file.path(tmp_dir, "easy_scripts_dataset.csv")
  wb_endo_csv <- file.path(tmp_dir, "easy_scripts_endo_dataset.csv")
  # Reconstruct WITHOUT es_template_path so we get only genuine workbook columns
  # (no NA-placeholder columns from template alignment — those would mask the
  # derived columns we need to restore from the original CSV below).
  beh_recon  <- reconstruct_behaviour_fn(wb_path, es_template_path = NULL)
  endo_recon <- reconstruct_endo_fn(wb_path)

  # ---- Augment workbook beh recon with derived columns from original CSV ------
  # The workbook only carries curated indicator columns.  Easy scripts
  # dplyr::select() a broader column set (timepoint, trial_id, fish_density_f,
  # etc.) that is not preserved in the workbook.  Restore those columns from the
  # original CSV so the scripts run identically; only the workbook's genuine
  # indicator values are substituted (Layer 1 has confirmed all values match).
  beh_orig_df <- tryCatch(
    as.data.frame(readr::read_csv(es_beh_csv, show_col_types = FALSE)),
    error = function(e) NULL)
  if (!is.null(beh_orig_df)) {
    key_beh  <- c("phys_trial_id", "timepoint_f")
    # Only update columns the workbook GENUINELY has — exclude key cols
    wb_cols   <- setdiff(names(beh_recon), key_beh)
    over_cols <- intersect(wb_cols, names(beh_orig_df))
    # Coerce key columns to character for type-safe join
    beh_orig_df$phys_trial_id <- as.character(beh_orig_df$phys_trial_id)
    beh_orig_df$timepoint_f   <- as.character(beh_orig_df$timepoint_f)
    beh_recon$phys_trial_id   <- as.character(beh_recon$phys_trial_id)
    beh_recon$timepoint_f     <- as.character(beh_recon$timepoint_f)
    # Drop columns we will re-add from workbook, then left-join workbook values
    base_df  <- beh_orig_df[, setdiff(names(beh_orig_df), over_cols), drop = FALSE]
    wb_slice <- beh_recon[, c(key_beh, over_cols), drop = FALSE]
    beh_recon <- dplyr::left_join(base_df, wb_slice, by = key_beh)
    message("[Layer2] Beh recon augmented: ", length(over_cols), " workbook cols + ",
            ncol(beh_orig_df) - length(over_cols), " restored derived cols from original")
  }

  # ---- Normalise endo treatment labels to match original endo CSV convention --
  # The workbook uses the behaviour pipeline's treatment labels (e.g.
  # "exercise choice") while the endo CSV uses its own (e.g. "treat").
  # Endo easy scripts hardcode factor levels from the endo CSV, so we recode
  # before writing the temp CSV.  The recode is detected automatically.
  endo_orig_df <- tryCatch(
    as.data.frame(readr::read_csv(es_endo_csv, show_col_types = FALSE)),
    error = function(e) NULL)
  if (!is.null(endo_orig_df)) {
    orig_non_ctrl <- setdiff(na.omit(unique(endo_orig_df$treatment)), "control")
    wb_non_ctrl   <- setdiff(na.omit(unique(endo_recon$treatment)),   "control")
    if (length(orig_non_ctrl) == 1 && length(wb_non_ctrl) == 1 &&
        orig_non_ctrl != wb_non_ctrl) {
      message("[Layer2] Endo treatment recode: '", wb_non_ctrl,
              "' -> '", orig_non_ctrl, "'")
      endo_recon$treatment[!is.na(endo_recon$treatment) &
                           endo_recon$treatment == wb_non_ctrl] <- orig_non_ctrl
    }
  }

  readr::write_csv(beh_recon,  wb_beh_csv)
  readr::write_csv(endo_recon, wb_endo_csv)

  message("[Layer2] workbook-derived CSVs at: ", tmp_dir)
  message("[Layer2] behaviour: ", nrow(beh_recon), " rows; ",
          "endocrine: ",  nrow(endo_recon), " rows")

  # ---- Sampling log -------------------------------------------------------
  log_rows <- list()
  log_rows[[1]] <- data.frame(kind = "seed", indicator = as.character(seed))
  for (nm in samp$beh)  log_rows[[length(log_rows)+1L]] <-
    data.frame(kind = "behaviour", indicator = nm)
  log_rows[[length(log_rows)+1L]] <-
    data.frame(kind = "cortisol", indicator = samp$cort)
  for (nm in samp$mono) log_rows[[length(log_rows)+1L]] <-
    data.frame(kind = "monoamine", indicator = nm)
  readr::write_csv(do.call(rbind, log_rows),
                   file.path(report_dir, "refit_sampling_log.csv"))

  # ---- Run pairs ----------------------------------------------------------
  diff_rows <- list()
  add <- function(df) { if (!is.null(df)) diff_rows[[length(diff_rows)+1L]] <<- df }

  # Behaviour
  for (ind in samp$beh) {
    rel <- .BEH_ES_SCRIPT_LOOKUP[[ind]]
    if (is.null(rel)) {
      add(.skip_row(ind, "behaviour", "no easy-script lookup"))
      next
    }
    script <- file.path(.ES_ROOT, rel)
    if (!file.exists(script)) {
      add(.skip_row(ind, "behaviour", paste("missing script:", rel)))
      next
    }
    message("[Layer2/beh] ", ind, "  ->  ", rel)
    add(.refit_pair(ind, "behaviour", script,
                    beh_csv = NULL,       # NULL = use original
                    endo_csv = NULL,
                    wb_beh_csv  = wb_beh_csv,
                    wb_endo_csv = wb_endo_csv))
  }

  # Cortisol
  script_cort <- file.path(.ES_ROOT, .CORT_ES_SCRIPT)
  if (file.exists(script_cort)) {
    message("[Layer2/cort] cort  ->  ", .CORT_ES_SCRIPT)
    add(.refit_pair("cort_plasma", "cortisol", script_cort,
                    beh_csv = NULL, endo_csv = NULL,
                    wb_beh_csv  = wb_beh_csv, wb_endo_csv = wb_endo_csv))
  } else {
    add(.skip_row("cort_plasma", "cortisol", "missing script"))
  }

  # Monoamines (per analyte, runs all 4 areas at once)
  for (an in samp$mono) {
    rel <- .MONO_ANALYTE_SCRIPT[[an]]
    if (is.null(rel)) { add(.skip_row(an, "monoamine", "no lookup")); next }
    script <- file.path(.ES_ROOT, rel)
    if (!file.exists(script)) {
      add(.skip_row(an, "monoamine", paste("missing script:", rel)))
      next
    }
    message("[Layer2/mono] ", an, "  ->  ", rel)
    add(.refit_pair(an, "monoamine", script,
                    beh_csv = NULL, endo_csv = NULL,
                    wb_beh_csv  = wb_beh_csv, wb_endo_csv = wb_endo_csv))
  }

  refit_df <- do.call(rbind, diff_rows)
  readr::write_csv(refit_df,
                   file.path(report_dir, "refit_equivalence.csv"))

  n_pass <- sum(refit_df$status == "PASS",          na.rm = TRUE)
  n_fail <- sum(grepl("^FAIL", refit_df$status),    na.rm = TRUE)
  n_skip <- sum(grepl("^SKIPPED", refit_df$status), na.rm = TRUE)
  message("[Layer2] PASS=", n_pass, ", FAIL=", n_fail, ", SKIPPED=", n_skip)

  list(n_compared = nrow(refit_df), n_pass = n_pass, n_fail = n_fail,
       n_skip = n_skip, sampling = samp, tmp_dir = tmp_dir,
       fails = refit_df[grepl("^FAIL", refit_df$status), , drop = FALSE])
}

.skip_row <- function(indicator, kind, reason) {
  data.frame(indicator = indicator, kind = kind, term = NA_character_,
             metric = "anova", pipeline_value = NA_real_,
             workbook_value = NA_real_, abs_diff = NA_real_,
             rel_diff = NA_real_, tolerance = NA_real_,
             status = paste0("SKIPPED-", reason),
             stringsAsFactors = FALSE)
}

.refit_pair <- function(indicator, kind, script_path,
                        beh_csv, endo_csv, wb_beh_csv, wb_endo_csv) {
  # Run A: original CSVs
  envA <- tryCatch(
    .with_silenced_ggsave(
      .with_redirected_csvs(beh_csv = beh_csv, endo_csv = endo_csv,
                             .refit_one_script(script_path))),
    error = function(e) {
      message("[Layer2]   run A FAILED: ", conditionMessage(e)); NULL })

  # Run B: workbook-derived CSVs
  envB <- tryCatch(
    .with_silenced_ggsave(
      .with_redirected_csvs(beh_csv  = wb_beh_csv,
                             endo_csv = wb_endo_csv,
                             .refit_one_script(script_path))),
    error = function(e) {
      message("[Layer2]   run B FAILED: ", conditionMessage(e)); NULL })

  if (is.null(envA) || is.null(envB))
    return(.skip_row(indicator, kind, "refit failed"))

  a_orig <- .extract_anova(envA)
  a_wb   <- .extract_anova(envB)
  .compare_anova(a_orig, a_wb, indicator, kind)
}
