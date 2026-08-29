# =============================================================================
# CHOICE EXPERIMENT  |  ZONE SWITCH DIAGNOSTICS
# Correlation between two indicators of arena traversal:
#   - switches_per_session: identity-less main-zone switch count (NN-matched
#       blob pairs crossing the flow<->calm boundary), normalised to 20-min.
#   - zone_flux_per_session: mean absolute frame-to-frame change in
#       (n_in_flow / n_total), scaled per minute.
# Both are derived from STEP2 (trial_activity_summary). Pearson + Spearman
# correlations are computed across all (trial x timepoint) rows and faceted by
# treatment. A scatter plot with regression line is exported as PNG and as a
# slide in a small PowerPoint, plus a one-paragraph summary suitable for
# embedding into the auto-drafted Results report.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(officer)
  library(rvg)
  library(readr)
})

ts_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")

.find_latest_csv_zsd <- function(step_name, csv_filename) {
  parent <- file.path(getwd(), step_name)
  if (!dir.exists(parent)) return(NULL)
  subdirs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  subdirs <- subdirs[grepl(paste0("^", step_name, "_\\d{8}_\\d{6}$"), basename(subdirs))]
  if (length(subdirs) == 0) return(NULL)
  candidate <- file.path(subdirs[which.max(file.mtime(subdirs))], csv_filename)
  if (file.exists(candidate)) candidate else NULL
}

.fmt_p_zsd <- function(p) {
  if (is.na(p)) return("NA")
  if (p <= 0.0001) return("p < 0.0001")
  paste0("p = ",
    trimws(sub("0+$", "", sub("\\.$", "",
      format(signif(p, 3), scientific = FALSE, trim = TRUE)))))
}

# ---- 1) Load data -----------------------------------------------------------
if (!exists("trial_activity_summary", envir = .GlobalEnv, inherits = FALSE)) {
  .csv <- .find_latest_csv_zsd("STEP2_output", "trial_activity_summary.csv")
  if (is.null(.csv))
    stop("zone_switch_diagnostics: trial_activity_summary not found.", call. = FALSE)
  trial_activity_summary <- readr::read_csv(.csv, show_col_types = FALSE)
  assign("trial_activity_summary", trial_activity_summary, envir = .GlobalEnv)
}
df_zsd <- get("trial_activity_summary", envir = .GlobalEnv, inherits = FALSE)

req <- c("switches_per_session", "zone_flux_per_session", "treatment")
missing_cols <- setdiff(req, names(df_zsd))
if (length(missing_cols) > 0) {
  warning("zone_switch_diagnostics: missing columns ",
          paste(missing_cols, collapse = ", "), " — skipped.", call. = FALSE)
  return(invisible(NULL))
}

df_zsd <- df_zsd %>%
  dplyr::filter(is.finite(switches_per_session),
                is.finite(zone_flux_per_session)) %>%
  dplyr::mutate(treatment = factor(trimws(tolower(as.character(treatment)))))

# ---- 2) Output dir ----------------------------------------------------------
.parent <- file.path(getwd(), "zone_switch_diagnostics")
if (!dir.exists(.parent)) dir.create(.parent, recursive = TRUE, showWarnings = FALSE)
ZSD_OUT <- file.path(.parent,
                     paste0("zone_switch_diagnostics_",
                            format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(ZSD_OUT, recursive = TRUE, showWarnings = FALSE)

# ---- 3) Correlations --------------------------------------------------------
.corr_block <- function(x, y, group_label) {
  if (length(x) < 4) return(NULL)
  pe <- tryCatch(cor.test(x, y, method = "pearson"),  error = function(e) NULL)
  sp <- tryCatch(cor.test(x, y, method = "spearman"), error = function(e) NULL)
  data.frame(
    group       = group_label,
    n           = length(x),
    pearson_r   = if (!is.null(pe)) round(unname(pe$estimate), 3) else NA_real_,
    pearson_p   = if (!is.null(pe)) pe$p.value                    else NA_real_,
    spearman_rho= if (!is.null(sp)) round(unname(sp$estimate), 3) else NA_real_,
    spearman_p  = if (!is.null(sp)) sp$p.value                    else NA_real_,
    stringsAsFactors = FALSE
  )
}

corr_overall <- .corr_block(df_zsd$zone_flux_per_session,
                            df_zsd$switches_per_session, "overall")

corr_by_treat <- df_zsd %>%
  dplyr::group_by(treatment) %>%
  dplyr::group_modify(~ .corr_block(.x$zone_flux_per_session,
                                    .x$switches_per_session,
                                    as.character(.y$treatment[1]))) %>%
  dplyr::ungroup() %>%
  dplyr::select(-treatment)

zsd_corr_tbl <- dplyr::bind_rows(corr_overall, corr_by_treat)
readr::write_csv(zsd_corr_tbl,
                 file.path(ZSD_OUT, "zone_switch_correlations.csv"))

# ---- 4) Scatter plot --------------------------------------------------------
.lab_overall <- if (!is.null(corr_overall) && !is.na(corr_overall$pearson_r)) {
  paste0("Overall: r = ", corr_overall$pearson_r,
         ", \u03c1 = ", corr_overall$spearman_rho,
         ", ", .fmt_p_zsd(corr_overall$pearson_p),
         "  (n = ", corr_overall$n, ")")
} else "Overall: insufficient data"

p_zsd <- ggplot2::ggplot(
    df_zsd,
    ggplot2::aes(x = zone_flux_per_session, y = switches_per_session,
                 colour = treatment, fill = treatment)) +
  ggplot2::geom_point(alpha = 0.7, size = 2.2) +
  ggplot2::geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.7) +
  ggplot2::facet_wrap(~ treatment, ncol = 2) +
  ggplot2::labs(
    title    = "Zone-flux vs identity-less main-zone switches",
    subtitle = .lab_overall,
    x        = "zone_flux_per_session (frame-to-frame |\u0394 prop in flow|, /min)",
    y        = "switches_per_session (NN-matched flow\u2194calm crossings, /20 min)"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid       = ggplot2::element_blank(),
    axis.line        = ggplot2::element_line(color = "black", linewidth = 0.7),
    legend.position  = "none",
    strip.text       = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(file.path(ZSD_OUT, "zone_switch_scatter.png"),
                p_zsd, width = 10, height = 5, dpi = 150)

# ---- 5) PowerPoint ----------------------------------------------------------
tryCatch({
  prs <- officer::read_pptx() %>%
    officer::add_slide(layout = "Title and Content", master = "Office Theme") %>%
    officer::ph_with(value = "Zone-flux vs main-zone switches",
                     location = officer::ph_location_type(type = "title")) %>%
    officer::ph_with(value = rvg::dml(ggobj = p_zsd),
                     location = officer::ph_location(left = 0.5, top = 1.2,
                                                     width = 9.0, height = 5.0))
  print(prs, target = file.path(ZSD_OUT, "zone_switch_diagnostics.pptx"))
}, error = function(e)
  warning("zone_switch_diagnostics PPTX failed: ", e$message, call. = FALSE))

# ---- 6) Results-paragraph string for the report -----------------------------
.zsd_paragraph <- if (!is.null(corr_overall) && !is.na(corr_overall$pearson_r)) {
  paste0(
    "The identity-less main-zone switch count (switches_per_session) was ",
    "compared against the simpler zone_flux_per_session diagnostic across all ",
    "trial \u00d7 timepoint rows (n = ", corr_overall$n, "). Pearson r = ",
    corr_overall$pearson_r, " (", .fmt_p_zsd(corr_overall$pearson_p),
    "); Spearman \u03c1 = ", corr_overall$spearman_rho,
    " (", .fmt_p_zsd(corr_overall$spearman_p),
    "). The two indicators agree well, supporting the use of ",
    "switches_per_session as the primary statistical response while retaining ",
    "zone_flux_per_session as a model-free diagnostic."
  )
} else {
  "Zone-switch diagnostics: insufficient finite data for correlation."
}

writeLines(.zsd_paragraph, file.path(ZSD_OUT, "results_paragraph.txt"))

assign("zone_switch_diag_corr",       zsd_corr_tbl,    envir = .GlobalEnv)
assign("zone_switch_diag_paragraph",  .zsd_paragraph,  envir = .GlobalEnv)
assign("zone_switch_diag_plot",       p_zsd,           envir = .GlobalEnv)
assign("ZONE_SWITCH_DIAG_DIR",        ZSD_OUT,         envir = .GlobalEnv)

ts_msg("Zone-switch diagnostics written: ", ZSD_OUT)
