# easy_scripts/_re_candidates.R
# -----------------------------------------------------------------------------
# Frozen copies of the random-effect candidate sets and global constants
# from activity_analysis_STATS_choice_exp.R (§3 + top-of-file constants).
# Mini-scripts source this file directly so they can be inspected in
# isolation without sourcing the full STATS file.
# -----------------------------------------------------------------------------

# Factor levels (mirror lines 180-186 of STATS file)
TREATMENT_LEVELS_g  <- c("control", "exercise choice")
DENSITY_LEVELS_g    <- c(4L, 8L, 12L, 16L)
TIMEPOINT_LEVELS_g  <- 1:3
PIPELINE_FPS        <- 25

# Zone area constants (cm² per main + sub zone — used by alr / CLR cell-mean)
ZONE_AREA_UNITS <- c(high = 10, medium = 18, low = 13, calm = 42)

# ---- Random-effect candidate sets (mirror STATS §3, lines 586-617) ----------

# Trial × timepoint grain (N = 48 sessions, 3 per physical trial)
RE_WIDE <- c(
  phys_trial   = "(1 | phys_trial_id)",
  tank         = "(1 | tank)",
  trial_tank   = "(1 | phys_trial_id) + (1 | tank)",
  date         = "(1 | trial_date)",
  trial_date_r = "(1 | phys_trial_id) + (1 | trial_date)",
  density      = "(1 | fish_density_f)"
)

# Long-format cell-mean (zone rows linked within physical trial)
RE_LONG <- c(
  phys_trial   = "(1 | phys_trial_id)",
  tank         = "(1 | tank)",
  trial_tank   = "(1 | phys_trial_id) + (1 | tank)",
  nested_tank  = "(1 | tank / phys_trial_id)"
)

# Aggregated grain (N = 16 physical trials)
RE_WIDE_AGG <- c(
  tank    = "(1 | tank)"
)
RE_LONG_AGG <- c(
  tank    = "(1 | tank)"
)

# ---- Endocrine RE candidates (mirror analysis_b3.R §4) ----------------------
RE_CORT_CAND <- list(
  "(1|tank)"            = "(1|tank)",
  "(1|trial)"           = "(1|trial)",
  "(1|plate)"           = "(1|plate)",
  "(1|tank)+(1|trial)"  = "(1|tank) + (1|trial)",
  "(1|tank)+(1|plate)"  = "(1|tank) + (1|plate)",
  "(1|trial)+(1|plate)" = "(1|trial) + (1|plate)"
)

# ---- Plot constants (frozen copies from STATS §6) ---------------------------
# These were originally computed via intermediate vars (e.g. `.pal_tmp` is
# `rm()`ed before the easy_scripts bridge parses the file), so we hard-code the
# resolved values here for robustness.

# Okabe palette pre-resolved (STATS line 3740-3743)
.pal_okabe <- c("#E69F00","#56B4E9","#009E73","#F0E442",
                "#0072B2","#D55E00","#CC79A7","#000000")
TREATMENT_COLORS_PAL <- c("control" = .pal_okabe[6],
                          "exercise choice" = .pal_okabe[5])

TREATMENT_COLORS <- c("control" = "#2166AC",
                      "exercise choice" = "#D6604D")
TREATMENT_SHAPES <- c("control" = 17, "exercise choice" = 16)

color_map        <- c(flow = "red4", calm = "lightskyblue3",
                      high = "tomato3", medium = "mediumturquoise",
                      low  = "goldenrod3")
color_map_broad  <- c(flow = "red4", calm = "lightskyblue3")

TIMEPOINT_LABELS <- c("1" = "5–25 min",
                       "2" = "45–65 min",
                       "3" = "85–105 min")
LINE_TP_LABELS   <- TIMEPOINT_LABELS

invisible(NULL)
