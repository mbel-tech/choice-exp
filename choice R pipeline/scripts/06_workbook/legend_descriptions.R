# =============================================================================
# legend_descriptions.R — hard-coded column descriptions
# =============================================================================
# Lookups consulted in this order:
#   1) METADATA_DESCRIPTIONS  — exact column-name match for IDs / metadata
#   2) INDICATOR_MAP_BEH / INDICATOR_MAP_ENDO — for indicator columns
#   3) Heuristic patterns                     — for source-file columns
# Anything not resolved gets description = "Needs review" — never guess.
# =============================================================================

# ---- Identifier / metadata columns -----------------------------------------
METADATA_DESCRIPTIONS <- list(
  # behaviour_wide ---------
  phys_trial      = list(type = "identifier", unit = "integer (1-16)",
    description = "Trial identifier (one holding tank x one date x one treatment). 16 trials; this is the numbering used throughout the manuscript and its figures. NOT chronological -- see trial_seq."),
  trial_seq       = list(type = "metadata",   unit = "integer (1-16)",
    description = "The same 16 trials in chronological running order, rebuilt from date and, within a day, from holding density (the fuller tank is always emptied first). It differs from phys_trial: trials 1 and 2, 13 and 14, and 15 and 16 are swapped between the two numberings."),
  interval        = list(type = "metadata",   unit = "factor (3 levels)",
    description = "Within-trial observation interval: 1, 2 or 3, in order. Each trial contributes three rows."),
  trial_id        = list(type = "identifier", unit = "text",
    description = "Original trial identifier from the source spreadsheet (concatenates video ID and metadata)."),
  # 2026-08-29: the separate timepoint_f / timepoint / timepoint_num entries
  # were removed. The source frames carry all three spellings of one variable;
  # the workbook and the deposit carry one column, `interval`, described above.
  obs_seconds     = list(type = "qc",         unit = "s",
    description = "Length of the scored interval. Denominator for every rate in this sheet."),
  n_frames        = list(type = "qc",         unit = "count",
    description = "Number of video frames contributing to the interval, after the tracking-quality filter."),
  trial_date      = list(type = "metadata",   unit = "date (DD.MM.YYYY)",
    description = "Date the trial was recorded. 8 distinct dates across the experiment."),
  treatment       = list(type = "metadata",   unit = "factor (2 levels)",
    description = "Experimental condition: 'control' or 'exercise choice'."),
  tank            = list(type = "metadata",   unit = "factor (4 levels)",
    description = "Holding tank ID: 27, 28, 29, or 31. Each tank held one group of fish; each date is fully nested in one tank."),
  fish_density    = list(type = "metadata",   unit = "count (fish)",
    description = "Number of fish in the school for that trial: 4, 8, 12, or 16."),
  motor_side      = list(type = "metadata",   unit = "factor",
    description = "Side of the tank where the motor (flow generator) was positioned during the trial."),
  fish_density_f  = list(type = "metadata",   unit = "factor (4 levels)",
    description = "Fish density as a factor (4 / 8 / 12 / 16)."),

  # endo_wide --------------
  sample_id       = list(type = "identifier", unit = "text",
    description = "Sample identifier, e.g. 'B3_01'. One value per individual fish sample (cortisol + monoamine assays share this ID)."),
  trial_key       = list(type = "metadata",   unit = "text",
    description = "Experimental unit this fish came from, coded '<holding tank>_<number of fish>', e.g. '27_16'. This pair is the one key every dataset in the study shares, and it maps one-to-one onto the 16 behavioural trials. Use phys_trial to join to the behaviour sheet."),
  trial           = list(type = "metadata",   unit = "text",
    description = "Superseded 2026-08-29. Earlier versions of this workbook carried an INTEGER trial column running 1-19 in the two endocrine sheets. It was not a trial identifier -- it was the N-index of the sampling video -- and joining behaviour to physiology on it produced wrong answers. Replaced by trial_key and phys_trial."),
  sex             = list(type = "metadata",   unit = "factor (2 levels)",
    description = "Fish sex: 'F' or 'M'."),
  plate           = list(type = "metadata",   unit = "factor",
    description = "Cortisol assay plate ID (5, 6, 7, 8, Evg). NA for monoamine samples."),
  plate_date      = list(type = "metadata",   unit = "date",
    description = "Date the cortisol assay plate was run. NA for monoamine samples.")
)

# ---- Occupancy context columns ----------------------------------------------
# Model INPUTS rather than model outputs, so they have no manifest row. They are
# deposited so a reader can recompute every log-ratio, CLR and Jacobs' D from
# the workbook instead of having to trust ours.
#
# The area correction is applied ONCE, in STEP2 section 6b, where each sub-zone
# proportion is divided by its area constant (high/10, medium/18, low/13,
# calm/42) and the four parts rescaled to sum to 1. Until 2026-08-18 the stats
# script divided by those constants a second time, which put a zone-specific
# offset on every sub-zone log-ratio; values in this workbook are post-fix.
.OCCUPANCY_DESCRIPTIONS <- list(
  prop_flow      = list(type = "indicator", unit = "proportion (0-1)",
    description = "Proportion of the interval the school spent in the flow zone, uncorrected for zone area. prop_flow + prop_calm does not equal 1: the remainder is time not assignable to either zone (fish untracked, or between polygons)."),
  prop_calm      = list(type = "indicator", unit = "proportion (0-1)",
    description = "Proportion of the interval spent in the calm zone, uncorrected for zone area. See prop_flow on why the two do not sum to 1."),
  prop_high_ac   = list(type = "indicator", unit = "proportion (0-1)",
    description = "Area-corrected occupancy of the high-velocity sub-zone. The four prop_*_ac columns are a closed composition summing to 1, and are the basis of alr_high/medium/low and clr_high/medium/low/calm."),
  prop_medium_ac = list(type = "indicator", unit = "proportion (0-1)",
    description = "Area-corrected occupancy of the medium-velocity sub-zone. Part of the four-part composition summing to 1."),
  prop_low_ac    = list(type = "indicator", unit = "proportion (0-1)",
    description = "Area-corrected occupancy of the low-velocity sub-zone. Part of the four-part composition summing to 1."),
  prop_calm_ac   = list(type = "indicator", unit = "proportion (0-1)",
    description = "Area-corrected occupancy of the calm sub-zone. The reference part of the additive log-ratios: alr_high = log(prop_high_ac / prop_calm_ac), with a Haldane floor of 1e-4 applied to both parts before the ratio is taken."),
  clr_high       = list(type = "derived",   unit = "CLR",
    description = "Centred log-ratio of prop_high_ac against the geometric mean of all four area-corrected sub-zone proportions."),
  clr_medium     = list(type = "derived",   unit = "CLR",
    description = "Centred log-ratio of prop_medium_ac against the geometric mean of the four-part composition."),
  clr_low        = list(type = "derived",   unit = "CLR",
    description = "Centred log-ratio of prop_low_ac against the geometric mean of the four-part composition."),
  clr_calm       = list(type = "derived",   unit = "CLR",
    description = "Centred log-ratio of prop_calm_ac against the geometric mean of the four-part composition.")
)
METADATA_DESCRIPTIONS <- c(METADATA_DESCRIPTIONS, .OCCUPANCY_DESCRIPTIONS)

# ---- Overrides for indicator columns whose auto-built text is too thin ------
# build_indicator_descriptions() writes "Indicator '<label>'." for every mapped
# column. For a deposited dataset that is not documentation. These replace it
# for the columns where the name alone is genuinely ambiguous or misleading.
.INDICATOR_OVERRIDES <- list(
  logit_flow            = "Additive log-ratio of flow-zone against calm-zone occupancy, ALR(flow vs calm). Reported outcome 1. Equal to the logit of prop_flow. Shipped as 'logit_flow' before 2026-08-29, which invited readers to treat it as a plain logit rather than as the compositional coordinate it is.",
  lr_high               = "Additive log-ratio of high-velocity sub-zone against calm sub-zone occupancy, ALR(high vs calm) = log(prop_high_ac / prop_calm_ac) with a Haldane floor of 1e-4. Reported outcome 2.",
  lr_medium             = "Additive log-ratio of the medium-velocity sub-zone against the calm sub-zone, ALR(medium vs calm). Reported outcome 3.",
  lr_low                = "Additive log-ratio of the low-velocity sub-zone against the calm sub-zone, ALR(low vs calm). Reported outcome 4.",
  zone_flux_per_session = "Flow-calm crossings of the school during the interval. THIS is the crossings measure reported in the Results (F(1,14) = 45.00, p < 0.001), not switches_per_session.",
  switches_per_session  = "School transitions between the two main zones during the interval, from the zone-switch detector. A near-duplicate of crossings_per_session and NOT the measure reported in the Results; kept as a descriptive cross-check. Rounded to an integer: the source carries float noise from a rate x duration round-trip.",
  bout_rate_flow        = "Flow bouts begun per minute of the interval. Reported outcome 6. A bout is a run in which more than half the tracked fish are in the flow zone, on a 1 s binned collective state sequence.",
  max_flow_bout_s       = "Duration of the longest single flow bout in the interval. Reported outcome 7. Taken over the full bout inventory, including bouts truncated by the interval boundary.",
  max_calm_bout_s       = "Duration of the longest single calm bout in the interval. Reported outcome 8.",
  dwell_Flow            = "Mean duration of a flow bout in the interval. Reported outcome 9; it is the mean of the same distribution whose maximum longest_flow_bout_s reports. Empty where the school never completed a flow bout within the interval.",
  dwell_Calm            = "Mean duration of a calm bout in the interval. Reported outcome 10. Empty where the school never completed a calm bout within the interval.",
  mean_nnd_cm           = "Mean nearest-neighbour distance across the interval. Reported outcome 11.",
  mean_iid_cm           = "Mean inter-individual distance across the interval. Reported outcome 12.",
  mean_hull_area_cm2    = "Mean area of the convex hull of the school across the interval. Reported outcome 13.",
  mean_centroid_spd_cm  = "Mean speed of the school centroid across the interval. Reported outcome 14.",
  D_main_flow           = "Jacobs' D selectivity index for the flow zone, comparing occupancy against the fraction of arena area the zone occupies. Bounded -1 (complete avoidance) to +1 (complete selection). Supplementary validation path, not a primary outcome.",
  D_sec_high            = "Jacobs' D selectivity index for the high-velocity sub-zone. Supplementary validation path.",
  D_sec_medium          = "Jacobs' D selectivity index for the medium-velocity sub-zone. Supplementary validation path.",
  D_sec_low             = "Jacobs' D selectivity index for the low-velocity sub-zone. Supplementary validation path."
)

# ---- Indicator column descriptions (auto-built from INDICATOR_MAP_*) -------
build_indicator_descriptions <- function(map_beh, map_endo) {
  if (is.null(map_beh) || is.null(map_endo))
    stop("Both INDICATOR_MAP_BEH and INDICATOR_MAP_ENDO must be available", call. = FALSE)
  out_beh <- vector("list", nrow(map_beh))
  names(out_beh) <- map_beh$wb_col
  for (i in seq_len(nrow(map_beh))) {
    .lbl <- map_beh$manifest_label[i]
    .dv  <- map_beh$derive[i]
    .nt  <- switch(.dv,
      "interval"       = "Per-interval value (one row per phys_trial x interval).",
      "broadcast_mean" = "Trial-level aggregate (mean across the 3 intervals). Value is broadcast across the 3 interval rows of the same phys_trial.",
      "representative" = "Representative cell (e.g. flow / high) from a multi-cell zone-occupancy model. The underlying inferential model fits ALL cells jointly; this column shows the headline cell.",
      "direct"         = "Trial-level value as provided by the source assembly.",
      "missing"        = "Source data not available in the current pipeline run; column is all-NA.",
      "")
    .es <- map_beh$es_col[i]
    .ov <- if (!is.na(.es) && nzchar(.es) && !is.null(.INDICATOR_OVERRIDES[[.es]]))
             .INDICATOR_OVERRIDES[[.es]] else NULL
    out_beh[[i]] <- list(
      type           = "indicator",
      unit           = map_beh$unit[i],
      indicator_name = .lbl,
      description    = if (is.null(.ov)) sprintf("Indicator '%s'.", .lbl) else .ov,
      notes          = .nt
    )
  }
  out_endo <- vector("list", nrow(map_endo))
  names(out_endo) <- map_endo$wb_col
  for (i in seq_len(nrow(map_endo))) {
    .lbl <- map_endo$manifest_label[i]
    out_endo[[i]] <- list(
      type           = "indicator",
      unit           = map_endo$unit[i],
      indicator_name = .lbl,
      description    = sprintf("Indicator '%s'.", .lbl),
      notes          = if (map_endo$wb_col[i] == "cort_plasma")
                          "Plasma cortisol, ng/mL, one value per sampled fish." else
                          "Tissue monoamine concentration per brain region, ng/mg tissue. Noradrenaline is quantified but not deposited: it is not entered into any reported model, so the analysed grid is 24 analyte x region cells, not 28."
    )
  }
  list(beh = out_beh, endo = out_endo)
}

# ---- Source-file column descriptions (heuristic) ---------------------------
# Used for src_* sheets. Returns "Needs review" when no rule matches.
SRC_COLUMN_HEURISTICS <- list(
  # exact-name → description
  exact = list(
    video_ID        = "Video identifier; key to behavioural recording.",
    Trial           = "Trial number.",
    trial           = "Trial number.",
    tank            = "Holding tank ID.",
    side            = "Tank side where stimulus/motor was placed.",
    `motor side`    = "Tank side where the flow motor was placed.",
    motor_side      = "Tank side where the flow motor was placed.",
    density         = "Number of fish in the trial.",
    fish_density    = "Number of fish in the trial.",
    time            = "Trial start time (clock time of recording).",
    interval        = "Within-trial interval (1, 2, 3).",
    treatment       = "Experimental condition.",
    condition       = "Experimental condition (alternative naming for 'treatment').",
    experience      = "Fish prior-experience flag.",
    date            = "Trial date.",
    trial_date      = "Trial date.",
    SAMPLE_ID       = "Sample identifier (one per fish).",
    Sample_ID       = "Sample identifier (one per fish).",
    sample_id       = "Sample identifier (one per fish).",
    sample          = "Sample identifier (one per fish).",
    Sample          = "Sample identifier (one per fish).",
    FISH_ID         = "Fish identifier.",
    sex             = "Fish sex (F/M).",
    weight          = "Fish wet weight.",
    `wt (g)`        = "Fish wet weight (grams).",
    wt              = "Fish wet weight.",
    lenght          = "Fish standard length (note: source uses misspelling 'lenght').",
    `lenght (cm)`   = "Fish standard length, cm (note: source uses misspelling 'lenght').",
    length          = "Fish standard length.",
    batch           = "Experimental batch (B1 / B2 / B3).",
    Sample          = "Sample identifier.",
    Dilution        = "Assay dilution factor.",
    mucus_c         = "Mucus cortisol concentration.",
    `m-plate`       = "Mucus assay plate ID.",
    plate           = "Assay plate ID.",
    plate_date      = "Assay plate run date.",
    area            = "Brain area: DM / POA / VV / VD.",
    Sample_ID       = "Sample identifier.",
    `5-HT`          = "Serotonin concentration (ng/mg tissue).",
    `5-HIAA`        = "5-Hydroxyindoleacetic acid concentration (ng/mg tissue).",
    `5-HIAA/5-HT`   = "5-HIAA to 5-HT ratio (serotonin turnover index).",
    DA              = "Dopamine concentration (ng/mg tissue).",
    DOPAC           = "3,4-Dihydroxyphenylacetic acid concentration (ng/mg tissue).",
    `DOPAC/DA`      = "DOPAC to DA ratio (dopamine turnover index).",
    NE              = "Norepinephrine concentration (ng/mg tissue).",
    notes           = "Free-text notes from the source spreadsheet.",
    ROI_index       = "Region-of-interest index in the zone reference layout.",
    ROI_name        = "Region-of-interest name (e.g. flow / calm / high / medium / low).",
    Vertex_index    = "Vertex index within the ROI polygon.",
    X               = "X coordinate of polygon vertex (pixels).",
    Y               = "Y coordinate of polygon vertex (pixels).",
    high            = "Time spent in high-flow zone (snapshot fraction).",
    medium          = "Time spent in medium-flow zone (snapshot fraction).",
    low             = "Time spent in low-flow zone (snapshot fraction).",
    calm            = "Time spent in calm zone (snapshot fraction).",
    flow            = "Time spent in flow zone (snapshot fraction).",
    tot             = "Total observed time (snapshot scaffold).",
    timestamp       = "Frame timestamp.",
    VIDEO_ID        = "Video identifier (uppercase variant)."
  ),
  # prefix → description (matched if exact missed)
  prefix = list(
    `mucus_`       = "Mucus assay scaffolding / measurement.",
    `cortisol_`    = "Cortisol assay scaffolding / measurement.",
    `cort_`        = "Cortisol-related quantity.",
    `mono_`        = "Monoamine-related quantity.",
    `tank_`        = "Tank metadata."
  )
)

describe_src_column <- function(colname) {
  if (is.null(colname) || is.na(colname) || !nzchar(colname))
    return(list(type = "unknown", description = "Needs review", notes = ""))
  # Unnamed Excel placeholders like "...18"
  if (grepl("^\\.\\.\\.\\d+$", colname))
    return(list(type = "assay_scaffolding",
                description = "Unnamed source column (empty header in original spreadsheet).",
                notes = "Needs review"))
  ex <- SRC_COLUMN_HEURISTICS$exact[[colname]]
  if (!is.null(ex))
    return(list(type = "metadata", description = ex, notes = ""))
  for (pfx in names(SRC_COLUMN_HEURISTICS$prefix)) {
    if (startsWith(colname, pfx))
      return(list(type = "assay_scaffolding",
                  description = SRC_COLUMN_HEURISTICS$prefix[[pfx]],
                  notes = ""))
  }
  list(type = "unknown", description = "Needs review", notes = "")
}
