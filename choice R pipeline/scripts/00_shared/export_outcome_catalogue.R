# =============================================================================
# export_outcome_catalogue.R
# =============================================================================
# Catalogue of EVERY behavioural outcome analysed at the trial x interval level:
# what it measures, how it is computed, and its Type III ANOVA for the treatment
# main effect and the treatment x interval interaction.
#
# Interval is a 3-level factor in every model, so interval and interaction terms
# are 2-df tests throughout (METHODS_CHANGES.md section 2). Where df1 = 1 appears
# it is rank deficiency from empty treatment x interval cells, not a
# parameterisation difference; those rows are flagged.
#
# Markdown is the source; the .docx is produced from it via Pandoc so the two
# cannot diverge.
# =============================================================================

PIPE  <- file.path(PROJECT_ROOT, "choice R pipeline")
OUT   <- file.path(PROJECT_ROOT, "interval_contrasts_export")
if (!exists("PANDOC")) source(file.path(PROJECT_ROOT, "config.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Latest run that is actually COMPLETE. `sentinel` is a file written late in the
# run; a directory lacking it is still being written, and picking it up silently
# yields a catalogue full of "—" for everything the run had not yet reached.
.latest <- function(step, sentinel = NULL) {
  p <- file.path(PIPE, "output", step)
  d <- list.dirs(p, full.names = TRUE, recursive = FALSE)
  d <- d[grepl(paste0("^", step, "_\\d{8}_\\d{6}$"), basename(d))]
  if (!is.null(sentinel)) {
    ok <- vapply(d, function(x) file.exists(file.path(x, sentinel)), logical(1))
    if (any(!ok))
      message("skipping ", sum(!ok), " incomplete ", step, " run(s): ",
              paste(basename(d[!ok]), collapse = ", "))
    d <- d[ok]
  }
  if (!length(d)) stop("No complete ", step, " run found under ", p)
  d[which.max(file.mtime(d))]
}
STEP5 <- .latest("STEP5_stats",  "centroid_speed_timepoint/anova.csv")
SEQ   <- .latest("SEQ_output",   "seq_anova_trial_x_timepoint.csv")
BOUT  <- .latest("BOUT_output",  "bout_anova_trial_x_timepoint.csv")
cat("STEP5:", STEP5, "\nSEQ  :", SEQ, "\nBOUT :", BOUT, "\n")

# ---- registry: what each outcome measures and how it is computed ------------
# key = engine::identifier
# `interp` = how to read the number (which direction means what), so the table
# is usable without going back to the code.
# `final` = panel in the final figure set, or "" if not reported in the paper.
R <- function(key, family, label, measures, computed, interp = "", final = "")
  data.frame(key = key, family = family, label = label,
             measures = measures, computed = computed,
             interp = interp, final = final, stringsAsFactors = FALSE)

REG <- do.call(rbind, list(
# ---------------- zone use / preference -------------------------------------
R("STEP5::zone_main_timepoint", "Zone use", "ALR(flow vs calm)",
  "Preference for flowing over still water. The primary endpoint.",
  "log(p_flow / p_calm) on AREA-NORMALISED occupancy proportions (flow/10, calm/42 of arena area). An additive log-ratio, which respects the constraint that time in one zone denies time to the other.",
  interp = "0 = flow used exactly in proportion to its area; >0 = flow over-used relative to calm. Higher means stronger preference for flowing water.",
  final = "Fig A-A"),
R("STEP5::zone_sec_high_timepoint", "Zone use", "ALR(high vs calm)",
  "Whether the preference is graded — use of the FASTEST sub-zone specifically.",
  "log(p_high / p_calm), area-normalised (high/10, calm/42). Sub-zone resolution of the same log-ratio.",
  interp = "Higher = the fastest band is over-used relative to calm. Distinguishes a graded preference from a generic one.",
  final = "Fig A-B"),
R("STEP5::zone_sec_medium_timepoint", "Zone use", "ALR(medium vs calm)",
  "Use of the intermediate-velocity sub-zone.",
  "log(p_medium / p_calm), area-normalised (medium/18, calm/42).",
  interp = "Higher = mid-velocity band over-used. A null here with a positive high-band result localises the preference at the top of the gradient.",
  final = "Fig A-C"),
R("STEP5::zone_sec_low_timepoint", "Zone use", "ALR(low vs calm)",
  "Use of the slowest flowing sub-zone.",
  "log(p_low / p_calm), area-normalised (low/13, calm/42).",
  interp = "Higher = slowest flowing band over-used. Expected near zero if fish select for velocity rather than for 'not-calm'.",
  final = "Fig A-D"),
R("STEP5::zone_main_cells_timepoint", "Zone use", "Main-zone cell means",
  "Occupancy of each main zone as a proportion, tested as a treatment x zone x interval grid.",
  "Beta GLMM on raw proportion of tracked time per zone; three-way design, so the reported interaction is treatment x interval averaged over zones.",
  interp = "Raw occupancy per zone. Reported as a cross-check on the log-ratio; not compositionally correct on its own.",
  final = ""),
R("STEP5::zone_sec_cells_timepoint", "Zone use", "Sub-zone cell means",
  "As above at sub-zone resolution (high/medium/low/calm).",
  "Beta GLMM on raw proportions across the four sub-zones.",
  interp = "As above at sub-zone resolution. Cross-check only.",
  final = ""),
R("STEP5::zone_sec_lr_pairs_timepoint", "Zone use", "Sub-zone log-ratio pairs",
  "All pairwise sub-zone log-ratios in one model.",
  "Long-format stack of the three ALR contrasts with an added zone_lr factor.",
  interp = "Tests whether the three sub-zone log-ratios differ from each other. Cross-check on Fig A panels B-D.",
  final = ""),
R("STEP5::flux_timepoint", "Zone use", "Zone flux per session",
  "How often the school crosses between zones — traffic across the boundary.",
  "Count of zone crossings per session; GLMM (Poisson/NB).",
  interp = "Higher = more traffic across the zone boundary, i.e. more fragmented use of the arena. Lower means the school settles.",
  final = "Fig B-A"),
R("STEP5::switches_timepoint", "Zone use", "Switches per session",
  "Count of main-zone changes.",
  "Count of flow<->calm transitions per session; GLMM.",
  interp = "Higher = more zone changes. Near-duplicate of zone flux; flux is the reported version.",
  final = ""),
R("STEP5::active_timepoint", "Zone use", "prop_active",
  "Fraction of time swimming above 1 body-length/s.",
  "Beta GLMM on the proportion of bins above threshold. MODEL DID NOT CONVERGE — no result exists.",
  interp = "No result exists (model did not converge). Report as non-convergence, not as absence of effect.",
  final = ""),
# ---------------- collective movement ---------------------------------------
R("STEP5::nnd_timepoint", "Collective", "Mean nearest-neighbour distance",
  "Local packing — how close each fish is to its closest neighbour.",
  "Per frame, mean over fish of the distance to the nearest other tracked fish; averaged over the interval. Identity-free (a function of the unlabelled point set).",
  interp = "Lower = fish packed more closely to their nearest neighbour, i.e. a tighter school.",
  final = "Fig C-A"),
R("STEP5::iid_timepoint", "Collective", "Mean inter-individual distance",
  "Global spread — mean separation of all pairs.",
  "Per frame, mean of all pairwise distances; averaged over the interval. Identity-free.",
  interp = "Lower = the group is less spread overall. With NND, separates local packing from global extent.",
  final = "Fig C-B"),
R("STEP5::hull_area_timepoint", "Collective", "School area",
  "Footprint — the area the group physically covers.",
  "Per frame, area of the convex hull of tracked positions; averaged over the interval. Identity-free.",
  interp = "Lower = the school covers less water. The most directly interpretable measure for holding-space arguments.",
  final = "Fig C-C"),
R("STEP5::centroid_speed_timepoint", "Collective", "School speed",
  "How fast the group as a whole moves, as distinct from where it is.",
  "Frame-to-frame displacement of the centroid of tracked positions / dt. Identity-free (the centroid does not depend on labels).",
  interp = "Higher = the group translates faster. A null here means any zone preference is not simply faster fish being carried into the current.",
  final = "Fig C-D"),
# ---------------- state sequence (binary alphabet) ---------------------------
R("SEQ::binary switch_rate", "Sequence (binary)", "Switch rate",
  "How often the school changes state per minute — fragmentation of behaviour.",
  "n_changes / obs_min, where n_changes is off-diagonal mass of the 1-s-bin transition-count matrix.",
  interp = "Higher = more state changes per minute, i.e. more fragmented behaviour. Collinear with entropy rate.",
  final = ""),
R("SEQ::binary entropy_rate", "Sequence (binary)", "Entropy rate",
  "Predictability of the next state given the current one.",
  "Normalised Shannon entropy rate of the fitted first-order Markov transition matrix (0 = perfectly predictable, 1 = maximally unpredictable).",
  interp = "0 = the next second is perfectly predictable from the current one; 1 = maximally unpredictable. Lower means the school has settled into one state.",
  final = "Fig B-B"),
R("SEQ::binary occ_Flow", "Sequence (binary)", "Occupancy, Flow",
  "Raw fraction of time in the flow state.",
  "Proportion of 1-s bins classified Flow.",
  interp = "Raw fraction of time in flow. Same quantity as Fig A-A but without the compositional correction.",
  final = ""),
R("SEQ::binary dwell_Flow", "Sequence (binary)", "Dwell time, Flow",
  "Mean length of an uninterrupted stay in flow.",
  "Mean run length (s) of consecutive Flow bins.",
  interp = "Higher = each visit to the flow lasts longer. Read with flow bouts per minute: few entries + long dwells = sustained engagement; many entries + short dwells = fragmented.",
  final = "Fig B-C"),
R("SEQ::binary dwell_Calm", "Sequence (binary)", "Dwell time, Calm",
  "Mean length of an uninterrupted stay in calm — whether rest is retained and how it is structured.",
  "Mean run length (s) of consecutive Calm bins. Undefined for a session that never enters Calm (n = 40/48).",
  interp = "Higher = each rest lasts longer. Tests whether rest is retained and how its structure changes.",
  final = "Fig B-D"),
R("SEQ::binary p_stay_Flow", "Sequence (binary)", "P(stay in Flow)",
  "Persistence: probability of remaining in flow from one second to the next.",
  "Self-transition probability tp[Flow,Flow] from the per-session transition matrix.",
  interp = "Probability of still being in flow one second later. Higher = stickier engagement. Collinear with dwell time, which is the reported version.",
  final = ""),
R("SEQ::binary p_stay_Calm", "Sequence (binary)", "P(stay in Calm)",
  "Persistence of the calm state.",
  "Self-transition probability tp[Calm,Calm].",
  interp = "Probability of still being in calm one second later. Higher = stickier rest.",
  final = ""),
R("SEQ::binary obs_min", "Sequence (binary)", "Observation minutes",
  "Tracked duration — a data-quality covariate, not an outcome.",
  "n_bins * bin width / 60.",
  interp = "Data-quality covariate. A treatment difference here would indicate unequal tracking, not behaviour.",
  final = ""),
# ---------------- state sequence (graded alphabet) ---------------------------
R("SEQ::graded switch_rate", "Sequence (graded)", "Switch rate (graded)",
  "State-change rate over the Low/Medium/High engagement alphabet.",
  "As binary switch rate but on the 3-level graded state.",
  interp = "Higher = more transitions among Low/Medium/High engagement bands.",
  final = ""),
R("SEQ::graded entropy_rate", "Sequence (graded)", "Entropy rate (graded)",
  "Predictability over the graded alphabet.",
  "Normalised entropy rate of the graded transition matrix.",
  interp = "As binary entropy rate, over three bands.",
  final = ""),
R("SEQ::graded occ_High", "Sequence (graded)", "Occupancy, High",
  "Time in the highest engagement band.",
  "Proportion of bins in the High state.",
  interp = "Higher = more time in the hardest-effort band.",
  final = ""),
R("SEQ::graded occ_Low", "Sequence (graded)", "Occupancy, Low",
  "Time in the lowest engagement band.",
  "Proportion of bins in the Low state.",
  interp = "Higher = more time in the easiest band.",
  final = ""),
R("SEQ::graded dwell_High", "Sequence (graded)", "Dwell time, High",
  "Mean uninterrupted stay in the highest band — sustained hard effort.",
  "Mean run length of consecutive High bins.",
  interp = "Higher = sustained hard effort lasts longer per visit.",
  final = ""),
R("SEQ::graded dwell_Low", "Sequence (graded)", "Dwell time, Low",
  "Mean uninterrupted stay in the lowest band.",
  "Mean run length of consecutive Low bins.",
  interp = "Higher = low-effort periods last longer per visit.",
  final = ""),
R("SEQ::graded p_stay_High", "Sequence (graded)", "P(stay in High)",
  "Persistence of the highest band.",
  "Self-transition probability tp[High,High].",
  interp = "Stickiness of the hardest band.",
  final = ""),
R("SEQ::graded p_stay_Low", "Sequence (graded)", "P(stay in Low)",
  "Persistence of the lowest band.",
  "Self-transition probability tp[Low,Low].",
  interp = "Stickiness of the easiest band.",
  final = ""),
R("SEQ::graded obs_min", "Sequence (graded)", "Observation minutes (graded)",
  "Tracked duration; covariate, not an outcome.",
  "As above.",
  interp = "Covariate, as above.",
  final = ""),
# ---------------- bout structure -- Tier 1 (defined for ~every session) ------
R("BOUT::commitment_index", "Bout (Tier 1)", "Commitment index",
  "Behavioural stickiness: how much of the interval is taken by ONE unbroken bout. NOTE it is state-agnostic — the longest bout of EITHER state.",
  "max(longest Flow bout, longest Calm bout) / interval duration. Because controls switch constantly, in practice it tracks flow commitment here, but the metric itself does not distinguish the two states.",
  interp = "Fraction of the interval taken by ONE unbroken bout, of EITHER state. Near 1 = the school essentially did one thing all interval. Does NOT say which state.",
  final = "Fig B-F"),
R("BOUT::max_flow_bout_s", "Bout (Tier 1)", "Longest flow bout",
  "The single longest uninterrupted engagement with the flow.",
  "max(duration_s) over Flow bouts, in seconds. Defined for all 48 sessions (unlike the bout-shape metrics, a maximum needs only one bout).",
  interp = "Higher = a longer single uninterrupted engagement with the flow. The cleanest 'long bouts' measure, and the natural partner to flow bouts per minute: few entries + a long maximum = sustained engagement.",
  final = "analysed, not figured"),
R("BOUT::max_calm_bout_s", "Bout (Tier 1)", "Longest calm bout",
  "The single longest uninterrupted stay in calm.",
  "max(duration_s) over Calm bouts.",
  interp = "Higher = a longer single uninterrupted rest.",
  final = ""),
R("BOUT::n_flow_ge30", "Bout (Tier 1)", "Sustained flow bouts (>=30 s)",
  "How many genuinely sustained engagements occur, ignoring brief incursions.",
  "Count of Flow bouts with duration >= 30 s.",
  interp = "Count of engagements lasting at least 30 s. Higher = more genuinely sustained visits, ignoring brief incursions.",
  final = ""),
R("BOUT::long_flow_frac", "Bout (Tier 1)", "Fraction of interval in long flow bouts",
  "Share of the whole interval spent inside long flow bouts.",
  "sum(duration of Flow bouts >= LONG_S) / interval duration.",
  interp = "Share of the whole interval spent inside long flow bouts. Combines how many and how long.",
  final = ""),
R("BOUT::flow_centroid", "Bout (Tier 1)", "Flow centroid (timing)",
  "WHEN within the interval flow is used: 0.5 = uniform, >0.5 back-loaded, <0.5 front-loaded. A pure timing measure, independent of amount.",
  "Time-weighted mean position of Flow bins within the interval, scaled to [0,1]. NA if fewer than 30 Flow bins.",
  interp = "0.5 = flow use spread evenly through the interval; >0.5 = back-loaded (builds up); <0.5 = front-loaded (fades). A pure timing measure, independent of amount.",
  final = ""),
R("BOUT::t_first_sustained_s", "Bout (Tier 1)", "Latency to first sustained bout",
  "How quickly the school commits to a sustained engagement.",
  "Start time of the first Flow bout lasting >= LATENCY_MIN_S; set to full interval length and flagged censored if none occurs.",
  interp = "Lower = commits to a sustained engagement sooner. Censored at interval length when no sustained bout occurs.",
  final = ""),
R("BOUT::engage_depth", "Bout (Tier 1)", "Engagement depth",
  "Average engagement band occupied, on a ranked scale.",
  "Time-weighted mean rank of the graded engagement state (1 = lowest .. K = highest).",
  interp = "Higher = occupies higher engagement bands on average.",
  final = ""),
R("BOUT::mean_intensity_sel", "Bout (Tier 1)", "Mean selected intensity",
  "Average flow intensity the school selects, on a ranked scale.",
  "Time-weighted mean rank of the intensity state.",
  interp = "Higher = selects faster water on average.",
  final = ""),
# ---------------- bout structure -- Tier 2 (CONDITIONAL, see caveat) ---------
R("BOUT::med_bout_flow", "Bout (Tier 2)", "Median flow bout",
  "Typical length of a flow engagement.",
  "Median of UNCENSORED Flow bout durations; NA below 5 uncensored bouts.",
  interp = "Typical (median) flow engagement length. Robust to the one long bout that drives the maximum.",
  final = ""),
R("BOUT::med_bout_calm", "Bout (Tier 2)", "Median calm bout",
  "Typical length of a calm stay.",
  "Median of uncensored Calm bout durations; NA below 5.",
  interp = "Typical calm stay length.",
  final = ""),
R("BOUT::burst_flow", "Bout (Tier 2)", "Burstiness, flow bouts",
  "Whether flow engagement is bursty (a few long bouts among many short) or regular.",
  "Kim & Jo (2016) finite-size-corrected burstiness of UNCENSORED Flow bout durations; NA below 5 uncensored bouts. UNDEFINED FOR THE WHOLE EXERCISE-CHOICE ARM AT INTERVAL 3 (mean 1.75 total bouts).",
  interp = "+1 = extremely bursty (a few long bouts among many short); 0 = random; -1 = regular. NOT REPORTED: undefined for the whole exercise-choice arm at interval 3.",
  final = ""),
R("BOUT::burst_calm", "Bout (Tier 2)", "Burstiness, calm bouts",
  "As above for calm.",
  "Same computation on Calm bouts.",
  interp = "As above for calm bouts.",
  final = ""),
R("BOUT::memory_flow", "Bout (Tier 2)", "Memory coefficient, flow",
  "REINFORCEMENT PROXY: does the length of one flow bout predict the length of the next?",
  "Lag-1 Pearson correlation of consecutive uncensored Flow bout durations within a session. Requires >= 5 uncensored bouts, so it is undefined for the exercise-choice arm exactly where commitment is strongest (defined 1/6/0 of 8).",
  interp = "Positive = a long bout tends to be followed by another long one (the reinforcement signature); negative = alternation. NOT REPORTED: defined for only 1/6/0 of 8 exercise-choice sessions.",
  final = ""),
R("BOUT::cv_bout_flow", "Bout (Tier 2)", "CV of flow bout duration",
  "Dispersion of flow bout lengths; a monotone transform of burstiness.",
  "(1 + B) / (1 - B) from the burstiness coefficient. Same definedness limits.",
  interp = "Higher = more variable flow bout lengths. Monotone transform of burstiness, same definedness limits.",
  final = ""),
R("BOUT::cv_bout_calm", "Bout (Tier 2)", "CV of calm bout duration",
  "Dispersion of calm bout lengths.",
  "As above on Calm.",
  interp = "Higher = more variable calm bout lengths.",
  final = ""),
R("BOUT::log2_med_ratio", "Bout (Tier 2)", "log2(median flow / median calm)",
  "Whether flow bouts are typically longer or shorter than calm bouts.",
  "log2(med_bout_flow / med_bout_calm); needs both medians defined.",
  interp = "0 = typical flow and calm bouts are equally long; >0 = flow bouts longer.",
  final = ""),
# ---------------- bout structure -- descriptive -----------------------------
R("BOUT::bout_rate_flow", "Bout (descriptive)", "Flow bouts per minute",
  "How frequently the school ENTERS the flow — fragmented vs sustained engagement.",
  "n_uncens_flow / obs_min. Defined for all 48 sessions, because a count of zero is defined where a dispersion of zero items is not.",
  interp = "Higher = enters the flow more often. Read WITH longest flow bout: few entries + long bouts = sustained engagement; many entries + short bouts = fragmented.",
  final = "Fig B-E"),
R("BOUT::mean_bout_flow", "Bout (descriptive)", "Mean flow bout",
  "Average flow engagement length.",
  "Mean of uncensored Flow bout durations.",
  interp = "Mean flow engagement length.",
  final = ""),
R("BOUT::mean_bout_calm", "Bout (descriptive)", "Mean calm bout",
  "Average calm stay length. Numerically identical to SEQ dwell_Calm.",
  "Mean of uncensored Calm bout durations.",
  interp = "Mean calm stay length. Numerically identical to SEQ dwell_Calm; report once.",
  final = ""),
R("BOUT::onset_hazard_flow", "Bout (descriptive)", "Flow onset hazard",
  "Instantaneous rate of leaving calm for flow.",
  "60 * (1 - P(stay in Calm)), recomputed locally. Numerically identical to SEQ p_stay_Calm.",
  interp = "Higher = leaves calm for flow more readily. Numerically identical to SEQ p_stay_Calm; report once.",
  final = ""),
R("BOUT::si_markov_max", "Bout (descriptive)", "Markov surrogate index",
  "How far the observed sequence departs from a first-order Markov surrogate.",
  "Maximum standardised deviation of observed vs Markov-shuffled statistics.",
  interp = "Higher = the sequence departs further from a first-order Markov process, i.e. has longer-range structure.",
  final = ""),
NULL
))

# ---- ANOVA harvest ----------------------------------------------------------
grab <- function(key) {
  parts <- strsplit(key, "::", fixed = TRUE)[[1]]
  eng <- parts[1]; id <- parts[2]
  if (eng == "STEP5") {
    f <- file.path(STEP5, id, "anova.csv")
    if (!file.exists(f)) return(NULL)
    a <- read.csv(f, stringsAsFactors = FALSE)
    tr <- a[a$term == "treatment", ]; it <- a[a$term == "treatment:timepoint_f", ]
    if (!nrow(tr)) return(NULL)
    list(tF = tr$chisq[1], td1 = tr$df[1], td2 = tr$df_denom[1], tp = tr$p_value[1],
         iF = if (nrow(it)) it$chisq[1] else NA, id1 = if (nrow(it)) it$df[1] else NA,
         id2 = if (nrow(it)) it$df_denom[1] else NA, ip = if (nrow(it)) it$p_value[1] else NA)
  } else {
    f <- if (eng == "SEQ") file.path(SEQ, "seq_anova_trial_x_timepoint.csv")
         else file.path(BOUT, "bout_anova_trial_x_timepoint.csv")
    if (!file.exists(f)) return(NULL)
    a <- read.csv(f, stringsAsFactors = FALSE)
    if (eng == "SEQ") {
      sp <- strsplit(id, " ", fixed = TRUE)[[1]]
      a <- a[a$alphabet == sp[1] & a$metric == sp[2], ]
    } else a <- a[a$metric == id, ]
    tr <- a[a$term == "treatment", ]; it <- a[a$term == "treatment:timepoint_f", ]
    if (!nrow(tr)) return(NULL)
    list(tF = tr$F[1], td1 = tr$df1[1], td2 = tr$df2[1], tp = tr$p[1],
         iF = if (nrow(it)) it$F[1] else NA, id1 = if (nrow(it)) it$df1[1] else NA,
         id2 = if (nrow(it)) it$df2[1] else NA, ip = if (nrow(it)) it$p[1] else NA)
  }
}

fp <- function(p) if (!is.finite(p)) "—" else if (p < 0.001) "< 0.001" else sprintf("%.3f", p)
st <- function(p) if (!is.finite(p)) "" else if (p < 0.001) " \\*\\*\\*" else if (p < 0.01) " \\*\\*" else if (p < 0.05) " \\*" else ""
ftxt <- function(F, d1, d2, p) {
  if (!is.finite(F)) return("—")
  # df as integers; fractional KR/Satterthwaite values stay in anova.csv.
  d2s <- if (is.finite(d2)) format(round(d2)) else "NA"
  sprintf("F(%s, %s) = %.2f, p = %s%s", format(round(d1)), d2s, F, fp(p), st(p))
}

L <- c()
add <- function(...) L <<- c(L, paste0(...))
add("# Catalogue of behavioural outcomes at the trial x interval level")
add("")
add("Every outcome analysed at the trial x interval level (N = 48 sessions; 16 trials x 3 intervals of 20 min, beginning 5, 45 and 85 min after onset), with what it measures, how it is computed, and its Type III ANOVA.")
add("")
add("Interval is fitted as a **three-level factor** in every model, so interval and treatment x interval terms are 2-df tests throughout. Where `df1 = 1` appears in an interaction it is **rank deficiency from empty treatment x interval cells**, not a different parameterisation.")
add("")
add("Significance: \\* p < 0.05, \\*\\* p < 0.01, \\*\\*\\* p < 0.001. \"Qualifies\" = treatment main effect **or** interaction significant.")
add("")

add("> **Reporting level.** The paper reports the **trial × interval analysis ",
    "only**. The aggregated / trial-level module is still produced by the ",
    "pipeline and retained as a cross-check, but it is not reported.")
add("")
add("> **Degrees of freedom** are shown as integers. Kenward-Roger and ",
    "Satterthwaite return fractional denominator df; the fractional values stay ",
    "in each engine's `anova.csv`, which is the analytic record.")
add("")
add("## Final indicator set")
add("")
add("Fifteen outcomes across three figures. Figure letters are PLACEHOLDERS — ",
    "manuscript numbering is assigned when the figures are placed.")
add("")
add("| Figure | Panels | Question it answers |")
add("| :--- | :--- | :--- |")
add("| **A** — zone preference | 2 × 2: ALR flow / high / medium / low, each vs calm | Do fish spend time in the flow rather than the calm, and is the preference graded across the velocity gradient? |")
add("| **B** — engagement pattern | 2 wide × 3 tall: flow↔calm crossings, entropy rate, dwell time Flow, dwell time Calm, flow bouts per minute, commitment index | Long uninterrupted bouts, or short bouts spaced by returns to calm? |")
add("| **C** — collective movement | 2 × 2: NND, IID, school area, school speed | Does exercise choice change how the group is organised in space? |")
add("")
add("**Longest flow bout** is analysed and reported in the contrast tables but ",
    "is not given a panel, so the final set is 15 outcomes in 14 panels.")
add("")
add("The **Final** column marks each outcome's place in the reported analysis:")
add("")
add("- a panel reference (e.g. **Fig B-E**) = reported and figured;")
add("- **analysed, not figured** = retained in the analysis and quotable in the ",
    "  text, with full contrasts in the companion Tukey document, but not given ",
    "  a panel;")
add("- **—** = not reported, either as a cross-check duplicate or excluded for ",
    "  the reason given in its interpretation.")
add("")
add("**Longest flow bout** is in the second category: it is kept in the analysis ",
    "and reported in the contrast tables, and was removed from Figure B only to ",
    "hold that figure at six panels.")
add("")

for (fam in unique(REG$family)) {
  add("## ", fam); add("")
  sub <- REG[REG$family == fam, ]
  add("| Outcome | Final | What it measures | How it is computed | How to read it | Treatment | Treatment × Interval | Qualifies |")
  add("| :--- | :---: | :--- | :--- | :--- | :--- | :--- | :---: |")
  for (i in seq_len(nrow(sub))) {
    g <- grab(sub$key[i])
    fin <- if (nzchar(sub$final[i])) paste0("**", sub$final[i], "**") else "—"
    if (is.null(g)) {
      add(sprintf("| **%s** | %s | %s | %s | %s | — | — | — |",
                  sub$label[i], fin, sub$measures[i], sub$computed[i], sub$interp[i]))
      next
    }
    q <- (is.finite(g$tp) && g$tp < 0.05) || (is.finite(g$ip) && g$ip < 0.05)
    add(sprintf("| **%s** | %s | %s | %s | %s | %s | %s | %s |",
                sub$label[i], fin, sub$measures[i], sub$computed[i], sub$interp[i],
                ftxt(g$tF, g$td1, g$td2, g$tp),
                ftxt(g$iF, g$id1, g$id2, g$ip),
                if (q) "**yes**" else "no"))
  }
  add("")
}

md <- file.path(OUT, "Behavioural_outcome_catalogue.md")
writeLines(L, md, useBytes = TRUE)
cat("markdown:", md, "\n")
dx <- file.path(OUT, "Behavioural_outcome_catalogue.docx")
if (file.exists(PANDOC)) {
  system2(PANDOC, c(shQuote(md), "-o", shQuote(dx)), stdout = TRUE, stderr = TRUE)
  cat("word:", dx, "exists:", file.exists(dx), "\n")
}
