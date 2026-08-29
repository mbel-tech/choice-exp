# Methods changes log — interval-level behavioural analysis

Running record of every analytical decision taken in the August 2026 revision,
written for direct transfer into Materials & Methods. Each entry gives **what
changed**, **why**, **what it supersedes**, and **draft M&M text** that can be
lifted with minimal editing.

Superseded statements are named explicitly so the M&M edit is a diff rather than
a rewrite, and so any sentence in the current Methods that is now false can be
found and corrected rather than left standing.

---

## 1. Analysis level: trial × interval only

**What changed.** Behavioural outcomes are reported at the trial × interval level
only (N = 48 sessions: 16 physical trials × 3 tracked intervals). The trial-level
(aggregated) module is no longer the basis of any main-text behavioural claim.

**Why.** The two levels answer different questions, and the interval level answers
the one the design was built for. Aggregating across intervals discards the
within-trial time course, which is where the distinctive finding lives — the
preference is not a constant offset but a pattern that changes across the trial.
Reporting both levels for every outcome also doubled the number of tests without
adding a construct.

**Trial key.** `phys_trial_id` (16 levels), not `trial_id` (48 levels, which is
session-level). This matters for the random effect: `(1 | phys_trial_id)` is what
makes the three intervals of a trial repeated measures on the same school.

**Supersedes.** Any Methods sentence presenting the aggregated analysis as the
primary behavioural module. The aggregated fits are still produced by the
pipeline and remain available for the Supplementary.

> **Draft M&M.** Behavioural outcomes were analysed at the level of the tracked
> interval (N = 48 sessions; 16 trials × 3 intervals of 20 min, beginning 5, 45
> and 85 min after trial onset). Trials, not sessions, were the unit of
> replication: every model included a random intercept for physical trial
> (`1 | phys_trial_id`, 16 levels), so the three intervals of a trial were
> treated as repeated measures on the same school.

---

## 2. Uniform interval parameterisation — the 2-df rule

**What changed.** Every interval model now fits interval as a **three-level
factor** (`timepoint_f`). The continuous (linear-trend, 1-df) alternative has
been removed from the AICc candidate set entirely, so the interval and
treatment × interval terms are **2-df tests in every outcome** and are directly
comparable across outcomes and figures.

**Implementation.** `.ALLOW_CONTINUOUS_TP <- FALSE` in
`scripts/01_pipeline_analysis/activity_analysis_STATS_choice_exp.R`; each of the
three model runners forces `cont_fixed_str <- NULL` at entry. Call sites still
pass `cont_fixed_str = "treatment * timepoint"` and were deliberately left
untouched, so the switch is a single override rather than an edit scattered
across ~15 call sites.

**Why this, and why not the reverse.** Five models previously fitted a continuous
interval (`nnd`, `iid`, `hull_area`, `zone_sec_medium`, `zone_sec_low`), which
made their interaction a 1-df treatment × linear-trend test rather than the 2-df
factor test used elsewhere — so interaction terms were not comparable across
outcomes. The obvious alternative fix, making interval continuous *everywhere*,
was rejected: the two parameterisations fail asymmetrically, and the AICc tables
quantify by how much.

| outcome | penalty paid by the losing parameterisation | shape |
|---|---|---|
| **ALR(flow)** — primary endpoint | continuous is **18.19 AICc worse** | strongly non-linear |
| ALR(high) | continuous is **15.13 worse** | strongly non-linear |
| zone flux | continuous is **15.63 worse** | strongly non-linear |
| polarisation | continuous is **8.77 worse** | non-linear |
| mean NND | factor is 4.95 worse | purely linear |
| mean IID | factor is 4.96 worse | purely linear |
| school area | factor is 4.79 worse | purely linear |
| ALR(medium) | factor is 5.27 worse | purely linear |
| ALR(low) | factor is 2.55 worse | purely linear |

Where the **factor** loses, it loses only its parameter penalty: `timepoint_f`
costs exactly two extra fixed parameters, worth ≈ 4–5 AICc units on their own,
and the observed gaps of 4.79–5.27 are almost entirely that penalty. The extra
parameters buy essentially no log-likelihood, so the deviation-from-linear
component is nil and **nothing is lost** by using the factor — the pattern is
still fully represented, merely tested with one spare degree of freedom.

Where the **continuous** form loses, it loses the signal itself. The factor model
wins ALR(flow) by 18.19 AICc *despite* paying that same penalty, implying a
non-linear component worth ≈ 23 units of 2·logLik. The treatment × interval
pattern is non-monotonic — a mid-trial dip and recovery — and a straight line
cannot represent it. Fitted to the commitment-index cell means (exercise choice
0.770 → 0.372 → 0.960), a linear contrast has a slope near zero and would report
"no change across the trial", which is false.

The factor is the general model; the linear model is a special case of it, and
the wrong special case for the headline outcomes. Forcing the factor costs power
on genuinely linear outcomes only, and §4 recovers that power exactly.

**Consequence, now confirmed by the refit (run `STEP5_stats_20260816_131610`).**
NND, IID and school area were predicted to fall from p ≈ 0.019–0.031 (1 df) to
p ≈ 0.06–0.09 (2 df) and thus to fail the omnibus test. They did:

| outcome | before (1 df) | after (2 df) |
|---|---|---|
| mean NND | F₁,₃₀ = 5.14, p = 0.031 | F₂,₂₈ = 2.47, **p = 0.103** |
| mean IID | F₁,₃₀ = 6.10, p = 0.019 | F₂,₂₈ = 3.21, **p = 0.056** |
| school area | F₁,₃₀ = 5.24, p = 0.029 | F₂,₂₈ = 2.66, **p = 0.088** |

Every model that already used the factor is numerically unchanged, confirming the
switch touched only the five continuous models. ALR(medium) also lost its
treatment main effect (F₁,₄₄ = 4.69, p = 0.036 → F₁,₁₄ = 1.19, p = 0.295), because
its denominator degrees of freedom corrected from 44 to 14.

If an outcome drops out under the 2-df rule, **it drops out**; reverting it to the
1-df fit because the 2-df result is unfavourable would be exactly the analytical
flexibility disclaimed in response to Reviewer 1.

**Incidental defect fixed.** `hull_area_timepoint` previously reported its ANOVA
from the continuous fit while writing its compact letter display from the factor
fit (`cld_treatmentxtimepoint_f.csv`), so the published figure's letters and its
quoted F-test came from two different models. One parameterisation makes them one
model.

**Supersedes.** (a) Any Methods statement that the fixed-effect parameterisation
of interval was selected per outcome by AICc. (b) The rationale block at the head
of `scripts/00_shared/rebuild_figure_caption_table_interval.R`, which documents
the old per-outcome policy and is now obsolete. (c) The comment in
`.mg_read_cld_tp()` about preferring the continuous CLD file — updated.

> **Draft M&M.** Observation interval was modelled as a three-level factor in
> every model, so that the interval and treatment × interval terms were tested on
> two degrees of freedom throughout and were directly comparable across outcomes.
> A continuous (linear-trend) parameterisation of interval was not used, because
> the treatment × interval pattern was strongly non-monotonic for the primary
> occupancy outcomes — for ALR(flow), the categorical parameterisation was
> favoured by ΔAICc = 18.2 over the linear one — so a linear trend would have
> mis-specified precisely the effect of interest.

---

### 2.3 What 1 df and 2 df actually test — for the Methods paragraph

Worth stating explicitly, because it is the first thing a statistical reviewer
will probe.

**1 df (interval numeric).** `treatment * timepoint` fits a single slope. The
interval term is one parameter; the interaction is the *difference in slopes*
between arms. It asks whether the outcome moves steadily across the trial, and
at a different rate in the two arms.

**2 df (interval a 3-level factor).** `treatment * timepoint_f` fits three free
cell means per arm. It asks whether the three interval means differ *in any way*,
and whether that pattern differs between arms.

The linear model is **nested inside** the factor model. With equally spaced
intervals the factor's 2 df partition orthogonally into linear + quadratic, so
the 2-df omnibus *is* the 1-df linear test plus a test of departure from
linearity. With only three timepoints the factor model is **saturated in time** —
three points determine three parameters — so there is no more flexible option
available and the choice is genuinely binary: constrain the shape to a straight
line, or do not constrain it.

Constraining buys power when the constraint is true and costs validity when it is
not. **In this dataset the two outcome families have different temporal shapes**,
which is why no single 1-df model could serve both:

| outcome | omnibus (2 df) | linear | quadratic | carried by |
|---|---|---|---|---|
| ALR(flow) — primary | p = 0.0023 | p = 0.051 | **p = 0.0024** | quadratic |
| ALR(high) | p = 0.0073 | p = 0.033 | **p = 0.0147** | quadratic |
| mean NND | p = 0.103 | **p = 0.035** | p = 0.835 | linear |
| mean IID | p = 0.056 | **p = 0.023** | p = 0.445 | linear |
| school area | p = 0.088 | **p = 0.033** | p = 0.590 | linear |
| school speed | p = 0.260 | p = 0.104 | p = 0.967 | linear |

Zone use is non-monotonic (dip at interval 2, recovery by interval 3); collective
structure is monotone drift. Forcing the linear form everywhere would have
reported the **primary endpoint as non-significant** (p = 0.051 rather than
0.002). Choosing the form per outcome is data-driven and is exactly the
flexibility disclaimed in response to Reviewer 1. Fitting the factor everywhere
costs nothing that cannot be recovered, because the linear contrast inside the
factor fit *is* the 1-df test.

**The honest counter-argument, for the record.** For an ordered factor with
equally spaced levels, a pre-specified linear trend test is a legitimate planned
contrast and is standard in time-course designs. Had linear-as-primary been
pre-specified, NND, IID and school area would all be significant (p = 0.023–0.035)
and ALR(flow) would not. That is a real trade. The omnibus was chosen because it
does not require committing in advance to a shape that turns out to be wrong for
the outcomes the study is about, and because sacrificing the primary endpoint to
rescue three secondary ones is the wrong direction. What matters most is that one
rule applies to all outcomes and is stated in advance.

**Rejected variant:** a hierarchical rule interpreting the linear contrast only
when the omnibus passes. It sounds principled but would gate out precisely the
collective outcomes whose linear component is the real signal, while adding
nothing for the quadratic-carried ones.

### 2.4 A second, larger consequence of the refit — Type III with a numeric covariate

The continuous parameterisation was not only costing comparability; it was making
the **treatment main effect uninterpretable** for the affected outcomes.

With `contr.sum` coding and a Type III ANOVA, the main effect of `treatment` is
evaluated where the interacting covariate equals zero. With interval as a factor
that is the average across intervals — the quantity of interest. With interval
numeric, `timepoint` ∈ {1, 2, 3}, so the treatment effect was being evaluated at
`timepoint = 0`: **extrapolated outside the observed range**. That inflates the
standard error and drives p towards 1.

The refit corrects this, and the change is large:

| outcome | treatment effect, continuous fit | treatment effect, factor fit |
|---|---|---|
| mean NND | F₁,₄₂ = 0.15, p = 0.699 | **F₁,₁₄ = 10.70, p = 0.006** |
| school area | F₁,₄₁ = 0.00, p = 0.996 | **F₁,₁₄ = 6.48, p = 0.023** |
| mean IID | F₁,₃₅ = 0.11, p = 0.740 | F₁,₁₄ = 3.19, p = 0.096 |

So the collective family did not disappear under the 2-df rule — it **changed
character**, from interactions that were partly an artefact of the
parameterisation to genuine treatment main effects. NND and school area now
qualify on the treatment main effect. This should be stated plainly in Methods
rather than left as an unexplained difference from any earlier draft.

---

## 3. Random-effect specification — unchanged

**What changed.** Nothing. Recorded here so the revision does not appear to
reopen a question already settled in the previous response to reviewers.

`(1 | phys_trial_id)` is fixed **by design**, not selected. AICc is retained only
to decide whether optional nuisance random effects (housing tank, trial date)
earn their place. Removing the continuous fixed-effect candidates (§2) does not
touch this: the candidate set still spans the same three random-effect
structures, now under a single fixed-effect formula.

> **Draft M&M.** The trial-level random intercept was fixed by design and was not
> subject to selection. Information criteria were used only to decide whether
> optional nuisance random effects (housing tank, trial date) were retained.

---

## 4. Orthogonal polynomial decomposition of the interval effect

**What changed.** New output, `poly_contrasts.csv` per model directory (STEP5),
plus `seq_poly_trial_x_timepoint.csv` and `bout_poly_trial_x_timepoint.csv` for
the sequence and bout engines. Each splits the interval effect into two 1-df
components.

**Why.** The three tracked intervals are **equally spaced** — midpoints 15, 55 and
95 min, 40 min apart — so the 2-df interval effect decomposes *orthogonally* into:

- **linear** — a monotone trend across the trial. Numerically what the withdrawn
  continuous parameterisation tested, so it recovers essentially the power that
  parameterisation had (for NND, ≈ F₁,₃₀ = 5.1 rather than the diluted
  F₂,₃₀ ≈ 2.6).
- **quadratic** — deviation from linear, i.e. the mid-trial dip.

This is what makes §2 costless: the uniform 2-df rule gives comparability, and
the decomposition gives back the sensitivity it spends. It also converts the
mid-trial dip from a description into a **formally tested claim**, which the
manuscript previously had no test for.

Equal spacing is what makes the `contr.poly` weights correct. If the tracking
schedule ever changes, the weights must be revisited.

**Reported exhaustively — this is the load-bearing detail.** The decomposition is
computed and written for **every** interval model, including those whose omnibus
is null, and a file is written even when nothing is significant. Emitting it only
where it happened to be significant would turn it into a second, unadjusted
significance filter. Note the deliberate contrast with the existing simple-effects
post-hoc blocks in the SEQ and BOUT engines, which *are* gated on a significant
omnibus — that gating is correct there, because simple effects are only
interpretable once an interaction is established.

**Multiplicity.** `adjust = "none"` within the decomposition is correct, not an
oversight: the two components are mutually orthogonal, pre-specified, and reported
in full rather than selected post hoc.

**Interpretive rule.** The **omnibus interaction alone** governs whether an
outcome qualifies for reporting. The decomposition characterises *shape*; it never
promotes a null omnibus to a positive result. An outcome whose omnibus fails is
reported as null, with any surviving linear component described as a suggestive
monotone trend.

**Known gap.** Several BOUT metrics have rank-deficient treatment × interval cells
(effective n = 30–33), which is why their omnibus interaction carries df1 = 1
rather than 2. For those the quadratic component is not estimable and is written
as `NA` rather than dropped, so the gap is visible rather than looking like an
untested metric.

> **Draft M&M.** Because the three intervals were equally spaced, the two-degree
> -of-freedom interval effect was additionally decomposed into orthogonal linear
> (monotone trend) and quadratic (mid-trial deviation) components. These contrasts
> were pre-specified and are reported for every outcome, significant or not
> (Supplementary Table SX). The omnibus interaction remained the primary test and
> alone determined which outcomes were reported in the main text; the
> decomposition was used only to characterise the shape of an effect, and no
> outcome was declared significant on the basis of a component contrast alone.

---

## 5. Outcome selection: construct coverage, not significance

**What changed.** The outcomes reported in the main text are a **pre-specified set
chosen on construct coverage**, each reported whatever it shows — not the subset
that reached significance.

**Why.** Selecting which outcomes to display by their p-values would reintroduce,
at the reporting stage, exactly the analytical flexibility disclaimed at the
modelling stage in response to Reviewer 1. A construct-based set applies the same
commitment consistently. The complete interval-level ANOVA for every outcome goes
to the Supplementary, so nothing is hidden and a reader can verify that the
reported set was not cherry-picked.

**The set (9 outcomes, 3 figures).**

| # | outcome | construct | why this one |
|---|---|---|---|
| 1 | ALR(flow) | preference for flowing water | Primary endpoint. The log-ratio respects the compositional constraint that time in one zone denies time to the other; every other zone-use measure is downstream of it. |
| 2 | ALR(high) | is the preference *graded*? | The arena offers a velocity gradient; whether fish concentrate in the fastest water is a claim a uniform-current design cannot make. |
| 3 | Dwell time, Calm | is rest retained? | The only outcome speaking to the calm side; underpins the claim that the calm zone was not redundant. Clears the collinearity screen. |
| 4 | Commitment to the flow | sustained or fragmented engagement? | Longest unbroken flow bout as a fraction of the interval. Collinear with #1, so reported as a time-domain re-expression. |
| 5 | Burstiness of flow bouts | is the *pattern* restructured? | Clears the screen against every ALR outcome — the one temporal-structure result that is not a re-description of #1. |
| 6 | Mean NND | local packing | Tightest local cohesion measure; the one the robustness checks were run on. |
| 7 | Mean IID | global spread | A school can tighten locally while spreading overall; #6 and #7 separate packing from extent. |
| 8 | School area | footprint | Space physically covered, independent of internal arrangement. Directly interpretable for holding-space arguments. |
| 9 | School speed | collective displacement | Separates *where* the group is from *how fast it moves*; guards against reading the preference as faster fish drifting into the current. |

**Excluded as redundant.** Switch rate, entropy rate and P(stay in flow) are
re-expressions of #1 and #4 (|r| ≈ 0.86–0.98). Longest flow bout is a rescaling
of #4. CV of flow-bout duration measures the same dispersion construct as #5.
Graded occupancy and dwell duplicate #2.

**Property worth stating.** Five of the nine qualify on significance anyway; the
four at risk are precisely those the 2-df rule affects. The two selection rules
agree on everything that is not parameterisation-sensitive.

> **Draft M&M.** Outcomes reported in the main text were specified in advance on
> the basis of construct coverage — one measure per question the study asks —
> rather than selected on their results, and each is reported irrespective of
> whether it reached significance. Interval-level results for all analysed
> outcomes are given in Supplementary Table SX.

---

## 6. The collective family is reported entire

**What changed.** All four identity-invariant collective metrics (NND, IID,
school area, school speed) are reported together, including nulls, as one figure.

**Why.** They measure related but non-identical aspects of group structure — local
packing, global spread, footprint, displacement — and reporting the family entire
removes any question of which member was chosen. Three of the four are expected to
fail the 2-df omnibus while retaining a monotone linear component; individually
those read as failed tests, but together they are one coherent bounded statement:
*the group tightens gradually and covers less water as the trial proceeds, at a
pace that does not separate the treatments cleanly, and without swimming faster.*
The school-speed null does specific work — it forecloses the reading that the flow
preference is merely faster-swimming fish being carried into the current.

Had the outcomes been selected on significance instead, the collective construct
would have disappeared from the main text entirely, with no statement either way.

**Figure captioning requirement.** Figure INT3's caption must be written around
the family and quote both the omnibus and the linear contrast per panel. A caption
reporting only the non-significant omnibus makes the panels look broken.

---

## 7. Polarisation excluded — identity invariance

**What changed.** Polarisation is excluded from the collective family and from
Figure INT3. This is a decision, not an omission, and it must be **stated
explicitly** rather than left silent.

**Why.** NND, IID, hull area and centroid speed are all computable from an
**unlabelled point cloud**: they are permutation-invariant functions of the set of
per-frame positions, so an identity swap between two fish cannot affect them.
Polarisation is not. It is an alignment order parameter over per-individual
headings, and a heading is the frame-to-frame displacement of *a given tracked
fish*. An identity swap corrupts both headings involved and hence the alignment
estimate. idtracker.ai cannot guarantee persistent identity in the choice arena —
the same limitation that removed fish-level random effects from every model — so
polarisation is not a trustworthy measurement here.

**Implementation (2026-08-16).** Removed from the pipeline entirely, not merely
hidden: the A9/A10 model fits, both BH family rows, the skinny-graph and
`gd_specs` plot entries, the sensitivity and equivalence lists, the two Word-report
sections, the descriptive column vectors, and all generated figure-caption text.
The analysis numbering A7–A16 is left unchanged so archived output stays
comparable. One reference is retained **deliberately**: the comment documenting
why Benjamini–Hochberg was withdrawn, because that reasoning depends on
polarisation having been a family member and deleting it would destroy the
justification for a decision the manuscript relies on.

Nothing was removed from `group_dynamics_STEP2b_choice_exp.R`, which still
computes the raw metric, nor from the workbook/reporting modules that reference it
historically — those are outside the analysis path and other code depends on them.

**A caveat worth recording.** Polarisation would have been the **strongest**
result in this family (treatment F₁,₁₄ = 12.26, p = 0.004; treatment × interval
F₂,₂₈ = 5.03, p = 0.014). Excluding the strongest result is defensible only if the
ground is stated; excluded silently it invites the suspicion that it was dropped
for the opposite reason. The exclusion no longer appears in any figure caption, so
**the Methods text must carry the identity-invariance argument explicitly** — it
is now the only place the reader can learn why an alignment measure is absent.

**Supersedes.** The brief exclusion note in the present Methods, and the figure
captions that previously stated the exclusion.

> **Draft M&M.** Group polarisation was measured but not analysed. Unlike
> nearest-neighbour distance, inter-individual distance, convex-hull area and
> centroid speed — each a permutation-invariant function of the per-frame set of
> positions, and therefore unaffected by identity errors — polarisation is
> computed from per-individual headings, i.e. the frame-to-frame displacement of a
> specific tracked fish. Because identity could not be maintained reliably across
> frames in the choice arena, polarisation was excluded on measurement-validity
> grounds. This exclusion was decided on those grounds alone: the treatment effect
> on polarisation was in fact the largest among the collective indicators.

---

## 8. Exclusions, failures and reduced sample sizes

State each of these rather than letting them appear as silent gaps.

- **`prop_active` cannot be reported at interval level.** The beta GLMM did not
  converge; `active_timepoint/` is an empty directory and row A5 of
  `bh_exploratory.csv` is `NA` throughout. Report the non-convergence rather than
  omitting the outcome.
- **Noradrenaline excluded throughout** — 24 analyte × region cells, not 28.
- **BOUT metrics duplicating SEQ metrics are reported once**, from SEQ:
  `mean_bout_calm` ≡ `dwell_Calm`, `mean_bout_flow` ≡ `dwell_Flow`,
  `onset_hazard_flow` ≡ `p_stay_Calm`.
- **Effective n below 48**, to be stated wherever it occurs: `dwell_Calm` n = 40
  (undefined for a session that never entered the Calm state); `burst_flow` n = 31
  (undefined without enough flow bouts to estimate a dispersion — also why its
  interaction carries df1 = 1).
- **Zone-assignment correction.** The earlier zone-assignment error and its
  correction are documented separately; the corrected values are the ones used
  throughout.

---

## 9. Multiple comparisons — unchanged, with one correction

**What changed.** The multiple-comparison structure is unchanged: dopaminergic
outcomes form a flat primary family across all areas, serotonergic outcomes are
exploratory. No adjustment is applied to the behavioural p-values (author
decision, 2026-08-07).

**The correction.** No BH-adjusted p-values may be quoted anywhere. `p_BH` and
`sig_BH` are `NA` in every `bh_*.csv` in the current output — the columns exist
but were never populated. Any text implying that BH-adjusted values were computed
or consulted is wrong and must be removed.

**Compact letter displays.** Where CLD letters appear, note that emmeans
substitutes Šidák for Tukey on the six-cell treatment × interval family (Tukey is
exact only for a single set of pairwise comparisons). Panels whose letters come
from a 6-cell family are therefore Šidák-adjusted; this must not be described as
Tukey.

---

## 10. Figures — three construct families

**What changed.** The interval-level results are presented as three figures, one
per construct family, replacing the previous single 4-panel interval figure.

| figure | panels | outcomes |
|---|---|---|
| INT1 — preference and rest | 3 (1 × 3) | ALR(flow), ALR(high), dwell time in Calm |
| INT2 — engagement structure | 2 (1 × 2) | commitment to the flow, burstiness of flow bouts |
| INT3 — collective movement | 4 (2 × 2) | NND, IID, school area, school speed |

**Implementation note worth keeping.** Every panel — including those sourced from
the SEQ and BOUT engines — is drawn by the single STEP5 panel builder
(`.mg_make_line_with_tukey()`). Using each engine's own builder makes the two
modules' point sizes, line widths and axis themes visible as a mismatch *inside a
single figure*. Routing every panel through one function removes the problem at
source rather than patching it with per-panel scale factors.

**Panels sourced from SEQ/BOUT carry no compact letter display**, because those
engines produce no CLD for their metrics. The captions say so; it is an omission
that is stated rather than silent.

**A panel must report the term that carries its claim.** Burstiness qualifies on
the **treatment main effect** (F₁,₂₂.₈ = 5.70, p = 0.026), not the interaction
(p = 0.176). Its panel reports the main effect first and the interaction beneath
it, so the null is visible but the panel does not appear to rest on it.

**Caption geometry.** Panels are sized at ~130 mm because that is the width at
which an 11 pt one-line ANOVA statement is confirmed to render without horizontal
clipping. At narrower widths the statement is silently cut mid-number — the canvas
is therefore sized from the caption, not from the plot.

---

## 11. Compact letter display — two defects corrected

Both were found while building the contrast export, and both had been changing
what the published figures asserted. The letter rule now lives in one file,
`scripts/00_shared/cld_letter_policy.R`, sourced by the figures and by the export,
so a letter in the document is by construction the letter in the panel. This is
verified, not assumed: a cross-check confirms the two agree for every outcome.

**Defect 1 — letters were truncated to two characters before relabelling.**
`.mg_read_cld_tp()` cut `.group` to `substr(..., 1, 2)`. With six cells,
three-letter groups are common, and the cut silently changed which cells appeared
to differ. IID at interval level is the worked example:

```
true groups     a | b | c | abc | abc | abc
after the cut   a | b | c | ab  | ab  | ab      <- "c" destroyed
```

The three `abc` cells genuinely do not differ from the `c` cell; after truncation
they no longer shared a letter with it. NND was affected identically
(`a | ab | abc | bc | bc | c`). The truncation also broke label merging, which
compares letter strings — and the `collapse_tp` workaround on the IID panel had
been written to hide that symptom. Letters are now kept in full and the
workaround is removed.

**Defect 2 — two panels displayed hand-typed letters.** `.fig9_A` (NND) and
`.fig9_C` (school area) passed a `label_override` column that replaced the
computed letters with a literal string, so those panels had stopped tracking their
own compact letter display. After the refit the NND override still happened to
agree with the model; the school-area one did not, and the disagreement was
substantive:

| cell pair, 85–105 min | override displayed | actual CLD | family-adjusted p |
|---|---|---|---|
| control vs exercise choice | `a` vs `b` → differ | `ab` and `a` → share a letter | **0.053** |

The override was showing the unadjusted simple-effect conclusion (p = 0.005)
beneath a caption stating that letters are the compact letter display and that
cells sharing a letter do not differ at p < 0.05. Both overrides are removed;
letters are read from `cld_treatmentxtimepoint_f.csv`.

**The rule itself.** Cells are walked control-first, then interval ascending, and
letters renamed in order of first appearance, so control at interval 1 always
begins with "a". This is a relabelling only — which cells share a letter is
untouched — and it exists so letters are comparable across panels.
`.cld_check_control_first()` warns if the invariant is ever violated.

**Two house rules registered 2026-08-17.**

1. *Control at interval 1 begins with "a"* (above). It cannot be forced to be "a"
   **alone**: the number of letters on a cell equals the number of maximal
   non-difference groups it belongs to, which is a property of the significance
   pattern rather than of the labelling. In school area, for instance, the only
   pair that differs is control/interval-2 vs exercise-choice/interval-3, so every
   other cell — control/interval-1 included — must carry both letters. Printing a
   single letter there would assert a difference the model does not support. Four
   panels are in this state and each emits a build-time note.
2. *No letters when nothing differs.* If every cell shares one group, the display
   is omitted rather than printed as six identical marks; six panels are affected.

> **Draft M&M.** Compact letter displays denote pairwise groupings across the six
> treatment × interval cells: cells sharing a letter did not differ at p < 0.05.
> Letters were relabelled so that the control group at the first interval always
> begins with "a", a relabelling that leaves the grouping itself unchanged. Panels
> without letters had no significant pairwise differences among any cells.

---

## Verification checklist

1. Every interval `anova.csv` reports `df = 2` for the interval and interaction
   terms; no `timepoint` (continuous) rows remain in any model directory.
2. All numeric statements in §2 and §5 re-derived from the regenerated
   `anova.csv` files — no value carried forward from a 1-df fit.
3. Decomposition checked against its parent test: NND's linear contrast should
   land near F₁,₃₀ = 5.1, the value the discarded continuous fit produced. If it
   does not, the contrast weights are wrong.
4. `poly_contrasts.csv` present for **every** interval model, including the nulls
   — a missing file is a silent selection.
5. Plotted cell means checked directly against source CSVs.
6. Every caption statistic matches its `anova.csv` / `poly_contrasts.csv` row.
7. `anova_statement_table_interval.csv` rebuilt from the new run before the
   figures are regenerated, or captions will silently show stale 1-df values.
8. Manuscript PDF re-exported and line numbers re-derived before any are cited in
   the response to reviewers — printed numbers differ from Word's internal count.

---

## 11. Second round (2026-08-17) — outcome set and the mean/maximum rationale

### 11.1 "Bout inventory" rewritten as "Bout-based outcomes" (as applied)

> Bout-based outcomes. Each 1 s bin was assigned to a state by majority vote:
> bins in which more than half of the zoned detections lay in the flow zone were
> labelled flow, and the remainder calm. A bout was then defined as a maximal
> run of consecutive bins in the same state. Because tracking is not continuous
> — each interval is a separate observation window, and bins carrying no valid
> zoned detection are discarded — bouts were never allowed to bridge a break in
> the sequence, and a bout lying at the start or end of a contiguous stretch of
> tracked bins is window-censored, so its true duration may exceed the observed
> value. The reported bout outcomes therefore draw on two inventories. The
> longest flow bout and the longest calm bout are taken over the full inventory,
> censored bouts included, since excluding the longest bout of an interval is
> precisely what would defeat them. Mean bout duration and the number of flow
> bouts per minute are computed from interior, uncensored bouts only, so that a
> bout of unknown true length neither biases the average nor is counted as
> though it had been observed to completion. Mean bout duration obtained this
> way is identical by construction to the dwell time reported from the
> state-sequence analysis, and the two agree exactly.

**Three corrections against the earlier draft, each traced to the script**
(`bout_structure_analysis_choice_exp.R`):

| Draft said | Script does | Line |
|---|---|---|
| full inventory = "maximum bout duration, count of sustained bouts, long-bout time fraction" | only the maximum is a reported outcome; `flow_bouts` carries no censoring filter | 390, 392 |
| uncensored = "mean and median duration, lag-1 duration correlation" | only the mean is reported | 416, 428 |
| — (omitted) | **flow bouts per minute uses the uncensored inventory** and is reported; it is a COUNT, so "location or shape of the duration distribution" does not describe it | 427 |
| "a bout touching either edge of its interval" | censoring is per contiguous SEGMENT: breaks come from discarded bins as well as from the gap between intervals | 288, 326 |

The majority-vote threshold (`prop_flow > 0.5`, line 235) was not stated anywhere
in the manuscript and is now given once, here. The 1 s binning and the flow
fraction are already defined in the preceding "Collective state sequence"
paragraph and are not restated.

Metrics dropped from the paragraph because they are not reported: count of
sustained bouts, long-bout time fraction, median duration, lag-1 duration
correlation, and the graded engagement / self-selected intensity alphabets. All
five terms occurred **only** in this paragraph, so nothing else in the
manuscript depends on them.

### 11.2 Mean-vs-maximum rationale (new paragraph, as applied)

> For each of the two states, three complementary statistics of the same
> bout-duration distribution are reported: how often the state is entered (flow
> bouts per minute), how long a typical completed visit lasts (dwell time, i.e.
> mean bout duration), and how long the single most sustained visit lasts
> (longest bout). The mean and the maximum are not interchangeable, and they do
> not see the same behaviour. Because the mean and the per-minute count are
> computed on interior, uncensored bouts only whereas the maximum is taken over
> the full inventory, a visit still running when the interval ends is invisible
> to the mean by construction and detectable only by the maximum. Where the two
> agree, the whole duration distribution has shifted and the typical unit of
> behaviour has changed. Where they disagree, the change is confined to part of
> that distribution: a longer maximum at an unchanged mean indicates that
> sustained episodes have been added to an unchanged population of brief visits,
> whereas a changed mean at an unchanged maximum indicates that time has been
> redistributed among visits without any change in how long the state is ever
> held continuously. Reporting the number and the duration of bouts together,
> and reading their joint pattern rather than either alone, follows established
> practice in the analysis of behavioural state architecture in fish (Yokogawa
> et al., 2007).

### 11.3 Also removed from Methods

Burstiness, from the shape-metric list, and its full specification (the
finite-size-corrected Kim & Jo 2016 estimator). Results already stated that
burstiness is not reported, so the two sections contradicted each other.
**Kim & Jo was never in the reference list**, so this also closes a dangling
citation.

### 11.4 Reference added / removed

Added — Yokogawa, T., Marin, W., Faraco, J., Pézeron, G., Appelbaum, L., Zhang,
J., Rosa, F., Mourrain, P., Mignot, E., 2007. Characterization of sleep in
zebrafish and insomnia in hypocretin receptor mutants. PLoS Biol. 5, e277.
https://doi.org/10.1371/journal.pbio.0050277

Removed — Trepka et al. (2021), cited only by the normalised entropy rate, which
leaves the paper.
