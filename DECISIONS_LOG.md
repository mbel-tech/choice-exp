# Decisions log — behavioural analysis revision

Running record of **what changed, when, and why**, across four categories:
modelling choices, indicator selection, outcome tables, and results. Companion to
`METHODS_CHANGES.md`, which holds the drafted Materials & Methods text; this file
is the chronology and the rationale.

Entries are newest-last. Every claim here is traceable to a script or an output
file, named inline.

---

## Standing decisions

| # | Decision | Status |
|---|---|---|
| D1 | Behavioural outcomes reported at the **trial × interval level only** (N = 48 sessions). The aggregated / trial-level module is still produced and retained as a cross-check, but is not reported. | fixed |
| D2 | Interval is a **three-level factor** in every model; the continuous (1-df) alternative is withdrawn from the candidate set. | fixed |
| D3 | Random effect `(1 \| phys_trial_id)` is **fixed by design**; information criteria decide only optional nuisance terms. | unchanged from previous revision |
| D4 | Orthogonal **linear + quadratic** decomposition reported for every interval outcome, exhaustively. The omnibus remains the primary test. | fixed |
| D5 | Indicators selected on **construct coverage**, not on significance. | fixed |
| D6 | **Polarisation excluded** on identity-invariance grounds and removed from the pipeline. | fixed |
| D7 | Denominator **degrees of freedom printed as integers**; fractional values retained in `anova.csv`. | fixed |
| D8 | Compact letter display on **every panel of every figure**, under one shared letter policy. | fixed |
| D9 | **Control at interval 1 begins with "a"** in every panel (relabelling only). | fixed |
| D10 | **No letters at all when every cell shares one group** — nothing differs, so letters are omitted rather than printed as six identical marks. | fixed |

---

## 1. Modelling choices

### 1.1 Uniform 2-df interval parameterisation
**Changed:** every interval model now fits `treatment * timepoint_f`. The
continuous candidates (`*_cont`) were removed from the AICc candidate set via
`.ALLOW_CONTINUOUS_TP <- FALSE` in `activity_analysis_STATS_choice_exp.R`.

**Why:** five models previously fitted a continuous interval, making their
interaction a 1-df trend test not comparable with the 2-df factor test used
elsewhere. The reverse fix (continuous everywhere) was rejected because the two
parameterisations fail asymmetrically — where the factor loses on AICc it loses
only its parameter penalty (~4–5 units, observed gaps 4.79–5.27), whereas the
continuous form loses the signal itself (ΔAICc = 18.19 on the primary endpoint).

**Consequence:** confirmed by refit — see §4.1.

### 1.2 Type III main effects were previously untestable for five outcomes
**Found:** with a numeric `timepoint`, the Type III treatment main effect is
evaluated at `timepoint = 0`, i.e. extrapolated outside the observed range
{1, 2, 3}. The affected outcomes' treatment tests were therefore meaningless
(p = 0.699, 0.740, 0.996). The factor parameterisation evaluates the main effect
averaged across intervals, which is the quantity of interest.

### 1.3 Orthogonal polynomial decomposition added
**Added:** `.write_poly_contrasts()` (STEP5) plus equivalent blocks in the SEQ and
BOUT engines, writing `poly_contrasts.csv` / `*_poly_trial_x_timepoint.csv`.

**Why:** recovers the power the 2-df rule costs on genuinely linear outcomes (the
linear contrast is numerically the withdrawn continuous test), and turns the
mid-trial dip into a formally tested claim. Reported **exhaustively** — a missing
file would itself be a silent selection.

### 1.4 Degrees of freedom displayed as integers
**Changed:** caption builders in `rebuild_figure_caption_table_interval.R`,
`.mg_seqbout_cap()`, `.mg_step5_cap()` and both export documents.
**Why:** Kenward-Roger/Satterthwaite return fractional df; "F(2, 24.2)" invites
the reader to wonder what a fifth of a degree of freedom is. Display only — the
fractional value stays in `anova.csv`.

---

## 2. Indicator selection

### 2.1 Final set — 15 outcomes across three figures

Figure letters are **placeholders**; manuscript numbering is assigned when the
figures are placed.

**Figure A — zone preference (2 × 2).** Does the school prefer flowing water, and
is the preference graded across the velocity gradient?

| Panel | Outcome | How to read it |
|---|---|---|
| A | ALR(flow vs calm) | 0 = flow used exactly in proportion to its area; >0 = over-used relative to calm. Higher = stronger preference. |
| B | ALR(high vs calm) | Higher = the fastest band is over-used. Distinguishes a graded preference from a generic one. |
| C | ALR(medium vs calm) | Higher = mid-velocity band over-used. A null here alongside a positive high-band result localises the preference at the top of the gradient. |
| D | ALR(low vs calm) | Higher = slowest flowing band over-used. Expected near zero if fish select for velocity rather than for "not-calm". |

**Figure B — engagement pattern (2 wide × 3 tall).** Long uninterrupted bouts, or
short bouts spaced by returns to calm?

| Panel | Outcome | How to read it |
|---|---|---|
| A | Flow↔calm crossings | Higher = more traffic across the boundary, i.e. more fragmented use of the arena. Lower = the school settles. |
| B | Entropy rate | 0 = the next second is perfectly predictable from the current one; 1 = maximally unpredictable. Lower = settled into one state. |
| C | Dwell time, Flow | Higher = each visit to the flow lasts longer. |
| D | Dwell time, Calm | Higher = each rest lasts longer. Defined for 40 of 48 sessions. |
| E | Flow bouts per minute | Higher = enters the flow more often. Read WITH C and with longest flow bout: few entries + long dwells = sustained engagement; many entries + short dwells = fragmented. |
| F | Commitment index | Fraction of the interval taken by ONE unbroken bout, of **either** state. Near 1 = the school essentially did one thing all interval. Does not say which state. |
| — | **Longest flow bout** | **Analysed, not figured.** Higher = a longer single uninterrupted engagement. The natural partner to panel E. Removed from the figure only to hold it at six panels; kept in the analysis and in the contrast tables. |

**Figure C — collective movement (2 × 2).** Does exercise choice change how the
group is organised in space?

| Panel | Outcome | How to read it |
|---|---|---|
| A | Mean NND | Lower = fish packed more closely to their nearest neighbour, i.e. a tighter school. |
| B | Mean IID | Lower = the group is less spread overall. With NND, separates local packing from global extent. |
| C | School area | Lower = the school covers less water. Most directly interpretable for holding-space arguments. |
| D | School speed | Higher = the group translates faster. A null here means any zone preference is not simply faster fish being carried into the current. |

### 2.2 Indicators removed, and why

| Indicator | Removed because |
|---|---|
| **Polarisation** | Requires persistent per-fish identity (heading = a given fish's frame-to-frame displacement); every other collective indicator is a permutation-invariant function of the unlabelled point set. Removed from the pipeline entirely, not merely hidden. Would have been the strongest collective result — so the Methods text must state the ground explicitly. |
| **Burstiness, CV of bout duration, memory coefficient** | Undefined for the whole exercise-choice arm at interval 3. Not a tunable artefact: censoring costs only ~1.3–1.8 bouts/session, and the arm averages **1.75 total flow bouts** at interval 3, so even with no censoring and a threshold of 2 only 3 of 8 sessions would qualify. Their treatment effects rest on the minority of sessions that kept switching — selection on the outcome, a caveat the bout engine's own header flags. |
| ~~**Longest flow bout**~~ | **NOT removed — kept in the analysis.** Dropped from Figure B (2026-08-17) only to hold that figure at six panels. It remains in the outcome catalogue, in the Tukey contrast tables, and is quotable in the text, marked "analysed, not figured". |
| **Switches per session** | Near-duplicate of zone flux; flux is the reported version. |
| **`prop_active`** | Beta GLMM did not converge. Must be reported as non-convergence, not as absence of effect. |

### 2.3 Combined two-state dwell panel — built and withdrawn
A single panel showing Flow and Calm dwell together (colour = state) was built on
2026-08-17 and removed the same day: it forced colour to mean "state" in that
panel while meaning "treatment" in every other panel of the same figure. Dwell
time is now two panels on the standard treatment palette.

---

## 3. Outcome tables and documents

| Document | Contents | Generator |
|---|---|---|
| `Behavioural_outcome_catalogue.{md,docx}` | **Every** interval-level outcome: what it measures, how it is computed, **how to read it**, both ANOVA terms, whether it qualifies, and which final panel it occupies | `export_outcome_catalogue.R` |
| `Interval_treatment_Tukey_contrasts.{md,docx}` | For the final set: omnibus, all 15 pairwise cell contrasts, both simple-effect directions, and the CLD | `export_tukey_interval_module.R` |
| `METHODS_CHANGES.md` | Drafted Materials & Methods text, per decision | maintained by hand |
| `DECISIONS_LOG.md` | This file | maintained by hand |

**Interpretation column added** (2026-08-17) so each table is usable without
returning to the code — it states which direction of the number means what.

**Both documents are Markdown-first**; the `.docx` is produced from the `.md` via
Pandoc, so the two cannot drift apart.

**Guard added:** both generators now refuse to read an incomplete pipeline run,
requiring a late-stage sentinel file. Without it the catalogue silently read a
run that was still being written and filled the collective rows with em-dashes.

---

## 4. Changes in results

### 4.1 Collective outcomes — interactions lost, main effects gained
The 2-df refit moved these in **both** directions. Predicted (F ≈ 2.6 / 3.0 / 2.6)
and observed:

| outcome | interaction before (1 df) | interaction after (2 df) | treatment before | treatment after |
|---|---|---|---|---|
| mean NND | F₁,₃₀ = 5.14, p = 0.031 | F₂,₂₈ = 2.47, p = 0.103 | p = 0.699 | **p = 0.006** |
| mean IID | F₁,₃₀ = 6.10, p = 0.019 | F₂,₂₈ = 3.21, p = 0.056 | p = 0.740 | p = 0.096 |
| school area | F₁,₃₀ = 5.24, p = 0.029 | F₂,₂₈ = 2.66, p = 0.088 | p = 0.996 | **p = 0.023** |

The collective family did not disappear; it **changed character**, from partly
parameterisation-dependent interactions to genuine treatment main effects.

### 4.2 ALR(medium) lost its treatment effect
F₁,₄₄ = 4.69, p = 0.036 → F₁,₁₄ = 1.19, p = 0.295, because its denominator df
corrected from 44 to 14. It no longer qualifies, and is retained in Figure A only
to show where on the gradient the preference sits.

### 4.3 Every already-factor model is numerically unchanged
Confirming the switch touched only the five continuous models.

### 4.4 The primary interaction is quadratic, not linear
ALR(flow): treatment × linear **p = 0.051**, treatment × quadratic **p = 0.0024**.
A linear-only parameterisation would have reported the paper's headline
interaction as null. Within arms, control is flat (linear p = 0.51, quadratic
p = 0.44) while exercise choice is strongly non-monotonic (p = 0.0014 and
p = 6.9 × 10⁻⁶).

### 4.5 Reinforcement: no first-vs-later correlation, but a within-session trend
Computed directly from `bout_inventory.csv`, because the intended metric
(`memory_flow`) is defined for only 1/6/0 of 8 exercise-choice sessions.

- **First vs later bout duration:** r = −0.04 (control), −0.05 (exercise choice),
  both p > 0.84. Largely untestable — 30 of 38 sessions have a left-censored
  first bout, and restricting to uncensored leaves n = 6 and n = 2.
- **Within-session trend** (log duration ~ treatment × bout index, 668 uncensored
  bouts, 38 sessions): interaction **F₁,₃₃₉ = 4.71, p = 0.031**. Control slope
  −0.018 log-units per bout (CI −0.032 to −0.004) — bouts shorten through a
  session. Exercise choice +0.033 (CI −0.012 to +0.078) — no shortening.

This is a new analysis, not an existing pipeline metric; it would need adding to
the engine to appear in the paper.

---

## 5. Defects found and corrected

| Defect | Effect | Fix |
|---|---|---|
| CLD letters truncated to 2 characters **before** relabelling | Destroyed the third letter of any 3-group cell, so cells that do **not** differ stopped sharing a letter. IID's true groups are `a\|b\|c\|abc\|abc\|abc`; the cut turned the three `abc` cells into `ab`. | Truncation removed; the `collapse_tp` workaround written to hide the symptom removed with it |
| Two panels displayed **hand-typed** letters (`label_override`) | Those panels had stopped tracking their own CLD. School area displayed `a` vs `b` at 85–105 min — a difference — where the CLD gives `ab` and `a`, which share a letter (family-adjusted p = **0.053**). It was showing the unadjusted simple-effect result under a caption stating the letters are the CLD. | Both overrides removed; letters read from `cld_treatmentxtimepoint_f.csv` |
| `commitment_index` mislabelled "commitment to the flow" | The metric is `max(longest Flow, longest Calm) / interval` — **state-agnostic**. | Relabelled "commitment index (longest bout, either state)" |
| Caption strings silently clipped mid-number | `.mg_wrap_caption` is a no-op, so long statements were cut ("p = 0.0"). | Panels sized from the caption (~130–140 mm each) |
| Figure modules auto-detected an **in-progress** pipeline run | Panels and tables rendered from half-written output. | Sentinel-file guard in the figures module and both exporters |
| Polarisation removal left orphaned references | `.p2_pol` / `.p3_pol` used but no longer defined → run failed after ~50 min of fitting | 2×2 diagnostic grids rebuilt as 3-panel rows; static orphan check added |

---

## 5b. Compact-letter-display house rules (registered 2026-08-17)

Both rules live in `scripts/00_shared/cld_letter_policy.R` and are applied by the
figures, the outcome catalogue and the Tukey export from that single definition.

### Rule 1 — control at interval 1 begins with "a"
Cells are walked control-first, then interval ascending, and letters renamed in
order of first appearance. Because control/interval-1 is first in that order, its
letters are the first seen, so its group always begins with "a". This is a
**relabelling only**: which cells share a letter is untouched.

**Limit of the rule, and why it is not forced further.** The rule cannot make that
group be "a" *alone*, and no relabelling can:

> the number of letters on a cell = the number of maximal non-difference groups it
> belongs to, which is a property of the significance pattern, not of the labels.

Worked example, **school area**: the only pair that genuinely differs is
control/interval-2 vs exercise-choice/interval-3. With exactly one differing pair
the display is forced — those two cells take distinct single letters, and every
other cell, control/interval-1 included, must carry **both**, because it differs
from neither. Printing "a" there instead of "ab" would assert that
control/interval-1 differs from control/interval-2, which the model does not
support.

Four panels are in this state (control/interval-1 = "ab"): **school area, flow
bouts per minute, commitment index, longest flow bout**. Each emits a build-time
note rather than being silently rewritten.

### Rule 2 — no letters when nothing differs
If every cell carries the same single group, no pair differs and the letters
convey nothing; six identical "a"s make a reader hunt for a contrast that is not
there. Those panels are drawn without a compact letter display, and the caption
states that a panel without letters has no significant pairwise differences.

Six panels are affected: **ALR(medium), ALR(low), dwell time Flow, dwell time
Calm, mean IID, school speed**.

### Where the convention is stated to the reader
One shared caption paragraph (`.CLD_CAPTION_NOTE`) is spliced into all three
figure captions, so they cannot describe the letters differently. It states the
sharing convention, the relabelling, the meaning of a two-letter cell, and the
no-letters case.

---

## 6. Open items

- **Verification run** of the polarisation-removed pipeline has not yet completed
  cleanly end to end. The statistics are unaffected (all reported numbers come
  from `STEP5_stats_20260816_131610`, which completed before the removal); what
  needs confirming is only that the script runs to completion.
- **`(Table S__)` placeholders** in the reviewer response, row 14.
- **Manuscript defects** at LN 191 ("Section XX" placeholder, double period) and
  LN 930 ("A-B min, C-D min, E-F min") — flagged, not fixed, because fixing
  reflows all line numbers.
- **Treatment × zone CLD reader** has the same truncation defect as the interval
  one (`figures_step5_module.R`, `.mg_read_cld_treatxzone`). Left alone because it
  feeds already-approved Figure 8 panels C/D, where longer letters could overflow
  the layout. Should be checked before those figures are final.
- **Within-session bout-duration slope** (§4.5) is not a pipeline metric.

---

## 7. Second figure/manuscript round (2026-08-17)

### 7a. Outcome set

| Outcome | Change | Reason |
|---|---|---|
| Normalised entropy rate | **removed** from Table 1, Results, figures | leaves the paper; retained in Table S10 |
| Commitment index | **removed** from Table 1, Results, Discussion, figures | state-agnostic by construction (longest bout of EITHER state), so it could not support a claim about the flow specifically; retained in Table S10 |
| Longest calm bout (s) | **added** to Table 1, Results, Figure B panel D | calm-state counterpart to the longest flow bout; its null result is what makes the flow effect flow-specific rather than a general lengthening of bouts |
| Dwell time, flow / calm | **retained** in Table 1 and Results, no longer plotted | they are the MEAN of the same distribution whose MAXIMUM is plotted |

### 7b. D11 — mean and maximum are both reported

Established from the code rather than asserted:

* Dwell time **is** mean bout duration. `bout_structure_analysis_choice_exp.R:23`
  asserts `mean_bout_flow == dwell_Flow` to 1e-8 and warns if it breaks
  (lines 512-524).
* The mean and the per-minute count use **interior, uncensored bouts only**
  (`behavioural_sequence_analysis_choice_exp.R:188-190`;
  `bout_structure_analysis_choice_exp.R:427-429`); the maximum uses the **full
  inventory**. A bout still running at the interval edge is therefore invisible
  to the mean *by construction* and visible only to the maximum.

Reading of accordance vs discordance, now stated in Methods and Discussion:

| | maximum changes | maximum unchanged |
|---|---|---|
| **mean changes** | whole distribution shifts; the unit of behaviour changed | redistribution within an unmoved ceiling |
| **mean unchanged** | tail-only: sustained episodes added to unchanged brief traffic (**consolidation**) | bout structure unchanged; occupancy change carried by bout NUMBER |

Observed, and dissociating in **opposite directions for the two states**:

| | mean (dwell) | maximum (longest bout) | reading |
|---|---|---|---|
| Flow | F(1,13) = 0.59, p = 0.458 | F(1,13) = 86.68, p < 0.001 | tail-only -> consolidation |
| Calm | T x I F(2,24) = 4.01, p = 0.031 | T x I F(2,28) = 2.20, p = 0.129 | centre-only -> restructuring |

Either statistic alone returns a null on one of the two states, which is the
justification for reporting both. Citation added: **Yokogawa et al. (2007)**,
*PLoS Biology* 5(10) e277 — zebrafish, reports bout number and bout duration as
separate architecture measures and reads their joint pattern as fragmentation;
our flow result is its mirror image. This is the only reference added.

### 7c. Figures

* **Figure B rebuilt 2x3 -> 2x2** (`figures_step5_module.R`): A crossings,
  B flow bouts per minute, C longest flow bout, D longest calm bout. Saved at
  280 x 240 mm, matching Figures A and C.
* `build_seqbout_cld.R` gained `max_calm_bout_s`, without which panel D would
  have been the only one of four without letters.
* `export_tukey_interval_module.R` Figure B block remapped; entropy, both dwell
  times and commitment moved to `panel = "not figured"` (they stay quotable).
  A `stopifnot(!anyDuplicated(...))` guard was added on the `no` field after a
  renumbering collision.
* **CLD label positions** (all three requested, all verified in the render):
  Figure A panel A control/interval 2 shifted left (x_adj = -0.22); Figure A
  panel B interval 1 split with `force_separate` on both rows before control
  could be placed left of its mean; Figure C panel A the **exercise-choice**
  "ab" moved, not control's "a" — geometry showed control's letter was already
  clear and "ab" was the one sitting on the control marker.

### 7d. Corrections made in passing

* Methods specified burstiness in full (Kim & Jo 2016) while Results said it was
  not reported. Specification removed. **Kim & Jo was never in the reference
  list**, so this also closes a dangling citation.
* Discussion said "the single longest unbroken bout of flow use (i.e.
  commitment)". Commitment was the longest bout of *either* state, so the
  sentence misdescribed its own measure. Restated as the longest **flow** bout,
  which keeps the claim and its significance (T x I F(2,28) = 9.78, p < 0.001)
  and makes it literally true.
* Table 1 row 2 said "per interval" in Explanation and "per session" in
  Computation. Harmonised; the Explanation now explains rather than restating
  the computation.
* Table 1 typo "AThe single longest..." fixed.
* Subject/verb agreement in the Discussion rest-structure sentence.

### 7e. Tooling

`redline.py` gained `to_del_text()`: inside `<w:del>`, `<w:instrText>` must
become `<w:delInstrText>`, not only `<w:t>` -> `<w:delText>`. Field
instructions carry no visible text, so this was invisible until `validate.py`
rejected the deleted Trepka reference, which is a Zotero field.

### 7f. Open item raised, NOT actioned

**Section 3.5 still reports a "switch-count main effect and interaction" as
non-significant.** The plan was to rename it to the flow<->calm crossings
outcome, but that is **wrong**: crossings are highly significant
(F(1,14) = 45.00, p < 0.001), as is interval-level `switch_rate`
(F(1,14) = 62.05, p < 0.001, Table S10). Renaming would have created a direct
contradiction with 3.1.2. What §3.5 actually tested cannot be recovered from the
current outputs, so the paragraph was left untouched and needs an author
decision.

* The Methods bout paragraph ("Bout inventory") listed metrics that are not
  reported (count of sustained bouts, long-bout time fraction, median duration,
  lag-1 duration correlation) and **omitted flow bouts per minute**, which is
  reported and which uses the uncensored inventory. Rewritten as "Bout-based
  outcomes", stating the split by outcome rather than by metric class, adding
  the majority-vote threshold (`prop_flow > 0.5`), and correcting censoring from
  "either edge of its interval" to the edge of any contiguous segment. See
  METHODS_CHANGES.md §11.1 for the line-by-line trace.
