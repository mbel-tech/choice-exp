# HANDOFF

How this project treats trajectory, monoamine, and stress-marker (cortisol)
data — written so someone running a **new** experiment with the **same data
types** can reproduce the analytical treatment without reproducing the design.
Experimental design varies; data type does not. Read this document front to
back once; use [`CONVENTIONS.md`](CONVENTIONS.md) and
[`RESULTS_STYLE.md`](RESULTS_STYLE.md) as lookup references afterward.

**Tags used throughout:** `[INV]` invariant — carries unchanged to any
experiment with these data types; breaking it is a defect. `[EX]`
this-experiment example — illustrative only, delete it and the rule must
still stand. `[DR]` decision rule — *if your data look like X, do Y* — this
is what survives a design change.

---

## §0 — How to use this bundle

| File | Read it when |
|---|---|
| `HANDOFF.md` (this file) | day one, front to back |
| [`CONVENTIONS.md`](CONVENTIONS.md) | while building a figure — numbered rules `C1`–`C12` |
| [`RESULTS_STYLE.md`](RESULTS_STYLE.md) | while writing Results prose — numbered rules `R1`–`R12` |
| [`house_style.R`](house_style.R) | sourced by any script that builds a figure or caption |
| [`QUICKCARD.md`](QUICKCARD.md) | pinned reference, one page |

**Authority order** — when two sources disagree, the earlier one wins:
`house_style.R` > `CONVENTIONS.md` > `RESULTS_STYLE.md` > `HANDOFF.md` >
anything in the repository predating this bundle (including
`design_memo.md`, `results_template.md`, and `CHANGELOG_stats.md`, all now
archived in `_superseded/` — see §12).

**If you read nothing else:** the unit of replication is the level at which
treatment was applied (§2); an observation window is a factor, never a
covariate (§5b); every p-value carries an effect size (§6); every result is
conditional on choices the data didn't force, and those choices get checked
(§7); no correction is applied to behavioural p-values, and that must be
stated, not assumed (§8).

---

## §1 — What this pipeline is for

Three data streams, one question shape: does a treatment change **how fish
move as a group**, **what their brains are doing chemically**, and **how
stressed they are physiologically**. Trajectory data (idtracker.ai output)
answers the first at the level of the group's spatial and temporal
organisation, not any one fish's identity. Brain monoamines and their
metabolites/turnover ratios (dopaminergic and serotonergic systems, by
region) answer the second. Plasma cortisol answers the third. All three are
measured on the same experimental units and are expected to be read together,
not as three separate studies bolted side by side.

---

## §2 — The unit of replication

`[INV]` **The unit of replication is the level at which treatment was
applied — not the level at which the response was measured.** Every model
carries a random intercept for it, **forced into every candidate random-effect
structure, never selected by model comparison**. Every figure plots one point
per unit (§4 for the raw-point-per-unit rule). Every caption *N* counts units,
not measurements.

`[EX]` `(1 | phys_trial_id)`; 8 trials per arm × 5 fish per trial; 48 trial ×
interval rows collapse to 16 physical trials.

`[DR]`
- Response measured per-fish but treatment applied per-tank → aggregate to
  the trial level for the primary analysis; keep the per-fish model as a
  **labelled sensitivity check**, never as an alternative primary.
- A unit contributing exactly one observation → its random effect is
  statistically unidentifiable. Drop to a fixed effect and say so in the
  model-selection record, rather than fitting a random intercept that AICc
  cannot meaningfully evaluate.

---

## §3 — The pipeline in seven stages

| Stage | What it does |
|---|---|
| STEP1 | idtracker.ai ingestion — raw per-frame trajectories, jump/identity-switch filtering, zone assignment |
| STEP2 | per-trial indicator computation — speed, occupancy, activity, zone metrics |
| STEP2b | group dynamics — NND, IID, polarisation (computed, excluded from analysis), hull area, centroid speed |
| STEP4 | exploratory graphs — not manuscript figures, a first look |
| STEP5 | statistics + manuscript figures — model fitting, selection, inference, plotting |
| validation | confirms the statistic reported matches the test actually run |
| reporting | builds caption tables, supplementary exports, the manuscript-ready text |

`[INV]` **Stage-boundary law:** each stage writes a timestamped, immutable
directory and never edits an upstream stage's output. Downstream stages
**auto-discover** the newest matching directory rather than being told which
run to use. An indicator is **defined in exactly one file** — if two stages
need the same transform or formula, it moves to a shared module both source,
never pasted twice.

`[DR]` A new indicator needing frame-level data → it belongs in the ingestion
stage even if only the statistics stage consumes it, because that is where
the raw data still exists. Two stages independently needing the same
computation → extract it to a shared module immediately; the moment it is
pasted twice is the moment it starts drifting (see the `number_formatting.R`
seven-copies history in §11).

---

## §4 — Data intake and QC, per data type

### §4.1 Trajectories

`[INV]` QC gates are fixed **before** any indicator is computed: an
identity-swap check, a missing-frame ceiling, and an arena-geometry
reference. Any metric that requires **persistent per-fish identity** across
frames is either validated for identity tracking quality, or excluded from
analysis outright — and the exclusion is stated as a **design ground**, never
silently omitted.

`[EX]` Polarisation is computed but excluded from analysis on
identity-invariance grounds — the arena's design does not support persistent
per-fish identity, so any metric requiring it (polarisation, individual
leadership/consistency measures) is not analysable here without validating
identity separately.

`[DR]`
- Identity validated for the new experiment → identity-dependent metrics
  become available as a **declared family**, analysed and reported as such.
- Identity not validated → restrict to permutation-invariant functions of the
  unlabelled point set — nearest-neighbour distance, inter-individual
  distance, convex hull area, centroid speed, zone occupancy. These do not
  require knowing *which* fish is where, only *how many* and *how far apart*.

`[INV]` An observation window (time bin, interval) is a **factor, never a
covariate**, unless a monotone trend across it is the actual hypothesis being
tested. Window boundaries are declared before analysis, not chosen after
looking at the data.

`[EX]` Three 20-minute intervals, boundaries at 5–25 / 45–65 / 85–105 minutes
→ a 2-degrees-of-freedom factor, tested as an omnibus plus an exhaustive
orthogonal polynomial decomposition (linear + quadratic).

`[DR]` *k* intervals → a (*k*−1)-degrees-of-freedom omnibus test, plus an
**exhaustive** orthogonal polynomial decomposition to order *k*−1 — reported
regardless of whether the omnibus is significant, because the decomposition
characterises shape, not presence of an effect. *k* = 2 → the omnibus test
**is** the linear contrast; report one, not both, under two different names.
Unequally spaced intervals → the contrasts use actual elapsed time, not rank
position — orthogonal polynomials on ranks assume equal spacing they don't
have.

### §4.2 Monoamines

`[INV]` The analyte and region set is declared **before** unblinding, with
turnover ratios named explicitly as derived (metabolite/parent), not as
independent measurements. An analyte excluded from analysis is excluded from
**analyses, tables, and figures alike**, with a dated decision recorded.
Analyte × region defines a **cell**, and the total cell count is stated
wherever the family structure is discussed.

`[EX]` 6 analytes (DA, DOPAC, 5-HT, 5-HIAA, DOPAC/DA, 5-HIAA/5-HT) × 4 regions
= 24 cells. Noradrenaline (NE) excluded throughout — it belongs to neither
the dopaminergic nor serotonergic family and its measurement was judged
unreliable for this experiment; had it been included the cell count would be
28, not 24.

`[DR]` A region is added to the design → the cell count and any declared
multiple-comparison family's *k* both change — restate both, don't let a
stale *k* survive the change. An analyte falls below the limit of detection
in more than a small fraction of samples → report the detection frequency
directly rather than a mean concentration built on imputed zeros.

### §4.3 Cortisol

`[INV]` Terminal blood sampling makes plasma cortisol a **single-timepoint
endpoint with no within-animal baseline** — the design structurally cannot
separate a trait-like baseline difference from an acute response to the
treatment or to sampling itself. This is a **mandatory limitation sentence**
in Methods and Results, not an optional caveat (see `RESULTS_STYLE.md` R9).

`[EX]` Plasma cortisol, ng/mL, one sample per fish, trial-clustered.

`[DR]` Repeated non-terminal sampling becomes available in a future design →
cortisol becomes a genuine repeated-measures endpoint under the same
interval-as-factor rule as §4.1's trajectory windows, and directional/causal
language becomes licensable in a way it currently is not.

**The endocrine↔behaviour join key.** A behavioural trial and an endocrine
sample are linked by **`(tank, testing date, treatment)`**, not by the raw
numeric `trial` field alone — the same integer label is independently
assigned within each data stream and does not denote the same physical trial
across streams. Joining on the bare integer silently mis-pairs a fraction of
the dataset. This warning previously lived only inside an archived changelog
(`CHANGELOG_stats.md` M5); it is restated here because it is intake-level, not
a modelling detail, and belongs wherever the two streams first meet.

---

## §5 — The modelling ladder

`[INV]` A **fixed order**, walked for every outcome, not a menu to pick from:

1. choose the response scale (§ transformations, below)
2. build the random-effect candidate set, with the design-imposed random
   effect **forced into every candidate**
3. select among the remaining candidates by AICc, fit with **ML** (not REML —
   REML likelihoods are not comparable across differing random-effect
   structures)
4. apply the **singularity guard**: if some but not all candidates are
   singular, singular ones are made ineligible; if every candidate is
   singular, the guard keeps the full set rather than discarding everything
5. **refit the selected structure with REML** before running inference
6. fit with **Type III sums of squares** under sum-to-zero (`contr.sum`)
   contrast coding, so main effects are marginal
7. run the inference engine appropriate to the model family (below)

**Design random effects are forced; AICc is only ever allowed to decide
nuisance terms.** Model selection never touches the fixed-effect structure —
the fixed effects are the hypothesis being tested, not a free parameter.

`[DR]` By outcome type:

| Response shape | Family / test |
|---|---|
| Gaussian LMM | F-test, **Kenward-Roger** denominator df primary, **Satterthwaite** as a stated sensitivity |
| bounded (0,1) proportion | beta GLMM, logit link; treatment term tested by **parametric-bootstrap LRT** (see §7 — KR/Satterthwaite are undefined here) |
| counts | Poisson GLMM by default; **dispersion ratio > 1.5** → switch to negative binomial, and name the family actually used in the caption |
| singular fit | the singularity guard drops it to the next eligible candidate, logged in `aicc_selection.csv`, never silently dropped |
| model does not converge | reported **as non-convergence**, never as an absence of effect (`RESULTS_STYLE.md` R7) |

`[EX]` The singularity guard is not cosmetic: on one outcome it moved the
reported test from `F(1,14) = 5.28, p = 0.038` to `F(1,39) = 8.93, p = 0.005`
— the naive AICc winner had a design random-effect variance estimated at
exactly zero, which the guard correctly refused to accept as "the design
random effect doesn't matter here."

---

## §6 — Effect sizes

`[INV]` **No p-value is reported without an effect size, and no effect size
without knowing which standardiser it uses.** An effect size is emitted for
**every** term with a valid p-value, regardless of significance — gating the
effect size on p < 0.05 discards exactly the information that makes a null
result interpretable (Nakagawa & Cuthill 2007).

`[INV]` Three tiers, reported together wherever available:

| Tier | Quantity | Role |
|---|---|---|
| Model term | **partial η²** (`η²p = F·df1 / (F·df1 + df2)`), with a noncentral-F 95% CI | the in-caption, in-text effect size |
| Model term | **partial ω²** (bias-corrected) | **the value to carry into any future power calculation** — not η²p, which is upward-biased at small df2 |
| Contrast | **Hedges' g**, signed, noncentral-t inversion (not Wald), J-corrected | the 1-degree-of-freedom pairwise effect |
| Variance | R²m / R²c / adjusted ICC | how much the model explains, and what clustering costs |
| Biological | **% change** vs. the named reference group | the reader's actual effect size — mandatory in text whenever means are given |

`[INV]` **Two Hedges' g standardisers, reported as a pair and never pooled**:
`g` (total-variance, between-unit standardiser) and `g_resid`
(within-cluster, δ/σ_resid). Averaging or substituting one for the other is a
category error, not a simplification. Likewise a cluster-level *d* and a
fish-level *d* are different quantities and must carry a `unit` label.

`[INV]` **A Wald χ² carries no η²p** — this is enforced structurally in
`house_style.R`'s `effect_size()`, which returns `NA` for `stat_type =
"chisq"`, so a caption builder cannot attach one even by accident. Where a
beta-GLMM treatment term is tested by parametric bootstrap, the reported
effect is the **back-transformed contrast** (odds ratio, ratio of means)
with its CI — not a pseudo-η².

`[EX]` `η²p = 0.738` on the primary flow-preference term. Planning effect
size `d = 0.70` at `ICC = 0.20`, 80% power, α = 0.05 — reused as the TOST
smallest-effect-size-of-interest (see §7). Design-realised power:
`DE = 1 + (m−1)·ICC`, `λ = d / √(2·DE / (m·k))`, `df = 2(k−1)`, at *m* = 5
fish per trial, *k* = 8 trials per arm.

`[DR]`
- Gaussian LMM with an F-test → η²p (+ CI) in the caption, ω²p in the
  supplementary table.
- GLMM / beta GLMM with a Wald χ² → **no η²p**; report the back-transformed
  contrast and its CI instead.
- A 1-degree-of-freedom contrast → signed Hedges' g with a noncentral-t CI,
  both standardisers named.
- Planning the *next* experiment's sample size → use **ω²p** and the
  observed ICC, never the η²p printed in a figure caption.
- An effect size below the display floor prints `< 0.001`, never `0.00`
  (`CONVENTIONS.md` C7).

---

## §7 — Controlling for modelling choices

The premise, `[INV]`: **every result is conditional on choices the data did
not force** — the response scale, the aggregation level, the random-effect
structure, the denominator-df approximation. A result that survives only one
combination of those choices is a property of the analysis pipeline, not of
the animals. So each choice gets a **named check**, the check runs **whether
or not it is expected to matter**, and its output ships as a data artifact
alongside the manuscript (see `CONVENTIONS.md` C12 for the exact filenames).

`[INV]` The mandatory set:

| Choice being probed | Check | Rule |
|---|---|---|
| Denominator-df approximation | Satterthwaite refit alongside every Kenward-Roger F | KR primary, Satterthwaite the sensitivity; disagreement across the α boundary must be reported, not buried |
| Test statistic where KR is undefined | **parametric-bootstrap LRT** for every beta-GLMM treatment term | simulate under the reduced model, refit full and reduced per replicate; `p_pb = (1 + #{LRT_sim ≥ LRT_obs}) / (n_ok + 1)`; **abort if fewer than 50% of replicates converge** |
| Bootstrap resolution | N is a **reporting constraint**, not an implementation detail | smallest reportable p is `1/(N+1)` — never quote a p finer than that. `N = 200` → floor ≈ 0.005; anything landing near a decision threshold is **re-run at N = 1000** before being reported |
| Random-effect structure | **forced-RE sensitivity** — where AICc dropped a design random effect, add it back to the same winning structure | never replace, never touch the fixed effects; compare F, p, ICC, singularity side by side |
| Singular fits | singularity guard inside AICc selection (§5) | excluded candidates stay visible in `aicc_selection.csv`, flagged, never silently dropped |
| Aggregation level / pseudoreplication | refit at the unit level **and** the level above it | report both β, both p, the β ratio, and a direction-agree flag; a conclusion that flips with aggregation level is **not reportable** as a treatment effect |
| Influential units | leave-one-unit-out refit across every unit | report the range of the term's estimate and p across all leave-one-out fits |
| The whole choice set at once | **multiverse / specification-curve** sweep — response scale × aggregation × RE structure × transform, all combinations | the primary endpoint gets one; the curve ships as a supplementary table — this is the honest, checkable version of "we checked" |
| Distinguishing a null from a non-detection | **TOST equivalence** against a declared smallest-effect-size-of-interest | the SESOI is declared **from the power analysis, before looking at the data**; an untested null is "inconclusive", never "confirmed absence" |
| Redundant outcomes masquerading as replication | **collinearity gate** — every secondary metric screened against the primary endpoint | `\|r\| ≥ 0.70` → demoted to descriptive, described as a re-expression of the primary result, never as independent confirmation |
| Design confounds | design-balance check, run **before** any model is fitted | flags factor↔factor 1:1 mapping, `\|cor\| > 0.9`, unintended nesting |
| Family / dispersion | Pearson dispersion ratio | ratio > 1.5 → Poisson gives way to negative binomial, and the caption names the family actually fitted |
| The whole reported set | `validate_stat_types.R` | confirms the string printed in the manuscript matches the test actually run |

`[EX]` Order, density, and tank are perfectly aliased in the current
design — order effects are structurally non-identifiable and must be
described as "order-or-density", never as chronological drift, because the
data cannot distinguish the two explanations.

`[INV]` **`stat_type` is authoritative; `df_denom` is not.** `df_denom` is
overloaded to carry the bootstrap replicate count on parametric-bootstrap-LRT
rows — this once caused nine beta-GLMMs to be mis-reported as `F(1, 1000)`
because a reader trusted `df_denom` over `stat_type`. Read `stat_type` first,
always.

`[INV]` **Reproducibility of every stochastic check**: seed recorded, N
recorded, both stated in Methods. A bootstrap without a stated N and seed is
not a check — it's an unverifiable claim.

`[DR]` Scaling the set to a new experiment:
- One primary endpoint → it gets the full set above. Secondary endpoints get,
  at minimum, the dispersion check, the singularity guard, and the
  aggregation-level check.
- No clustering in the design (e.g. one fish per unit, no shared tank) → the
  pseudoreplication and forced-RE checks are dropped, and the document says
  *why* they were dropped — because the design has one level, not because
  they were skipped.
- A check that changes a conclusion → it is not a footnote. It moves into the
  main text and the primary analysis is reconsidered in light of it.
- A check that is computationally expensive → run it at reduced N first to
  locate borderline cases, then re-run **only those** at full N. Never report
  the reduced-N p-value as if it were final.

---

## §8 — Families, multiplicity, what α means here

`[INV]` Multiple-comparison families are **declared before looking at the
data**, with the family size *k* written down. Whatever regime is chosen —
correction applied, or explicitly none — is stated in Methods, and **every
p-value in the manuscript is labelled adjusted or raw**, never left
ambiguous. A correction declared in Methods must be **computed by a module
the reporting engine actually sources** — a declaration with no code behind
it is worse than no declaration, because it reads as done when it isn't (see
A-2 in §12).

`[EX]` Dopaminergic outcomes (DA, DOPAC, DOPAC/DA) form a **flat primary
family across all four brain regions**, *k* = 12. Serotonergic outcomes
(5-HT, 5-HIAA, 5-HIAA/5-HT) are exploratory, also *k* = 12, flat across
regions. NE is excluded from both families (§4.2). **No adjustment is
applied to behavioural p-values** — raw p at α = 0.05, stated explicitly in
Methods, not left to be inferred from silence.

`[DR]` A CLD family spans more than 2 levels on one factor → `emmeans`
substitutes **Šidák** for the CI adjustment while the letters remain
Tukey-derived internally; the caption and the prose must say Šidák where that
applies, not default to naming Tukey out of habit (`CONVENTIONS.md` C5,
`RESULTS_STYLE.md` R6).

---

## §9 — What is reported, what is withheld

`[INV]` Outcomes reported in the main text are selected on **construct
coverage** — one measure per question the study asks — **never on
significance**. A decomposition, once begun (e.g. the orthogonal polynomial
breakdown of an interval effect), is reported **exhaustively** — a missing
row in that decomposition is itself a selection, even if no one intended it
as one. A panel or sentence must report the **term that actually carries its
claim**: if the claim rests on a main effect and not an interaction, the main
effect is stated first, with the interaction reported beneath it so the null
there is visible rather than implied away.

---

## §10 — Run hygiene and reproducibility

`[INV]` Every stochastic step records its seed. Every pipeline stage writes a
new timestamped directory rather than overwriting the last run. Downstream
consumers **auto-discover** the newest matching directory and **refuse an
in-progress run** (a completed-run sentinel file, checked before use). A
derived cache records which run it was built from, and the consumer halts on
a mismatch rather than silently reading stale numbers. Running
`validate_stat_types.R` against a finished analysis run is **mandatory**
before any figure or caption built from it is treated as final.

`[DR]` Auto-discovery must occasionally be overridden (e.g. rebuilding a
withdrawn run) → a single dated override at the top of one file, with a
stated reason in a comment — never a buried, unexplained pin mid-script (see
§11 and A-12 for what happens when this discipline lapses).

---

## §11 — Known fragilities of THIS repo

Scoped explicitly: these are **repo-specific**, not properties of the method.
A new project built from this handoff should not inherit them.

- **Hardcoded absolute paths** (`D:/…`) scattered across dozens of files
  rather than confined to one config block.
- **`setwd()`** called directly inside several pipeline scripts, making
  output location depend on the working directory a script happens to be
  launched from.
- **Frozen run-directory pins** to a pre-reorganisation path that no longer
  exists — see A-12 in §12; this is the most consequential fragility in the
  repo, because it affects exactly the sensitivity artifacts §7 mandates.
- **Two parallel output trees** for the same pipeline stage, with no single
  source of truth between them.
- **`.BAK_`/`.BACKUP_` sibling files and folders living *inside* directories
  that mtime-based auto-discovery scans** — meaning a stale backup can be
  silently picked up as "the latest run" (§10's auto-discovery guard exists
  specifically because of this).
- The **`number_formatting.R` seven-copies problem** (§3, §11): a formatting
  rule pasted into seven engines instead of sourced once, which is exactly
  the failure mode `house_style.R` in this bundle exists to close off.

---

## §12 — Adjudication register

An unadjudicated contradiction is worse than an absent document — the next
person finds the losing side and follows it. Every row below names what
disagrees, which side wins, and why.

| # | Contradiction | Ruling |
|---|---|---|
| **A-1** | Three incompatible multiplicity regimes: `CHANGELOG_stats.md` M1 (BH-FDR, per brain region, 4 families of 7 analytes, q-values surfaced into figure significance) vs. `METHODS_CHANGES.md` §9 (no adjustment; `p_BH`/`sig_BH` never populated) vs. `correction_methods.R` (flat dopa/sero, *k* = 12) | **`METHODS_CHANGES.md` wins.** No adjustment is in force on behavioural p-values; raw p at α = 0.05, stated in Methods. The dopa/sero family structure is a **declaration only** (see A-2). `CHANGELOG_stats.md` is archived — note its "7 monoamines per region" predates the NE exclusion, so its family size is wrong on two counts, not one |
| **A-2** | `correction_methods.R` declares the dopa/sero family structure but is not sourced by either analysis engine | `[INV]` *a correction declared in Methods is computed by a module the reporting engine actually sources, or the declaration is removed.* Currently the declaration stands with no code behind it — stated plainly rather than left implicit |
| **A-3** | Two live constants both named `TREATMENT_COLORS` hold different palettes — one dead non-Okabe pair, one live Okabe pair, same name; same pattern for `ZONE_COLORS_MAIN` | **Highest-severity drift trap in the repo.** Okabe-only is `[INV]` (`CONVENTIONS.md` C1); the dead constants are named explicitly so they cannot be copied by mistake; the correct values live in `house_style.R` under **different, unambiguous names** |
| **A-4** | `CLD_SIZE <- 4.0` is defined but never read by a draw call; every live `geom_text()` hardcodes `size = 6` | 6 is the convention. `house_style.R` defines exactly one `CLD_SIZE` constant, set to 6 |
| **A-5** | Two theme dialects coexist: `BASE_THEME` + crossbar idiom (behaviour) vs. `theme_legacy` (endocrine — different base size, axis linewidth, jitter width, and mean-mark geometry) | `BASE_THEME` + the crossbar idiom is canonical for any new work; `theme_legacy` is retained only to avoid re-rendering already-approved figures. `[INV]` one theme object, one mark idiom, per manuscript |
| **A-6** | "No box plots" as an absolute rule vs. two live `geom_boxplot`/violin calls in the actual codebase | The rule needs a scope or it gets broken then ignored. Binds anything reaching a manuscript or supplement; assay-QC and pipeline-diagnostic figures are exempt **and must live in a `diagnostics/` tree no manuscript figure assembly reads from** |
| **A-7** | Divergent `.BAK_`/`.BACKUP_` file siblings live **inside** the shared-code directory itself, not in a separate archive location, and differ from the current file | `[INV]` *a backup is a commit, never a sibling file in a directory that gets scanned.* This bundle's own `_superseded/` folder is deliberately **outside** any directory a build script auto-discovers from |
| **A-8** | CLD-letter truncation (`substr(.group, 1, 2)`) — the exact defect the CLD policy document exists to prevent — is still live at five call sites in the current codebase | Never-truncate is `[INV]`, **and** the live violations are named explicitly rather than presenting the rule as universally enforced today |
| **A-9** | Two p-value floors in circulation: `p < 0.001` (current code, dated decision) vs. `p < 0.0001` with appended significance stars (an older template and an older formatter) | **`p < 0.001` wins**, no stars, in prose and captions alike. The older floor and star-appending are superseded and not reproduced in `house_style.R` |
| **A-10** | Descriptive-value precision: an older template says "3 decimals for everything, including means and percentages"; the number-formatting rule explicitly carves descriptives **out of scope** and says to keep each quantity's established precision | The number-formatting rule wins. This is the rule most likely to be misapplied going forward, because "3 decimals for everything" is easier to remember than the correct, narrower rule — stated loudly here for that reason |
| **A-11** | An older figure-design memo states three rules the code has since superseded: CLD letters truncated to 2 characters; letters remapped so "the leftmost cell is always `a`"; treatment "double-encoded" | All three rewritten from the actual code in `CONVENTIONS.md` C2/C5; the memo is archived precisely so the stale text cannot be found and followed by a future reader who doesn't know it's outdated |
| **A-12** | Stale run-directory pins point at a pre-reorganisation path that **no longer exists** — including the validation script and the sensitivity analyses that back the robustness claims in §7 | These are exactly the artifacts §7 mandates, pinned to a directory nobody can regenerate from today. §10 states one-tree + auto-discovery + sentinel-guard as `[INV]`; the current violations are named here as `[EX]`, not silently inherited. `[DR]` *a pin is a single dated override at the top of one file, or it is a bug* |
| **A-13** | Five separate files disagree about the definition of at least one outcome (a flux metric described as a Poisson/negative-binomial count in one place when it is actually a continuous Gaussian-LMM quantity elsewhere) | `CONVENTIONS.md` C11 `[INV]`: *one definition of record per outcome, everywhere.* The five disagreeing sources are named there so a future reader knows where to look and reconcile |
| **A-14** | Bootstrap replicate count was lowered from 1000 to 200 for runtime reasons, moving the smallest reportable p from ≈0.001 to ≈0.005 | `[INV]` N is a **reporting constraint**: the floor is stated in Methods, and any p landing near a decision threshold is re-run at the higher N before being reported as final (§7) |
| **A-15** | **OPEN.** A design-balance diagnostic reports no tank–density confound, contradicting the in-code justification used elsewhere to drop density as a random-effect candidate; compounded by the order/density/tank aliasing noted in §7 | Carried here as **explicitly unresolved**. A register containing only resolved items would not be credible — this one is left open on purpose |

---

## Appendix A — the assertion-selection criterion

Every runnable check in `house_style.R`'s `preflight_check()`, and every item
in `CONVENTIONS.md` Appendix C, corresponds either to a defect that **actually
happened in this project's history** or to a rule this bundle states as
invariant. Nothing is included for a hypothetical failure mode. This is
deliberate: a checklist that checks for imagined problems trains the reader
to skim past it; one that checks for problems that actually occurred earns
attention.
