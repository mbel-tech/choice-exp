# RESULTS_STYLE

How a statistical output becomes a sentence in the Results section. Written
for whoever drafts the manuscript prose — this reader may not run R. For the
exact figures and captions these sentences describe, see
[`CONVENTIONS.md`](CONVENTIONS.md). For why the pipeline is shaped the way it
is, see [`HANDOFF.md`](HANDOFF.md).

**Tags:** `[INV]` invariant, carries to any experiment with these data types.
`[EX]` this-experiment example, illustrative only. `[DR]` decision rule — *if
your data look like X, do Y*.

Every rule is numbered `R1`–`R12`, permanently.

---

## R1 — Anatomy of the Results section

`[INV]` Behaviour before physiology. Within behaviour: trial-level first,
then interval-resolved. Within physiology: **system → analyte → metabolite →
turnover ratio → region**, then circulating endpoints (cortisol). Robustness
comes **last**, as its own subsection — not folded into the primary claim.
Each block opens with **one** orientation sentence naming the analysis unit
and the predictors tested.

`[EX]` "Behavioural responses were analysed at the trial level to test the
effects of treatment, interval, and their interaction on flow-zone
preference."

`[DR]` No time structure in the design → delete the interval-resolved
subsection; don't leave it as an empty heading. A third data type is added →
it becomes a sibling subsection ordered after cortisol, following the same
system → analyte → region logic if it applies.

---

## R2 — The canonical sentence

`[INV]` Every result paragraph instantiates one skeleton, slots in this
order:

> **indicator → model term → statistic → p → effect size → raw mean ± SE by
> group → % change → post hoc → figure reference → cautious interpretation**

A paragraph missing a slot is **incomplete, not concise** — omit only a slot
that genuinely does not apply (e.g. no post hoc for a 2-level factor), never
one that is inconvenient.

**Filled skeleton — significant result:**
> "[Indicator] differed between [groups] ([statistic], p = [value], [effect
> size]). Raw mean ± SE values were [value] in [reference group] and [value]
> in [comparison group], corresponding to a [XX%] [increase/reduction] in
> [comparison group]. [Post hoc contrast if relevant]. These results
> [suggest/are consistent with] [cautious interpretation] (Figure NX)."

**Filled skeleton — null:**
> "[Indicator] did not differ between [groups] ([statistic], p = [value],
> [effect size]). Raw mean ± SE values were [value] and [value], corresponding
> to a [XX%] numerical [increase/reduction] that was not statistically
> supported."

**Filled skeleton — interaction:**
> "[Indicator] showed a significant [factor A × factor B] interaction
> ([statistic], p = [value], [effect size]). [Post hoc contrasts describing
> the shape of the interaction]. [Simple main effect, if it carries the
> claim, stated separately]."

**Filled skeleton — non-convergence:**
> "The [family] model for [indicator] did not converge; this is reported as
> non-convergence, not as an absence of effect."

---

## R3 — Ordering laws

`[INV]` **Statistic → p → effect size, always** — in running text exactly as
in figure captions (`CONVENTIONS.md` C6). The reference group is named
**first**, every time, in the same order the figure's x-axis and legend use.

`[DR]` If the figure and the text disagree on group order, **the figure is
right** — fix the text, not the figure.

---

## R4 — Descriptives

`[INV]` **Mean ± SE, never median [IQR]** — the same summary statistic the
figure draws, so text and figure cannot silently disagree. A back-transformed
value is flagged inline as such, e.g. "(back-transformed EMM)".

`[EX]` "Raw mean ± SE values were 0.661 ± 0.022 in control and 0.868 ± 0.038
in exercise choice."

`[DR]` A genuinely skewed distribution where a median is wanted → it goes in
a supplementary table, and the *figure's* mark changes accordingly
(`CONVENTIONS.md` C4 `[DR]`); the main text still reports mean ± SE.

---

## R5 — Effect sizes in prose

`[INV]` Every statistical claim carries **two** effect sizes, not one:

1. a **statistical** effect size — η²p for an F-test, or the back-transformed
   contrast (odds ratio, ratio of means) where η²p is undefined;
2. a **biological** one — **% change against the named reference group**,
   mandatory whenever means are given.

`[INV]` A Wald χ² sentence carries **no η²p** — it carries the odds ratio or
ratio of means, with direction stated in brackets, instead. `[INV]` A
contrast reported as Hedges' g **names its standardiser** — `g` (total
variance, between-unit) and `g_resid` (within-cluster) are different
quantities and are never pooled or presented as interchangeable.

`[EX]`
- Behavioural % change to 1 decimal place: "corresponding to a 31.4%
  increase."
- Neurochemical % change as a whole number, hedged: "corresponding to an
  approximately 154% increase in exercise choice relative to control."
- Odds ratio with direction: "odds ratio [control/exercise choice] = 0.293, p
  < 0.001."
- Variance-explained trailing sentence, where a mixed model carries the
  claim: "R²m = 0.101, R²c = 0.261, adjusted ICC(trial) = 0.178."

`[DR]` An effect size without a confidence interval is acceptable in a figure
caption (space-constrained); it is **never** acceptable in the text of a
primary claim — the CI belongs in the sentence or the adjacent table.

---

## R6 — Post hoc and CLD in prose

`[INV]` A CLD letter is a **non-difference statement**, not a ranking. Cells
sharing a letter "do not differ" — never "are similar" or "are equivalent",
which claim more than the test supports. A two-letter cell differs from
**neither** group and must be reported that way, not simplified to one side.
A letterless panel means **no pairwise difference at all**, and the text
**says that explicitly** rather than going silent about it.

`[EX]` "Compact letter display: cells sharing a letter do not differ at p <
0.05" (verbatim from `CLD_CAPTION_NOTE`, `house_style.R`).

`[DR]`
- The adjustment applied is Šidák, not Tukey → the text says Šidák.
- Joint Tukey and simple-effects contrasts disagree → report the one the
  **figure displays**, and say which; never quote whichever is more
  favourable to the claim.

---

## R7 — Nulls, non-convergence, exclusions

`[INV]` **Four distinct forms, never collapsed into one register:**

| Form | Written as |
|---|---|
| (a) tested, not significant | statistic, p, effect size, means, % change, "not statistically supported" |
| (b) did not converge | reported as non-convergence, explicitly **not** as absence of effect |
| (c) undefined / excluded | the definitional ground and the number of affected units, stated |
| (d) equivalent | TOST against the declared SESOI passed — a **confirmed** absence, distinct from (a) |

`[INV]` (a) and (d) must never be written the same way. An untested null is
"inconclusive rather than a confirmed absence" — never described as evidence
of specificity or absence without the equivalence test behind it.

---

## R8 — Robustness in prose

`[INV]` The robustness subsection states, for the primary endpoint, **which
choices were varied and whether the conclusion held** — never the bare claim
"results were robust" with no named check behind it. Each check is named with
its artifact (see `CONVENTIONS.md` C12).

`[INV]` A **bootstrap-tested term always names N and the method** in text. A
p sitting at the bootstrap floor is written `p ≤ 1/(N+1)` (e.g. `p ≤ 0.005`
at N = 200), never as a spuriously precise exact value.

`[INV]` A metric that fails the collinearity gate (`|r| ≥ 0.70` against the
primary endpoint) is introduced as **a re-expression of the primary result,
not independent confirmation** — in the running text *and* in its own figure
caption, not one without the other.

`[INV]` A non-identifiable factor (perfectly aliased with another in the
design) is named as such — "order-or-density" — and never given a causal or
directional reading.

`[EX]` "Most of what follows re-expresses the ALR(flow) result in the time
domain rather than confirming it independently (Table S13)."

`[DR]` A check that flips a conclusion → it leads the robustness subsection,
and the primary analysis is reconsidered, not footnoted around. A check that
was not run for this endpoint → say so and why; never leave the omission to
the reader's inference.

---

## R9 — Endocrine specifics

`[INV]` Declared multiple-comparison families and the adjustment status are
stated **once**, in Methods, and Results never contradicts that statement.
Excluded analytes are named as excluded, with the ground stated. Terminal
single-timepoint sampling carries a **mandatory limitation sentence** — the
design cannot separate trait from response. Cautious verbs only: "consistent
with", "may indicate", "suggests", "could reflect" — never causal language
("increased", "caused", "drove") without a design that licenses it.

`[DR]` A directional claim is wanted → it needs a design that supports
direction (e.g. repeated non-terminal sampling); otherwise, rewrite as an
association.

---

## R10 — Cross-referencing

`[INV]` The figure reference sits at the **end of the statistical clause**,
before the interpretation sentence. A panel and its sentence report the
**same term** — if the claim rests on a main effect, the sentence leads with
the main effect too, matching the caption (`CONVENTIONS.md` C6 `[DR]`).

`[EX]` "…; Figure 8A", "(Table S13)", "(3.2.1)". Combined form,
figure-before-table: "(F1,10 = 85.12, p < 0.001, η²p = 0.895; Figure 10D;
Table S19)".

---

## R11 — Forbidden phrasings

- No causal verbs without a licensing design.
- No "trend towards significance" / "marginally significant".
- No star-only reporting — a p-value is always given.
- No `p = 0.000`.
- No unqualified "Tukey" where the adjustment applied was Šidák.
- No adjusted p quoted when none was computed.
- No bare "robust" / "results were consistent" without naming the check
  (R8).
- No η²p attached to a Wald χ² result (R5).
- No bootstrap p quoted finer than the replicate count resolves (R8).

---

## R12 — Worked examples

**EXAMPLE — this experiment only.**

> "Flow-zone preference differed significantly between treatments (F1, 8 =
> 22.58, p = 0.001, η²p = 0.738). The mean ± SE proportion of trial time in
> the flow zone was 0.661 ± 0.022 in control, whereas it was 0.868 ± 0.038 in
> exercise choice, corresponding to a 31.4% increase (Figure 8A)."

**EXAMPLE — this experiment only.**

> "The treatment × zone interaction was highly significant (χ²1 = 49.84, p <
> 0.001); the treatment main effect within zone was not (χ²1 = 0.095, p =
> 0.757). Post hoc contrasts showed exercise choice schools used the flow
> zone significantly more than control (t = −3.83, p < 0.001) and the calm
> zone significantly less (t = 2.37, p = 0.021); the medium and low sub-zones
> did not differ between treatments (both p > 0.05; Figure 8)."

**EXAMPLE — bootstrap-tested term.**

> "The treatment effect on [proportion outcome] was tested by
> parametric-bootstrap likelihood-ratio test (N = 200 replicates, seed
> 20260808), since Kenward-Roger and Satterthwaite are undefined for a beta
> generalised linear mixed model (LRT = 8.41, p ≤ 0.005). Because this result
> sits at the bootstrap floor, it was re-run at N = 1000 to confirm."

**EXAMPLE — equivalence result.**

> "Equivalence testing at a smallest-effect-size-of-interest of d = 0.70
> confirmed that the absence of a treatment effect on dopamine concentration
> reflects genuine equivalence, not merely non-detection, in the Vv, Vd and
> POA regions; the corresponding nulls for 5-HT and for the turnover ratios in
> most regions remain inconclusive rather than confirmed absences, and are
> described as such."

**EXAMPLE — non-convergence.**

> "The beta generalised linear mixed model for `prop_active` did not
> converge at the interval level. This is reported as non-convergence, not as
> an absence of a treatment effect."

---

## Appendix — skeletons and checklist

Lift-and-fill from R2. Before finalising each result paragraph, confirm it
includes:

- clear indicator name
- analysis unit (the unit of replication, not the measurement level)
- model term tested
- correctly formatted statistic (`CONVENTIONS.md` C7)
- exact p-value, or the stated floor (`p < 0.001` / bootstrap floor)
- statistical effect size, with its standardiser named if it is Hedges' g
- raw mean ± SE by group
- biological effect size as % increase/reduction
- consistent group order, matching the figure
- post hoc result if relevant, with the correct adjustment named
- figure/table reference, positioned at the end of the statistical clause
- cautious interpretation, only if the design supports it
- no `p = 0.000`, no forbidden phrasing from R11
