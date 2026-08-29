# Instructions for Claude — Statistical revisions (Hormones & Behavior rejection)

**Manuscript:** *Behavioural and neuroendocrine correlates of swimming exercise choice in
juvenile Atlantic salmon* (Bellio et al.), HB-D-26-00125.

**Purpose of this file.** A future Claude session must implement the **statistical**
modifications demanded by the editor and reviewers, and **propagate every change end-to-end**:
re-fit models → regenerate tables/figures → update every reported number → update Abstract,
Results, Discussion, Conclusions, and Supplementary Materials. No statistic may be updated in
one place and left stale in another.

> **Golden rule — single source of truth.** Every F/χ²/t, df, p-value, %-change, N, asterisk,
> and significance letter in the manuscript must be *regenerated from code* and must match the
> tables and figures exactly. If a number changes anywhere, trace it to every location it
> appears and update all of them. Never hand-edit a statistic that a script produces.

Reviewer codes below match the comments in the annotated manuscript
(`…RECORDS - ANNOTATED (Word).docx`): **R1-4, R1-6, R1-7, R1-8, R2-5**.

**Scope note.** This file covers only the **statistical** modifications (R1-4, R1-6, R1-7,
R1-8, R2-5). The remaining reviewer/editor points — reward-language reframing (ED, R1-1, R1-2,
R1-3, R1-5), and structural/methods items (R2-1, R2-2, R2-3, R2-4, R2-6, R2-7) — are tracked in
`Revision plan - to-do and roadmap.docx` in this same folder. The two documents intersect at
**M1** and **M5** above: once FDR-adjusted q-values (M1) and the neurochemistry↔behaviour
limitation (M5) are settled here, feed the outcome back into that roadmap's Phase 1 language
pass (R1-2) and Phase 4 discussion items (R1-4), since the wording those items need depends on
which effects survive these statistical revisions. Do the statistical work in this file first;
it constrains what language changes are truthful to make in the roadmap.

---

## 0. Before changing anything — build the map

1. **Locate the current analysis assets** (confirm exact paths; these are the likely ones):
   - Manuscript text: `Behavioural and neuroendocrine correlates RECORDS.docx` (and the `.md`).
   - Data: `Behavioural and neuroendocrine correlates datasets FINAL.xlsx`
     (verify vs `FINAL Behavioural and neuroendocrine correlates datasets.xlsx` — use the one
     the manuscript's Supplementary references; confirm with the user if ambiguous).
   - R scripts: `~/OneDrive - Norwegian University of Life Sciences/Documents/EXPERIMENTS/choice exp/`
     and `…/Documents/choice exp for claude/`. The latest-dated / `def` / `review` variants are
     the current ones (e.g. `choice statistics MD may 2025.R`, `cort stats graph choice nov 2025.R`,
     `monoamine B3 tests and graphs NOV 2025.R`, `statistics_bathtub_review.R`,
     `choice exp occupancy through time UPDATED.R`, `trial date ato assignment mono choice exp.R`).
     **Do not guess — inventory them and confirm which script produces which manuscript result.**
2. **Create `stat_map.md`**: a table mapping every reported statistic to
   (manuscript location → Table/Figure → R script + code block that generates it). This map is
   the checklist you will walk for propagation. Nothing gets marked "done" until every row that
   an edit touches is regenerated.
3. Work on **copies** of scripts and data; keep the originals untouched. Set a fixed random seed
   where any procedure is stochastic, and record package versions (`sessionInfo()`).

---

## Required modifications

### M1 — Multiple-comparison control (R1-8)  ← highest priority

**Problem.** The neurochemistry conclusions rest on many analyte × brain-region tests with no
Type-I-error control. The headline effects (Dm: DA p=0.007, DOPAC p=0.022, 5-HT p=0.025) may not
survive adjustment — and the paper's central "reward" claim depends on them.

**Do this.**
1. Define the **test families explicitly** and document the choice:
   - *Preference experiment (Table 2):* treatment effect for each monoamine measure × brain
     region = the primary family (≈ 4 regions × 6 measures: 5-HT, 5-HIAA, 5-HIAA/5-HT, DA,
     DOPAC, DOPAC/DA). Apply **Benjamini–Hochberg FDR** across this family. Report **both raw p
     and adjusted q** (`p.adjust(p, method = "BH")`).
   - Decide and justify whether to pool all 24 tests into one family or split into biologically
     motivated sub-families (e.g. by neurotransmitter system, or a *primary* family = the four
     Dm effects vs a *secondary* family = other regions). State the rule **a priori** in Methods.
   - *Observational study (Table 1):* region contrasts are already Tukey-controlled *within* each
     measure. State this; additionally consider BH across the 6 omnibus ANOVAs and report.
2. Recompute significance using **q-values**. Effects with q ≥ 0.05 are **not** significant —
   report them as such (trends at most).

**Propagate M1 to:**
- **Table 2** — add an adjusted-p (q) column; re-mark bold/significant rows by q.
- **Table 1** — add adjustment note / q column as decided.
- **Figures 11 & 12** (Dm 5-HT/5-HIAA/ratio; DA/DOPAC/ratio) — asterisks must reflect **q**, not
  raw p; update the panels that lose/keep significance.
- **Results 3.2.3.1 / 3.2.3.2** — rewrite each sentence's p to include q; change "significant"
  wording where q ≥ 0.05.
- **Abstract** — the claim "higher dopamine and serotonin signalling in the dorsomedial pallium"
  must match what survives adjustment; soften or qualify if needed.
- **Discussion** (Dm reward-salience paragraphs) and **Conclusions** ("dopaminergic signalling in
  the Dm is particularly relevant") — align with adjusted results.
- **Methods 2.8.2** — add the FDR procedure and family definition.
- **Supplementary** — full table of raw p and q for every test.

### M2 — Experimental unit / pseudoreplication (R1-6)

**Problem.** Fish were tested in groups; occupancy is a group/trial-level measure. The
independent replicate must be the **trial/testing group**, not the individual fish, and this must
be consistent across all analyses.

**Do this.**
1. Audit **every** model. For each, state the replicate and confirm the random-effect structure
   prevents pseudoreplication:
   - *Behavioural occupancy / transitions / collective metrics:* unit = **trial (group)**.
     Confirm scan data are aggregated to trial level (or trial is a random intercept) and that
     N reflects groups (e.g. N = 8 groups/treatment), not individuals or frames.
   - *Cortisol & monoamines:* unit = **individual fish, nested in trial and tank**. Ensure a
     random intercept for **trial (and/or tank/date)** is present so multiple fish from one
     trial are not treated as independent. Refit any model missing this nesting.
2. Wherever a model changes, its **df, F/χ², p all change** — recompute and propagate.

**Propagate M2 to:**
- Every affected model's numbers in **Results** (behaviour, cortisol, monoamines).
- **Methods 2.8** — add one explicit sentence per analysis stating the unit of replication and
  the random structure that enforces it.
- **All Table/Figure captions** — state the replicate and N at the correct level (e.g. "N = 8
  groups of 5 fish"). Figures 6–13 captions and asterisks/letters as needed.

### M3 — Pre-specification & robustness of the analytical strategy (R1-7)

**Problem.** Many model-selection steps (AICc random-effect selection, Box–Cox, ALR/CLR, Jacobs,
beta GLMMs, transformation switches) raise "researcher-degrees-of-freedom" concerns.

**Do this.**
1. Build a **Statistical Decisions Log** (new Supplementary table): one row per decision, marking
   **a priori vs data-driven**, the rule used, and the justification.
2. Add **sensitivity analyses** demonstrating the primary treatment effects are robust to:
   - alternative random-effect structures (not only the AICc-minimal one);
   - with vs without the chosen transformation;
   - ALR vs CLR vs Jacobs' D for zone preference (already partly done — consolidate and report
     side-by-side).
   Summarise as a supplementary robustness table (effect direction, p/q across specifications).
3. If any headline conclusion is **not** robust, soften it in text and propagate.

**Propagate M3 to:** Methods 2.8 (add rationale + reference to the log), new Supplementary tables,
and any main-text claim whose robustness is now qualified.

### M4 — Treatment-order / between-trial-variability confound (R2-5)

**Problem.** In the preference experiment, order was (choice, control, control, choice) over 8
days; Figure 9 shows high between-trial variability. Order/day may confound treatment.

**Do this.**
1. Add **trial order (or trial date)** as a covariate — or as an additional random effect and a
   fixed-effect test — to **all preference-experiment models** (occupancy, transitions, NND, IID,
   polarisation, school area, school speed, cortisol, monoamines). Formally test the order/day
   effect and its interaction with treatment.
2. Report the outcome: if order is non-significant, present it as evidence against the confound;
   if significant, re-estimate treatment effects adjusting for order and re-interpret.

**Propagate M4 to:** Results 3.2 (behaviour, cortisol, monoamines) with the new terms; the
Discussion variability paragraph (currently argues variability away qualitatively — replace with
the quantitative order test); Figure 9/10 interpretation; Methods 2.8.

### M5 — No directional inference for neurochemistry↔behaviour (R1-4)

**Problem.** Neurochemistry was measured only at session end and never linked to the magnitude of
behavioural preference, yet the text implies a directional link.

**Do this.**
1. If feasible, add a **trial-level association** analysis between behavioural flow-preference and
   trial-aggregated monoamine levels (only where both exist for the same trial). Report r and
   p/q. Note the limitation that terminal, individual-level sampling vs group-level behaviour
   constrains this.
2. If not feasible, **state so explicitly** and strip directional wording throughout.

**Propagate M5 to:** Results (new association or explicit statement), Discussion (remove
"underlie/drive" phrasing), Abstract/Conclusions framing.

---

## Global propagation matrix — walk this for EVERY change

| Layer | Action |
|-------|--------|
| 1. R scripts | Refit model with the modification (M1–M5). Keep old + new side by side in the change log. |
| 2. Model output | Regenerate ANOVA / emmeans / `p.adjust` objects. Export to a tidy results file. |
| 3. Tables 1 & 2 | Regenerate from the results file (add q columns, update bold/significant rows, N, notes). |
| 4. Figures 6–13 | Regenerate; update asterisks (use q), significance letters, error bars, and N labels. |
| 5. Results prose | Update every F/χ²/df/p/q/%-change sentence to match layers 3–4. |
| 6. Abstract | Update every quantitative or significance claim that changed. |
| 7. Discussion | Update every sentence that cites a specific effect/p; re-argue points now handled quantitatively (esp. M4 variability). |
| 8. Conclusions | Align headline claims (esp. Dm DA/5-HT "reward") with adjusted results. |
| 9. Methods 2.8 | Document FDR family (M1), replicate/random structure (M2), decisions log + sensitivity (M3), order term (M4), association/limitation (M5). |
| 10. Supplementary | Add: full raw-p/q table, decisions log, sensitivity table, order-effect results; update S5–S7 as needed. |
| 11. Captions | Every Table/Figure caption states replicate N and the correction/adjustment applied. |

---

## Consistency & reproducibility rules

- **One results object, many views.** Tables, figures, and in-text numbers must all read from the
  same regenerated results file — never typed independently.
- **Significance markers follow q, not raw p**, wherever M1 applies. Asterisk legends and letter
  groupings must be regenerated, not edited by hand.
- **Honest reporting.** If an effect loses significance after adjustment (M1) or after adding the
  order term (M4), say so plainly and follow it through to the Abstract/Discussion/Conclusions.
  This directly addresses the editor's core objection that the data overstate "reward."
- **Reproducible.** Fixed seed, recorded `sessionInfo()`, scripts runnable top-to-bottom on the
  `FINAL` dataset. Save a `CHANGELOG_stats.md` (old value → new value → files touched) per result.

---

## Verification / acceptance criteria (do all before declaring done)

1. **Traceability:** every statistic in the manuscript maps to a line in a regenerated script
   (`stat_map.md` fully populated, no orphans).
2. **Cross-document match:** extract all p/q-values and Ns from the manuscript text and diff them
   against the regenerated Tables/Figures — zero mismatches.
3. **Adjustment present:** Methods state the FDR family and correction; Table 2 shows q-values;
   Figures 11–12 asterisks reflect q.
4. **Replicate stated:** each analysis in Methods and each caption names its replicate/N.
5. **Order tested:** every preference-experiment model reports the order/day test (M4).
6. **Robustness reported:** sensitivity + decisions log in Supplementary (M3).
7. **Claims aligned:** Abstract, Discussion, Conclusions contain no significance claim that failed
   adjustment; no directional neurochemistry↔behaviour claim without support (M5).
8. **Reviewer coverage:** produce a short table mapping R1-4, R1-6, R1-7, R1-8, R2-5 → what was
   changed → where, for the response-to-reviewers letter.

---

## Deliverables

- Updated R scripts (copies), runnable end-to-end, seeded.
- Regenerated Tables 1–2 and Figures 6–13.
- Updated manuscript with all propagated numbers (track changes if editing the `.docx`).
- New Supplementary: raw-p/q table, statistical decisions log, sensitivity/robustness table,
  order-effect results.
- `stat_map.md`, `CHANGELOG_stats.md`, and the reviewer-coverage table.

**Do not** start rewriting prose before the models are refit and the results file regenerated —
the numbers drive the text, not the other way around.
