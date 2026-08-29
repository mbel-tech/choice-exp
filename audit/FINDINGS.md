# Audit findings — Pass 1 and Pass 3 (Results outcomes)

Run 2026-08-18 against `STEP1_output_20260807_162823`, `STEP2_output_20260807_183736`,
`STEP2b_output_20260807_214103`, `BOUT_output_20260816_134054`.

Rule applied throughout: **find, don't fix.** Nothing below has been changed.

Evidence files: `audit_pins.csv`, `audit_chain.csv`, `audit_recompute.csv`,
`audit_recompute_v2.csv`, `audit_recompute_detail.csv`.

---

## What passed

Worth stating plainly, because it is most of the surface area.

**All six bout outcomes reproduce exactly** from an independent implementation written
from the Table 1 definitions — 1 s binning, majority vote, segment splitting at
untracked gaps, edge-censoring, and the non-edge filter. `n_bins`, `max_flow_bout_s` and
`max_calm_bout_s` are bit-exact; the means and the rate agree to 5e-14 relative. The
whole engagement-pattern chain is sound.

**All zone-use outcomes reproduce exactly** once the denominator is read from the code
rather than from Table 1 (see F3). Agreement is 1e-14 to 1e-16 relative across all 48
sessions for `prop_time_in_flow`, `prop_time_in_calm`, the four area-corrected sub-zone
proportions, and `zone_flux_per_session`.

**14 of 15 derived-layer transform checks pass**: the area correction is applied exactly
once, the sub-zone composition closes to 1, the flux and switch formulas reproduce from
stored intermediates, bout rates use non-edge counts, and all four collective metrics
divide by the same trial `length_unit`.

**One suspicion tested and dismissed.** The bout module filters on `main_zone` only and
never applies STEP2's `valid_pos` filter, which looked like two different denominators
for "the detected school". Checked directly: of the 19,248 detections lacking a finite
interpolated position, **zero** carry a flow/calm label, so they cannot enter the
majority vote. Harmless. Not a finding.

---

## F1 — Five scripts pin a STEP5 run that is not on disk

**Class C. Pass 1.** `audit_pins.csv`

All five point at `STEP5_stats_20260511_175522` (11 May) *and* use the
pre-reorganisation path `choice R pipeline/STEP5_stats/` rather than
`choice R pipeline/output/STEP5_stats/`, so the directory does not exist:

| file:line | role |
|---|---|
| `05_stats_validation/_run_validate.R:1` | the validation entry point |
| `04_reporting/build_reviewer_response_report.R:36` | reviewer response report |
| `00_main/export_indicators_excel.R:10` | indicator workbook export |
| `00_shared/effect_size_and_sample_size.R:75` | effect sizes / sample size |
| `02_utilities/diag_cells.R:3` | cell diagnostics |

Severity: these fail loudly rather than silently, but the validation entry point being
among them means the audit trail has not run against a current model set since May.

## F2 — Four sensitivity analyses pin the 8 August run

**Class C. Pass 1.** `audit_pins.csv`

`sensitivity_forced_trial_re.R:34`, `_v2.R:29`, `sensitivity_part2_loo_icc_equiv.R:40`,
`sensitivity_multiverse.R:39` pin `STEP5_stats_20260808_123821`; three also pin
`STEP2b_output_20260807_214103`.

These paths resolve, so they run and produce output — silently, against a model set that
predates the 17 August re-run and the 18 August area-normalisation fix. Given that
today's open question is itself random-effect sensitivity on ALR(medium), the
forced-trial-RE and multiverse scripts are the ones that would normally settle it, and
they are the ones reading stale input.

## F3 — Table 1's denominator is not the denominator in the code

**Class E. Pass 3, raw leg.** `audit_recompute_v2.csv`

Table 1 defines *r* as "the proportion of detections (i.e. fish) present in a given zone
per frame". STEP2 §5 restricts to detections with a finite interpolated position
(`valid_pos <- is.finite(x_interp) & is.finite(y_interp)`); 19,248 rows (0.269% of
7,156,680) fail that test and are excluded from the denominator.

Recomputing under both candidates settles it — the pipeline reproduces the as-coded
denominator, not the documented one:

| outcome | vs documented denominator | vs as-coded denominator |
|---|---|---|
| `prop_time_in_flow` | 35/48 disagree, worst 6.5% | **0/48, worst 3e-16** |
| `prop_time_ac_low` | 35/48 disagree, worst 6.9% | **0/48, worst 2e-14** |
| `zone_flux_per_session` | 35/48 disagree, worst 75% | **0/48, worst 5e-16** |

The code is right; the sentence in Table 1 is incomplete. Fix is one clause, but the
definition as written is not reconstructible from the manuscript, which is the point of
the finding.

## F4 — `zone_flux_per_session` is fragile to denominator noise

**Robustness, not a defect. Pass 3, raw leg.**

A 0.269% change in which rows count shifts flux by up to **75% relative** on a single
session, against ~6% for the occupancy proportions. Flux is a sum of frame-to-frame
absolute differences, so any per-frame denominator perturbation adds a spurious blip at
that frame *and* the next, and in a committed session with few real transitions those
blips dominate the sum.

This is the reported crossings outcome (manuscript §3.1.2, Figure 5A). Nothing is
computed wrongly, but the metric amplifies small upstream changes by more than an order
of magnitude relative to the other zone outcomes, and it is worth knowing that before
any future change to tracking or interpolation.

## F5 — ALR(flow)'s calm reference is `1 − p_flow`, and the gap is treatment-correlated

**Class B. Pass 3, derived layer.** `audit_chain.csv`

The code computes `log(p_flow / (1 − p_flow))`. But `p_flow + p_calm` runs 0.886–1.000:
up to 11.4% of valid-position detections carry no flow/calm label, so `1 − p_flow` is not
`p_calm`. Those unlabelled detections are not evenly distributed:

```
control          mean 3.89%   max 11.4%
exercise choice  mean 0.78%   max  6.7%
treatment effect on the unlabelled fraction:  p = 0.0016
```

The reference is therefore inflated more for controls, in the direction that widens the
treatment gap — the same one-sided shape as the `.eps_lr` floor bug fixed earlier today.

**It does not change the conclusion.** Refitting with the Haldane-floored
`log(p_flow/p_calm)`, consistent with how the sub-zones are handled:

| | treatment | interval | treatment × interval |
|---|---|---|---|
| as coded, `1 − p_flow` | F=17.15, p=0.0010 | F=12.92, p<0.001 | F=7.02, p=0.0034 |
| floored `log(flow/calm)` | F=15.75, p=0.0014 | F=13.28, p<0.001 | F=7.63, p=0.0023 |

And the code's choice is defensible: one session has `p_calm = 0` exactly (a school that
never entered calm), which makes the literally-documented ratio unfittable — hence the
`prop_flow > 0 & prop_flow < 1` filter alongside it. This is a guard whose documented
meaning drifted from its behaviour, not a miscalculation.

Note the manuscript is now internally inconsistent on this point: Table 1 was corrected
today to `log(r flow / (1 − r flow))`, which is accurate to the code, while the Figure 4
caption still reads "the log of time in the named zone over time in the calm zone".

---

## Not covered by this run

- **Collective movement (NND, IID, school area, school speed)** were checked only at the
  derived layer (F-pass: consistent `length_unit` across all four). The raw leg needs
  every position of a frame in memory at once and is queued behind the running STEP5 job.
- **Pass 3 takes `master_fish_by_frame.csv` as ground truth.** Nothing here validates
  STEP1's tracking, interpolation, identity-switch correction or zone geometry.
- Passes 2 (duplicated-block drift) and 4 (documentation reconciliation) were not run.
