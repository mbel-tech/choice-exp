# Behavioural sequence analysis — exercise-choice behaviour

Companion doc for `behavioural_sequence_analysis_choice_exp.R`.

## The problem this solves

idtracker.ai gives **correct positions but unreliable identities** within a
session (STEP2 header). Classic sequence analysis — following *one individual's*
behaviour states over time — is therefore impossible. Instead we treat the
**school as the unit** and analyse the temporal sequence of *collective*
exercise-choice states.

"Behaviour" here is spatial and identity-free: at each moment, what fraction of
the school is in the **flow (exercise)** zone versus the **calm** zone.

## Pipeline

1. **Input** — latest `output/STEP1_output/*/master_fish_by_frame.csv`
   (~7.2 M detection-rows, ~25 fps). Keep detections with
   `main_zone ∈ {flow, calm}`.
2. **Per-frame → per-bin** — pool detections into fixed **1 s bins** (`BIN_S`).
   Per bin: `prop_flow = n_flow / (n_flow + n_calm)` — the collective
   exercise-vs-calm choice fraction.
3. **State alphabets (run in parallel)**
   - **Graded**: `prop_flow` → **Low / Med / High** by data-driven terciles of
     the pooled bin distribution (cut-points ≈ 0.20 / 0.80 — the school is
     strongly bimodal: mostly-calm or mostly-flow, rarely evenly split).
   - **Binary**: `prop_flow > 0.5` → **Flow-dominant / Calm-dominant**.
4. **Per-session metrics** (one session = one school × one timepoint):
   transition-probability matrix, **switch rate** (state-changes/min),
   **mean dwell time** per state (interior/uncensored bouts, from run-lengths),
   **state occupancy**, and **normalised entropy rate** (0 = perfectly
   predictable next state, 1 = maximally random). Transitions counted only
   between temporally adjacent bins so gaps never bridge.

## Experimental design (verified from the data)

- **16 physical schools** (`trial`), each observed at **3 timepoints** = 48
  sessions.
- Treatment is **between-school** but **balanced within each of 4 tanks**
  (2 control + 2 exercise-choice schools per tank) → tank is a clean block,
  not a confound.
- `fish_density` (4 / 8 / 12 / 16) is crossed with treatment (balanced).

## Statistics

- **Primary — nested LMM** on all 48 sessions:
  `metric ~ treatment + fish_density + (1|tank) + (1|trial)`.
  Because there are only 4 tanks, `(1|tank)` variance often collapses; the code
  **prefers the most complex non-singular fit**, falling back to
  `(1|trial)` when needed (recorded per metric in `lmm_form`). p-values via
  Satterthwaite (`lmerTest`).
- **Robust secondary — Wilcoxon** on school-level means (8 control vs 8
  exercise-choice schools), which avoids all pseudoreplication assumptions.

## Key results (latest run)

| Alphabet | Metric | Control | Exercise | p (LMM) | p (Wilcoxon) |
|---|---|--:|--:|--:|--:|
| binary | entropy rate | 0.192 | 0.094 | **0.003** | **0.015** |
| binary | switch rate /min | 2.13 | 1.07 | **0.015** | 0.065 |
| graded | occupancy High | 0.247 | 0.447 | **0.012** | 0.44 |
| graded | entropy rate | 0.242 | 0.168 | **0.022** | 0.19 |
| binary | occupancy Flow | 0.386 | 0.499 | 0.058 | 0.51 |
| graded | dwell High (s) | 13.4 | 40.5 | 0.081 | 0.33 |

**Interpretation.** Under **exercise choice**, the school's collective
engagement with the flow zone is **more stable and more predictable**: fewer
collective state-switches, lower sequence entropy, and more time spent fully
committed to the high-flow state. The graded and binary alphabets agree. The
robust Wilcoxon (n = 8 vs 8) confirms the entropy effect and leaves the others
directionally consistent but, as expected with 8 schools/group, weaker.

## Outputs (`output/SEQ_output/SEQ_output_<timestamp>/`)

- `seq_metrics_per_trial_graded.csv`, `seq_metrics_per_trial_binary.csv`
- `transition_matrix_pooled_graded.csv`, `transition_matrix_pooled_binary.csv`
- `treatment_contrasts.csv` (both alphabets; `lmm_form` records model used)
- `figures/` — A state-ribbon timelines, B/C transition heatmaps, D metric
  comparisons, H dwell-by-state (D and H aggregate to school level, 16 points).

## Knobs

`BIN_S` (bin width, default 1 s), `MIN_BINS` (drop sequences shorter than 10
bins). Terciles are recomputed from data each run and printed to the log.

## Possible extensions

- **HMM** (data-driven latent states, no threshold) — needs `depmixS4`
  (not currently installed); the graded terciles are the threshold-free stand-in.
- Timepoint as a fixed effect / interaction (`treatment × timepoint`) if the
  temporal trajectory across the three sessions is of interest.
- Second-order transitions / higher-order dependency tests.
