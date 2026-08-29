# QUICKCARD

One page. Full detail in `HANDOFF.md` / `CONVENTIONS.md` / `RESULTS_STYLE.md`.
Values below are pointers — copy from `house_style.R`, never retype.

**Unit of replication** = level treatment was applied at, not measured at.
Forced random intercept, never selected. One plotted point per unit.

**Palette** (Okabe–Ito): control `#D55E00` △ dashed · exercise choice
`#0072B2` ● solid. Zones: flow `red4` · calm `lightskyblue3` · high
`tomato3` · medium `mediumturquoise` · low `goldenrod3`. **Statistical ink is
always black.**

**Theme**: `theme_minimal(base_size=13)`, no grid, black axis lines 0.85,
bold axis text, x-title blanked when self-describing.

**Marks**: jittered points → black mean crossbar → black ±1 SE bar. **No box
plots.** Raw data only, never EMMs.

**CLD** (in order, never reorder): `cld()` → drop redundant letters → remap
control-first → guard. **Never truncate.** Omit entirely if uniform.

**Caption**: `statistic → p → effect size`, always that order.
`F1,8 = 22.58, p = 0.001, η²p = 0.738`. **Wald χ² gets no η²p — structural,
not optional.**

**Numbers**: `fmt_F` = 4 sig figs. `fmt3` = 3 dp, one trailing zero trimmed.
Floor `< 0.001` for stats/ES (df exempt). **`p = 0.000` never appears.**
Bootstrap p never finer than `1/(N+1)`.

**Effect sizes**: every p gets one. η²p (F-test) / no ES (χ², back-transform
the contrast instead) / signed Hedges' g (contrast, name the standardiser) +
mandatory % change vs. reference group. ω²p, not η²p, feeds future power
calcs.

**Robustness**: every modelling choice gets a named check with a named CSV —
KR vs. Satterthwaite, singularity guard, forced-RE, pseudoreplication,
leave-one-out, multiverse, TOST equivalence, collinearity gate, design
balance, `validate_stat_types.R`. Bootstrap: state N + seed; abort <50%
convergence.

**Multiplicity**: family declared before looking, *k* stated, adjustment
status (or "none") stated in Methods — never left implicit.

**Time structure**: interval is a factor, never a covariate. *k* intervals →
(k−1)-df omnibus + exhaustive orthogonal polynomial decomposition.

**Reporting order** (`RESULTS_STYLE.md` R2): indicator → term → statistic →
p → effect size → mean±SE by group → %change → post hoc → figure ref →
cautious interpretation.

**Nulls — 4 distinct forms, never collapsed**: (a) tested n.s. (b)
non-convergent — say so, not "no effect" (c) excluded — state the ground (d)
equivalent — TOST-confirmed absence, not just untested.

**Never**: causal verbs without a licensing design · "trend towards
significance" · star-only reporting · unqualified "Tukey" where it was
Šidák · bare "robust" without naming the check.

**Sensitivity artifacts** (must exist per model dir, or a stated reason for
absence): `anova_satterthwaite.csv` · `pb_lrt_<term>.csv` ·
`aicc_selection.csv` · `sensitivity_forced_re.csv` · `sensitivity_loo.csv` ·
`multiverse_<endpoint>.csv` · `equivalence_tost.csv` ·
`pseudoreplication_check.csv` · `metric_collinearity.csv` ·
`design_balance.csv` · `stat_type_validation.csv`.

**Authority order**: `house_style.R` > `CONVENTIONS.md` > `RESULTS_STYLE.md`
> `HANDOFF.md` > anything predating this bundle.
