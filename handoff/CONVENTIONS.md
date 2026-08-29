# CONVENTIONS

Numbered lookup rules for figures and their statistical annotation. This is a
**reference to consult while building a figure, never to read front to back**.
For the narrative — why the pipeline is shaped this way, the modelling ladder,
robustness checks — see [`HANDOFF.md`](HANDOFF.md). For how a statistic
becomes a sentence in the Results section, see [`RESULTS_STYLE.md`](RESULTS_STYLE.md).

Every rule number is permanent. A retired rule keeps its number, struck
through, rather than being renumbered out from under anyone who cited it.

**Tags used throughout:** `[INV]` invariant — carries to any experiment with
these data types, breaking it is a defect. `[EX]` this-experiment example —
illustrative only, delete it and the rule must still stand. `[DR]` decision
rule — *if your data look like X, do Y*.

**Authority.** [`house_style.R`](house_style.R) is canonical for every value
that is a number or a string. This document explains and cross-references it;
it never restates a value that could drift out of sync. Where a fenced block
below is marked `DO NOT RETYPE — COPY`, copy it from `house_style.R`, not from
this file.

---

## C1 — Palette

`[INV]` Okabe–Ito, colourblind-safe. Treatment levels take Okabe indices in a
fixed order. Zone/state colours are a **separate, non-overlapping family** so
a zone-coloured panel can never be mistaken for a treatment-coloured one.
**Statistical ink is always black** — means, error bars, CLD letters,
significance brackets. Colour carries group identity; black carries
inference; the two channels never mix.

`[EX]`

| Role | Value | Okabe index |
|---|---|---|
| control | `#D55E00` (vermillion) | 6 |
| exercise choice | `#0072B2` (blue) | 5 |
| flow | `red4` | — (zone family) |
| calm | `lightskyblue3` | — |
| high | `tomato3` | — |
| medium | `mediumturquoise` | — |
| low | `goldenrod3` | — |
| sequential ramp | `white → #0072B2` | — |
| diverging ramp | `#B2432F → white → #2E5FA3`, centred on 0 | — |

`[DR]`
- 2 levels → Okabe 6, then 5, control (or reference group) first.
- 3–4 levels → extend along Okabe 6, 5, 3, 1. Never index 4 (`#F0E442`, yellow
  — invisible on white) or index 8 (`#000000` — reserved for statistical ink).
- \>4 levels → not a colour variable; facet it instead.
- A new categorical family (a third treatment axis, sex, region) gets its
  **own** C1 entry and never borrows the treatment pair.

**ADJUDICATION A-3.** `figures_step5_module.R:563-565` defines two dead
constants — `TREATMENT_COLORS <- c(control = "#2166AC", "exercise choice" =
"#D6604D")` and `ZONE_COLORS_MAIN <- c(flow = "#4DAC26", calm = "#8073AC")` —
that are **not** Okabe-Ito and have zero live draw calls in that file. They
share a name with the live, correct objects used elsewhere in the same
codebase. **Do not copy either constant from that file.** `house_style.R`
deliberately does not reproduce them, under any name.

---

## C2 — Encoding channels

`[INV]` Treatment is **triple-encoded**: colour **and** shape **and** (on line
plots) linetype. A reader in greyscale must still be able to tell the groups
apart. One channel carries one variable per figure — if colour means
treatment in one panel of a figure, it may not mean state or zone in another
panel of the same figure.

`[EX]` control = shape 17 (filled triangle), linetype dashed. exercise choice
= shape 16 (filled circle), linetype solid.

This corrects `design_memo.md`'s "double-encoded" — the linetype channel is
live on every line plot (`house_line()` in `house_style.R`), making the
encoding triple, not double.

---

## C3 — Theme

`[INV]` `theme_minimal(base_size = 13)`, then everything decorative stripped:
**no grid at all**, no panel background, no border, no strip background. Two
black axis rules; bold text on the categorical axis. **The x-axis title is
blanked whenever the axis is self-describing** — do not print "Treatment"
under labels that already read "Control" / "Exercise choice". Two-word
category labels **wrap, never shrink or rotate**. En-dashes in numeric ranges
(`5–25`, not `5-25`).

`[EX]` The full element table — copy from `house_style.R`, block `A2`:

```r
# DO NOT RETYPE — COPY from house_style.R
BASE_THEME <- ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(
  panel.grid = element_blank(), panel.background = element_blank(),
  axis.line.x = element_line("black", linewidth = 0.85),
  axis.line.y = element_line("black", linewidth = 0.85),
  axis.ticks  = element_line("black", linewidth = 0.7),
  axis.text.x = element_text(size = 13, face = "bold", color = "black"),
  axis.text.y = element_text(size = 11, color = "black"),
  axis.title  = element_text(size = 13, face = "bold", color = "black"),
  axis.title.x = element_blank(),
  legend.title = element_text(size = 11, face = "bold"),
  legend.text  = element_text(size = 10),
  strip.text   = element_text(size = 12, face = "bold"),
  strip.background = element_blank(),
  plot.margin  = margin(8, 10, 16, 10),
  plot.caption = ggtext::element_markdown(size = 14, hjust = 0, face = "italic",
                                          lineheight = 1.3, margin = margin(t = 10)),
  plot.subtitle = element_text(size = 9, hjust = 0, face = "italic", colour = "grey40")
)
```

Wrapped label example: `"Exercise\nchoice"`, `lineheight = 0.9`, centred.

`[DR]` A label exceeds ~14 characters even wrapped → shorten the label; never
rotate it and never reduce `base_size` to make it fit.

---

## C4 — Marks and geometry

`[INV]` House scatter idiom, in draw order:

1. raw points, jittered horizontally only (`height = 0`), coloured by group;
2. group mean, black `geom_crossbar`;
3. dispersion, black `geom_errorbar` at **mean ± 1 SE**.

**No box plots.** **Raw data is plotted, never EMMs** — EMMs feed CLD letters
and contrast tables only, they are never drawn as points or bars. Nothing is
plotted that hides *n*.

`[EX]` Geometry constants (`house_style.R` block `A3/A4`):

| Constant | Value |
|---|---|
| `JITTER_W` | 0.12 |
| `PT_SIZE` | 2.4 |
| `PT_ALPHA` | 0.80 |
| `MEAN_W` (crossbar width) | 0.40 |
| `LW_MEAN` | 0.85 |
| `ERR_W` | 0.13 |
| `LW_ERR` | 0.65 |
| line plots | `linewidth 0.9`, point size ≈ 3.8, error bar width 0.10 |

`[DR]`
- *n* per cell exceeds ~40 so points overplot → the mark changes to an
  **explicit distributional** one and is **labelled as such** in the caption.
  It never silently becomes a box plot of units.
- Model fit on a transformed scale → **plot the raw scale**; state the model
  scale in the caption. Never substitute back-transformed EMMs for the raw
  data layer.
- Assay-QC or pipeline-diagnostic figures fall outside this rule entirely and
  belong in a `diagnostics/` tree no manuscript figure assembly reads from
  (see A-6).

**ADJUDICATION A-6.** "No box plots" needs this scope or it gets broken then
silently ignored: two `geom_boxplot`/violin calls exist in the repo (a
cortisol violin+box kept for a Word appendix, an assay intra-day CV QC
figure). The rule binds anything that can reach a manuscript or supplement;
assay-QC and pipeline diagnostics are the only exemption, and only inside
`diagnostics/`.

---

## C5 — Compact letter display (CLD)

`[INV]` Fixed, ordered pipeline — never reorder these four steps:

1. `multcomp::cld(em, Letters = letters, adjust = "tukey")`
2. **drop redundant letters** — a letter that conveys no non-difference claim
   not already carried by another letter on the same cell is removed
3. **remap** — control-first (default) or `lone_first` (per-panel opt-in)
4. **post-condition guard** — verifies control/interval-1 begins with `a`
   under `control_first` mode, warns rather than silently passes if not

Letters are **never truncated** — truncation destroys sharing relations and
changes what the panel asserts (worked example: `abc | abc | abc | c` becomes
`ab | ab | ab | c` under a 2-character cut, which silently claims a difference
the model does not support). If every cell shares one group, the CLD is
**omitted entirely** and the caption states that explicitly. Remapping is
**relabelling only** — it permutes letter identities, never which cells share
one. A cell may legitimately carry 2–3 letters and must not be "cleaned up"
to a single letter. Letters are black, bold, centred on their own cell,
positioned `CLD_LETTER_FRAC = 0.15` of the panel's data range above that
cell's highest plotted point; the panel expands by `CLD_EXPAND_FRAC = 0.26` so
nothing clips. A colliding letter moves **vertically** into the gap between
the two groups' error bars — never sideways off its own cell.

`[EX]` Letter size **6** in every live `geom_text()` call.

`[DR]`
- Supplied letters inconsistent with the p-matrix → the reduction step is
  **skipped with a warning**, never forced through.
- A panel whose engine produces no CLD at all → the caption states the
  omission explicitly; it is never silently absent.
- A family with >2 levels on one factor → `emmeans` substitutes **Šidák** for
  the CI adjustment while the letters remain Tukey-derived; state which
  applies (see C6, R6).

```r
# DO NOT RETYPE — COPY from house_style.R block A8
cld_result <- cld_pipeline(cld_df, cell_ids, pmat, alpha = 0.05,
                           mode = "control_first", label = "my_outcome")
```

**ADJUDICATION A-4.** `CLD_SIZE <- 4.0` in `figures_step5_module.R:575` is
dead — no draw call reads it. The three live `geom_text()` calls hardcode
`size = 6`. `house_style.R` defines exactly one `CLD_SIZE` constant, set to
6, and any geom in new work reads it — never redefine a second one.

**ADJUDICATION A-8.** CLD truncation (`substr(.group, 1, 2)`) — the exact
defect this rule exists to prevent — is **still live at five call sites**:
`figures_step5_module.R:1029` and `activity_analysis_STATS_choice_exp.R:5531,
5559, 6109, 6535`. State never-truncate as `[INV]`, and do not present it as
uniformly enforced in the current repo — Figure 8 C/D read through the
truncating path. Any new work must not add a sixth site.

---

## C6 — The statistics caption line

`[INV]` One line, bottom-left (`hjust = 0`), italic, black, rendered with
`ggtext::element_markdown()`. **Order is law: statistic → p → effect size.**
Denominator df are displayed as **integers** in the caption (fractional
Kenward-Roger/Satterthwaite values stay in `anova.csv`). The caption breaks to
a second line only after the p-value, and only where the panel is too narrow
— measure the longest caption against the *rendered* image, don't predict it.

`[EX]` `F<sub>1,8</sub> = 22.58, p = 0.001, η<sup>2</sup><sub>p</sub> = 0.738;
Figure 8A` — verified byte-for-byte reproducible by `anova_caption()` in
`house_style.R` against the manuscript's own example.

`[DR]` Which effect size the caption carries, by test type:

| Test | Caption carries |
|---|---|
| Gaussian LMM, F-test | `η²p` |
| GLMM / beta-GLMM, Wald χ² | **none** — the caption must not invent one |
| Parametric-bootstrap LRT (beta-GLMM treatment term) | the bootstrap p and N (`fmt_p_boot`), plus the back-transformed contrast in the caption text |

Claim resting on a main effect rather than the interaction → the caption
reports the **main effect first**, the interaction beneath it, so the null is
visible but the panel does not appear to rest on it.

```r
# DO NOT RETYPE — COPY from house_style.R block A7
anova_caption("Treatment", "F", stat = 22.58, df1 = 1, df2 = 8, p = 0.001)
# -> "Treatment: F<sub>1,8</sub> = 22.58, p = 0.001, η<sup>2</sup><sub>p</sub> = 0.738"

anova_caption("Treatment x Zone", "chisq", stat = 235.7, df1 = 1, p = 0.0001)
# -> "Treatment x Zone: χ<sup>2</sup><sub>1</sub> = 235.7, p < 0.001"   (no eta2p — structural, not a rule to remember)
```

---

## C7 — Number formatting

`[INV]` Exactly **two** base formatters, in one file, sourced not copied.

- `fmt_F` — F-statistics and Wald χ²: **4 significant figures**, trailing
  zeros trimmed **only when a decimal point is present** (so an integer
  denominator df of `30` never becomes `3`).
- `fmt3` — p, η²p, ω²p, Cohen's d/Hedges' g, R²m/R²c, ICC: **3 decimal
  places**, with a **single** trailing-zero trim (`0.500 → "0.50"`, never
  `"0.5"`).

Floor: any statistic or effect size with magnitude **< 0.001 prints `<
0.001`** (`fmt_Fstat`/`fmt_es3`), never a long decimal or an uninformative
`0.00`. **Denominator df are exempt** — a df is never below 1. **`p = 0.000`
may never appear** anywhere — the p floor is `p < 0.001` (`fmt_p_house`).

`[INV]` **A bootstrap p is never printed finer than `1/(N+1)`** — the
replicate count the check actually ran. `fmt_p_boot(p, N)` enforces this;
at the floor it prints `p ≤ <value>`, not an exact decimal (see C12).

Descriptives — concentrations, cm/s, masses, and other quantities that
predate this formatting pass — are **explicitly out of scope** and keep
their established precision.

**ADJUDICATION A-9.** `results_template.md` states a `p < 0.0001` floor and
`analysis_b3_REVISED.R`'s `fmt_p()` implements it, appending significance
stars. Both are **superseded**. `p < 0.001`, no stars, is the current and
sole rule (dated author decision, 2026-08-07). `fmt_p_house()` in
`house_style.R` implements the winner; `fmt_p()` is not reproduced.

**ADJUDICATION A-10.** `results_template.md` also states "3 decimals for
everything, including means, SEs, percentages" — this **directly
contradicts** `number_formatting.R`'s explicit descriptives-out-of-scope
carve-out. `number_formatting.R` wins. State this scope boundary loudly: it
is the rule most likely to be misapplied, because "3 decimals for everything"
is easier to remember than the correct, narrower rule.

```r
# DO NOT RETYPE — COPY from house_style.R block A5
fmt_F(30)         # -> "30"       (never "3")
fmt3(0.5)         # -> "0.50"     (never "0.5")
fmt_Fstat(0.0009) # -> "< 0.001"
fmt_p_house(0.0004) # -> "p < 0.001"
fmt_p_boot(0.003, 200) # -> "p ≤ 0.005"   (floor is 1/201)
```

---

## C8 — Panel assembly and legends

`[INV]` Numbered manuscript panels carry **no `plot.title`** — the title
lives in the sibling caption text file (C10). Tags are `A, B, C…`, bold, 18
pt, top-left, in the margin. **One legend per figure**, horizontal, centred
below the bottom row of panels; individual panels set `legend.position =
"none"` and a legend is never repeated between rows. Grids are built with
`patchwork::wrap_elements(full = …)` so each panel keeps its own margins and
captions align across a row. **Every panel in one figure is drawn by one
builder function** — not each analysis engine's own — so sizes and widths
cannot mismatch inside a figure.

`[DR]` Standalone/supplementary figures may relax exactly three things, and
only these three: keep a `plot.title`; keep their own individual legend;
use black (not coloured) cell outlines on grid-summary heat maps (they carry
no group channel to spend colour on).

**ADJUDICATION A-5.** Two theme dialects exist in the repo: `BASE_THEME` +
the crossbar idiom (behaviour) vs `theme_legacy` (endocrine — base 14, axis
line 1.0, jitter 0.30, mean drawn via `geom_segment` + `geom_point` instead
of a crossbar). `BASE_THEME` + the crossbar idiom is canonical for **any new
work**; `theme_legacy` is retained only to avoid re-rendering already-approved
endocrine figures. `[INV]` *one theme object and one mark idiom per
manuscript, sourced by every engine* — do not present the two dialects as a
single style.

---

## C9 — Export

`[INV]` **PNG and PDF**, 300 dpi, dimensions in **mm**, `bg = "white"`. PDFs
go through `cairo_pdf` so fonts embed — standard PDF fonts render the caption's
sub/superscripts as blank boxes otherwise. **One output tree per run.** A
filename encodes figure identity, never a date.

`[EX]` Working sizes: multi-panel grids 260–280 × 220 mm; single panels 140 ×
130 mm; standalone supplementary panels 120–200 mm wide. `ggsave` with
`cairo_pdf` takes inches, hence the `mm / 25.4` conversion in `house_save()`.

`[DR]` A figure is rebuilt → overwrite in place. **Never** create a
`.BACKUP_`/`.BAK_` sibling in a directory that auto-discovery scans (see A-7).

```r
# DO NOT RETYPE — COPY from house_style.R block A10
house_save(p, "Figure_8", out_dir, width_mm = 280, height_mm = 220)
```

---

## C10 — Figure caption files

`[INV]` A plain-text file sits next to every figure, containing, in order:
figure number and title; N per cell and what a unit is; panel key (A/B/C…);
what the metric means and why it is parameterised that way; the model term
and df convention; what the error bars represent; the CLD note (spliced from
`CLD_CAPTION_NOTE`, identical across every figure by construction); **any
robustness caveat that qualifies the panel** (e.g. a collinearity
re-expression note); and any stated omission (missing CLD, non-convergence).

`[DR]` A figure with no CLD → the CLD paragraph is **replaced** by an
explicit statement of why, never deleted outright.

---

## C11 — Naming

`[INV]` Run directories: `<STAGE>_YYYYMMDD_HHMMSS`. **One name per outcome,
everywhere** — the indicator label equals the directory name equals the CSV
prefix equals the catalogue string.

`[DR]` An outcome is renamed → it is renamed in the catalogue, the column
map, the legend descriptions, the manifest, **and** the manuscript table, in
the same commit.

**ADJUDICATION A-13.** Five sources currently disagree about outcome
definitions (manuscript Table 1, `export_outcome_catalogue.R`,
`legend_descriptions.R`, `indicator_column_map.R`, `manifest_indicators.R`) —
e.g. `flux_timepoint` is described in the catalogue as "Count of zone
crossings per session; GLMM (Poisson/NB)" when it is actually a continuous
flux fit as a Gaussian LMM on `log1p`. There is no single definition of
record today; this rule is what fixes that going forward.

---

## C12 — Sensitivity outputs

`[INV]` Every robustness check named in [`HANDOFF.md §7`](HANDOFF.md) writes a
**named CSV in the run directory**, whether or not it changed anything — a
check with no artifact did not happen. Fixed filenames, so a reader can look
for them by name without hunting:

| File | Check |
|---|---|
| `anova_satterthwaite.csv` | Satterthwaite df sensitivity alongside every KR fit |
| `pb_lrt_<term>.csv` | parametric-bootstrap LRT, per beta-GLMM treatment term |
| `aicc_selection.csv` | RE-candidate selection, with `singular` / `eligible` columns |
| `sensitivity_forced_re.csv` | forced-RE sensitivity |
| `sensitivity_loo.csv` | leave-one-unit-out |
| `multiverse_<endpoint>.csv` | specification-curve sweep for the primary endpoint |
| `equivalence_tost.csv` | TOST equivalence against the declared SESOI |
| `pseudoreplication_check.csv` | aggregation-level agreement |
| `metric_collinearity.csv` | collinearity gate vs the primary endpoint |
| `design_balance.csv` | design-confound screen |
| `stat_type_validation.csv` | the reported string matches the test actually run |

`[INV]` Each stochastic check's file records the **seed and N** it ran with.

---

## Appendix A — `house_style.R`

The full file lives at [`house_style.R`](house_style.R). It is verified to
run cleanly under R 4.6.0 with `ggplot2`, `ggtext` and `effectsize` installed
— `preflight_check()` passes, and `anova_caption()` reproduces the
manuscript's own worked examples (`F1, 8 = 22.58, p = 0.001, η²p = 0.738`;
`χ²1 = 235.7, p < 0.001` with no η²p) byte-for-byte. Source it; do not retype
any block from it.

---

## Appendix B — the forbidden list

Negative rules are cheaper to check than positive ones.

- No box plots outside `diagnostics/` (C4, A-6).
- No CLD letter truncation (C5, A-8).
- No `TREATMENT_COLORS` / `ZONE_COLORS_MAIN` copied from
  `figures_step5_module.R` (C1, A-3).
- No `CLD_SIZE <- 4.0` (C5, A-4).
- No `theme_legacy` in new work (C8, A-5).
- No `p = 0.000`, no `p < 0.0001` floor, no significance stars in prose or
  captions (C7, A-9).
- No η²p attached to a Wald χ² caption (C6).
- No bootstrap p printed finer than `1/(N+1)` (C7, C12).
- No `.BACKUP_`/`.BAK_` sibling inside a directory that auto-discovery scans
  (C9, A-7).
- No EMMs plotted as the raw data layer (C4).
- No second output tree per run (C9).

---

## Appendix C — pre-flight checklist

Before any figure set is called final:

1. Palette matches `PALETTE_TREATMENT` / `PALETTE_ZONE` in `house_style.R` exactly.
2. No `TREATMENT_COLORS` or `ZONE_COLORS_MAIN` object in scope.
3. Theme is `BASE_THEME`, sourced, not retyped.
4. No `geom_boxplot` outside `diagnostics/`.
5. Every CLD frame has gone through `cld_pipeline()`.
6. No `.group` string shortened relative to its source CSV.
7. Caption order is statistic → p → effect size, every panel.
8. No caption contains both `chi` and `eta`.
9. Every F-test caption carries an η²p.
10. No rendered number is `0.000`; no `p = 0.000`.
11. Every bootstrap p obeys the `1/(N+1)` floor; every PB-LRT row records N and seed.
12. Every model directory carries its full C12 sensitivity-artifact set, or a stated reason for each absence.
13. `stat_type` is read before `df_denom` anywhere the latter is consumed.
14. The discovered run directory has a completed-sentinel and does not match `BACKUP|BAK|PRE_|_old`.
15. Every saved figure has both `.png` and `.pdf`, 300 dpi, mm dimensions.
16. `fmt_F(30) == "30"`, `fmt3(0.5) == "0.50"`, `fmt3(0.77) == "0.77"` — run `preflight_check()`.
17. Every panel's N counts the unit of replication, not the fish.
18. No plotted point is an EMM.
19. Legend appears once per figure, below the bottom row.
20. Every figure has a sibling caption `.txt` file with the CLD note spliced in verbatim.
21. Any collinear metric (|r| ≥ 0.70 vs the primary endpoint) is captioned as a re-expression, not independent confirmation.
22. Filenames encode figure identity, not a build date.
23. No path in the build script is hardcoded outside a single config block.
24. Superseded conventions are not being read from an un-adjudicated source file.
