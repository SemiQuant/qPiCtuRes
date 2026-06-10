# qPiCtuRes

[![PyPI - Version](https://img.shields.io/pypi/v/qpictures.svg)](https://pypi.org/project/qpictures/)
[![Python](https://img.shields.io/pypi/pyversions/qpictures.svg)](https://pypi.org/project/qpictures/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Interactive plots and analyses for qPCR data exported from Thermo Fisher
QuantStudio Design & Analysis.

The first module turns a *Well Results* CSV into an interactive Plotly bar plot
of **relative quantities** (RQ = 2<sup>−ΔCq</sup>) against a chosen reference
sample, with log-10 y-axis and asymmetric error bars propagated from the Cq
SDs (Livak & Schmittgen, 2008).

## Install

### From PyPI (end users)

```bash
pip install qpictures
# with optional static-image export (PDF / PNG / SVG via kaleido)
pip install "qpictures[export]"
```

### From source (developers)

The repo ships with an `environment.yml` (Python 3.11) for `micromamba`:

```bash
micromamba env create -f environment.yml
micromamba activate qpictures
# the env file already runs `pip install -e .`; if you skipped that:
pip install -e ".[dev]"
```

## CLI

```text
qpictures inspect    <input.csv>                      # per-sample Cq summary
qpictures relative   <input.csv> -r BB_PBS [-t atpE]  # RQ vs one reference sample
qpictures paired     <input.csv> -b BB    [-t atpE]   # per-condition paired RQ
qpictures standcurve <input.csv>          [-t atpE]   # absolute Q from a std curve
```

### `relative` — single reference sample

```bash
qpictures relative data/ex1/WellResult-2026-06-09-144852.csv \
  --reference BB_PBS \
  --out results/ex1_relative.html \
  --csv-out results/ex1_relative.csv
```

Useful flags:

| Flag | Default | What it does |
| --- | --- | --- |
| `-r, --reference` | *required* | Reference sample (case-insensitive). |
| `-t, --target` | auto | Which assay/target to use (auto if only one). |
| `-o, --out` | `<input>.relative.html` | Output Plotly HTML. |
| `--csv-out` | *(off)* | Also write the RQ table as CSV. |
| `--pdf [PATH]` | *(off)* | Also write a PDF (auto-named from `--out` if no path). |
| `--png [PATH]` | *(off)* | Also write a PNG. |
| `--svg [PATH]` | *(off)* | Also write an SVG. |
| `--test` | `welch` | Per-sample significance test vs reference: `welch`, `student`, `mannwhitney`, `none`. |
| `--padjust` | `none` | Multiple-testing correction: `none`, `bh` (FDR), `holm`, `bonferroni`. |
| `--no-signif` | off | Hide significance stars on the plot. |
| `--background-rem [NAME]` | off | Use a named control sample as the background threshold. As a bare flag it looks up a sample called `control` (case-insensitive); pass an explicit name to override, e.g. `--background-rem PBS`. Any sample whose mean Cq exceeds the control's mean Cq is flagged with a `below_background` column and rendered as a dimmed grey bar. The control is **independent of `--reference`**, so you can normalize against one sample and threshold against another. |
| `--task TEXT` (repeatable) | `Unknown` | Which `Task` values to keep. |
| `--all-tasks` | off | Keep every `Task` (NTC / standards / etc.). |
| `--no-log` | log10 on | Linear y-axis instead of log10. |
| `--no-open` | opens | Do not auto-open the resulting HTML. |

The bars show the significance star (e.g. `**`) above the RQ value when a
test is run. See **Significance testing** below for the star thresholds and
caveats — the same rules apply in both `relative` and `paired` modes.

### `paired` — per-condition paired comparisons

For sample names of the form `<base><sep><condition>` (e.g. `BB_pbs`,
`TL_culture_swab`), normalize every base within each condition against a
**reference base** (e.g. `BB`). Each condition gets a group on the x-axis,
bars within a group are colored by base, and the reference base sits at
RQ = 1 (with its own SD propagated).

```bash
qpictures paired data/ex1/WellResult-2026-06-09-144852.csv \
  --reference-base BB \
  --out results/ex1_paired.html \
  --csv-out results/ex1_paired.csv
```

| Flag | Default | What it does |
| --- | --- | --- |
| `-b, --reference-base` | *required* | Reference base (case-insensitive), e.g. `BB`. |
| `-s, --separator` | `_` | Split `base<sep>condition` on this. |
| `-t, --target` | auto | Which assay/target to use. |
| `-o, --out` | `<input>.paired.html` | Output Plotly HTML. |
| `--csv-out` | *(off)* | Also write the paired RQ table as CSV. |
| `--pdf [PATH]` | *(off)* | Also write a PDF (auto-named from `--out` if no path). |
| `--png [PATH]` | *(off)* | Also write a PNG. |
| `--svg [PATH]` | *(off)* | Also write an SVG. |
| `--test` | `welch` | Per-pair significance test: `welch`, `student`, `mannwhitney`, `none`. |
| `--padjust` | `none` | Multiple-testing correction: `none`, `bh` (FDR), `holm`, `bonferroni`. |
| `--no-signif` | off | Hide significance stars on the plot. |
| `--background-rem [NAME]` | off | Use a named control sample as a single global background threshold. Bare flag defaults to a sample called `control`; pass an explicit name (e.g. `--background-rem PBS`) to override. Lookup happens *before* the separator split, so a control sample without `_` (e.g. `Control`) is still found. Any (base, condition) row with a mean Cq above the control's mean Cq is flagged in `below_background` and dimmed in the plot. |
| `--no-split` | off | Suppress splitting on `--separator`. Every sample is treated as its own base under one synthetic condition (`all`), and `--reference-base` must match a *full* sample name (e.g. `BB_pbs`). Equivalent to a single-reference comparison rendered through the paired layout. |
| `--task TEXT` (repeatable) | `Unknown` | Which `Task` values to keep. |
| `--all-tasks` | off | Keep every `Task`. |
| `--no-log` | log10 on | Linear y-axis. |
| `--no-open` | opens | Don't auto-open the HTML. |

Samples lacking the separator (`Input`, etc.) are skipped with a notice,
as are any conditions where the reference base is missing.

The paired plot adds three visual cues so the comparisons inside each group
are unambiguous:

- **Alternating background bands** behind each condition box each comparison
  unit as a single visual group.
- **Diagonal hatching** on the reference-base bars marks the baseline within
  every group at a glance.
- A small **`vs <ref_sample>` label under each x-tick** (e.g. `vs BB_pbs`)
  spells out the specific reference sample that group is being compared
  against.

The Python API exposes these as `group_shading`, `mark_reference_pattern`
and `show_per_group_reference` on `paired_relative_quantity_barplot` (all
default `True`).

#### Example: `--background-rem` and `--no-split`

```bash
# Bare flag: look up a sample called "control" and use its mean Cq as the threshold
qpictures paired data/ex1/WellResult-2026-06-09-14570.csv \
  --reference-base BB --background-rem

# Explicit control sample (case-insensitive)
qpictures paired data/ex1/WellResult-2026-06-09-14570.csv \
  --reference-base BB --background-rem PBS

# Don't split: use BB_pbs as a single reference, every other sample compared to it
qpictures paired data/ex1/WellResult-2026-06-09-14570.csv \
  --reference-base BB_pbs --no-split --background-rem Control
```

#### Significance testing

For each non-reference (base, condition) pair we run a two-sided test on the
**raw `Cq` values** (Cq is already log-scaled, so a t-test on Cq is the
conventional choice). Multiple-testing correction across all non-reference
pairs is **off by default**; opt in with `--padjust bh` (FDR), `holm`, or
`bonferroni`:

| Threshold (adjusted p) | Star |
| --- | --- |
| p < 1e-4 | `****` |
| p < 1e-3 | `***`  |
| p < 1e-2 | `**`   |
| p < 0.05 | `*`    |
| else     | `ns`   |

Output table columns: `p_value`, `p_adj`, `signif`. The star is drawn above
the RQ value inside each non-reference bar.

**Caveat:** typical qPCR designs have very few replicates (often `n = 3–4`),
so these p-values have low statistical power. Treat them as a guide, not a
verdict. With `n=2`, no test can run and `p_value` will be `NaN`.

### `standcurve` — absolute quantification from a standard curve

For runs that include a serial dilution of a known quantity (e.g. `DNASD_4`,
`DNASD_6`, … with `Task == Standard` and a `Quantity` column), fit

```text
Cq = slope · log10(Quantity) + intercept
```

and back-calculate the **absolute Quantity** for every `Unknown` sample.

```bash
qpictures standcurve data/sc_well_table_results-2026-06-10-092019.csv \
  --out results/std.html \
  --csv-out results/std.csv \
  --curve-csv results/std_curve.csv \
  --unit "copies/µL"
```

The output HTML stacks two interactive panels:

1. **Standard curve**: scatter of the standard replicates plus the OLS fit,
   a ± 1 residual-SD band, and an annotation box summarising
   `slope`, `intercept`, `R²`, `efficiency = 10^(-1/slope) - 1`,
   `slope SE`, `residual SD`, `n` and number of levels. Unknown samples are
   overlaid at their `(log10 Q_hat, Cq_mean)` so you can see whether each one
   falls inside the calibrated range.
2. **Absolute-quantity bar plot**: one bar per sample, sorted by quantity,
   with asymmetric error bars from `SD(Cq) / |slope|` propagated through the
   curve. Bars whose mean Cq lies *outside* the standard range are recolored
   (`extrapolated == True`); dotted horizontal guide lines mark the
   calibrated min/max.

| Flag | Default | What it does |
| --- | --- | --- |
| `-t, --target` | auto | Which assay/target to fit (auto if only one is present). |
| `-o, --out` | `<input>.standcurve.html` | Combined interactive HTML (curve + bars). |
| `--curve-html FILE` | *(off)* | Also write a standalone HTML containing only the standard curve. |
| `--bars-html FILE` | *(off)* | Also write a standalone HTML containing only the bar plot. |
| `--csv-out FILE` | *(off)* | Per-sample absolute-quantity table as CSV. |
| `--curve-csv FILE` | *(off)* | Per-replicate standard wells used for the fit. |
| `--pdf [PATH]` | *(off)* | Also write the combined figure to PDF. |
| `--png [PATH]` | *(off)* | Also write PNG. |
| `--svg [PATH]` | *(off)* | Also write SVG. |
| `--unit TEXT` | *(none)* | Unit string appended to the y-axis title (e.g. `copies/µL`). |
| `--standard-task TEXT` | `Standard` | `Task` value identifying the standards. |
| `--unknown-task TEXT` | `Unknown` | `Task` value identifying the samples. |
| `--background-rem [NAME]` | *(off)* | Use a named control sample (e.g. `Water`) as the background threshold; samples with mean Cq above it are flagged in `below_background` and dimmed. |
| `--task TEXT` (repeatable) | *(none)* | Extra `Task` values to keep when reading the CSV (in addition to the standard and unknown tasks, which are always kept). |
| `--no-log` | log10 on | Linear y-axis for the bar plot. |
| `--no-open` | opens | Do not auto-open the HTML. |

**Output table columns** (`--csv-out`):

`Target, Sample, n, Cq_mean, Cq_sd, log10_Q, log10_Q_sd, Quantity,
Quantity_low, Quantity_high, Q_err_minus, Q_err_plus, extrapolated,
Background_control, below_background`.

**Quality cues to watch for**:

- `Efficiency` should sit around 90–110 % (slope ≈ −3.32 for 100 %). The
  reported value reflects the *slope*, not the absolute accuracy of any
  single quantity.
- `R²` should be ≥ 0.98 for a well-behaved dilution series.
- Any row with `extrapolated == True` had a mean Cq outside the calibrated
  log10-Quantity range — its back-calculated Quantity is an extrapolation
  and should be interpreted with care (the bar is recolored).

### Filtering defaults

`read_thermo_well_results()` (and every CLI command) by default keeps **only
`Omit == false`** wells and **only `Task == "Unknown"`** rows, so NTC /
standard / endogenous-control wells are excluded automatically. Override
with `--all-tasks` or `--task <value>` (repeatable) on the CLI, or
`tasks=None` / `tasks=("Unknown", "Standard")` from Python.

`qpictures standcurve` automatically keeps both `Standard` and `Unknown`
wells (configurable via `--standard-task` / `--unknown-task`); use
`--task <value>` to also keep additional task categories.

## Python API

```python
from qpictures import (
    read_thermo_well_results,
    relative_quantification,
    relative_quantity_barplot,
    paired_relative_quantification,
    paired_relative_quantity_barplot,
    fit_standard_curve,
    absolute_quantification,
    standard_curve_plot,
    absolute_quantity_barplot,
    standard_curve_with_quantities_figure,
)

df = read_thermo_well_results("data/ex1/WellResult-2026-06-09-144852.csv")

# 1) Single reference sample
rq, ref = relative_quantification(df, reference="BB_PBS")  # case-insensitive
relative_quantity_barplot(rq, log_y=True).write_html("rq.html")

# 2) Per-condition paired (e.g. TL_pbs vs BB_pbs, TL_sputum vs BB_sputum, ...)
#    Welch's t-test on Cq by default; multiple-testing correction is off
#    by default — pass padjust="bh" for Benjamini–Hochberg FDR.
pairs, summary = paired_relative_quantification(
    df, reference_base="BB", test="welch", padjust="none"
)
paired_relative_quantity_barplot(
    pairs, reference_base=summary.reference_base_resolved, log_y=True
).write_html("paired.html")

# 3) Absolute quantification from a standard curve. The CSV here includes
#    Task == Standard wells (DNASD_*) with known Quantity values plus the
#    usual Unknown wells; we must read both tasks to fit a curve.
df_sc = read_thermo_well_results(
    "data/sc_well_table_results-2026-06-10-092019.csv",
    tasks=("Unknown", "Standard"),
)
curve = fit_standard_curve(df_sc)                          # auto-picks target
print(f"slope={curve.slope:.3f}  R^2={curve.r_squared:.4f}"
      f"  efficiency={curve.efficiency*100:.1f}%")
abs_table, _ = absolute_quantification(df_sc, curve=curve)
standard_curve_with_quantities_figure(
    curve, abs_table, quantity_unit="copies/µL"
).write_html("standcurve.html")
```

## Math

For a single target (no endogenous control):

- `Cq_mean(sample)` and `Cq_sd(sample)` are computed across technical replicates.
- `ΔCq      = Cq_mean(sample) - Cq_mean(reference)`
- `SD(ΔCq)  = sqrt( SD(sample)^2 + SD(reference)^2 )`
- `RQ       = 2^(-ΔCq)`
- `RQ_high  = 2^( -(ΔCq - SD(ΔCq)) )`
- `RQ_low   = 2^( -(ΔCq + SD(ΔCq)) )`

The bar plot uses asymmetric error bars `[RQ - RQ_low, RQ_high - RQ]` so that
the visual error is symmetric on the log axis.

### Standard curve (absolute quantification)

For runs that include serial dilutions of a known input (`Task == Standard`,
`Quantity` known), OLS on individual replicates yields

- `Cq        = slope · log10(Q) + intercept`
- `efficiency = 10^(-1/slope) - 1` (1.0 == 100 %; slope ≈ −3.32 ideal)

For each unknown sample with mean Cq `Cq_mean` and SD `Cq_sd`:

- `log10(Q_hat)  = (Cq_mean - intercept) / slope`
- `Q_hat         = 10 ^ log10(Q_hat)`
- `SD(log10 Q)   = Cq_sd / |slope|`     (replicate scatter only)
- `Q_low         = 10 ^ (log10(Q_hat) - SD(log10 Q))`
- `Q_high        = 10 ^ (log10(Q_hat) + SD(log10 Q))`

Uncertainty in `slope` / `intercept` is *not* propagated into per-sample
bounds — they reflect technical-replicate scatter only. The curve's
`slope_se`, `intercept_se` and `residual_sd` are still surfaced as part of
the fit (and the ± 1 residual-SD band on the standard-curve plot).

## Notes on the Thermo export

The parser is tolerant of two common quirks:

1. A metadata header block above the actual table is skipped automatically.
2. Sample names containing **unquoted commas** (e.g. `beads,3mm,WS, Slow`) are
   detected (row has more fields than the header) and re-glued into the
   `Sample` column.
