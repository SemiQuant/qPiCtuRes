# Changelog

All notable changes to **qPiCtuRes** are documented here. This project
follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html) and the
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format.

## [Unreleased]

### Added

- Bio-Rad CFX Maestro exports (CFX Opus, CFX96, ...) are accepted wherever a
  Well Results CSV is: the *Quantification Cq Results* CSV, the folder from
  *Export All Data Sheets*, or a zip of it. Handles semicolon-separated files
  with decimal commas. `Content` codes map to the usual task labels and
  standard quantities come from the `SQ` column. New `read_bio_rad_cfx` in the
  Python API. A `.pcrd` file gets an error explaining how to export.

## [0.1.3] - 2026-11-24

### Added

- QuantStudio Design & Analysis `.eds` run files are accepted wherever a Well
  Results CSV is (every CLI command and `read_thermo_well_results`). Sample
  names, tasks and standard quantities come from the plate setup, Cq values
  from the analysis stored in the file. New `read_thermo_eds` in the Python
  API. A file with no stored analysis raises an error asking for it to be
  analysed and saved; an out-of-date analysis raises a warning.
- An offline web app: `qpictures.html`, attached to each GitHub release, runs
  the relative, paired and standard-curve analyses in the browser from a CSV or
  `.eds` file, with a plate map to leave wells out. No install, no network.

### Fixed

- `.eds` files from 48-well blocks get 6 x 8 well names (they were laid out as
  96-well plates).

## [0.1.2] - 2026-07-10

### Added

- `qpictures paired --base-order` and `--condition-order` to control bar /
  legend ordering and x-axis condition groups (comma-separated and/or repeated
  flags; case-insensitive).
- Biological-replicate count `N` annotated on every bar plot (`relative`,
  `paired`, `standcurve`). qPCR is assumed to run in technical triplicate, so
  `N = ceil(n_wells / tech_reps)` (default `tech_reps = 3`); a biological
  replicate still counts as long as not all of its technical wells were
  omitted. Configurable via `--tech-reps` and suppressible with `--no-n`
  (Python API: `show_n` / `tech_replicates` on the bar-plot builders). `N` is
  also shown alongside the raw well count `n` in each bar's hover tooltip.

### Changed

- CI publish workflow verifies that a release tag matches `pyproject.toml`
  version before building.

## [0.1.1] - 2026-06-16


## [0.1.0] - 2026-06-10

Initial public release.

### Added

- `qpictures inspect` — quick per-sample Cq summary of a Thermo QuantStudio
  Design & Analysis "Well Results" CSV.
- `qpictures relative` — relative quantification (Livak ΔCq) against a
  single reference sample, with asymmetric error bars on a log y-axis,
  optional Welch / Student / Mann–Whitney significance testing, optional
  multiple-testing correction (BH / Holm / Bonferroni), and an optional
  background-control threshold.
- `qpictures paired` — per-condition paired RQ plot. Splits sample names
  on a separator (e.g. `BB_pbs` → base=`BB`, condition=`pbs`) and
  normalizes every base in each condition against a reference base in the
  same condition. Adds visual grouping cues (alternating shading,
  reference-bar hatching, per-group `vs <ref_sample>` x-tick labels).
- `qpictures standcurve` — absolute quantification via an OLS standard
  curve (`Cq = slope · log10(Q) + intercept`). Reports slope, intercept,
  R², PCR efficiency, residual SD; back-calculates per-sample quantities
  with error bars propagated from Cq SD; flags extrapolations outside the
  calibrated standard range. Combined HTML stacks the curve diagnostic
  with the per-sample bar plot.
- Python API: `read_thermo_well_results`, `relative_quantification`,
  `paired_relative_quantification`, `fit_standard_curve`,
  `absolute_quantification`, plus matching plot builders
  (`relative_quantity_barplot`, `paired_relative_quantity_barplot`,
  `standard_curve_plot`, `absolute_quantity_barplot`,
  `standard_curve_with_quantities_figure`).
- Optional static-image export (PDF / PNG / SVG) for every plotting
  command via Plotly + `kaleido` (installed by `pip install qpictures[export]`).
- Tolerant CSV parser that auto-skips the Design & Analysis metadata
  preamble and re-glues sample names containing unquoted commas.

[Unreleased]: https://github.com/SemiQuant/qPiCtuRes/compare/v0.1.3...HEAD
[0.1.3]: https://github.com/SemiQuant/qPiCtuRes/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/SemiQuant/qPiCtuRes/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/SemiQuant/qPiCtuRes/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/SemiQuant/qPiCtuRes/releases/tag/v0.1.0
