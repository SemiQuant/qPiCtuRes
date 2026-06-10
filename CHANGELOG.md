# Changelog

All notable changes to **qPiCtuRes** are documented here. This project
follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html) and the
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format.

## [Unreleased]

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

[Unreleased]: https://github.com/SemiQuant/qPiCtuRes/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/SemiQuant/qPiCtuRes/releases/tag/v0.1.0
