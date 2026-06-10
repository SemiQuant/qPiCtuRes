"""Relative quantification (Livak / ΔCq) for qPCR well-result tables.

This module covers three workflows:

- :func:`relative_quantification` — single reference sample (Livak ΔCq).
- :func:`paired_relative_quantification` — per-condition paired comparisons.
- :func:`fit_standard_curve` + :func:`absolute_quantification` — absolute
  quantification from a Cq vs log10(Quantity) standard curve.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import pandas as pd
from scipy import stats as _sps

# Supported significance tests (per base, within each condition).
TestName = Literal["welch", "student", "mannwhitney", "none"]
# Supported multiple-testing correction methods.
AdjustName = Literal["bh", "fdr_bh", "bonferroni", "holm", "none"]


def _pairwise_pvalue(
    sample_cq: np.ndarray, ref_cq: np.ndarray, test: TestName
) -> float:
    """Two-sided p-value for the given test.

    qPCR Cq values are already on a log scale, so a *t*-test on Cq is the
    conventional choice. Returns ``NaN`` when the test can't run (e.g. a
    group has < 2 finite replicates).
    """
    if test in (None, "none"):
        return float("nan")
    s = np.asarray(sample_cq, dtype=float)
    r = np.asarray(ref_cq, dtype=float)
    s = s[np.isfinite(s)]
    r = r[np.isfinite(r)]
    if len(s) < 2 or len(r) < 2:
        return float("nan")
    if test == "welch":
        result = _sps.ttest_ind(s, r, equal_var=False)
    elif test == "student":
        result = _sps.ttest_ind(s, r, equal_var=True)
    elif test == "mannwhitney":
        result = _sps.mannwhitneyu(s, r, alternative="two-sided")
    else:
        raise ValueError(f"Unknown test {test!r}")
    return float(result.pvalue)


def _adjust_pvalues(pvals: np.ndarray, method: AdjustName) -> np.ndarray:
    """Multiple-testing correction. NaN p-values pass through unchanged."""
    p = np.asarray(pvals, dtype=float)
    if method in (None, "none"):
        return p.copy()

    mask = np.isfinite(p)
    if not mask.any():
        return p.copy()

    valid = p[mask]
    m = len(valid)
    order = np.argsort(valid)

    if method == "bonferroni":
        adj_unsorted = np.clip(valid * m, 0, 1)
        adj = adj_unsorted
    elif method in ("bh", "fdr_bh"):
        sorted_p = valid[order]
        ranks = np.arange(1, m + 1, dtype=float)
        # adjusted_(i) = min over k >= i of (m/k) * p_(k)
        raw = sorted_p * m / ranks
        adj_sorted = np.minimum.accumulate(raw[::-1])[::-1]
        adj_sorted = np.clip(adj_sorted, 0, 1)
        adj = np.empty_like(valid)
        adj[order] = adj_sorted
    elif method == "holm":
        sorted_p = valid[order]
        ranks = np.arange(1, m + 1, dtype=float)
        raw = sorted_p * (m - ranks + 1)
        adj_sorted = np.maximum.accumulate(raw)
        adj_sorted = np.clip(adj_sorted, 0, 1)
        adj = np.empty_like(valid)
        adj[order] = adj_sorted
    else:
        raise ValueError(f"Unknown padjust method {method!r}")

    out = np.full_like(p, np.nan)
    out[mask] = adj
    return out


def _significance_stars(pval: float) -> str:
    """Map a (possibly adjusted) p-value to a conventional star string."""
    if pval is None or pd.isna(pval):
        return ""
    p = float(pval)
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 5e-2:
        return "*"
    return "ns"


@dataclass(frozen=True)
class ReferenceMatch:
    """Result of resolving a user-supplied reference name to the data."""

    requested: str
    resolved: str
    n_replicates: int


def _match_reference(
    sample_values: pd.Series, requested: str, *, kind: str = "Reference"
) -> ReferenceMatch:
    """Case-insensitive lookup of a sample name (reference or control)."""
    unique = sample_values.dropna().unique().tolist()
    lowered = {s.lower(): s for s in unique}
    key = requested.strip().lower()
    if key not in lowered:
        raise ValueError(
            f"{kind} sample {requested!r} not found. "
            f"Available samples: {sorted(unique)}"
        )
    resolved = lowered[key]
    n_rep = int((sample_values == resolved).sum())
    return ReferenceMatch(requested=requested, resolved=resolved, n_replicates=n_rep)


def _sample_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Per (Sample, Target) mean/SD/n of Cq."""
    grouped = (
        df.groupby(["Target", "Sample"], dropna=False)["Cq"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"mean": "Cq_mean", "std": "Cq_sd", "count": "n"})
    )
    # std of a single replicate is NaN; treat as 0 so error bars don't disappear.
    grouped["Cq_sd"] = grouped["Cq_sd"].fillna(0.0)
    return grouped


def relative_quantification(
    df: pd.DataFrame,
    reference: str,
    *,
    target: str | None = None,
    test: TestName = "welch",
    padjust: AdjustName = "none",
    background_control: str | None = None,
) -> tuple[pd.DataFrame, ReferenceMatch]:
    """Compute relative quantities (RQ = 2^-ΔCq) against a reference sample.

    Error propagation follows Livak & Schmittgen (2008):
        SD(ΔCq) = sqrt(SD_sample^2 + SD_reference^2)
        RQ_high = 2 ^ -(mean(ΔCq) - SD(ΔCq))
        RQ_low  = 2 ^ -(mean(ΔCq) + SD(ΔCq))

    Significance testing
    --------------------
    For each non-reference sample we run a two-sided test on the *raw Cq
    values* against the reference sample's replicates (Cq is already
    log-scaled, so a t-test on Cq is conventional). Supported ``test``
    values: ``"welch"`` (default), ``"student"``, ``"mannwhitney"``,
    ``"none"``. Multiple-testing correction (``padjust``) is off by
    default; choose from ``"bh"`` / ``"fdr_bh"``, ``"holm"``,
    ``"bonferroni"``.

    With small replicate counts (typical qPCR ``n=3-4``) these p-values
    have low statistical power; treat them as a guide.

    Parameters
    ----------
    df:
        Tidy DataFrame as returned by :func:`qpictures.io.read_thermo_well_results`.
        Must contain ``Sample``, ``Target`` and ``Cq``.
    reference:
        Name of the reference sample (case-insensitive).
    target:
        If given, restrict the analysis to a single target / assay. If omitted
        and only one target exists in ``df``, that target is used; otherwise an
        error is raised.
    test, padjust:
        See "Significance testing" above.
    background_control:
        Optional name of a *control* sample (case-insensitive lookup against
        all samples in ``df``, independent of the reference). When given,
        every sample with a mean Cq greater than the control's mean Cq is
        flagged in the ``below_background`` column. The resolved control name
        is also stored in ``Background_control``. When ``None`` (default), the
        ``below_background`` column is ``False`` for all rows.

    Returns
    -------
    (pd.DataFrame, ReferenceMatch)
        A per-sample table with mean Cq, ΔCq, RQ, asymmetric error bounds and
        (when ``test != "none"``) ``p_value``, ``p_adj`` and ``signif`` columns
        (reference row gets NaN / empty), plus information about the resolved
        reference.
    """
    required = {"Sample", "Target", "Cq"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Input is missing required columns: {sorted(missing)}")

    work = df.dropna(subset=["Cq"]).copy()
    if work.empty:
        raise ValueError("No usable Cq values in the input table.")

    targets = sorted(t for t in work["Target"].dropna().unique() if t != "")
    if target is None:
        if len(targets) != 1:
            raise ValueError(
                "Multiple targets present; pass `target=` to select one. "
                f"Targets: {targets}"
            )
        target = targets[0]
    elif target not in targets:
        raise ValueError(f"Target {target!r} not found. Available: {targets}")

    work = work.loc[work["Target"] == target].copy()
    ref_match = _match_reference(work["Sample"], reference)

    stats = _sample_stats(work)
    stats = stats.loc[stats["Target"] == target].copy()

    ref_row = stats.loc[stats["Sample"] == ref_match.resolved].iloc[0]
    ref_mean = float(ref_row["Cq_mean"])
    ref_sd = float(ref_row["Cq_sd"])

    stats["Target"] = target
    stats["Reference"] = ref_match.resolved
    stats["dCq"] = stats["Cq_mean"] - ref_mean
    stats["dCq_sd"] = np.sqrt(stats["Cq_sd"] ** 2 + ref_sd**2)
    stats["RQ"] = np.power(2.0, -stats["dCq"])
    stats["RQ_high"] = np.power(2.0, -(stats["dCq"] - stats["dCq_sd"]))
    stats["RQ_low"] = np.power(2.0, -(stats["dCq"] + stats["dCq_sd"]))
    stats["RQ_err_plus"] = stats["RQ_high"] - stats["RQ"]
    stats["RQ_err_minus"] = stats["RQ"] - stats["RQ_low"]
    stats["is_reference"] = stats["Sample"] == ref_match.resolved

    # Below-background: a sample's mean Cq strictly greater than the *control*
    # sample's mean Cq (i.e. less signal than the control). The control is a
    # separate sample from the reference; if omitted, no sample is flagged.
    if background_control is not None:
        ctrl_match = _match_reference(
            work["Sample"], background_control, kind="Background control"
        )
        ctrl_cq = work.loc[work["Sample"] == ctrl_match.resolved, "Cq"]
        if ctrl_cq.empty:
            raise ValueError(
                f"Background control sample {ctrl_match.resolved!r} has no "
                "Cq values."
            )
        ctrl_mean = float(ctrl_cq.mean())
        stats["Background_control"] = ctrl_match.resolved
        stats["below_background"] = (
            (stats["Cq_mean"] > ctrl_mean)
            & (stats["Sample"] != ctrl_match.resolved)
        )
    else:
        stats["Background_control"] = ""
        stats["below_background"] = False

    # --- per-sample significance testing (raw Cq vs reference) ---
    ref_cq_values = work.loc[work["Sample"] == ref_match.resolved, "Cq"].to_numpy()
    pvals: list[float] = []
    for _, row in stats.iterrows():
        if row["is_reference"] or test in (None, "none"):
            pvals.append(float("nan"))
            continue
        sample_cq = work.loc[work["Sample"] == row["Sample"], "Cq"].to_numpy()
        pvals.append(_pairwise_pvalue(sample_cq, ref_cq_values, test))
    stats["p_value"] = pvals
    stats["p_adj"] = _adjust_pvalues(np.asarray(pvals, dtype=float), padjust)
    stats["signif"] = [
        "" if bool(row["is_reference"]) else _significance_stars(row["p_adj"])
        for _, row in stats.iterrows()
    ]

    ordered = stats[
        [
            "Target",
            "Sample",
            "n",
            "Cq_mean",
            "Cq_sd",
            "Reference",
            "dCq",
            "dCq_sd",
            "RQ",
            "RQ_low",
            "RQ_high",
            "RQ_err_minus",
            "RQ_err_plus",
            "p_value",
            "p_adj",
            "signif",
            "is_reference",
            "Background_control",
            "below_background",
        ]
    ].sort_values(["is_reference", "RQ"], ascending=[False, False])

    ordered.reset_index(drop=True, inplace=True)
    return ordered, ref_match


# --------------------------------------------------------------------------- #
# Paired / per-condition relative quantification
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class PairedSummary:
    """Bookkeeping for :func:`paired_relative_quantification`."""

    reference_base_requested: str
    reference_base_resolved: str
    separator: str
    bases: tuple[str, ...]
    conditions: tuple[str, ...]
    skipped_samples: tuple[str, ...] = field(default_factory=tuple)
    skipped_conditions: tuple[str, ...] = field(default_factory=tuple)
    test: str = "none"
    padjust: str = "none"


def _split_sample(sample: str, separator: str) -> tuple[str | None, str | None]:
    """Split ``sample`` into ``(base, condition)`` using the *first* separator.

    Returns ``(None, None)`` when the sample doesn't contain ``separator``.
    """
    if sample is None or pd.isna(sample) or separator not in sample:
        return None, None
    base, condition = sample.split(separator, 1)
    return base.strip() or None, condition.strip() or None


def _match_base(bases: list[str], requested: str) -> str:
    """Case-insensitive lookup of a base name."""
    lowered = {b.lower(): b for b in bases}
    key = requested.strip().lower()
    if key not in lowered:
        raise ValueError(
            f"Reference base {requested!r} not found. "
            f"Available bases: {sorted(bases)}"
        )
    return lowered[key]


def paired_relative_quantification(
    df: pd.DataFrame,
    reference_base: str,
    *,
    separator: str = "_",
    target: str | None = None,
    test: TestName = "welch",
    padjust: AdjustName = "none",
    split_samples: bool = True,
    background_control: str | None = None,
) -> tuple[pd.DataFrame, PairedSummary]:
    """Within each condition, normalize each base against a reference base.

    Sample names are split on the *first* occurrence of ``separator``:

        ``BB_pbs``           -> base=``BB``, condition=``pbs``
        ``TL_culture_swab``  -> base=``TL``, condition=``culture_swab``
        ``Input``            -> skipped (no separator)

    For every condition where the reference base is present, the RQ of each
    non-reference base is computed against the reference base in the *same*
    condition (Livak ΔCq with SD propagation, same as
    :func:`relative_quantification`).

    Significance testing
    --------------------
    For each non-reference (Base, Condition) pair we run a two-sided test on
    the *raw Cq values* (Cq is already log-scaled, so a t-test on Cq is
    conventional). Available ``test`` values:

    - ``"welch"`` (default): Welch's t-test (unequal variances).
    - ``"student"``: Student's t-test (equal variances).
    - ``"mannwhitney"``: Mann–Whitney U (non-parametric).
    - ``"none"``: skip; p-values will be NaN.

    The resulting p-values can be corrected for multiple comparisons across
    all non-reference (Base, Condition) pairs using ``padjust``:

    - ``"none"`` (default): no correction; ``p_adj`` equals ``p_value``.
    - ``"bh"`` / ``"fdr_bh"``: Benjamini–Hochberg FDR.
    - ``"holm"``: Holm step-down.
    - ``"bonferroni"``.

    Caveat: typical qPCR designs have very few replicates (often n = 3–4),
    so these p-values have low statistical power. Treat them as a guide,
    not a definitive verdict.

    Parameters
    ----------
    df:
        Tidy DataFrame as returned by :func:`qpictures.io.read_thermo_well_results`.
    reference_base:
        Base prefix to use as the per-condition reference (case-insensitive),
        e.g. ``"BB"``.
    separator:
        Sample-name separator. Default ``"_"``.
    target:
        Restrict to a single target/assay. Auto-detected if only one is present.
    test, padjust:
        See "Significance testing" above.
    split_samples:
        If ``True`` (default), split each sample name on ``separator`` into
        ``(base, condition)``. If ``False``, suppress splitting: every sample
        is treated as its own base and all samples share a single synthetic
        condition (``"all"``). In that mode ``reference_base`` must match a
        full sample name (e.g. ``"BB_pbs"``), and every other sample is
        compared directly to that one reference — equivalent to a single-
        reference comparison but rendered with the paired layout.
    background_control:
        Optional name of a *control* sample, looked up case-insensitively
        against the full target-filtered input (the lookup is done before
        the separator split, so a control like ``"Control"`` without a
        separator is still found). Every (base, condition) row whose mean Cq
        is greater than the control's mean Cq is flagged in
        ``below_background``. When ``None`` (default), the column is
        ``False`` for all rows.

    Returns
    -------
    (pd.DataFrame, PairedSummary)
        Per-(base, condition) RQ table plus a summary of bases/conditions and
        anything that was skipped. The table additionally contains
        ``p_value``, ``p_adj`` and ``signif`` columns (NaN/empty for the
        reference base rows).
    """
    required = {"Sample", "Target", "Cq"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Input is missing required columns: {sorted(missing)}")

    work = df.dropna(subset=["Cq"]).copy()
    if work.empty:
        raise ValueError("No usable Cq values in the input table.")

    targets = sorted(t for t in work["Target"].dropna().unique() if t != "")
    if target is None:
        if len(targets) != 1:
            raise ValueError(
                "Multiple targets present; pass `target=` to select one. "
                f"Targets: {targets}"
            )
        target = targets[0]
    elif target not in targets:
        raise ValueError(f"Target {target!r} not found. Available: {targets}")

    work = work.loc[work["Target"] == target].copy()

    # Resolve the (optional) background control against all target-filtered
    # samples *before* splitting, so a no-separator control sample like
    # "Control" is still found.
    if background_control is not None:
        ctrl_match = _match_reference(
            work["Sample"], background_control, kind="Background control"
        )
        ctrl_cq = work.loc[work["Sample"] == ctrl_match.resolved, "Cq"]
        if ctrl_cq.empty:
            raise ValueError(
                f"Background control sample {ctrl_match.resolved!r} has no "
                "Cq values."
            )
        ctrl_mean: float | None = float(ctrl_cq.mean())
        ctrl_resolved: str | None = ctrl_match.resolved
    else:
        ctrl_mean = None
        ctrl_resolved = None

    if split_samples:
        split = work["Sample"].apply(lambda s: _split_sample(s, separator))
        work["Base"] = [b for b, _ in split]
        work["Condition"] = [c for _, c in split]

        splittable = work.dropna(subset=["Base", "Condition"]).copy()
        skipped = sorted(set(work.loc[work["Base"].isna(), "Sample"]))
        if splittable.empty:
            raise ValueError(
                f"No samples contain the separator {separator!r}. "
                "Cannot do a paired analysis."
            )
    else:
        # No splitting: every sample is its own "base", all share one
        # synthetic condition. This is effectively a single-reference
        # comparison rendered through the paired pipeline.
        work["Base"] = work["Sample"]
        work["Condition"] = "all"
        splittable = work.dropna(subset=["Base", "Condition"]).copy()
        skipped = []
        if splittable.empty:
            raise ValueError("No usable samples after filtering.")

    bases = sorted(splittable["Base"].unique().tolist())
    ref_base = _match_base(bases, reference_base)

    stats = (
        splittable.groupby(["Target", "Base", "Condition"], dropna=False)["Cq"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"mean": "Cq_mean", "std": "Cq_sd", "count": "n"})
    )
    stats["Cq_sd"] = stats["Cq_sd"].fillna(0.0)

    ref_rows = stats.loc[stats["Base"] == ref_base, ["Condition", "Cq_mean", "Cq_sd"]]
    ref_rows = ref_rows.rename(columns={"Cq_mean": "Cq_ref_mean", "Cq_sd": "Cq_ref_sd"})

    conditions_with_ref = set(ref_rows["Condition"].unique())
    all_conditions = set(stats["Condition"].unique())
    skipped_conditions = sorted(all_conditions - conditions_with_ref)

    merged = stats.merge(ref_rows, on="Condition", how="inner")
    if merged.empty:
        raise ValueError(
            f"Reference base {ref_base!r} is not paired with any condition "
            f"that another base also covers."
        )

    if split_samples:
        merged["Sample"] = merged["Base"] + separator + merged["Condition"]
        merged["Reference_sample"] = ref_base + separator + merged["Condition"]
    else:
        # No splitting: Sample == Base, and the reference sample is the
        # full reference base name for every row.
        merged["Sample"] = merged["Base"]
        merged["Reference_sample"] = ref_base
    merged["Target"] = target
    merged["Reference_base"] = ref_base
    merged["dCq"] = merged["Cq_mean"] - merged["Cq_ref_mean"]
    merged["dCq_sd"] = np.sqrt(merged["Cq_sd"] ** 2 + merged["Cq_ref_sd"] ** 2)
    merged["RQ"] = np.power(2.0, -merged["dCq"])
    merged["RQ_high"] = np.power(2.0, -(merged["dCq"] - merged["dCq_sd"]))
    merged["RQ_low"] = np.power(2.0, -(merged["dCq"] + merged["dCq_sd"]))
    merged["RQ_err_plus"] = merged["RQ_high"] - merged["RQ"]
    merged["RQ_err_minus"] = merged["RQ"] - merged["RQ_low"]
    merged["is_reference"] = merged["Base"] == ref_base
    # Below-background flag: this row's mean Cq strictly greater than the
    # *control sample's* mean Cq (a single global threshold, not per-
    # condition). The control sample itself is excluded.
    if ctrl_mean is not None:
        merged["Background_control"] = ctrl_resolved
        merged["below_background"] = (
            (merged["Cq_mean"] > ctrl_mean)
            & (merged["Sample"] != ctrl_resolved)
        )
    else:
        merged["Background_control"] = ""
        merged["below_background"] = False

    # --- per-pair significance testing (raw Cq values) ---
    pvals: list[float] = []
    for _, row in merged.iterrows():
        if row["Base"] == ref_base or test in (None, "none"):
            pvals.append(float("nan"))
            continue
        sample_cq = splittable.loc[
            (splittable["Base"] == row["Base"])
            & (splittable["Condition"] == row["Condition"]),
            "Cq",
        ].to_numpy()
        ref_cq = splittable.loc[
            (splittable["Base"] == ref_base)
            & (splittable["Condition"] == row["Condition"]),
            "Cq",
        ].to_numpy()
        pvals.append(_pairwise_pvalue(sample_cq, ref_cq, test))
    merged["p_value"] = pvals
    merged["p_adj"] = _adjust_pvalues(np.asarray(pvals, dtype=float), padjust)
    merged["signif"] = [
        "" if row["Base"] == ref_base else _significance_stars(row["p_adj"])
        for _, row in merged.iterrows()
    ]

    ordered = merged[
        [
            "Target",
            "Condition",
            "Base",
            "Sample",
            "n",
            "Cq_mean",
            "Cq_sd",
            "Reference_base",
            "Reference_sample",
            "dCq",
            "dCq_sd",
            "RQ",
            "RQ_low",
            "RQ_high",
            "RQ_err_minus",
            "RQ_err_plus",
            "p_value",
            "p_adj",
            "signif",
            "is_reference",
            "Background_control",
            "below_background",
        ]
    ].sort_values(
        # reference base first within each condition, then by Base for stable order
        ["Condition", "is_reference", "Base"],
        ascending=[True, False, True],
    )
    ordered.reset_index(drop=True, inplace=True)

    summary = PairedSummary(
        reference_base_requested=reference_base,
        reference_base_resolved=ref_base,
        separator=separator,
        bases=tuple(bases),
        conditions=tuple(sorted(conditions_with_ref)),
        skipped_samples=tuple(skipped),
        skipped_conditions=tuple(skipped_conditions),
        test=str(test) if test is not None else "none",
        padjust=str(padjust) if padjust is not None else "none",
    )
    return ordered, summary


# --------------------------------------------------------------------------- #
# Absolute quantification via a standard curve
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class StandardCurve:
    """Linear fit of ``Cq = slope * log10(Quantity) + intercept``.

    Attributes
    ----------
    target:
        Target / assay the curve was fitted on.
    slope, intercept:
        Coefficients of the line ``Cq = slope * log10(Q) + intercept`` (so the
        slope is the *standard-curve slope*, conventionally negative — about
        ``-3.32`` for 100 % efficiency).
    slope_se, intercept_se:
        Standard errors of ``slope`` and ``intercept`` from OLS.
    r_squared:
        Coefficient of determination of the linear fit.
    efficiency:
        PCR efficiency derived from the slope: ``10**(-1/slope) - 1``. A
        perfectly efficient assay returns ``1.0`` (100 %). ``NaN`` when the
        slope is non-finite or zero.
    n_points:
        Number of standard replicates that went into the fit.
    n_levels:
        Number of distinct standard concentrations used.
    log10_q_min, log10_q_max:
        Range of ``log10(Quantity)`` covered by the standards (useful for
        flagging unknowns extrapolating outside the curve).
    residual_sd:
        Residual standard error (``sigma``) of the Cq ~ log10(Q) fit, in
        Cq units.
    standards:
        Per-well rows used for the fit (``Sample, Cq, Quantity, log10_Q``).
    """

    target: str
    slope: float
    intercept: float
    slope_se: float
    intercept_se: float
    r_squared: float
    efficiency: float
    n_points: int
    n_levels: int
    log10_q_min: float
    log10_q_max: float
    residual_sd: float
    standards: pd.DataFrame = field(repr=False)


def _resolve_target(values: pd.Series, target: str | None) -> str:
    targets = sorted(t for t in values.dropna().unique() if t != "")
    if target is None:
        if len(targets) != 1:
            raise ValueError(
                "Multiple targets present; pass `target=` to select one. "
                f"Targets: {targets}"
            )
        return targets[0]
    if target not in targets:
        raise ValueError(f"Target {target!r} not found. Available: {targets}")
    return target


def fit_standard_curve(
    df: pd.DataFrame,
    *,
    target: str | None = None,
    standard_task: str = "Standard",
) -> StandardCurve:
    """Fit a qPCR standard curve on wells with ``Task == standard_task``.

    The fit is ``Cq = slope * log10(Quantity) + intercept`` (ordinary least
    squares on individual replicates, not group means).

    Parameters
    ----------
    df:
        Tidy DataFrame produced by :func:`qpictures.io.read_thermo_well_results`.
        Must contain ``Sample``, ``Target``, ``Task``, ``Cq`` and ``Quantity``.
        Note that the default reader filters to ``Task == "Unknown"``; for
        absolute-quantification workflows pass ``tasks=("Unknown", "Standard")``
        (or ``tasks=None``) when reading.
    target:
        Target / assay to fit. Auto-detected when only one target is present.
    standard_task:
        Value of the ``Task`` column that identifies standards. Default
        ``"Standard"`` (matches the Design & Analysis export).

    Returns
    -------
    StandardCurve
        See :class:`StandardCurve` for the field list.
    """
    required = {"Sample", "Target", "Task", "Cq", "Quantity"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            f"Input is missing required columns: {sorted(missing)}. "
            "When reading with `read_thermo_well_results`, pass "
            "`tasks=('Unknown', 'Standard')` (or `tasks=None`) so the "
            "Standard wells survive task filtering."
        )

    work = df.copy()
    work = work.loc[work["Task"].astype(str).str.strip().str.lower()
                    == standard_task.strip().lower()].copy()
    if work.empty:
        raise ValueError(
            f"No wells with Task == {standard_task!r} found. Did you forget to "
            "pass `tasks=('Unknown', 'Standard')` when reading the CSV?"
        )

    target = _resolve_target(work["Target"], target)
    work = work.loc[work["Target"] == target].copy()

    work["Cq"] = pd.to_numeric(work["Cq"], errors="coerce")
    work["Quantity"] = pd.to_numeric(work["Quantity"], errors="coerce")
    work = work.loc[work["Cq"].notna() & work["Quantity"].notna()].copy()
    work = work.loc[work["Quantity"] > 0].copy()
    if work.empty:
        raise ValueError(
            f"No usable Cq / Quantity values among standards for target {target!r}."
        )

    work["log10_Q"] = np.log10(work["Quantity"].astype(float))

    x = work["log10_Q"].to_numpy(dtype=float)
    y = work["Cq"].to_numpy(dtype=float)
    n = int(x.size)
    n_levels = int(np.unique(np.round(x, 8)).size)
    if n < 2 or n_levels < 2:
        raise ValueError(
            "Need at least two distinct standard concentrations to fit a curve "
            f"(got n_points={n}, n_levels={n_levels})."
        )

    lr = _sps.linregress(x, y)
    slope = float(lr.slope)
    intercept = float(lr.intercept)
    r_squared = float(lr.rvalue) ** 2
    slope_se = float(lr.stderr)
    intercept_se = float(lr.intercept_stderr)

    residuals = y - (slope * x + intercept)
    dof = max(n - 2, 1)
    residual_sd = float(np.sqrt(np.sum(residuals**2) / dof))

    if np.isfinite(slope) and slope != 0:
        efficiency = float(10.0 ** (-1.0 / slope) - 1.0)
    else:
        efficiency = float("nan")

    standards = (
        work[["Sample", "Cq", "Quantity", "log10_Q"]]
        .sort_values(["log10_Q", "Sample"])
        .reset_index(drop=True)
    )

    return StandardCurve(
        target=target,
        slope=slope,
        intercept=intercept,
        slope_se=slope_se,
        intercept_se=intercept_se,
        r_squared=r_squared,
        efficiency=efficiency,
        n_points=n,
        n_levels=n_levels,
        log10_q_min=float(np.min(x)),
        log10_q_max=float(np.max(x)),
        residual_sd=residual_sd,
        standards=standards,
    )


def _quantity_from_cq(
    cq: np.ndarray | float, curve: StandardCurve
) -> np.ndarray:
    """Invert the standard curve: ``Q = 10 ** ((Cq - intercept) / slope)``."""
    cq_arr = np.asarray(cq, dtype=float)
    if not np.isfinite(curve.slope) or curve.slope == 0:
        return np.full_like(cq_arr, np.nan, dtype=float)
    return np.power(10.0, (cq_arr - curve.intercept) / curve.slope)


def absolute_quantification(
    df: pd.DataFrame,
    *,
    target: str | None = None,
    curve: StandardCurve | None = None,
    standard_task: str = "Standard",
    unknown_task: str = "Unknown",
    background_control: str | None = None,
) -> tuple[pd.DataFrame, StandardCurve]:
    """Estimate per-sample absolute quantities from a standard curve.

    For every ``unknown_task`` sample we compute mean ± SD of Cq across
    replicates and back-calculate ``Quantity`` from the fitted curve:

        log10(Q_hat) = (Cq_mean - intercept) / slope
        Q_hat        = 10 ** log10(Q_hat)
        SD(log10(Q)) = Cq_sd / |slope|
        Q_low        = 10 ** (log10(Q_hat) - SD(log10(Q)))
        Q_high       = 10 ** (log10(Q_hat) + SD(log10(Q)))

    The curve's uncertainty (``slope_se`` etc.) is *not* propagated here — the
    bars reflect technical-replicate scatter only.

    Parameters
    ----------
    df:
        Tidy DataFrame. Must contain at least ``Sample``, ``Target``, ``Task``
        and ``Cq``; ``Quantity`` is read for the standards if ``curve`` is
        ``None``.
    target:
        Target / assay to analyze. Auto-detected if only one is present.
    curve:
        Pre-fit :class:`StandardCurve`. When ``None`` (default), a curve is
        fitted on ``Task == standard_task`` rows in ``df``.
    standard_task, unknown_task:
        ``Task`` column values identifying standards vs unknowns. Defaults
        ``"Standard"`` / ``"Unknown"``.
    background_control:
        Optional case-insensitive name of a control sample (typically a NTC
        / water well). Any unknown whose mean Cq strictly exceeds the
        control's mean Cq is flagged in ``below_background``. When ``None``
        the column is ``False`` for every row.

    Returns
    -------
    (pd.DataFrame, StandardCurve)
        Per-sample table with columns: ``Target, Sample, n, Cq_mean, Cq_sd,
        log10_Q, log10_Q_sd, Quantity, Quantity_low, Quantity_high,
        Q_err_minus, Q_err_plus, extrapolated, Background_control,
        below_background`` — plus the curve that was used (fitted internally
        when ``curve=None``).
    """
    if curve is None:
        curve = fit_standard_curve(df, target=target, standard_task=standard_task)
        target = curve.target
    else:
        target = curve.target if target is None else target

    required = {"Sample", "Target", "Task", "Cq"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Input is missing required columns: {sorted(missing)}")

    work = df.copy()
    work["Cq"] = pd.to_numeric(work["Cq"], errors="coerce")
    work = work.loc[work["Target"] == target].copy()

    unknown_mask = (
        work["Task"].astype(str).str.strip().str.lower()
        == unknown_task.strip().lower()
    )
    unknowns = work.loc[unknown_mask & work["Cq"].notna()].copy()
    unknowns = unknowns.loc[unknowns["Sample"].astype(str).str.strip() != ""].copy()
    if unknowns.empty:
        raise ValueError(
            f"No usable Cq values found for Task == {unknown_task!r} and "
            f"Target == {target!r}."
        )

    grouped = (
        unknowns.groupby("Sample", dropna=False)["Cq"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"mean": "Cq_mean", "std": "Cq_sd", "count": "n"})
    )
    grouped["Cq_sd"] = grouped["Cq_sd"].fillna(0.0)

    abs_slope = abs(curve.slope) if np.isfinite(curve.slope) and curve.slope != 0 else float("nan")
    grouped["log10_Q"] = (grouped["Cq_mean"] - curve.intercept) / curve.slope
    grouped["log10_Q_sd"] = grouped["Cq_sd"] / abs_slope
    grouped["Quantity"] = np.power(10.0, grouped["log10_Q"])
    grouped["Quantity_low"] = np.power(10.0, grouped["log10_Q"] - grouped["log10_Q_sd"])
    grouped["Quantity_high"] = np.power(10.0, grouped["log10_Q"] + grouped["log10_Q_sd"])
    grouped["Q_err_minus"] = grouped["Quantity"] - grouped["Quantity_low"]
    grouped["Q_err_plus"] = grouped["Quantity_high"] - grouped["Quantity"]
    grouped["extrapolated"] = (
        (grouped["log10_Q"] < curve.log10_q_min) | (grouped["log10_Q"] > curve.log10_q_max)
    )
    grouped["Target"] = target

    if background_control is not None:
        ctrl_match = _match_reference(
            work["Sample"], background_control, kind="Background control"
        )
        ctrl_cq = work.loc[work["Sample"] == ctrl_match.resolved, "Cq"]
        if ctrl_cq.dropna().empty:
            raise ValueError(
                f"Background control sample {ctrl_match.resolved!r} has no "
                "finite Cq values (likely an undetermined NTC). Use a control "
                "with at least one detected Cq, or omit `background_control`."
            )
        ctrl_mean = float(ctrl_cq.dropna().mean())
        grouped["Background_control"] = ctrl_match.resolved
        grouped["below_background"] = (
            (grouped["Cq_mean"] > ctrl_mean)
            & (grouped["Sample"] != ctrl_match.resolved)
        )
    else:
        grouped["Background_control"] = ""
        grouped["below_background"] = False

    ordered = grouped[
        [
            "Target",
            "Sample",
            "n",
            "Cq_mean",
            "Cq_sd",
            "log10_Q",
            "log10_Q_sd",
            "Quantity",
            "Quantity_low",
            "Quantity_high",
            "Q_err_minus",
            "Q_err_plus",
            "extrapolated",
            "Background_control",
            "below_background",
        ]
    ].sort_values("Quantity", ascending=False)
    ordered.reset_index(drop=True, inplace=True)
    return ordered, curve
