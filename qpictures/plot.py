"""Interactive Plotly visualizations for qPCR relative quantification."""

from __future__ import annotations

from typing import TYPE_CHECKING, Iterable, Mapping

import numpy as np
import pandas as pd
import plotly.colors as pc
import plotly.graph_objects as go
from plotly.subplots import make_subplots

if TYPE_CHECKING:
    from .analysis import StandardCurve

# Default qualitative palette used to color "base" groups.
_DEFAULT_PALETTE: tuple[str, ...] = tuple(pc.qualitative.Set2)


def _build_color_map(
    bases: list[str],
    reference_base: str | None,
    *,
    reference_color: str,
    overrides: Mapping[str, str] | None = None,
) -> dict[str, str]:
    """Assign a stable color per base, with the reference base in grey.

    Non-reference bases get colors from a qualitative palette in the order they
    appear in ``bases``. ``overrides`` (if given) wins over both.
    """
    mapping: dict[str, str] = {}
    palette = list(_DEFAULT_PALETTE)
    palette_idx = 0
    for base in bases:
        if reference_base is not None and base == reference_base:
            mapping[base] = reference_color
        else:
            mapping[base] = palette[palette_idx % len(palette)]
            palette_idx += 1
    if overrides:
        mapping.update(overrides)
    return mapping


def _format_hover(row: pd.Series) -> str:
    stats_line = ""
    if "p_value" in row.index and pd.notna(row.get("p_value")):
        p = float(row["p_value"])
        if "p_adj" in row.index and pd.notna(row.get("p_adj")):
            p_adj = float(row["p_adj"])
            stats_line = (
                f"<br>p = {p:.3g}  ·  p<sub>adj</sub> = {p_adj:.3g}"
                f"  ·  {row.get('signif', '')}"
            )
        else:
            stats_line = f"<br>p = {p:.3g}  ·  {row.get('signif', '')}"
    bg_line = ""
    if bool(row.get("below_background", False)):
        ctrl = str(row.get("Background_control") or "control")
        bg_line = f"<br><i>below background (Cq &gt; {ctrl})</i>"
    return (
        f"<b>{row['Sample']}</b><br>"
        f"Target: {row['Target']}<br>"
        f"n = {int(row['n'])}<br>"
        f"Cq = {row['Cq_mean']:.2f} ± {row['Cq_sd']:.2f}<br>"
        f"ΔCq = {row['dCq']:.2f} ± {row['dCq_sd']:.2f}<br>"
        f"RQ = {row['RQ']:.3g}"
        f" (95%-style: {row['RQ_low']:.3g} – {row['RQ_high']:.3g})"
        f"{stats_line}{bg_line}"
    )


def _format_rq_label(value: float) -> str:
    """Compact, human-readable RQ label for plotting on bars."""
    if value is None or pd.isna(value):
        return ""
    v = float(value)
    if v == 0:
        return "0"
    abs_v = abs(v)
    if abs_v >= 100:
        return f"{v:.0f}"
    if abs_v >= 1:
        return f"{v:.2f}"
    if abs_v >= 0.01:
        return f"{v:.3f}"
    if abs_v >= 1e-4:
        return f"{v:.4f}"
    # Very small numbers: 2 sig figs in scientific notation.
    return f"{v:.1e}".replace("e-0", "e-").replace("e+0", "e+")


def relative_quantity_barplot(
    rq_table: pd.DataFrame,
    *,
    reference: str | None = None,
    target: str | None = None,
    title: str | None = None,
    sample_order: Iterable[str] | None = None,
    log_y: bool = True,
    show_values: bool = True,
    show_signif: bool = True,
    value_position: str = "inside",
    reference_color: str = "#7f7f7f",
    sample_color: str = "#1f77b4",
    dim_below_background: bool = False,
    below_background_color: str = "#cfcfcf",
) -> go.Figure:
    """Build an interactive bar plot of RQ values with asymmetric error bars.

    Parameters
    ----------
    rq_table:
        Output of :func:`qpictures.analysis.relative_quantification`.
    reference, target:
        Optional metadata for the chart title; auto-detected from the table when
        not supplied.
    title:
        Override the auto-generated title.
    sample_order:
        Explicit ordering of samples on the x-axis. Defaults to the order in
        ``rq_table`` (reference first, then descending RQ).
    log_y:
        Use a log-10 y-axis (equivalent to ``log_y=True`` in Plotly Express).
        Default ``True``.
    show_values:
        Draw the RQ value on each bar. Default ``True``.
    show_signif:
        If the input table has a ``signif`` column, draw the significance
        star above the value on each non-reference bar. Default ``True``.
    value_position:
        Where to draw value labels — ``"inside"`` (default, at the top inside
        the bar), ``"outside"`` (above the error bar) or ``"auto"``.
    dim_below_background:
        If ``True`` and the table has a ``below_background`` column, draw any
        sample below the control's mean Cq (RQ < 1) in ``below_background_color``
        instead of ``sample_color``.
    below_background_color:
        Color used for below-background bars when ``dim_below_background``
        is ``True``.
    """
    if rq_table.empty:
        raise ValueError("Cannot plot an empty RQ table.")

    table = rq_table.copy()
    if sample_order is not None:
        order = list(sample_order)
        table["__rank"] = table["Sample"].map({s: i for i, s in enumerate(order)})
        if table["__rank"].isna().any():
            missing = table.loc[table["__rank"].isna(), "Sample"].tolist()
            raise ValueError(f"sample_order is missing: {missing}")
        table = table.sort_values("__rank").drop(columns="__rank")

    if reference is None and "Reference" in table.columns:
        reference = str(table["Reference"].iloc[0])
    if target is None and "Target" in table.columns:
        target = str(table["Target"].iloc[0])

    has_bg = dim_below_background and "below_background" in table.columns
    colors: list[str] = []
    for _, row in table.iterrows():
        if bool(row.get("is_reference", False)):
            colors.append(reference_color)
        elif has_bg and bool(row.get("below_background", False)):
            colors.append(below_background_color)
        else:
            colors.append(sample_color)
    hover_text = [_format_hover(row) for _, row in table.iterrows()]

    # On a log axis, a zero/negative lower error bar is invalid. Clamp the
    # downward error so the bar can't drop to 0.
    if log_y:
        rq_safe_floor = np.clip(table["RQ"].to_numpy() * 1e-6, 1e-300, None)
        err_minus = np.minimum(
            table["RQ_err_minus"].to_numpy(),
            table["RQ"].to_numpy() - rq_safe_floor,
        )
    else:
        err_minus = table["RQ_err_minus"].to_numpy()

    # Compose per-bar text: significance star above the RQ value (multi-line),
    # so the bar reads e.g.  **<br>0.499  with `insidetextanchor="end"`.
    has_signif = show_signif and "signif" in table.columns
    if show_values or has_signif:
        labels: list[str] = []
        for _, row in table.iterrows():
            value_text = _format_rq_label(row["RQ"]) if show_values else ""
            star = ""
            if (
                has_signif
                and row.get("signif")
                and not bool(row.get("is_reference", False))
            ):
                star = str(row["signif"])
            if star and value_text:
                labels.append(f"{star}<br>{value_text}")
            else:
                labels.append(star or value_text)
        value_labels: list[str] | None = labels
    else:
        value_labels = None

    text_colors = [
        "#555" if c == below_background_color else "white" for c in colors
    ]
    fig = go.Figure(
        data=[
            go.Bar(
                x=table["Sample"],
                y=table["RQ"],
                marker=dict(color=colors, line=dict(color="#222", width=0.5)),
                error_y=dict(
                    type="data",
                    symmetric=False,
                    array=table["RQ_err_plus"],
                    arrayminus=err_minus,
                    thickness=1.2,
                    width=6,
                    color="#222",
                ),
                text=value_labels,
                textposition=value_position if value_labels is not None else "none",
                insidetextanchor="end",
                textangle=0,
                textfont=dict(size=11, color=text_colors),
                outsidetextfont=dict(size=11, color="#222"),
                cliponaxis=False,
                constraintext="none",
                hovertext=hover_text,
                hoverinfo="text",
                name="Relative quantity",
            )
        ]
    )

    bits = []
    if target:
        bits.append(f"target <b>{target}</b>")
    if reference:
        bits.append(f"vs reference <b>{reference}</b>")
    subtitle = " ".join(bits)
    fig.update_layout(
        title=dict(
            text=title or f"Relative quantity ({subtitle})" if subtitle else "Relative quantity",
            x=0.02,
            xanchor="left",
        ),
        template="plotly_white",
        xaxis=dict(title="Sample", tickangle=-30, categoryorder="array",
                   categoryarray=list(table["Sample"])),
        yaxis=dict(
            title="Relative quantity (2<sup>−ΔCq</sup>)",
            type="log" if log_y else "linear",
            zeroline=False,
        ),
        bargap=0.25,
        margin=dict(l=70, r=30, t=70, b=120),
        showlegend=False,
        hoverlabel=dict(bgcolor="white", font_size=12),
    )

    # Reference line at RQ = 1 (i.e., reference sample baseline).
    fig.add_hline(
        y=1.0,
        line=dict(color=reference_color, width=1, dash="dash"),
        annotation_text="reference (RQ = 1)",
        annotation_position="top right",
        annotation_font_color=reference_color,
    )
    return fig


def _format_hover_paired(row: pd.Series, condition: str) -> str:
    stats_line = ""
    if "p_value" in row.index and pd.notna(row.get("p_value")):
        p = float(row["p_value"])
        if "p_adj" in row.index and pd.notna(row.get("p_adj")):
            p_adj = float(row["p_adj"])
            stats_line = (
                f"<br>p = {p:.3g}  ·  p<sub>adj</sub> = {p_adj:.3g}"
                f"  ·  {row.get('signif', '')}"
            )
        else:
            stats_line = f"<br>p = {p:.3g}  ·  {row.get('signif', '')}"
    bg_line = ""
    if bool(row.get("below_background", False)):
        ctrl = str(row.get("Background_control") or "control")
        bg_line = f"<br><i>below background (Cq &gt; {ctrl})</i>"
    return (
        f"<b>{row['Sample']}</b><br>"
        f"Base: {row['Base']}  ·  Condition: {condition}<br>"
        f"Target: {row['Target']}<br>"
        f"n = {int(row['n'])}<br>"
        f"Cq = {row['Cq_mean']:.2f} ± {row['Cq_sd']:.2f}<br>"
        f"ΔCq = {row['dCq']:.2f} ± {row['dCq_sd']:.2f}<br>"
        f"vs reference: <b>{row['Reference_sample']}</b><br>"
        f"RQ = {row['RQ']:.3g}"
        f" (range: {row['RQ_low']:.3g} – {row['RQ_high']:.3g})"
        f"{stats_line}{bg_line}"
    )


def paired_relative_quantity_barplot(
    rq_table: pd.DataFrame,
    *,
    reference_base: str | None = None,
    target: str | None = None,
    title: str | None = None,
    log_y: bool = True,
    show_values: bool = True,
    show_signif: bool = True,
    value_position: str = "inside",
    condition_order: Iterable[str] | None = None,
    base_order: Iterable[str] | None = None,
    color_map: Mapping[str, str] | None = None,
    reference_color: str = "#7f7f7f",
    dim_below_background: bool = False,
    below_background_color: str = "#cfcfcf",
    group_shading: bool = False,
    mark_reference_pattern: bool = False,
    show_per_group_reference: bool = True,
) -> go.Figure:
    """Grouped bar plot of per-condition paired relative quantities.

    One x-axis group per ``Condition``, one bar per ``Base`` within each group,
    bars colored by ``Base``. The reference base sits at RQ = 1 in every group
    (with propagated SD from its own replicates).

    Expects the table emitted by
    :func:`qpictures.analysis.paired_relative_quantification`.

    Grouping cues (all on by default, can be turned off individually):

    - ``group_shading``: alternating subtle grey bands behind each condition
      group, so a comparison unit reads as a single visual box.
    - ``mark_reference_pattern``: the reference-base trace is rendered with
      diagonal hatching, so within each group the baseline bar is unmistakable.
    - ``show_per_group_reference``: each x-axis tick gets a small
      ``vs <ref_sample>`` line, so the specific sample each group is being
      compared against is spelled out.
    """
    if rq_table.empty:
        raise ValueError("Cannot plot an empty paired RQ table.")

    required = {"Base", "Condition", "RQ", "RQ_err_plus", "RQ_err_minus", "is_reference"}
    missing = required.difference(rq_table.columns)
    if missing:
        raise ValueError(f"Paired RQ table missing columns: {sorted(missing)}")

    table = rq_table.copy()

    if reference_base is None and "Reference_base" in table.columns:
        reference_base = str(table["Reference_base"].iloc[0])
    if target is None and "Target" in table.columns:
        target = str(table["Target"].iloc[0])

    conditions = (
        list(condition_order) if condition_order is not None
        else sorted(table["Condition"].dropna().unique().tolist())
    )
    bases_present = table["Base"].dropna().unique().tolist()
    if base_order is not None:
        bases = list(base_order)
    else:
        # reference base first, then alphabetical
        bases = ([reference_base] if reference_base in bases_present else []) + sorted(
            b for b in bases_present if b != reference_base
        )

    colors = _build_color_map(
        bases,
        reference_base,
        reference_color=reference_color,
        overrides=color_map,
    )

    # Pre-format bar labels: significance star above the RQ value (multi-line),
    # so the bar reads e.g.  **\n0.499  with `insidetextanchor="end"`.
    label_cache: dict[tuple[str, str], str] = {}
    for _, row in table.iterrows():
        key = (row["Base"], row["Condition"])
        value_text = _format_rq_label(row["RQ"]) if show_values else ""
        star = ""
        if (
            show_signif
            and "signif" in table.columns
            and row.get("signif")
            and not bool(row.get("is_reference", False))
        ):
            star = str(row["signif"])
        if star and value_text:
            label_cache[key] = f"{star}<br>{value_text}"
        else:
            label_cache[key] = star or value_text

    has_bg = dim_below_background and "below_background" in table.columns

    fig = go.Figure()
    for base in bases:
        sub = table.loc[table["Base"] == base].set_index("Condition")
        # Align this base's bars to the global condition axis (NaN-fill missing
        # so all traces have the same x categories).
        y = [float(sub.loc[c, "RQ"]) if c in sub.index else float("nan") for c in conditions]
        e_plus = [
            float(sub.loc[c, "RQ_err_plus"]) if c in sub.index else 0.0 for c in conditions
        ]
        e_minus_raw = [
            float(sub.loc[c, "RQ_err_minus"]) if c in sub.index else 0.0 for c in conditions
        ]
        if log_y:
            e_minus = []
            for yi, ei in zip(y, e_minus_raw):
                if np.isnan(yi):
                    e_minus.append(0.0)
                else:
                    floor = max(yi * 1e-6, 1e-300)
                    e_minus.append(min(ei, yi - floor))
        else:
            e_minus = e_minus_raw

        text = [
            label_cache.get((base, c), "")
            if (show_values or show_signif) and c in sub.index
            else ""
            for c in conditions
        ]
        hover = [
            _format_hover_paired(sub.loc[c], c) if c in sub.index else ""
            for c in conditions
        ]

        is_ref = base == reference_base
        # Per-bar color so we can dim individual (base, condition) cells
        # that fall below the per-condition control.
        bar_colors: list[str] = []
        bar_text_colors: list[str] = []
        for c in conditions:
            below = (
                has_bg
                and c in sub.index
                and bool(sub.loc[c].get("below_background", False))
            )
            if below and not is_ref:
                bar_colors.append(below_background_color)
                bar_text_colors.append("#555")
            else:
                bar_colors.append(colors[base])
                bar_text_colors.append("white")

        marker_kwargs: dict = dict(
            color=bar_colors,
            line=dict(color="#222", width=0.5),
        )
        if is_ref and mark_reference_pattern:
            # Diagonal hatching on every bar of the reference trace makes the
            # baseline visually obvious within each comparison group.
            marker_kwargs["pattern"] = dict(
                shape="/",
                solidity=0.35,
                fgcolor="rgba(255,255,255,0.75)",
            )

        fig.add_trace(
            go.Bar(
                name=f"{base} (ref)" if is_ref else base,
                x=conditions,
                y=y,
                marker=marker_kwargs,
                error_y=dict(
                    type="data",
                    symmetric=False,
                    array=e_plus,
                    arrayminus=e_minus,
                    thickness=1.2,
                    width=5,
                    color="#222",
                ),
                text=text,
                textposition=value_position if (show_values or show_signif) else "none",
                insidetextanchor="end",
                textangle=0,
                textfont=dict(size=10, color=bar_text_colors),
                outsidetextfont=dict(size=10, color="#222"),
                cliponaxis=False,
                constraintext="none",
                hovertext=hover,
                hoverinfo="text",
                legendgroup=base,
            )
        )

    # Alternating background bands per condition group make each comparison
    # unit read as a single visual box, even when groups have many bases.
    if group_shading and len(conditions) > 1:
        for i in range(len(conditions)):
            if i % 2 == 1:
                fig.add_shape(
                    type="rect",
                    xref="x", yref="paper",
                    x0=i - 0.5, x1=i + 0.5,
                    y0=0, y1=1,
                    fillcolor="rgba(0,0,0,0.045)",
                    line=dict(width=0),
                    layer="below",
                )

    # Per-group "vs <ref_sample>" tick label: spells out the actual reference
    # sample (e.g. BB_pbs) each group's bars are being compared against.
    xaxis_kwargs: dict = dict(
        title="Condition",
        tickangle=-20,
        categoryorder="array",
        categoryarray=conditions,
    )
    if show_per_group_reference and "Reference_sample" in table.columns:
        ref_sample_by_cond: dict[str, str] = {}
        for c in conditions:
            ref_rows = table.loc[
                (table["Condition"] == c) & table["is_reference"].astype(bool),
                "Reference_sample",
            ]
            if not ref_rows.empty:
                ref_sample_by_cond[c] = str(ref_rows.iloc[0])
        if ref_sample_by_cond:
            ticktext = []
            for c in conditions:
                ref = ref_sample_by_cond.get(c)
                if ref:
                    ticktext.append(
                        f"{c}<br>"
                        f"<span style='font-size:10px;color:#888'>vs {ref}</span>"
                    )
                else:
                    ticktext.append(c)
            xaxis_kwargs["tickmode"] = "array"
            xaxis_kwargs["tickvals"] = list(range(len(conditions)))
            xaxis_kwargs["ticktext"] = ticktext

    bits: list[str] = []
    if target:
        bits.append(f"target <b>{target}</b>")
    if reference_base:
        bits.append(f"reference base <b>{reference_base}</b>")
    subtitle = " · ".join(bits)
    fig.update_layout(
        title=dict(
            text=title or (f"Paired relative quantity ({subtitle})" if subtitle else "Paired relative quantity"),
            x=0.02,
            xanchor="left",
        ),
        template="plotly_white",
        barmode="group",
        bargap=0.25,
        bargroupgap=0.05,
        xaxis=xaxis_kwargs,
        yaxis=dict(
            title="Relative quantity (2<sup>−ΔCq</sup>)",
            type="log" if log_y else "linear",
            zeroline=False,
        ),
        legend=dict(title=dict(text="Base"), orientation="v"),
        margin=dict(l=70, r=30, t=70, b=130),
        hoverlabel=dict(bgcolor="white", font_size=12),
    )

    fig.add_hline(
        y=1.0,
        line=dict(color=reference_color, width=1, dash="dash"),
        annotation_text=(
            f"reference base = {reference_base} (RQ = 1)"
            if reference_base
            else "reference (RQ = 1)"
        ),
        annotation_position="top left",
        annotation_font_color=reference_color,
    )
    return fig


# --------------------------------------------------------------------------- #
# Standard curve / absolute quantification plots
# --------------------------------------------------------------------------- #


def _format_quantity_label(value: float) -> str:
    """Compact, human-readable quantity label for bar tops."""
    if value is None or pd.isna(value):
        return ""
    v = float(value)
    if v == 0:
        return "0"
    abs_v = abs(v)
    if abs_v >= 1000:
        return f"{v:.0f}"
    if abs_v >= 100:
        return f"{v:.0f}"
    if abs_v >= 10:
        return f"{v:.1f}"
    if abs_v >= 1:
        return f"{v:.2f}"
    if abs_v >= 0.01:
        return f"{v:.3f}"
    if abs_v >= 1e-4:
        return f"{v:.4f}"
    return f"{v:.1e}".replace("e-0", "e-").replace("e+0", "e+")


def _curve_annotation_text(curve: "StandardCurve") -> str:
    """One-block annotation summarising the fit (slope, R^2, efficiency, ...)."""
    eff_txt = "—" if pd.isna(curve.efficiency) else f"{curve.efficiency * 100:.1f}%"
    return (
        f"<b>Cq = {curve.slope:.3f} · log<sub>10</sub>(Q) + {curve.intercept:.3f}</b><br>"
        f"R² = {curve.r_squared:.4f}  ·  Efficiency = {eff_txt}<br>"
        f"slope SE = {curve.slope_se:.3f}  ·  residual SD = {curve.residual_sd:.3f}<br>"
        f"n = {curve.n_points} replicates across {curve.n_levels} levels"
    )


def standard_curve_plot(
    curve: "StandardCurve",
    *,
    title: str | None = None,
    unknowns_table: pd.DataFrame | None = None,
    standard_color: str = "#1f77b4",
    line_color: str = "#d62728",
    unknown_color: str = "#7f7f7f",
    show_band: bool = True,
    band_color: str = "rgba(214, 39, 40, 0.15)",
) -> go.Figure:
    """Plot the fitted standard curve (Cq vs log10 Quantity) plus diagnostics.

    Parameters
    ----------
    curve:
        Fitted :class:`qpictures.analysis.StandardCurve`.
    title:
        Override the auto-generated title.
    unknowns_table:
        Optional output of :func:`qpictures.analysis.absolute_quantification`.
        When given, each unknown sample is overlaid on the curve at
        ``(log10_Q, Cq_mean)`` so the data placement against the standards
        is visible at a glance.
    show_band:
        If ``True`` (default), shade the ±1 residual-SD band around the
        regression line.
    """
    from .analysis import StandardCurve  # local import to avoid cycle at import time

    if not isinstance(curve, StandardCurve):  # pragma: no cover - defensive
        raise TypeError(f"Expected StandardCurve, got {type(curve).__name__}.")

    standards = curve.standards
    if standards.empty:
        raise ValueError("Standard curve has no underlying standards to plot.")

    x_pad = max((curve.log10_q_max - curve.log10_q_min) * 0.05, 0.1)
    x_line = np.linspace(
        curve.log10_q_min - x_pad, curve.log10_q_max + x_pad, 200
    )
    y_line = curve.slope * x_line + curve.intercept

    fig = go.Figure()

    # Confidence/uncertainty band around the regression line (±1 residual SD).
    if show_band and np.isfinite(curve.residual_sd):
        fig.add_trace(
            go.Scatter(
                x=np.concatenate([x_line, x_line[::-1]]),
                y=np.concatenate([
                    y_line + curve.residual_sd,
                    (y_line - curve.residual_sd)[::-1],
                ]),
                fill="toself",
                fillcolor=band_color,
                line=dict(width=0),
                hoverinfo="skip",
                name="± 1 residual SD",
                showlegend=True,
            )
        )

    fig.add_trace(
        go.Scatter(
            x=x_line,
            y=y_line,
            mode="lines",
            line=dict(color=line_color, width=2),
            name="Fit",
            hoverinfo="skip",
        )
    )

    # Standards: hover shows the per-well sample name + concentration.
    hover_std = [
        f"<b>{row['Sample']}</b><br>"
        f"Quantity = {row['Quantity']:.4g}<br>"
        f"log<sub>10</sub>(Q) = {row['log10_Q']:.3f}<br>"
        f"Cq = {row['Cq']:.3f}"
        for _, row in standards.iterrows()
    ]
    fig.add_trace(
        go.Scatter(
            x=standards["log10_Q"],
            y=standards["Cq"],
            mode="markers",
            marker=dict(
                color=standard_color,
                size=9,
                line=dict(color="#222", width=0.5),
                symbol="circle",
            ),
            name="Standards",
            hovertext=hover_std,
            hoverinfo="text",
        )
    )

    if unknowns_table is not None and not unknowns_table.empty:
        unk = unknowns_table.copy()
        # Skip rows with no usable Cq (shouldn't happen, but be safe).
        unk = unk.loc[unk["Cq_mean"].notna() & unk["log10_Q"].notna()]
        if not unk.empty:
            hover_unk = [
                f"<b>{row['Sample']}</b><br>"
                f"n = {int(row['n'])}<br>"
                f"Cq = {row['Cq_mean']:.3f} ± {row['Cq_sd']:.3f}<br>"
                f"log<sub>10</sub>(Q) = {row['log10_Q']:.3f}<br>"
                f"Q = {row['Quantity']:.4g}"
                f"{' <i>(extrapolated)</i>' if bool(row.get('extrapolated', False)) else ''}"
                for _, row in unk.iterrows()
            ]
            symbols = [
                "x" if bool(row.get("extrapolated", False)) else "diamond"
                for _, row in unk.iterrows()
            ]
            fig.add_trace(
                go.Scatter(
                    x=unk["log10_Q"],
                    y=unk["Cq_mean"],
                    mode="markers",
                    marker=dict(
                        color=unknown_color,
                        size=8,
                        line=dict(color="#222", width=0.5),
                        symbol=symbols,
                        opacity=0.7,
                    ),
                    name="Unknowns (projected)",
                    hovertext=hover_unk,
                    hoverinfo="text",
                )
            )

    fig.add_annotation(
        xref="paper", yref="paper",
        x=0.02, y=0.02, xanchor="left", yanchor="bottom",
        showarrow=False,
        align="left",
        text=_curve_annotation_text(curve),
        bgcolor="rgba(255,255,255,0.85)",
        bordercolor="#888",
        borderwidth=1,
        borderpad=6,
        font=dict(size=11, color="#222"),
    )

    fig.update_layout(
        title=dict(
            text=title or f"Standard curve  (target <b>{curve.target}</b>)",
            x=0.02,
            xanchor="left",
        ),
        template="plotly_white",
        xaxis=dict(title="log<sub>10</sub>(Quantity)"),
        yaxis=dict(title="Cq", autorange="reversed", zeroline=False),
        margin=dict(l=70, r=30, t=70, b=70),
        hoverlabel=dict(bgcolor="white", font_size=12),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    return fig


def _format_absolute_hover(row: pd.Series, curve: "StandardCurve | None" = None) -> str:
    extrap = ""
    if bool(row.get("extrapolated", False)):
        extrap = "<br><i>extrapolated outside the standard range</i>"
    bg_line = ""
    if bool(row.get("below_background", False)):
        ctrl = str(row.get("Background_control") or "control")
        bg_line = f"<br><i>below background (Cq &gt; {ctrl})</i>"
    return (
        f"<b>{row['Sample']}</b><br>"
        f"Target: {row['Target']}<br>"
        f"n = {int(row['n'])}<br>"
        f"Cq = {row['Cq_mean']:.2f} ± {row['Cq_sd']:.2f}<br>"
        f"log<sub>10</sub>(Q) = {row['log10_Q']:.3f} ± {row['log10_Q_sd']:.3f}<br>"
        f"Q = {row['Quantity']:.4g} "
        f"(range: {row['Quantity_low']:.4g} – {row['Quantity_high']:.4g})"
        f"{extrap}{bg_line}"
    )


def absolute_quantity_barplot(
    abs_table: pd.DataFrame,
    *,
    target: str | None = None,
    title: str | None = None,
    sample_order: Iterable[str] | None = None,
    log_y: bool = True,
    show_values: bool = True,
    value_position: str = "inside",
    sample_color: str = "#1f77b4",
    extrapolated_color: str = "#ff7f0e",
    dim_below_background: bool = False,
    below_background_color: str = "#cfcfcf",
    quantity_unit: str | None = None,
    mark_extrapolated: bool = True,
    curve: "StandardCurve | None" = None,
) -> go.Figure:
    """Bar plot of absolute quantities (back-calculated from a standard curve).

    Parameters
    ----------
    abs_table:
        Output of :func:`qpictures.analysis.absolute_quantification`.
    target:
        Optional metadata for the chart title; auto-detected from the table.
    title:
        Override the auto-generated title.
    sample_order:
        Explicit ordering of samples on the x-axis. Default: descending
        ``Quantity``.
    log_y:
        Use a log-10 y-axis (typical for standard-curve outputs). Default ``True``.
    show_values:
        Draw the quantity value on each bar.
    sample_color:
        Default bar color.
    extrapolated_color:
        Bar color for samples whose mean Cq falls *outside* the calibrated
        standard range (``extrapolated == True``). Only applied when
        ``mark_extrapolated`` is ``True``.
    dim_below_background:
        If ``True`` and the table has a ``below_background`` column, samples
        whose Cq exceeds the control's are rendered in ``below_background_color``.
    quantity_unit:
        Unit string appended to the y-axis title (e.g. ``"copies/µL"``).
    mark_extrapolated:
        If ``True`` (default), recolor bars whose mean Cq is outside the
        standard range and add an annotation.
    curve:
        Optional :class:`StandardCurve` (used purely for the y-axis annotation
        of standard min/max guide lines, when log_y=True).
    """
    if abs_table.empty:
        raise ValueError("Cannot plot an empty absolute-quantity table.")

    table = abs_table.copy()
    if sample_order is not None:
        order = list(sample_order)
        table["__rank"] = table["Sample"].map({s: i for i, s in enumerate(order)})
        if table["__rank"].isna().any():
            missing = table.loc[table["__rank"].isna(), "Sample"].tolist()
            raise ValueError(f"sample_order is missing: {missing}")
        table = table.sort_values("__rank").drop(columns="__rank")

    if target is None and "Target" in table.columns:
        target = str(table["Target"].iloc[0])

    has_bg = dim_below_background and "below_background" in table.columns
    has_extrap = mark_extrapolated and "extrapolated" in table.columns

    colors: list[str] = []
    for _, row in table.iterrows():
        if has_bg and bool(row.get("below_background", False)):
            colors.append(below_background_color)
        elif has_extrap and bool(row.get("extrapolated", False)):
            colors.append(extrapolated_color)
        else:
            colors.append(sample_color)

    hover_text = [_format_absolute_hover(row, curve) for _, row in table.iterrows()]

    if log_y:
        q_safe_floor = np.clip(table["Quantity"].to_numpy() * 1e-6, 1e-300, None)
        err_minus = np.minimum(
            table["Q_err_minus"].to_numpy(),
            table["Quantity"].to_numpy() - q_safe_floor,
        )
    else:
        err_minus = table["Q_err_minus"].to_numpy()

    if show_values:
        labels = [_format_quantity_label(v) for v in table["Quantity"]]
    else:
        labels = None

    text_colors = []
    for c in colors:
        if c in (below_background_color, extrapolated_color):
            text_colors.append("#333")
        else:
            text_colors.append("white")

    fig = go.Figure(
        data=[
            go.Bar(
                x=table["Sample"],
                y=table["Quantity"],
                marker=dict(color=colors, line=dict(color="#222", width=0.5)),
                error_y=dict(
                    type="data",
                    symmetric=False,
                    array=table["Q_err_plus"],
                    arrayminus=err_minus,
                    thickness=1.2,
                    width=6,
                    color="#222",
                ),
                text=labels,
                textposition=value_position if labels is not None else "none",
                insidetextanchor="end",
                textangle=0,
                textfont=dict(size=11, color=text_colors),
                outsidetextfont=dict(size=11, color="#222"),
                cliponaxis=False,
                constraintext="none",
                hovertext=hover_text,
                hoverinfo="text",
                name="Absolute quantity",
            )
        ]
    )

    y_title = "Quantity"
    if quantity_unit:
        y_title = f"Quantity ({quantity_unit})"
    fig.update_layout(
        title=dict(
            text=title or (
                f"Absolute quantity (target <b>{target}</b>)" if target else "Absolute quantity"
            ),
            x=0.02,
            xanchor="left",
        ),
        template="plotly_white",
        xaxis=dict(
            title="Sample", tickangle=-30,
            categoryorder="array", categoryarray=list(table["Sample"]),
        ),
        yaxis=dict(
            title=y_title,
            type="log" if log_y else "linear",
            zeroline=False,
        ),
        bargap=0.25,
        margin=dict(l=70, r=30, t=70, b=120),
        showlegend=False,
        hoverlabel=dict(bgcolor="white", font_size=12),
    )

    if curve is not None and log_y:
        for q_log, label in (
            (curve.log10_q_min, f"std min (Q = {10**curve.log10_q_min:.3g})"),
            (curve.log10_q_max, f"std max (Q = {10**curve.log10_q_max:.3g})"),
        ):
            fig.add_hline(
                y=10**q_log,
                line=dict(color="#888", width=1, dash="dot"),
                annotation_text=label,
                annotation_position="top right",
                annotation_font_color="#666",
                annotation_font_size=10,
            )

    return fig


def standard_curve_with_quantities_figure(
    curve: "StandardCurve",
    abs_table: pd.DataFrame,
    *,
    title: str | None = None,
    log_y: bool = True,
    show_values: bool = True,
    quantity_unit: str | None = None,
    dim_below_background: bool = False,
    mark_extrapolated: bool = True,
) -> go.Figure:
    """Combined figure: standard curve (top) + absolute-quantity bars (bottom).

    A convenience wrapper that stacks :func:`standard_curve_plot` and
    :func:`absolute_quantity_barplot` into a single Plotly subplot figure,
    suitable for writing as one HTML.
    """
    curve_fig = standard_curve_plot(curve, unknowns_table=abs_table)
    bar_fig = absolute_quantity_barplot(
        abs_table,
        target=curve.target,
        log_y=log_y,
        show_values=show_values,
        quantity_unit=quantity_unit,
        dim_below_background=dim_below_background,
        mark_extrapolated=mark_extrapolated,
        curve=curve,
    )

    eff_txt = "—" if pd.isna(curve.efficiency) else f"{curve.efficiency * 100:.1f}%"
    fig = make_subplots(
        rows=2, cols=1,
        row_heights=[0.42, 0.58],
        vertical_spacing=0.22,
        subplot_titles=(
            f"Standard curve · R² = {curve.r_squared:.4f} · E = {eff_txt}",
            "Absolute quantity per sample",
        ),
    )

    for tr in curve_fig.data:
        fig.add_trace(tr, row=1, col=1)
    for tr in bar_fig.data:
        # Bars don't need to occupy a legend entry — the x-axis labels each
        # sample already, and the curve traces above are what benefit from
        # disambiguation in the shared legend.
        tr.showlegend = False
        fig.add_trace(tr, row=2, col=1)

    fig.update_xaxes(curve_fig.layout.xaxis.to_plotly_json(), row=1, col=1)
    fig.update_yaxes(curve_fig.layout.yaxis.to_plotly_json(), row=1, col=1)
    fig.update_xaxes(bar_fig.layout.xaxis.to_plotly_json(), row=2, col=1)
    fig.update_yaxes(bar_fig.layout.yaxis.to_plotly_json(), row=2, col=1)

    # Bring along the standard-curve annotations (stats box) but anchor to the
    # top subplot's axes instead of the standalone figure's paper coordinates.
    for ann in curve_fig.layout.annotations or ():
        a = ann.to_plotly_json()
        if a.get("xref") == "paper" and a.get("yref") == "paper":
            a["xref"] = "x domain"
            a["yref"] = "y domain"
        fig.add_annotation(a, row=1, col=1)

    # Standard-range guide lines on the bar subplot. Anchor labels on the left
    # so they don't collide with the right-side legend; tuck "std max" *inside*
    # so the label doesn't punch through the subplot title above it.
    if log_y:
        for q_log, label, anchor in (
            (curve.log10_q_min, f"std min (Q = {10**curve.log10_q_min:.3g})", "bottom left"),
            (curve.log10_q_max, f"std max (Q = {10**curve.log10_q_max:.3g})", "bottom left"),
        ):
            fig.add_hline(
                y=10**q_log,
                line=dict(color="#888", width=1, dash="dot"),
                annotation_text=label,
                annotation_position=anchor,
                annotation_font_color="#666",
                annotation_font_size=10,
                row=2, col=1,
            )

    fig.update_layout(
        title=dict(
            text=title or f"Absolute quantification (target <b>{curve.target}</b>)",
            x=0.02,
            xanchor="left",
            y=0.985,
            yanchor="top",
        ),
        template="plotly_white",
        showlegend=True,
        # Park the legend to the right of the top subplot so it neither
        # overlaps the figure title nor the bar subplot's title.
        legend=dict(
            orientation="v",
            yanchor="top", y=0.95,
            xanchor="left", x=1.02,
        ),
        margin=dict(l=70, r=160, t=100, b=140),
        hoverlabel=dict(bgcolor="white", font_size=12),
        bargap=0.25,
    )
    return fig
