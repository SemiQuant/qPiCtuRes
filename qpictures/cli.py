"""Command-line entrypoints for qPiCtuRes."""

from __future__ import annotations

import sys
from pathlib import Path

import click
import pandas as pd
from rich.console import Console
from rich.table import Table

from . import __version__
from .analysis import (
    absolute_quantification,
    fit_standard_curve,
    paired_relative_quantification,
    relative_quantification,
)
from .io import read_thermo_well_results
from .plot import (
    absolute_quantity_barplot,
    paired_relative_quantity_barplot,
    relative_quantity_barplot,
    standard_curve_plot,
    standard_curve_with_quantities_figure,
)

console = Console()

# Sentinel used by Click flag-with-optional-value support for --pdf / --png /
# --svg. When the flag is given on its own (no value), Click sets the option to
# this string, which we then expand to ``out_html.with_suffix(".pdf")`` etc.
_AUTO_PATH = "<auto>"


def _resolve_image_path(value: str | None, html_path: Path, ext: str) -> Path | None:
    """Convert a --pdf/--png/--svg option value to a concrete path or ``None``."""
    if value is None:
        return None
    if value == _AUTO_PATH or value == "":
        return html_path.with_suffix(f".{ext}")
    return Path(value)


def _tasks_option(default: tuple[str, ...] = ("Unknown",)):
    """Reusable decorator stack for the task filter options.

    - ``--task TEXT`` (repeatable): explicit list of Task values to keep.
    - ``--all-tasks``: keep every Task (controls, standards, ...).
    """

    def _decorator(func):
        func = click.option(
            "--all-tasks", "all_tasks", is_flag=True, default=False,
            help="Keep every Task (NTC / standards / etc.). "
                 f"Default keeps only {', '.join(default)}.",
        )(func)
        func = click.option(
            "--task", "task_filter", multiple=True, default=(),
            help="Task value(s) to keep, may be passed multiple times. "
                 f"Default {list(default)}.",
        )(func)
        return func

    return _decorator


def _resolve_tasks(task_filter: tuple[str, ...], all_tasks: bool,
                   default: tuple[str, ...] = ("Unknown",)) -> tuple[str, ...] | None:
    if all_tasks:
        return None
    return tuple(task_filter) if task_filter else default


def _image_export_options(func):
    """Reusable decorator stack for static image export options."""
    func = click.option(
        "--svg", "svg_out", default=None, is_flag=False, flag_value=_AUTO_PATH,
        help="Also write SVG. Use as flag (auto-derives from --out) or pass a path.",
    )(func)
    func = click.option(
        "--png", "png_out", default=None, is_flag=False, flag_value=_AUTO_PATH,
        help="Also write PNG. Use as flag (auto-derives from --out) or pass a path.",
    )(func)
    func = click.option(
        "--pdf", "pdf_out", default=None, is_flag=False, flag_value=_AUTO_PATH,
        help="Also write PDF. Use as flag (auto-derives from --out) or pass a path.",
    )(func)
    return func


def _export_static(fig, html_path: Path, *, pdf_out, png_out, svg_out) -> list[Path]:
    """Write any requested static image formats and return the produced paths."""
    written: list[Path] = []
    for value, ext in ((pdf_out, "pdf"), (png_out, "png"), (svg_out, "svg")):
        path = _resolve_image_path(value, html_path, ext)
        if path is None:
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        try:
            fig.write_image(path)
        except Exception as exc:  # pragma: no cover - best-effort surfacing
            console.print(
                f"[bold red]Could not write {ext.upper()}[/bold red] "
                f"({path}): {exc}\n"
                "  Hint: install the `kaleido` package "
                "(`pip install kaleido` / it's in environment.yml)."
            )
            continue
        written.append(path)
    return written


def _df_to_rich_table(df: pd.DataFrame, *, title: str | None = None,
                      float_fmt: str = "{:.3g}",
                      dim_column: str | None = None) -> Table:
    """Render a DataFrame as a rich Table.

    If ``dim_column`` is provided and the column exists, any row where that
    column is truthy is rendered with the ``dim`` style (visually greyed out).
    """
    table = Table(title=title, show_lines=False, header_style="bold cyan")
    for col in df.columns:
        justify = "right" if pd.api.types.is_numeric_dtype(df[col]) else "left"
        table.add_column(str(col), justify=justify, overflow="fold")
    for _, row in df.iterrows():
        cells = []
        for col in df.columns:
            val = row[col]
            if isinstance(val, float):
                cells.append("" if pd.isna(val) else float_fmt.format(val))
            elif isinstance(val, bool):
                cells.append("✓" if val else "")
            else:
                cells.append("" if val is None or (isinstance(val, float) and pd.isna(val)) else str(val))
        style = None
        if dim_column and dim_column in df.columns and bool(row.get(dim_column, False)):
            style = "dim"
        table.add_row(*cells, style=style)
    return table


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__, prog_name="qpictures")
def cli() -> None:
    """qPiCtuRes - interactive qPCR plots and analyses."""


@cli.command("inspect")
@click.argument("input_csv", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--keep-undetermined", is_flag=True,
              help="Keep wells with Undetermined Cq in the table.")
@_tasks_option()
def inspect_cmd(input_csv: Path, keep_undetermined: bool,
                task_filter: tuple[str, ...], all_tasks: bool) -> None:
    """Show a quick summary of a Thermo Well-Results CSV."""
    tasks = _resolve_tasks(task_filter, all_tasks)
    with click.progressbar(length=3, label="Reading", show_eta=False) as bar:
        df = read_thermo_well_results(
            input_csv, keep_undetermined=keep_undetermined, tasks=tasks
        )
        bar.update(1)
        targets = sorted(df["Target"].dropna().unique().tolist())
        bar.update(1)
        summary = (
            df.groupby(["Target", "Sample"], dropna=False)
            .agg(n=("Cq", "size"), n_amp=("Cq", lambda s: int(s.notna().sum())),
                 Cq_mean=("Cq", "mean"), Cq_sd=("Cq", "std"))
            .reset_index()
            .sort_values(["Target", "Cq_mean"], na_position="last")
        )
        bar.update(1)

    console.rule(f"[bold]{input_csv.name}[/bold]")
    console.print(f"Rows used: [bold]{len(df)}[/bold]   Targets: [bold]{', '.join(targets) or '<none>'}[/bold]")
    console.print(_df_to_rich_table(summary, title="Per-sample Cq summary"))


@cli.command("relative")
@click.argument("input_csv", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("-r", "--reference", required=True,
              help="Reference sample name (case-insensitive). Example: BB_PBS.")
@click.option("-t", "--target", default=None,
              help="Target / assay to analyze. Auto-detected if only one target is present.")
@click.option("-o", "--out", "out_html", type=click.Path(dir_okay=False, path_type=Path),
              default=None,
              help="Where to write the interactive HTML plot. "
                   "Defaults to <input>.relative.html next to the CSV.")
@click.option("--csv-out", "csv_out", type=click.Path(dir_okay=False, path_type=Path),
              default=None, help="Optionally write the RQ table to CSV.")
@click.option("--no-log", "log_y", flag_value=False, default=True,
              help="Use a linear y-axis instead of log10.")
@click.option("--open/--no-open", "open_browser", default=True,
              help="Open the resulting HTML in the default browser.")
@click.option("--title", default=None, help="Custom plot title.")
@click.option(
    "--test", "stat_test",
    type=click.Choice(["welch", "student", "mannwhitney", "none"], case_sensitive=False),
    default="welch", show_default=True,
    help="Per-sample significance test on raw Cq values vs the reference.",
)
@click.option(
    "--padjust", "padjust",
    type=click.Choice(["none", "bh", "fdr_bh", "holm", "bonferroni"], case_sensitive=False),
    default="none", show_default=True,
    help="Multiple-testing correction across all non-reference samples. "
         "Off by default; use e.g. `--padjust bh` for Benjamini–Hochberg FDR.",
)
@click.option("--no-signif", "show_signif", flag_value=False, default=True,
              help="Don't render significance stars on the plot.")
@click.option(
    "--background-rem", "background_rem",
    default=None, is_flag=False, flag_value="control",
    help="Use a named control sample as the background threshold (lookup is "
         "case-insensitive). Use as a flag to default to a sample called "
         "`control`, or pass an explicit name (e.g. `--background-rem PBS`). "
         "Any sample whose mean Cq exceeds the control's mean Cq is flagged "
         "with `below_background` and dimmed in the table and plot.",
)
@_tasks_option()
@_image_export_options
def relative_cmd(
    input_csv: Path,
    reference: str,
    target: str | None,
    out_html: Path | None,
    csv_out: Path | None,
    log_y: bool,
    open_browser: bool,
    title: str | None,
    stat_test: str,
    padjust: str,
    show_signif: bool,
    background_rem: str | None,
    task_filter: tuple[str, ...],
    all_tasks: bool,
    pdf_out: str | None,
    png_out: str | None,
    svg_out: str | None,
) -> None:
    """Build a relative-quantity bar plot from a Thermo Well-Results CSV."""
    out_html = out_html or input_csv.with_suffix(".relative.html")
    tasks = _resolve_tasks(task_filter, all_tasks)

    steps = ["read", "analyze", "plot", "write"]
    with click.progressbar(steps, label="qpictures relative",
                           item_show_func=lambda s: s or "") as bar:
        bar.update(0)
        df = read_thermo_well_results(input_csv, tasks=tasks)
        bar.update(1)

        try:
            rq, ref_match = relative_quantification(
                df, reference=reference, target=target,
                test=stat_test, padjust=padjust,
                background_control=background_rem,
            )
        except ValueError as exc:
            click.echo()
            console.print(f"[bold red]Error:[/bold red] {exc}")
            sys.exit(2)
        bar.update(1)

        fig = relative_quantity_barplot(
            rq,
            reference=ref_match.resolved,
            target=str(rq["Target"].iloc[0]),
            title=title,
            log_y=log_y,
            show_signif=show_signif,
            dim_below_background=background_rem is not None,
        )
        bar.update(1)

        out_html.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(out_html, include_plotlyjs="cdn", full_html=True)
        if csv_out:
            csv_out.parent.mkdir(parents=True, exist_ok=True)
            rq.to_csv(csv_out, index=False)
        image_paths = _export_static(
            fig, out_html, pdf_out=pdf_out, png_out=png_out, svg_out=svg_out
        )
        bar.update(1)

    if ref_match.resolved != ref_match.requested:
        console.print(
            f"[yellow]Reference matched case-insensitively:[/yellow] "
            f"{ref_match.requested!r} -> {ref_match.resolved!r} "
            f"(n={ref_match.n_replicates})"
        )

    display = rq.copy()
    for col in ("Cq_mean", "Cq_sd", "dCq", "dCq_sd", "RQ", "RQ_low", "RQ_high",
                "p_value", "p_adj"):
        if col in display.columns:
            display[col] = display[col].astype(float)
    stats_cols = [c for c in ("p_value", "p_adj", "signif") if c in display.columns]
    bg_enabled = background_rem is not None
    bg_cols = ["below_background"] if bg_enabled and "below_background" in display.columns else []
    display = display[
        ["Sample", "n", "Cq_mean", "Cq_sd", "dCq", "RQ", "RQ_low", "RQ_high",
         *stats_cols, "is_reference", *bg_cols]
    ]
    title_suffix_parts: list[str] = []
    if stat_test != "none":
        title_suffix_parts.append(f"test: {stat_test}, adj: {padjust}")
    if bg_enabled:
        n_below = int(rq["below_background"].sum()) if "below_background" in rq.columns else 0
        ctrl_resolved = (
            str(rq["Background_control"].iloc[0])
            if "Background_control" in rq.columns
            and len(rq) > 0
            and rq["Background_control"].iloc[0]
            else background_rem
        )
        title_suffix_parts.append(f"background vs {ctrl_resolved}: {n_below} below")
    title_suffix = f"  ({'; '.join(title_suffix_parts)})" if title_suffix_parts else ""
    console.print(
        _df_to_rich_table(
            display,
            title=f"Relative quantity vs {ref_match.resolved}{title_suffix}",
            dim_column="below_background" if bg_enabled else None,
        )
    )
    console.print(f"[green]Plot written to[/green] {out_html}")
    for p in image_paths:
        console.print(f"[green]Plot written to[/green] {p}")
    if csv_out:
        console.print(f"[green]Table written to[/green] {csv_out}")

    if open_browser:
        click.launch(str(out_html))


@cli.command("paired")
@click.argument("input_csv", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("-b", "--reference-base", required=True,
              help='Reference base prefix (case-insensitive). Example: "BB".')
@click.option("-s", "--separator", default="_", show_default=True,
              help="Sample-name separator splitting base from condition.")
@click.option("-t", "--target", default=None,
              help="Target / assay to analyze. Auto-detected if only one is present.")
@click.option("-o", "--out", "out_html", type=click.Path(dir_okay=False, path_type=Path),
              default=None,
              help="Where to write the interactive HTML plot. "
                   "Defaults to <input>.paired.html next to the CSV.")
@click.option("--csv-out", "csv_out", type=click.Path(dir_okay=False, path_type=Path),
              default=None, help="Optionally write the paired RQ table to CSV.")
@click.option("--no-log", "log_y", flag_value=False, default=True,
              help="Use a linear y-axis instead of log10.")
@click.option("--open/--no-open", "open_browser", default=True,
              help="Open the resulting HTML in the default browser.")
@click.option("--title", default=None, help="Custom plot title.")
@click.option(
    "--test", "stat_test",
    type=click.Choice(["welch", "student", "mannwhitney", "none"], case_sensitive=False),
    default="welch", show_default=True,
    help="Per-pair significance test on raw Cq values.",
)
@click.option(
    "--padjust", "padjust",
    type=click.Choice(["none", "bh", "fdr_bh", "holm", "bonferroni"], case_sensitive=False),
    default="none", show_default=True,
    help="Multiple-testing correction across all non-reference comparisons. "
         "Off by default; use e.g. `--padjust bh` for Benjamini–Hochberg FDR.",
)
@click.option("--no-signif", "show_signif", flag_value=False, default=True,
              help="Don't render significance stars on the plot.")
@click.option(
    "--background-rem", "background_rem",
    default=None, is_flag=False, flag_value="control",
    help="Use a named control sample as the background threshold (lookup is "
         "case-insensitive against the full target-filtered table, so a "
         "control without a separator like `Control` is still found). Use as "
         "a flag to default to `control`, or pass an explicit name "
         "(e.g. `--background-rem PBS`). Any (base, condition) row whose "
         "mean Cq exceeds the control's mean Cq is flagged with "
         "`below_background` and dimmed.",
)
@click.option(
    "--no-split", "no_split", is_flag=True, default=False,
    help="Suppress splitting sample names on the separator. Every sample is "
         "treated as its own base under a single synthetic condition, and "
         "`--reference-base` must match a *full* sample name "
         "(e.g. `BB_pbs`) — equivalent to a single-reference comparison.",
)
@_tasks_option()
@_image_export_options
def paired_cmd(
    input_csv: Path,
    reference_base: str,
    separator: str,
    target: str | None,
    out_html: Path | None,
    csv_out: Path | None,
    log_y: bool,
    open_browser: bool,
    title: str | None,
    stat_test: str,
    padjust: str,
    show_signif: bool,
    background_rem: str | None,
    no_split: bool,
    task_filter: tuple[str, ...],
    all_tasks: bool,
    pdf_out: str | None,
    png_out: str | None,
    svg_out: str | None,
) -> None:
    """Per-condition paired RQ plot.

    Splits each sample on the first separator (``base<sep>condition``) and
    normalizes every base in each condition against the reference base in the
    *same* condition (e.g. ``TL_pbs`` vs ``BB_pbs``).
    """
    out_html = out_html or input_csv.with_suffix(".paired.html")
    tasks = _resolve_tasks(task_filter, all_tasks)
    steps = ["read", "analyze", "plot", "write"]
    with click.progressbar(steps, label="qpictures paired",
                           item_show_func=lambda s: s or "") as bar:
        bar.update(0)
        df = read_thermo_well_results(input_csv, tasks=tasks)
        bar.update(1)

        try:
            rq, summary = paired_relative_quantification(
                df,
                reference_base=reference_base,
                separator=separator,
                target=target,
                test=stat_test,
                padjust=padjust,
                split_samples=not no_split,
                background_control=background_rem,
            )
        except ValueError as exc:
            click.echo()
            console.print(f"[bold red]Error:[/bold red] {exc}")
            sys.exit(2)
        bar.update(1)

        fig = paired_relative_quantity_barplot(
            rq,
            reference_base=summary.reference_base_resolved,
            target=str(rq["Target"].iloc[0]),
            log_y=log_y,
            title=title,
            show_signif=show_signif,
            dim_below_background=background_rem is not None,
        )
        bar.update(1)

        out_html.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(out_html, include_plotlyjs="cdn", full_html=True)
        if csv_out:
            csv_out.parent.mkdir(parents=True, exist_ok=True)
            rq.to_csv(csv_out, index=False)
        image_paths = _export_static(
            fig, out_html, pdf_out=pdf_out, png_out=png_out, svg_out=svg_out
        )
        bar.update(1)

    if summary.reference_base_resolved != summary.reference_base_requested:
        console.print(
            f"[yellow]Reference base matched case-insensitively:[/yellow] "
            f"{summary.reference_base_requested!r} -> {summary.reference_base_resolved!r}"
        )
    console.print(
        f"Bases:      [bold]{', '.join(summary.bases) or '<none>'}[/bold]"
    )
    console.print(
        f"Conditions: [bold]{', '.join(summary.conditions) or '<none>'}[/bold]"
    )
    if summary.skipped_samples:
        console.print(
            f"[yellow]Skipped (no '{summary.separator}'):[/yellow] "
            f"{', '.join(summary.skipped_samples)}"
        )
    if summary.skipped_conditions:
        console.print(
            f"[yellow]Skipped conditions (no reference base "
            f"{summary.reference_base_resolved!r}):[/yellow] "
            f"{', '.join(summary.skipped_conditions)}"
        )

    stats_cols = [c for c in ("p_value", "p_adj", "signif") if c in rq.columns]
    bg_enabled = background_rem is not None
    bg_cols = ["below_background"] if bg_enabled and "below_background" in rq.columns else []
    display = rq[
        [
            "Condition",
            "Base",
            "Sample",
            "n",
            "Cq_mean",
            "Cq_sd",
            "dCq",
            "RQ",
            "RQ_low",
            "RQ_high",
            *stats_cols,
            "is_reference",
            *bg_cols,
        ]
    ].copy()
    test_label = summary.test
    padj_label = summary.padjust
    title_suffix_parts: list[str] = []
    if test_label != "none":
        title_suffix_parts.append(f"test: {test_label}, adj: {padj_label}")
    if no_split:
        title_suffix_parts.append("no-split")
    if bg_enabled:
        n_below = int(rq["below_background"].sum()) if "below_background" in rq.columns else 0
        ctrl_resolved = (
            str(rq["Background_control"].iloc[0])
            if "Background_control" in rq.columns
            and len(rq) > 0
            and rq["Background_control"].iloc[0]
            else background_rem
        )
        title_suffix_parts.append(f"background vs {ctrl_resolved}: {n_below} below")
    title_suffix = f"  ({'; '.join(title_suffix_parts)})" if title_suffix_parts else ""
    console.print(
        _df_to_rich_table(
            display,
            title=f"Paired RQ vs base {summary.reference_base_resolved}{title_suffix}",
            dim_column="below_background" if bg_enabled else None,
        )
    )
    console.print(f"[green]Plot written to[/green] {out_html}")
    for p in image_paths:
        console.print(f"[green]Plot written to[/green] {p}")
    if csv_out:
        console.print(f"[green]Table written to[/green] {csv_out}")

    if open_browser:
        click.launch(str(out_html))


@cli.command("standcurve")
@click.argument("input_csv", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("-t", "--target", default=None,
              help="Target / assay to analyze. Auto-detected if only one is present.")
@click.option("-o", "--out", "out_html", type=click.Path(dir_okay=False, path_type=Path),
              default=None,
              help="Where to write the combined interactive HTML (standard curve "
                   "+ absolute-quantity bar plot). Defaults to "
                   "<input>.standcurve.html next to the CSV.")
@click.option("--curve-html", "curve_html", type=click.Path(dir_okay=False, path_type=Path),
              default=None,
              help="Also write a standalone HTML containing just the standard "
                   "curve plot (no bars).")
@click.option("--bars-html", "bars_html", type=click.Path(dir_okay=False, path_type=Path),
              default=None,
              help="Also write a standalone HTML containing just the "
                   "absolute-quantity bar plot.")
@click.option("--csv-out", "csv_out", type=click.Path(dir_okay=False, path_type=Path),
              default=None,
              help="Optionally write the per-sample absolute-quantity table to CSV.")
@click.option("--curve-csv", "curve_csv", type=click.Path(dir_okay=False, path_type=Path),
              default=None,
              help="Optionally write the per-replicate standard-curve table "
                   "(Sample, Cq, Quantity, log10_Q) to CSV.")
@click.option("--no-log", "log_y", flag_value=False, default=True,
              help="Use a linear y-axis on the bar plot instead of log10.")
@click.option("--open/--no-open", "open_browser", default=True,
              help="Open the resulting combined HTML in the default browser.")
@click.option("--title", default=None, help="Custom plot title.")
@click.option("--unit", "quantity_unit", default=None,
              help="Optional unit string for the y-axis (e.g. 'copies/µL').")
@click.option("--standard-task", "standard_task", default="Standard", show_default=True,
              help="Task value identifying standard wells.")
@click.option("--unknown-task", "unknown_task", default="Unknown", show_default=True,
              help="Task value identifying sample (unknown) wells.")
@click.option(
    "--background-rem", "background_rem",
    default=None, is_flag=False, flag_value="control",
    help="Use a named control sample as the background threshold "
         "(case-insensitive). Bare flag defaults to a sample called "
         "`control`; pass an explicit name (e.g. `--background-rem Water`). "
         "Any unknown whose mean Cq exceeds the control's mean Cq is "
         "flagged in `below_background` and dimmed in the bar plot.",
)
@click.option(
    "--task", "extra_tasks", multiple=True, default=(),
    help="Extra Task values to read from the CSV (on top of the standard / "
         "unknown tasks, which are always kept). Pass multiple times to "
         "include several.",
)
@_image_export_options
def standcurve_cmd(
    input_csv: Path,
    target: str | None,
    out_html: Path | None,
    curve_html: Path | None,
    bars_html: Path | None,
    csv_out: Path | None,
    curve_csv: Path | None,
    log_y: bool,
    open_browser: bool,
    title: str | None,
    quantity_unit: str | None,
    standard_task: str,
    unknown_task: str,
    background_rem: str | None,
    extra_tasks: tuple[str, ...],
    pdf_out: str | None,
    png_out: str | None,
    svg_out: str | None,
) -> None:
    """Absolute quantification via a standard curve.

    Reads a Well-Results CSV that includes ``Task == Standard`` wells with
    known ``Quantity`` values, fits ``Cq = slope · log10(Q) + intercept``,
    and back-calculates per-sample quantities for ``Task == Unknown`` wells.
    Emits a combined HTML with the standard curve on top and a bar plot of
    estimated absolute quantities (with propagated error bars) on the bottom.
    """
    out_html = out_html or input_csv.with_suffix(".standcurve.html")
    tasks = tuple(sorted({standard_task, unknown_task, *extra_tasks}))

    steps = ["read", "fit", "quantify", "plot", "write"]
    with click.progressbar(steps, label="qpictures standcurve",
                           item_show_func=lambda s: s or "") as bar:
        bar.update(0)
        df = read_thermo_well_results(input_csv, tasks=tasks)
        bar.update(1)

        try:
            curve = fit_standard_curve(
                df, target=target, standard_task=standard_task
            )
        except ValueError as exc:
            click.echo()
            console.print(f"[bold red]Error:[/bold red] {exc}")
            sys.exit(2)
        bar.update(1)

        try:
            abs_table, _ = absolute_quantification(
                df,
                target=curve.target,
                curve=curve,
                standard_task=standard_task,
                unknown_task=unknown_task,
                background_control=background_rem,
            )
        except ValueError as exc:
            click.echo()
            console.print(f"[bold red]Error:[/bold red] {exc}")
            sys.exit(2)
        bar.update(1)

        fig = standard_curve_with_quantities_figure(
            curve,
            abs_table,
            title=title,
            log_y=log_y,
            quantity_unit=quantity_unit,
            dim_below_background=background_rem is not None,
        )
        curve_fig = standard_curve_plot(curve, unknowns_table=abs_table)
        bar_fig = absolute_quantity_barplot(
            abs_table,
            target=curve.target,
            log_y=log_y,
            quantity_unit=quantity_unit,
            dim_below_background=background_rem is not None,
            curve=curve,
        )
        bar.update(1)

        out_html.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(out_html, include_plotlyjs="cdn", full_html=True)
        if curve_html:
            curve_html.parent.mkdir(parents=True, exist_ok=True)
            curve_fig.write_html(curve_html, include_plotlyjs="cdn", full_html=True)
        if bars_html:
            bars_html.parent.mkdir(parents=True, exist_ok=True)
            bar_fig.write_html(bars_html, include_plotlyjs="cdn", full_html=True)
        if csv_out:
            csv_out.parent.mkdir(parents=True, exist_ok=True)
            abs_table.to_csv(csv_out, index=False)
        if curve_csv:
            curve_csv.parent.mkdir(parents=True, exist_ok=True)
            curve.standards.to_csv(curve_csv, index=False)
        image_paths = _export_static(
            fig, out_html, pdf_out=pdf_out, png_out=png_out, svg_out=svg_out
        )
        bar.update(1)

    # --- Curve summary table (one row) ---
    eff_pct = "" if pd.isna(curve.efficiency) else f"{curve.efficiency * 100:.2f}%"
    curve_summary = pd.DataFrame(
        [{
            "Target": curve.target,
            "slope": curve.slope,
            "slope_SE": curve.slope_se,
            "intercept": curve.intercept,
            "intercept_SE": curve.intercept_se,
            "R^2": curve.r_squared,
            "Efficiency": eff_pct,
            "n_points": curve.n_points,
            "n_levels": curve.n_levels,
            "log10_Q_min": curve.log10_q_min,
            "log10_Q_max": curve.log10_q_max,
            "residual_SD": curve.residual_sd,
        }]
    )
    console.print(
        _df_to_rich_table(curve_summary, title=f"Standard curve fit · {curve.target}")
    )

    # --- Per-sample absolute quantity table ---
    bg_enabled = background_rem is not None
    extrap_count = int(abs_table["extrapolated"].sum())
    n_below = int(abs_table["below_background"].sum()) if bg_enabled else 0

    display_cols = [
        "Sample", "n", "Cq_mean", "Cq_sd", "log10_Q", "log10_Q_sd",
        "Quantity", "Quantity_low", "Quantity_high", "extrapolated",
    ]
    if bg_enabled:
        display_cols.append("below_background")
    display = abs_table[display_cols].copy()
    for col in ("Cq_mean", "Cq_sd", "log10_Q", "log10_Q_sd",
                "Quantity", "Quantity_low", "Quantity_high"):
        display[col] = display[col].astype(float)

    title_suffix_parts: list[str] = []
    if extrap_count:
        title_suffix_parts.append(f"{extrap_count} extrapolated")
    if bg_enabled:
        ctrl_resolved = (
            str(abs_table["Background_control"].iloc[0])
            if "Background_control" in abs_table.columns
            and len(abs_table)
            and abs_table["Background_control"].iloc[0]
            else background_rem
        )
        title_suffix_parts.append(f"background vs {ctrl_resolved}: {n_below} below")
    title_suffix = f"  ({'; '.join(title_suffix_parts)})" if title_suffix_parts else ""
    console.print(
        _df_to_rich_table(
            display,
            title=f"Absolute quantity per sample{title_suffix}",
            dim_column="below_background" if bg_enabled else None,
        )
    )

    console.print(f"[green]Plot written to[/green] {out_html}")
    if curve_html:
        console.print(f"[green]Curve plot written to[/green] {curve_html}")
    if bars_html:
        console.print(f"[green]Bar plot written to[/green] {bars_html}")
    for p in image_paths:
        console.print(f"[green]Plot written to[/green] {p}")
    if csv_out:
        console.print(f"[green]Quantity table written to[/green] {csv_out}")
    if curve_csv:
        console.print(f"[green]Standards table written to[/green] {curve_csv}")

    if open_browser:
        click.launch(str(out_html))


if __name__ == "__main__":  # pragma: no cover
    cli()
