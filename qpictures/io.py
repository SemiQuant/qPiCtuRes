"""Readers for Thermo Fisher (QuantStudio Design & Analysis) qPCR exports.

Two inputs are supported:

* The "Well Results" CSV that Design & Analysis produces. It is *mostly*
  well-behaved, but sample names are not always quoted, so a value like
  ``beads,3mm,WS, Slow`` silently inflates the column count. This module is
  intentionally tolerant of that quirk.
* The ``.eds`` run file itself (a zip of JSON). Sample names, tasks and
  standard quantities come from the plate setup; Cq values from the analysis
  stored in the file, so the run must have been analysed and saved in Design &
  Analysis.
"""

from __future__ import annotations

import csv
import io as _io
import json
import math
import re
import warnings
import zipfile
from pathlib import Path
from typing import Iterable

import pandas as pd

# Canonical column names we expect from a Design & Analysis "Well Results" export.
EXPECTED_COLUMNS: tuple[str, ...] = (
    "Well",
    "Omit",
    "Sample",
    "Target",
    "Task",
    "Dyes",
    "Cq",
    "Cq Conf",
    "Amp Score",
    "Amp Status",
    "Annotated",
    "Threshold",
    "Baseline Start",
    "Baseline End",
    "Curve Quality",
    "Result Quality Issues",
)

# Values used by Design & Analysis to indicate "no Cq".
_UNDETERMINED_TOKENS = {"undetermined", "no ct", "nan", ""}


def _find_header_row(lines: list[str]) -> int:
    """Return the 0-indexed line where the header begins.

    Design & Analysis often prepends a metadata block (``# Block Type:`` etc.)
    before the actual table. We scan for the first line that contains ``Well``
    and ``Cq``.
    """
    for idx, line in enumerate(lines):
        lowered = line.lower()
        if "well" in lowered and "cq" in lowered and "sample" in lowered:
            return idx
    return 0


def _coalesce_sample_columns(row: list[str], n_expected: int) -> list[str]:
    """Glue an over-long row back together by merging extra commas into ``Sample``.

    The Sample column sits at index 2 of the header. If we see ``len(row) > n``
    we assume the extra fields came from unquoted commas inside Sample, so we
    fold them back in.
    """
    extra = len(row) - n_expected
    if extra <= 0:
        return row
    sample_idx = EXPECTED_COLUMNS.index("Sample")
    merged_sample = ",".join(row[sample_idx : sample_idx + 1 + extra]).strip()
    return row[:sample_idx] + [merged_sample] + row[sample_idx + 1 + extra :]


def _to_float(value: str) -> float:
    """Convert a Cq-style cell to a float, mapping ``Undetermined`` to NaN."""
    if value is None:
        return float("nan")
    token = str(value).strip()
    if token.lower() in _UNDETERMINED_TOKENS:
        return float("nan")
    try:
        return float(token)
    except ValueError:
        return float("nan")


def _is_eds(path: Path) -> bool:
    return path.suffix.lower() == ".eds" or zipfile.is_zipfile(path)


def _filter_wells(
    df: pd.DataFrame,
    *,
    drop_omitted: bool,
    keep_undetermined: bool,
    tasks: Iterable[str] | None,
) -> pd.DataFrame:
    """Apply the Omit / Task / Undetermined filters shared by every reader."""
    if drop_omitted and "Omit" in df.columns:
        df = df.loc[~df["Omit"]].copy()

    if tasks is not None and "Task" in df.columns:
        wanted = {t.strip().lower() for t in tasks if t and t.strip()}
        if wanted:
            df = df.loc[df["Task"].str.lower().isin(wanted)].copy()

    if not keep_undetermined and "Cq" in df.columns:
        df = df.loc[df["Cq"].notna()].copy()

    df.reset_index(drop=True, inplace=True)
    return df


def read_thermo_well_results(
    path: str | Path,
    *,
    drop_omitted: bool = True,
    keep_undetermined: bool = False,
    tasks: Iterable[str] | None = ("Unknown",),
) -> pd.DataFrame:
    """Read a Thermo QuantStudio *Well Results* CSV into a tidy DataFrame.

    An ``.eds`` run file is also accepted and handed to
    :func:`read_thermo_eds`, so every caller works with either input.

    Parameters
    ----------
    path:
        File path to the exported ``.csv`` (or the ``.eds`` run file).
    drop_omitted:
        If ``True`` (default), rows with ``Omit == true`` are removed.
    keep_undetermined:
        If ``False`` (default), wells with no Cq (``Undetermined``) are dropped
        *after* being parsed. Set ``True`` to keep them as NaN rows.
    tasks:
        Iterable of ``Task`` values to keep (case-insensitive). Default
        ``("Unknown",)``, which drops NTC/Standard/etc. control wells. Pass
        ``None`` to disable this filter.

    Returns
    -------
    pandas.DataFrame
        Columns: all canonical Thermo columns plus typed ``Cq``, ``Cq Conf``,
        ``Amp Score``, ``Threshold`` (floats) and a boolean ``Omit``.
    """
    path = Path(path)
    if _is_eds(path):
        return read_thermo_eds(
            path,
            drop_omitted=drop_omitted,
            keep_undetermined=keep_undetermined,
            tasks=tasks,
        )

    text = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    if not text:
        raise ValueError(f"{path} is empty")

    header_idx = _find_header_row(text)
    payload = "\n".join(text[header_idx:])

    reader = csv.reader(_io.StringIO(payload))
    rows: Iterable[list[str]] = list(reader)
    header = [h.strip() for h in rows[0]]
    n_expected = len(header)

    cleaned: list[list[str]] = []
    for raw in rows[1:]:
        if not raw or all((cell or "").strip() == "" for cell in raw):
            continue
        cleaned.append(_coalesce_sample_columns(raw, n_expected))

    df = pd.DataFrame(cleaned, columns=header)

    # Normalize whitespace in string-like columns.
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].astype(str).str.strip()

    if "Omit" in df.columns:
        df["Omit"] = df["Omit"].str.lower().eq("true")

    for numeric in ("Cq", "Cq Conf", "Amp Score", "Threshold"):
        if numeric in df.columns:
            df[numeric] = df[numeric].map(_to_float)

    return _filter_wells(
        df, drop_omitted=drop_omitted, keep_undetermined=keep_undetermined, tasks=tasks
    )


# --------------------------------------------------------------------------- #
# .eds run files
# --------------------------------------------------------------------------- #

_PLATE_ROWS = "ABCDEFGHIJKLMNOP"

# .eds enums -> the labels Design & Analysis writes in its Well Results CSV, so
# `--task` filters and plots behave the same whichever file is read.
_EDS_TASKS = {
    "UNKNOWN": "Unknown",
    "STANDARD": "Standard",
    "NTC": "Negative Control (NTC)",
    "POSITIVE_CONTROL": "Positive Control",
}
_EDS_AMP_STATUS = {"AMP": "Amp", "NO_AMP": "No Amp", "INCONCLUSIVE": "Inconclusive", "NA": "N/A"}


def _eds_label(value: str | None, mapping: dict[str, str]) -> str:
    if not value:
        return ""
    return mapping.get(value, value.replace("_", " ").title())


def _eds_float(value) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return out if math.isfinite(out) else float("nan")


def read_thermo_eds(
    path: str | Path,
    *,
    drop_omitted: bool = True,
    keep_undetermined: bool = False,
    tasks: Iterable[str] | None = ("Unknown",),
) -> pd.DataFrame:
    """Read a QuantStudio Design & Analysis ``.eds`` run file into a tidy DataFrame.

    The result has the same columns as :func:`read_thermo_well_results` plus
    ``Quantity`` (the standard quantity from the plate setup), so it feeds the
    relative, paired and standard-curve analyses unchanged.

    Sample names, targets, tasks and quantities come from the current plate
    setup; Cq values and amplification calls from the analysis stored in the
    file. A well whose stored result has no matching target in the plate setup
    is skipped. A Cq of ``-1`` (the file's "no Cq") becomes NaN.

    Parameters
    ----------
    path:
        File path to the ``.eds`` run file.
    drop_omitted, keep_undetermined, tasks:
        As for :func:`read_thermo_well_results`.

    Raises
    ------
    ValueError
        If the file is not a Design & Analysis ``.eds`` or holds no stored
        analysis (the run was never analysed, or was saved by the older
        QuantStudio Real-Time PCR software).
    """
    path = Path(path)
    try:
        z = zipfile.ZipFile(path)
    except zipfile.BadZipFile as exc:
        raise ValueError(f"{path} is not a valid .eds file (not a zip archive)") from exc

    with z:
        names = set(z.namelist())
        if "setup/plate_setup.json" not in names:
            raise ValueError(
                f"{path} has no setup/plate_setup.json. Only .eds files saved by "
                "QuantStudio Design & Analysis are supported; export a Well "
                "Results CSV instead."
            )
        if "primary/analysis_result.json" not in names:
            raise ValueError(
                f"{path} holds no stored analysis. Open it in Design & Analysis, "
                "analyse and save, then try again (or export a Well Results CSV)."
            )

        def load(name: str) -> dict:
            return json.loads(z.read(name).decode("utf-8")) if name in names else {}

        summary = load("summary.json")
        setup = load("setup/plate_setup.json")
        result = load("primary/analysis_result.json")

    status = ((summary.get("analysis") or {}).get("primary") or {}).get("status", "")
    if status and status.upper() != "COMPLETED":
        warnings.warn(
            f"{path.name}: the stored analysis is out of date (status {status}). "
            "The plate setup changed after it was analysed; re-analyse and save in "
            "Design & Analysis for current Cq values.",
            stacklevel=2,
        )

    block = summary.get("blockType") or setup.get("blockType") or ""
    # 384-well blocks are 16 x 24, 48-well blocks 6 x 8, everything else 96-well 8 x 12.
    ncols = 24 if "384" in block else 8 if re.search(r"(^|[^0-9])48W", block, re.I) else 12
    dyes = {
        t.get("name", ""): "-".join(x for x in (t.get("reporter"), t.get("quencher")) if x)
        for t in setup.get("targets", [])
    }
    plate = {int(w["index"]): w for w in setup.get("wells", [])}

    rows: list[dict] = []
    for well in result.get("wellResults", []):
        idx = int(well["wellIndex"])
        setup_well = plate.get(idx)
        if setup_well is None:
            continue
        assigned = {t["targetName"]: t for t in setup_well.get("targetAssignments", [])}
        for reaction in well.get("reactionResults", []):
            target = reaction.get("targetName", "")
            if target not in assigned:
                continue
            amp = reaction.get("amplificationResult") or {}
            cq = _eds_float(amp.get("cq"))
            quantity = _eds_float(assigned[target].get("quantity"))
            rows.append({
                "Well": f"{_PLATE_ROWS[idx // ncols]}{idx % ncols + 1}",
                "Omit": bool(reaction.get("omitted", False)),
                "Sample": str(setup_well.get("sampleName") or well.get("sampleName") or "").strip(),
                "Target": target,
                "Task": _eds_label(assigned[target].get("task") or reaction.get("task"), _EDS_TASKS),
                "Dyes": dyes.get(target, ""),
                "Cq": cq if cq > 0 else float("nan"),
                "Cq Conf": _eds_float(amp.get("cqConf")),
                "Amp Score": _eds_float(amp.get("ampScore")),
                "Amp Status": _eds_label(amp.get("ampStatus"), _EDS_AMP_STATUS),
                "Threshold": _eds_float(amp.get("ctThreshold")),
                "Baseline Start": amp.get("ctBaselineStart"),
                "Baseline End": amp.get("ctBaselineEnd"),
                "Result Quality Issues": "; ".join(map(str, reaction.get("resultQCIssues") or [])),
                "Quantity": quantity if quantity > 0 else float("nan"),
                "_idx": idx,
            })

    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError(f"{path} has no analysed wells with a target in the plate setup")
    df = df.sort_values(["_idx", "Target"], kind="stable").drop(columns="_idx")

    return _filter_wells(
        df, drop_omitted=drop_omitted, keep_undetermined=keep_undetermined, tasks=tasks
    )
