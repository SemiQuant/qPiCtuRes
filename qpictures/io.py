"""Readers for Thermo Fisher (QuantStudio Design & Analysis) qPCR exports.

The "Well Results" CSV that Design & Analysis produces is *mostly* well-behaved,
but sample names are not always quoted, so a value like ``beads,3mm,WS, Slow``
silently inflates the column count. This module is intentionally tolerant of
that quirk.
"""

from __future__ import annotations

import csv
import io as _io
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


def read_thermo_well_results(
    path: str | Path,
    *,
    drop_omitted: bool = True,
    keep_undetermined: bool = False,
    tasks: Iterable[str] | None = ("Unknown",),
) -> pd.DataFrame:
    """Read a Thermo QuantStudio *Well Results* CSV into a tidy DataFrame.

    Parameters
    ----------
    path:
        File path to the exported ``.csv``.
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
        if drop_omitted:
            df = df.loc[~df["Omit"]].copy()

    if tasks is not None and "Task" in df.columns:
        wanted = {t.strip().lower() for t in tasks if t and t.strip()}
        if wanted:
            df = df.loc[df["Task"].str.lower().isin(wanted)].copy()

    for numeric in ("Cq", "Cq Conf", "Amp Score", "Threshold"):
        if numeric in df.columns:
            df[numeric] = df[numeric].map(_to_float)

    if not keep_undetermined and "Cq" in df.columns:
        df = df.loc[df["Cq"].notna()].copy()

    df.reset_index(drop=True, inplace=True)
    return df
