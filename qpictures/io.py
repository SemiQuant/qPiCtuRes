"""Readers for Thermo Fisher (QuantStudio Design & Analysis) and Bio-Rad CFX qPCR exports.

Three inputs are supported:

* The "Well Results" CSV that Design & Analysis produces. It is *mostly*
  well-behaved, but sample names are not always quoted, so a value like
  ``beads,3mm,WS, Slow`` silently inflates the column count. This module is
  intentionally tolerant of that quirk.
* The ``.eds`` run file itself (a zip of JSON). Sample names, tasks and
  standard quantities come from the plate setup; Cq values from the analysis
  stored in the file, so the run must have been analysed and saved in Design &
  Analysis.
* A Bio-Rad CFX Maestro export (CFX Opus, CFX96, ...): the *Quantification
  Cq Results* CSV, the folder from *Export All Data Sheets*, or a zip of it.
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


_NTC_TASK = r"ntc|negative control|no[ _-]?template"
_NTC_SAMPLE = r"^\s*(?:ntc|no[ _-]?template)"


def _ntc_mask(df: pd.DataFrame) -> pd.Series:
    """Rows that are no-template controls: an NTC Task, or a sample named NTC / no template."""
    mask = pd.Series(False, index=df.index)
    if "Task" in df.columns:
        mask |= df["Task"].astype(str).str.contains(_NTC_TASK, case=False, regex=True)
    if "Sample" in df.columns:
        mask |= df["Sample"].astype(str).str.contains(_NTC_SAMPLE, case=False, regex=True)
    return mask


def _cuttable_mask(df: pd.DataFrame) -> pd.Series:
    """Sample wells a Cq cutoff applies to (not no-template wells, not standards)."""
    mask = ~_ntc_mask(df)
    if "Task" in df.columns:
        mask &= ~df["Task"].astype(str).str.contains("standard", case=False, regex=False)
    return mask


def ntc_check(df: pd.DataFrame) -> pd.DataFrame:
    """Targets whose no-template wells gave a Cq (possible contamination or primer dimers).

    Pass a table read with ``tasks=None`` so the NTC wells are present. Returns one
    row per affected target with ``Target``, ``wells`` (NTC wells), ``amplified``
    (those with a Cq), ``min_cq`` (the earliest NTC Cq) and ``late`` (sample wells
    with a Cq at or after ``min_cq``). Empty when no NTC amplified.
    """
    cols = ["Target", "wells", "amplified", "min_cq", "late"]
    if "Cq" not in df.columns or "Target" not in df.columns:
        return pd.DataFrame(columns=cols)
    ntc = _ntc_mask(df)
    samples = _cuttable_mask(df)
    rows = []
    for target, grp in df.groupby("Target", dropna=False, sort=False):
        n = ntc.loc[grp.index]
        amp = grp.loc[n & grp["Cq"].notna(), "Cq"]
        if amp.empty:
            continue
        smp = grp.loc[samples.loc[grp.index] & grp["Cq"].notna(), "Cq"]
        rows.append({"Target": target, "wells": int(n.sum()), "amplified": int(len(amp)),
                     "min_cq": float(amp.min()), "late": int((smp >= amp.min()).sum())})
    return pd.DataFrame(rows, columns=cols)


def apply_cq_cutoff(df: pd.DataFrame, cutoff: float | None) -> pd.DataFrame:
    """Count sample wells with ``Cq >= cutoff`` as not amplified (Cq becomes NaN).

    Standards and no-template wells are left alone. ``None`` returns ``df`` unchanged.
    """
    if cutoff is None or "Cq" not in df.columns:
        return df
    out = df.copy()
    late = _cuttable_mask(out) & (out["Cq"] >= cutoff)
    out.loc[late, "Cq"] = float("nan")
    return out


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
    filters = dict(drop_omitted=drop_omitted, keep_undetermined=keep_undetermined, tasks=tasks)
    if path.suffix.lower() == ".pcrd":
        raise ValueError(
            f"{path.name} is a CFX Maestro run file, which is encrypted and cannot be "
            "read. In CFX Maestro use Export > Export All Data Sheets (CSV) and give "
            "qpictures the export folder, its zip, or the 'Quantification Cq Results' CSV."
        )
    if path.is_dir() or _is_cfx_zip(path):
        return read_bio_rad_cfx(path, **filters)
    if _is_eds(path):
        return read_thermo_eds(path, **filters)

    text = _decode(path.read_bytes()).splitlines()
    if not text:
        raise ValueError(f"{path} is empty")
    if _looks_like_cfx(text):
        return read_bio_rad_cfx(path, **filters)

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

    if "Cq" not in df.columns:
        raise ValueError(
            f"{path.name} has no Cq column. qpictures reads QuantStudio Well Results "
            "CSVs, .eds run files and Bio-Rad CFX 'Quantification Cq Results' CSVs."
        )
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


# --------------------------------------------------------------------------- #
# Bio-Rad CFX (CFX Maestro) exports
# --------------------------------------------------------------------------- #

# The CFX Maestro sheets we read, best first. "Quantification Cq Results" has
# sample names and starting quantities; "Quantification Summary" has no sample
# names, so it is only used when given directly.
_CFX_SHEET = "Quantification Cq Results"
_CFX_SHEETS = (_CFX_SHEET, "Quantification Summary")

# Header names for the Cq column across CFX Maestro and older CFX Manager.
_CFX_CQ_COLUMNS = ("Cq", "C(t)", "Ct")

# CFX "Content" codes -> the Task labels used elsewhere in qpictures.
_CFX_CONTENT = {
    "unkn": "Unknown",
    "std": "Standard",
    "ntc": "Negative Control (NTC)",
    "pos ctrl": "Positive Control",
    "neg ctrl": "Negative Control",
    "nrt": "No RT Control",
}

_DECIMAL_COMMA = re.compile(r"^-?\d+,\d+$")


def _decode(data: bytes) -> str:
    """Decode an export: UTF-16 or UTF-8 (with or without BOM), else Windows-1252."""
    if data[:2] in (b"\xff\xfe", b"\xfe\xff"):
        return data.decode("utf-16")
    try:
        return data.decode("utf-8-sig")
    except UnicodeDecodeError:
        return data.decode("cp1252", errors="replace")


def _is_cfx_zip(path: Path) -> bool:
    """True if ``path`` is a zip (not an ``.eds``) holding a CFX quantification CSV."""
    if path.suffix.lower() == ".eds" or not zipfile.is_zipfile(path):
        return False
    with zipfile.ZipFile(path) as z:
        return any(_CFX_SHEET.lower() in n.lower() and n.lower().endswith(".csv")
                   for n in z.namelist())


def _split_header(line: str) -> tuple[str, list[str]]:
    delim = ";" if line.count(";") > line.count(",") else "\t" if "\t" in line else ","
    return delim, [c.strip().strip('"').lower() for c in line.split(delim)]


def _is_cfx_header(cols: list[str]) -> bool:
    return {"well", "content"} <= set(cols) and any(c.lower() in cols for c in _CFX_CQ_COLUMNS)


def _looks_like_cfx(lines: list[str]) -> bool:
    """A CFX sheet has ``Well``, ``Content`` and ``Cq`` columns (QuantStudio has no ``Content``)."""
    return any(_is_cfx_header(_split_header(line)[1]) for line in lines[:20])


def _one_sheet(matches: list[str], where: str) -> str | None:
    """Pick the single matching sheet; refuse to choose silently between runs."""
    if len(matches) > 1:
        listed = "\n  ".join(matches)
        raise ValueError(
            f"{where} holds more than one CFX '{_CFX_SHEET}' CSV:\n  {listed}\n"
            "Pass the CSV of the run you want."
        )
    return matches[0] if matches else None


def _cfx_sheet_text(path: Path) -> tuple[str, str]:
    """Return ``(name, text)`` of the CFX quantification sheet at ``path``.

    ``path`` may be the CSV itself, the export folder, or a zip of it. A folder
    or zip must hold exactly one *Quantification Cq Results* CSV (a folder is
    searched at the top level and one level down).
    """
    missing = (
        f"No CFX '{_CFX_SHEET}' CSV found in {path.name}. Export it from CFX "
        "Maestro with Export > Export All Data Sheets (CSV)."
    )
    if path.is_dir():
        files = [p for pattern in ("*.csv", "*/*.csv") for p in path.glob(pattern)
                 if not p.name.startswith("._") and "__macosx" not in str(p).lower()]
        hit = _one_sheet(sorted(str(p) for p in files if _CFX_SHEET.lower() in p.name.lower()),
                         str(path))
        if hit is None:
            raise ValueError(missing)
        return Path(hit).name, _decode(Path(hit).read_bytes())
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as z:
            names = [n for n in z.namelist()
                     if n.lower().endswith(".csv") and not Path(n).name.startswith("._")
                     and "__macosx" not in n.lower()]
            hit = _one_sheet(sorted(n for n in names if _CFX_SHEET.lower() in Path(n).name.lower()),
                             path.name)
            if hit is None:
                raise ValueError(missing)
            return Path(hit).name, _decode(z.read(hit))
    return path.name, _decode(path.read_bytes())


def _cfx_number(value: str, decimal_comma: bool) -> float:
    token = (value or "").strip()
    if decimal_comma and "," in token:
        token = token.replace(".", "").replace(",", ".")
    return _to_float(token)


def _cfx_well(well: str) -> str:
    """``A01`` -> ``A1`` (the well names QuantStudio files use)."""
    m = re.fullmatch(r"\s*([A-Pa-p])0*(\d+)\s*", well or "")
    return f"{m.group(1).upper()}{m.group(2)}" if m else (well or "").strip()


def _cfx_task(content: str) -> str:
    base = re.sub(r"[-\s]*\d+$", "", (content or "").strip()).strip().lower()
    return _CFX_CONTENT.get(base, (content or "").strip())


def read_bio_rad_cfx(
    path: str | Path,
    *,
    drop_omitted: bool = True,
    keep_undetermined: bool = False,
    tasks: Iterable[str] | None = ("Unknown",),
) -> pd.DataFrame:
    """Read a Bio-Rad CFX (CFX Maestro, e.g. CFX Opus) export into a tidy DataFrame.

    ``path`` is the *Quantification Cq Results* CSV, the folder written by
    *Export All Data Sheets*, or a zip of that folder. A folder or zip must
    hold one run. Semicolon-delimited exports with decimal commas (European
    locale settings), UTF-16 and Windows-1252 files are handled. The
    *Quantification Summary* CSV can be read too, but it has no sample names.

    The result has the same columns as :func:`read_thermo_well_results` plus
    ``Quantity`` and ``Content``, so it feeds the relative, paired and
    standard-curve analyses unchanged:

    * ``Content`` codes map to ``Task``: ``Unkn`` -> Unknown, ``Std`` ->
      Standard, ``NTC`` -> Negative Control (NTC), ``Pos Ctrl``, ``Neg Ctrl``
      and ``NRT`` likewise.
    * A blank ``Target`` falls back to the ``Fluor`` (e.g. ``SYBR``); a blank
      ``Sample`` falls back to the ``Content`` code (e.g. ``Std-01``).
    * A non-blank ``Biological Set Name`` is appended to the sample
      (``Sample_Set``), so the same sample in two sets stays two groups and
      ``paired`` can use the set as the condition.
    * ``Quantity`` is the starting quantity CFX holds for standards; it is
      NaN for other wells (CFX's own back-calculated values are ignored).
    * Wells excluded in CFX Maestro are absent from the export, so ``Omit``
      is always ``False``.

    Parameters
    ----------
    path:
        CSV, export folder or zip as above.
    drop_omitted, keep_undetermined, tasks:
        As for :func:`read_thermo_well_results`.
    """
    path = Path(path)
    name, text = _cfx_sheet_text(path)
    lines = text.splitlines()
    header_idx = next(
        (i for i, line in enumerate(lines) if _is_cfx_header(_split_header(line)[1])), None
    )
    if header_idx is None:
        raise ValueError(f"{name} is not a CFX quantification sheet (no Well/Content/Cq header)")
    delim, _ = _split_header(lines[header_idx])
    rows = list(csv.reader(_io.StringIO("\n".join(lines[header_idx:])), delimiter=delim))
    header = [h.strip() for h in rows[0]]
    col = {h.lower(): i for i, h in enumerate(header) if h}
    sample_idx = col.get("sample")

    def get(row: list[str], *names: str) -> str:
        for n in names:
            i = col.get(n.lower())
            if i is not None and i < len(row):
                return row[i].strip()
        return ""

    body: list[list[str]] = []
    for raw in rows[1:]:
        if not raw or all(not c.strip() for c in raw):
            continue
        extra = len(raw) - len(header)
        if extra > 0:
            # An unquoted delimiter inside a sample name: glue it back together.
            if sample_idx is None:
                raise ValueError(f"{name}: a row has more fields than the header: {raw}")
            raw = (raw[:sample_idx] + [delim.join(raw[sample_idx:sample_idx + extra + 1])]
                   + raw[sample_idx + extra + 1:])
        if get(raw, "Well"):
            body.append(raw)

    # Decimal mark: commas when the separator is not a comma, or when a comma file
    # quotes its numbers ("15,44").
    cq_cells = [get(r, *_CFX_CQ_COLUMNS) for r in body]
    decimal_comma = delim != "," or any(_DECIMAL_COMMA.match(c) for c in cq_cells)

    records: list[dict] = []
    for raw in body:
        content = get(raw, "Content")
        task = _cfx_task(content)
        fluor = get(raw, "Fluor")
        sample = get(raw, "Sample") or content
        bioset = get(raw, "Biological Set Name")
        if bioset:
            sample = f"{sample}_{bioset}"
        sq = _cfx_number(get(raw, "Starting Quantity (SQ)", "SQ"), decimal_comma)
        records.append({
            "Well": _cfx_well(get(raw, "Well")),
            "Omit": False,
            "Sample": sample,
            "Target": get(raw, "Target") or fluor,
            "Task": task,
            "Dyes": fluor,
            "Cq": _cfx_number(get(raw, *_CFX_CQ_COLUMNS), decimal_comma),
            "Content": content,
            "Biological Set Name": bioset,
            "Quantity": sq if task == "Standard" and sq > 0 else float("nan"),
            "Well Note": get(raw, "Well Note"),
        })

    df = pd.DataFrame(records)
    if df.empty:
        raise ValueError(f"{name} has no wells")
    named = df["Task"] != "Standard"
    if sample_idx is None or (named.any() and (df.loc[named, "Sample"] == df.loc[named, "Content"]).all()):
        warnings.warn(
            f"{name} has no sample names, so wells are labelled by their CFX Content "
            f"code (e.g. Unkn-01). The '{_CFX_SHEET}' CSV carries the sample names.",
            stacklevel=2,
        )
    return _filter_wells(
        df, drop_omitted=drop_omitted, keep_undetermined=keep_undetermined, tasks=tasks
    )
