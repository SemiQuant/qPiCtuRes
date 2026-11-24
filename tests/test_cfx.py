"""Bio-Rad CFX (CFX Maestro) export reader."""

from __future__ import annotations

import math
import zipfile
from pathlib import Path

import pytest

import qpictures as q

DATA = Path(__file__).resolve().parents[1] / "data" / "cfx"
CQ = next(DATA.glob("*Quantification Cq Results.csv"))
SUMMARY = next(DATA.glob("*Quantification Summary.csv"))


def _all(path):
    return q.read_thermo_well_results(path, tasks=None, keep_undetermined=True)


def test_semicolon_decimal_comma_export():
    df = _all(CQ)
    assert len(df) == 54
    assert df["Task"].value_counts().to_dict() == {
        "Unknown": 30, "Standard": 20, "Negative Control (NTC)": 4,
    }
    a1 = df.iloc[0]
    assert (a1["Well"], a1["Sample"], a1["Target"]) == ("A1", "Std-01", "SYBR")
    assert a1["Cq"] == pytest.approx(15.44727431)
    assert a1["Quantity"] == pytest.approx(9.68)
    unk = df.loc[df["Task"] == "Unknown"]
    assert unk["Quantity"].isna().all()
    assert set(unk["Sample"]) == {f"Z{i}" for i in range(1, 6)} | {f"T{i}" for i in range(1, 6)}
    assert df.loc[df["Well"] == "B9", "Cq"].isna().all()
    assert df.loc[df["Well"] == "B12", "Cq"].iloc[0] == pytest.approx(39.46495284)


def test_default_filters_keep_unknowns_with_cq():
    df = q.read_thermo_well_results(CQ)
    assert len(df) == 30 and set(df["Task"]) == {"Unknown"}


def test_summary_sheet_matches_cq_results():
    a = _all(CQ)
    with pytest.warns(UserWarning, match="no sample names"):
        b = _all(SUMMARY)
    assert list(a["Well"]) == list(b["Well"])
    assert (a["Cq"] - b["Cq"]).abs().max() < 1e-7


def test_folder_and_zip(tmp_path):
    folder = tmp_path / "export"
    folder.mkdir()
    (folder / SUMMARY.name).write_bytes(SUMMARY.read_bytes())
    (folder / CQ.name).write_bytes(CQ.read_bytes())
    zpath = tmp_path / "export.zip"
    with zipfile.ZipFile(zpath, "w") as z:
        z.write(SUMMARY, f"export/{SUMMARY.name}")
        z.write(CQ, f"export/{CQ.name}")
        z.writestr(f"__MACOSX/export/._{CQ.name}", b"\x00\x01junk")
    ref = _all(CQ)
    for src in (folder, zpath):
        df = _all(src)
        assert list(df["Well"]) == list(ref["Well"])
        assert df["Cq"].equals(ref["Cq"])


def test_standard_curve_matches_cfx():
    # CFX Maestro's own "Standard Curve Results" for this run.
    curve = q.fit_standard_curve(_all(CQ))
    assert curve.target == "SYBR"
    assert curve.slope == pytest.approx(-3.47079129598117, rel=1e-9)
    assert curve.intercept == pytest.approx(18.4409238596741, rel=1e-9)
    assert curve.r_squared == pytest.approx(0.997992912489351, rel=1e-9)
    assert curve.efficiency * 100 == pytest.approx(94.1416586398971, rel=1e-9)


def test_comma_delimited_export(tmp_path):
    text = (
        "Well,Fluor,Target,Content,Sample,Biological Set Name,Cq,Starting Quantity (SQ)\n"
        "A01,FAM,atpE,Unkn-01,\"BB, slow\",,22.5,1.5\n"
        "A02,HEX,,Std-01,,,18.25,1000\n"
        "A03,FAM,atpE,Pos Ctrl,PC,,NaN,NaN\n"
    )
    p = tmp_path / "Quantification Cq Results.csv"
    p.write_text(text)
    df = _all(p)
    assert list(df["Sample"]) == ["BB, slow", "Std-01", "PC"]
    assert list(df["Target"]) == ["atpE", "HEX", "atpE"]
    assert list(df["Task"]) == ["Unknown", "Standard", "Positive Control"]
    assert df["Cq"].iloc[0] == 22.5 and math.isnan(df["Cq"].iloc[2])
    assert df["Quantity"].iloc[1] == 1000


def test_pcrd_gives_clear_error(tmp_path):
    p = tmp_path / "run.pcrd"
    p.write_bytes(b"PK")
    with pytest.raises(ValueError, match="Export All Data Sheets"):
        q.read_thermo_well_results(p)


def test_thermo_csv_still_reads():
    df = q.read_thermo_well_results(DATA.parent / "ex1" / "WellResult-2026-06-09-14570.csv")
    assert len(df) == 44 and "Content" not in df.columns


HEADER = "Well;Fluor;Target;Content;Sample;Biological Set Name;Cq;Starting Quantity (SQ)\n"


def _write(tmp_path, body, header=HEADER, name="Quantification Cq Results.csv", encoding="utf-8"):
    p = tmp_path / name
    p.write_bytes((header + body).encode(encoding))
    return p


def test_unquoted_delimiter_in_sample_name(tmp_path):
    p = _write(tmp_path, "A01;SYBR;;Unkn-01;beads;3mm;;15,44;NaN\n")
    df = _all(p)
    assert df["Sample"].iloc[0] == "beads;3mm" and df["Cq"].iloc[0] == pytest.approx(15.44)


def test_biological_set_kept_apart(tmp_path):
    p = _write(tmp_path, "A01;SYBR;;Unkn-01;P1;Day1;15;\nA02;SYBR;;Unkn-02;P1;Day2;25;\n")
    assert list(_all(p)["Sample"]) == ["P1_Day1", "P1_Day2"]


@pytest.mark.parametrize("encoding", ["cp1252", "utf-16", "utf-8-sig"])
def test_encodings(tmp_path, encoding):
    p = _write(tmp_path, "A01;SYBR;;Unkn-01;Probe \u00e4;;15;\nA02;SYBR;;Unkn-02;Probe \u00f6;;16;\n",
               encoding=encoding)
    assert list(_all(p)["Sample"]) == ["Probe \u00e4", "Probe \u00f6"]


def test_cfx_manager_ct_column(tmp_path):
    p = _write(tmp_path, "A01,FAM,gyrB,Unkn,S1,21.5\n", header="Well,Fluor,Target,Content,Sample,C(t)\n")
    df = _all(p)
    assert df["Cq"].iloc[0] == 21.5 and df["Task"].iloc[0] == "Unknown"


def test_comma_file_with_quoted_decimal_commas(tmp_path):
    p = _write(tmp_path, 'A01,SYBR,,Std-01,,,"15,5","9,68"\n', header=HEADER.replace(";", ","))
    df = _all(p)
    assert df["Cq"].iloc[0] == 15.5 and df["Quantity"].iloc[0] == pytest.approx(9.68)


def test_several_runs_in_a_folder_is_an_error(tmp_path):
    for run in ("runA", "runB"):
        _write(tmp_path, "A01;SYBR;;Unkn-01;S;;15;\n", name=f"{run} -  Quantification Cq Results.csv")
    with pytest.raises(ValueError, match="more than one"):
        _all(tmp_path)


def test_folder_without_cfx_sheet(tmp_path):
    (tmp_path / "notes.csv").write_text("a,b\n1,2\n")
    with pytest.raises(ValueError, match="Export All Data Sheets"):
        _all(tmp_path)


def test_default_output_names(tmp_path, monkeypatch):
    from qpictures.cli import _default_out

    run = tmp_path / "exp 1.5"
    run.mkdir()
    assert _default_out(run, "relative") == tmp_path / "exp 1.5.relative.html"
    monkeypatch.chdir(run)
    assert _default_out(Path("."), "paired") == tmp_path / "exp 1.5.paired.html"
    assert _default_out(Path("x.csv"), "standcurve").name == "x.standcurve.html"
