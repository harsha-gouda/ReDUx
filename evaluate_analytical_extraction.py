"""Score the analytical metadata columns of our ReDU contribution sheet against ReDU.

Compares plant_redu_contribution_full_<date>.tsv (our extraction) with the ReDU dump
all_sampleinformation.tsv (ReDU curated labels) for the analytical / sample-handling
columns only, file by file, and writes metric tables plus plots.

Join: each of our files is matched to a ReDU record by, in order,
  1. exact USI (case-insensitive),
  2. accession + relative path (ReDU filename without its "f." prefix),
  3. accession + filename stem, only where the stem is unique on both sides.

Per cell outcome (one file x one column):
  exact      both filled, identical after whitespace/case normalisation
  partial    both filled, differ, but agree at a coarser level (same instrument family,
             same ionisation source, same chromatography mode, same solvents, ...)
  conflict   both filled and disagree
  missed     ReDU filled, ours empty
  added      ours filled, ReDU empty (cannot be verified against ReDU)
  both_empty neither filled

Metrics per column:
  precision         exact / both filled
  lenient precision (exact + partial) / both filled
  recall            exact / ReDU filled
  F1                harmonic mean of precision and recall
  deposit macro     mean of per-deposit precision (files within a deposit are not independent)

Usage:
  python evaluate_analytical_extraction.py
  python evaluate_analytical_extraction.py --ours <tsv> --redu <tsv> --out <dir> --fields MassSpectrometer SampleType
"""

import argparse
import csv
import os
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_OURS = PROJECT / "Claude outputs" / "plant_redu_contribution_full_2026-10-21.tsv"
DEFAULT_REDU = PROJECT / "ReDU extraction" / "ReDU data" / "all_sampleinformation.tsv"
DEFAULT_OUT = PROJECT / "Evaluvation" / "analytical_accuracy_2026-10-21"

ANALYTICAL_FIELDS = [
    "SampleType",
    "SampleTypeSub1",
    "SampleCollectionMethod",
    "SampleExtractionMethod",
    "InternalStandardsUsed",
    "MassSpectrometer",
    "IonizationSourceAndPolarity",
    "ChromatographyAndPhase",
]

EMPTY = {"", "missing value", "not applicable", "not specified", "nan", "none specified", "unknown", "na", "n/a"}

OUTCOMES = ["exact", "partial", "conflict", "missed", "added", "both_empty"]
OUTCOME_LABELS = {
    "exact": "Exact match",
    "partial": "Partial match",
    "conflict": "Conflict",
    "missed": "Missed (ReDU only)",
    "added": "Added (ours only)",
    "both_empty": "Both empty",
}
# Reference categorical palette (dataviz skill), fixed order; both_empty is neutral.
OUTCOME_COLORS = {
    "exact": "#2a78d6",
    "partial": "#1baf7a",
    "conflict": "#e34948",
    "missed": "#eb6834",
    "added": "#eda100",
    "both_empty": "#d6d5d0",
}
METRIC_COLORS = ["#2a78d6", "#1baf7a", "#eb6834", "#6250d6"]
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"


# ---------------------------------------------------------------- normalisation

def norm(v):
    v = "" if v is None else str(v)
    return re.sub(r"\s+", " ", v.strip()).lower()


def is_empty(v):
    return norm(v) in EMPTY


def _before_paren(v):
    return norm(v).split("(")[0].strip()


def _instrument(v):
    name = norm(v).split("|")[0].replace("+", "").strip()
    return name


GENERIC_INSTRUMENT_TOKENS = {"q", "xevo", "orbitrap", "maxis", "impact", "6500", "linear"}


def coarse_mass_spec(a, b):
    a, b = _instrument(a), _instrument(b)
    if a == b or a in b or b in a:
        return True
    ta, tb = a.split(), b.split()
    return bool(ta and tb and ta[0] == tb[0] and ta[0] not in GENERIC_INSTRUMENT_TOKENS)


def coarse_solvents(a, b):
    solv = lambda v: frozenset(s.strip() for s in _before_paren(v).split("-") if s.strip())
    return solv(a) == solv(b)


def coarse_collection(a, b):
    return norm(a).split(",")[0].strip() == norm(b).split(",")[0].strip()


PLANT_DERIVED = {"plant", "food", "beverage"}


def coarse_sample_type(a, b):
    # ReDU types most plant-derived foods as "plant"; our sheet uses "food"/"beverage".
    return norm(a) in PLANT_DERIVED and norm(b) in PLANT_DERIVED


def coarse_sub1(a, b):
    return norm(a).split("_")[0] == norm(b).split("_")[0]


def coarse_standards(a, b):
    sa = {s.strip() for s in re.split(r"[;,]", norm(a)) if s.strip()}
    sb = {s.strip() for s in re.split(r"[;,]", norm(b)) if s.strip()}
    return bool(sa & sb)


def coarse_ionisation(a, b):
    # Same source counts as partial only when polarity is not a hard positive/negative flip.
    if _before_paren(a) != _before_paren(b):
        return False
    pa, pb = norm(a), norm(b)
    flip = {"(positive)", "(negative)"}
    return not any(p in pa for p in flip) or not any(p in pb for p in flip)


COARSE = {
    "SampleType": coarse_sample_type,
    "SampleTypeSub1": coarse_sub1,
    "SampleCollectionMethod": coarse_collection,
    "SampleExtractionMethod": coarse_solvents,
    "InternalStandardsUsed": coarse_standards,
    "MassSpectrometer": coarse_mass_spec,
    "IonizationSourceAndPolarity": coarse_ionisation,
    "ChromatographyAndPhase": lambda a, b: _before_paren(a) == _before_paren(b),
}


# ---------------------------------------------------------------- adjusted tier
#
# For both-filled cells that are not an exact match, say *why* they differ. Classes:
#   convention     same thing under a different labelling convention
#   specificity    one side is the generic / NOS / family form of the other (not wrong, less precise)
#   ours_supported independent evidence in the record (filename, ReDU's own other columns) backs ours
#   conflict       a real disagreement
# Adjusted precision counts exact + convention + specificity + ours_supported as correct.

ADJ_CORRECT = ["exact", "convention", "specificity", "ours_supported"]
ADJ_CLASSES = ADJ_CORRECT + ["conflict"]
ADJ_LABELS = {
    "exact": "Exact match",
    "convention": "Labelling convention",
    "specificity": "Specificity difference",
    "ours_supported": "Ours backed by evidence",
    "conflict": "Real conflict",
}
ADJ_COLORS = {
    "exact": "#2a78d6",
    "convention": "#86b6ef",
    "specificity": "#1baf7a",
    "ours_supported": "#6250d6",
    "conflict": "#e34948",
}

# PSI-MS instrument classes and the model-name prefixes they cover.
GENERIC_INSTRUMENT_CLASSES = {
    "linear ion trap": ("ltq",),
    "triple quadrupole": ("tsq", "qtrap", "xevo tq", "triple quad"),
    "time-of-flight": ("tof", "tripletof", "xevo g2", "synapt"),
}


def _paren(v):
    m = re.search(r"\(([^)]*)\)", norm(v))
    return m.group(1).strip() if m else ""


def _filename_polarity(usi):
    """Explicit polarity words in the path only; single letters (N_, P_) are too ambiguous."""
    s = str(usi).lower()
    pos = re.search(r"(?<![a-z])(pos|positive)(?![a-z])", s)
    neg = re.search(r"(?<![a-z])(neg|negative)(?![a-z])", s)
    if pos and not neg:
        return "positive"
    if neg and not pos:
        return "negative"
    return ""


# Manually verified ground truth, (deposit, field) -> [(regex on the USI file path, value or
# {captured code: value}, source)]. Where it applies it overrides ReDU as the reference: a mismatch
# is "ours_supported" when the truth is ours, "conflict" when it is not, and an "added" cell (ReDU
# empty) whose value contradicts it is flagged as verified-wrong. All checked 2026-10-21.
REPOSITORY_TRUTH = {
    ("MSV000097043", "SampleExtractionMethod"): [
        # PlaE02_..., Ara50_01_..., LabW13_...: E = ethanol, 50 = 50% ethanol, W = water
        (r"(?:^|/)[A-Za-z]{3}(E|50|W)[^/]*$",
         {"E": "ethanol (100%)", "50": "ethanol-water (1:1)", "W": "water (100%) (deg_C_NOS)"},
         "filename solvent code"),
    ],
    ("MTBLS2382", "ChromatographyAndPhase"): [
        # Polar fraction on HSS T3 C18, lipid fraction on a C8 column
        (r"/Polar_(pos|neg)/", "reverse phase (C18)", "MetaboLights Chromatography protocol"),
        (r"/Lipid_pos/", "reverse phase (C8)", "MetaboLights Chromatography protocol"),
    ],
    ("ST002169", "IonizationSourceAndPolarity"): [
        # mwTab SUBJECT_SAMPLE_FACTORS: P_* files Treatment:Positive (AN003553), N_* Treatment:Negative (AN003554)
        (r"(?:^|/)P_[^/]*$", "electrospray ionization (positive)", "Workbench AN003553/AN003554 sample factors"),
        (r"(?:^|/)N_[^/]*$", "electrospray ionization (negative)", "Workbench AN003553/AN003554 sample factors"),
    ],
}


def repository_truth(usi, field):
    parts = str(usi).split(":", 2)
    if len(parts) < 3:
        return ""
    for pattern, value, _source in REPOSITORY_TRUTH.get((parts[1], field), []):
        m = re.search(pattern, parts[2], flags=re.IGNORECASE)
        if m:
            if isinstance(value, dict):
                code = m.group(1)
                return norm(value.get(code.upper() if code.isalpha() else code, ""))
            return norm(value)
    return ""


def adjudicate(field, ours, redu, usi="", redu_chrom=""):
    a, b = norm(ours), norm(redu)
    truth = repository_truth(usi, field)
    if truth:
        return "ours_supported" if truth == a and truth != b else "conflict"
    if field == "SampleType":
        if a in PLANT_DERIVED and b in PLANT_DERIVED:
            return "convention"
        return "conflict"
    if field == "SampleTypeSub1":
        pa, pb = a.rsplit("_", 1)[0], b.rsplit("_", 1)[0]
        if pa == pb and ("nos" in (a.rsplit("_", 1)[-1], b.rsplit("_", 1)[-1])):
            return "specificity"
        return "conflict"
    if field == "MassSpectrometer":
        ia, ib = _instrument(ours), _instrument(redu)
        # GC-Orbitrap vs its LC namesake: ReDU's own chromatography column decides.
        if re.search(r"\bgc\b", ia) and not re.search(r"\bgc\b", ib):
            return "ours_supported" if norm(redu_chrom).startswith("gas chromatography") else "conflict"
        if ia in ib or ib in ia:
            return "specificity"
        for cls, prefixes in GENERIC_INSTRUMENT_CLASSES.items():
            if (ia == cls and ib.startswith(prefixes)) or (ib == cls and ia.startswith(prefixes)):
                return "specificity"
        return "conflict"
    if field == "IonizationSourceAndPolarity":
        if _before_paren(a) != _before_paren(b):
            return "conflict"
        fp = _filename_polarity(usi)
        if fp and fp == _paren(a) and fp != _paren(b):
            return "ours_supported"
        if "alternating" in (_paren(a), _paren(b)) or "nos" in (_paren(a), _paren(b)):
            return "specificity"
        return "conflict"
    if field == "ChromatographyAndPhase":
        if _before_paren(a) == _before_paren(b) and "nos" in (_paren(a), _paren(b)):
            return "specificity"
        return "conflict"
    if field == "SampleExtractionMethod":
        if coarse_solvents(ours, redu) and (not _paren(a) or not _paren(b) or "nos" in (_paren(a) + _paren(b))):
            return "specificity"
        return "conflict"
    if field == "SampleCollectionMethod":
        # "liquid" vs "liquid, solid phase extraction (C18)": one is the parent of the other
        if a.startswith(b) or b.startswith(a):
            return "specificity"
        return "conflict"
    if field == "InternalStandardsUsed":
        return "specificity" if coarse_standards(ours, redu) else "conflict"
    return "conflict"


def outcome(field, ours, redu):
    eo, er = is_empty(ours), is_empty(redu)
    if eo and er:
        return "both_empty"
    if eo:
        return "missed"
    if er:
        return "added"
    if norm(ours) == norm(redu):
        return "exact"
    if COARSE.get(field, lambda a, b: False)(ours, redu):
        return "partial"
    return "conflict"


# ---------------------------------------------------------------- loading + join

def load_redu(path, accessions, fields):
    """Stream the ReDU dump (quoted fields contain newlines) keeping only our deposits."""
    csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
    keep = ["ATTRIBUTE_DatasetAccession", "filename", "USI"] + fields
    rows = []
    with open(path, encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        missing = [c for c in keep if c not in header]
        if missing:
            raise SystemExit(f"ReDU dump lacks columns: {missing}")
        idx = [header.index(c) for c in keep]
        acc_i = header.index("ATTRIBUTE_DatasetAccession")
        for row in reader:
            if len(row) == len(header) and row[acc_i] in accessions:
                rows.append([row[i] for i in idx])
    return pd.DataFrame(rows, columns=keep)


def _stem(p):
    return os.path.splitext(os.path.basename(str(p).replace("\\", "/")))[0].lower()


def _usi_path(usi):
    parts = str(usi).split(":", 2)
    return parts[2].lower() if len(parts) == 3 else ""


def join(ours, redu):
    ours = ours.copy()
    redu = redu.copy()
    ours["_row"] = np.arange(len(ours))
    redu["_rid"] = np.arange(len(redu))

    ours["k_usi"] = ours["USI"].str.lower()
    redu["k_usi"] = redu["USI"].str.lower()
    ours["k_path"] = ours["MassiveID"] + "|" + ours["USI"].map(_usi_path)
    redu["k_path"] = redu["ATTRIBUTE_DatasetAccession"] + "|" + redu["filename"].str.lower().str.replace(r"^f\.", "", regex=True)
    ours["k_stem"] = ours["MassiveID"] + "|" + ours["filename"].map(_stem)
    redu["k_stem"] = redu["ATTRIBUTE_DatasetAccession"] + "|" + redu["filename"].map(_stem)

    matched = []
    used_o, used_r = set(), set()
    for key, method in [("k_usi", "usi"), ("k_path", "path"), ("k_stem", "stem")]:
        o = ours[~ours["_row"].isin(used_o)]
        r = redu[~redu["_rid"].isin(used_r)]
        o = o[~o[key].duplicated(keep=False)]
        r = r[~r[key].duplicated(keep=False)]
        m = o[["_row", key]].merge(r[["_rid", key]], on=key)
        m["join_method"] = method
        matched.append(m[["_row", "_rid", "join_method"]])
        used_o.update(m["_row"])
        used_r.update(m["_rid"])
    pairs = pd.concat(matched, ignore_index=True)
    return pairs


# ---------------------------------------------------------------- metrics

def score(ours, redu, pairs, fields):
    o = ours.iloc[pairs["_row"].to_numpy()].reset_index(drop=True)
    r = redu.iloc[pairs["_rid"].to_numpy()].reset_index(drop=True)
    long = []
    redu_chrom = r["ChromatographyAndPhase"] if "ChromatographyAndPhase" in r else pd.Series("", index=r.index)
    for f in fields:
        oc = [outcome(f, a, b) for a, b in zip(o[f], r[f])]
        adj = [
            "exact" if c == "exact"
            else adjudicate(f, a, b, u, ch) if c in ("partial", "conflict")
            else ""
            for c, a, b, u, ch in zip(oc, o[f], r[f], o["USI"], redu_chrom)
        ]
        truth = [repository_truth(u, f) if not is_empty(a) else "" for u, a in zip(o["USI"], o[f])]
        verified = ["" if not t else "ours_right" if norm(a) == t else "ours_wrong" for t, a in zip(truth, o[f])]
        long.append(pd.DataFrame({
            "MassiveID": o["MassiveID"],
            "filename": o["filename"],
            "USI": o["USI"],
            "join_method": pairs["join_method"],
            "field": f,
            "ours": o[f],
            "redu": r[f],
            "outcome": oc,
            "adjusted": adj,
            "repository_verified": verified,
        }))
    return pd.concat(long, ignore_index=True)


def summarise(cells, group):
    counts = (cells.groupby(group + ["outcome"]).size().unstack("outcome", fill_value=0)
              .reindex(columns=OUTCOMES, fill_value=0))
    both = counts["exact"] + counts["partial"] + counts["conflict"]
    redu_filled = both + counts["missed"]
    ours_filled = both + counts["added"]
    out = counts.copy()
    out["n_files"] = counts.sum(axis=1)
    out["both_filled"] = both
    out["redu_filled"] = redu_filled
    out["ours_filled"] = ours_filled
    out["precision"] = counts["exact"] / both.replace(0, np.nan)
    out["lenient_precision"] = (counts["exact"] + counts["partial"]) / both.replace(0, np.nan)
    out["recall"] = counts["exact"] / redu_filled.replace(0, np.nan)
    p, rc = out["precision"], out["recall"]
    out["f1"] = 2 * p * rc / (p + rc)
    out["coverage_vs_redu"] = both / redu_filled.replace(0, np.nan)

    adj = (cells[cells["adjusted"] != ""].groupby(group + ["adjusted"]).size()
           .unstack("adjusted", fill_value=0).reindex(columns=ADJ_CLASSES, fill_value=0)
           .reindex(out.index, fill_value=0))
    adj_ok = adj[ADJ_CORRECT].sum(axis=1)
    for c in ADJ_CLASSES:
        out[f"adj_{c}"] = adj[c]
    out["adjusted_precision"] = adj_ok / both.replace(0, np.nan)
    out["adjusted_recall"] = adj_ok / redu_filled.replace(0, np.nan)
    p, rc = out["adjusted_precision"], out["adjusted_recall"]
    out["adjusted_f1"] = 2 * p * rc / (p + rc)
    # "added" cells (ReDU empty) that the repository record shows to be wrong: invisible to ReDU scoring
    out["added_verified_wrong"] = (cells[(cells["outcome"] == "added") & (cells["repository_verified"] == "ours_wrong")]
                                   .groupby(group).size().reindex(out.index, fill_value=0))
    return out.reset_index()


# ---------------------------------------------------------------- plots

def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=TEXT_2, labelsize=9)
    ax.xaxis.label.set_color(TEXT_2)
    ax.yaxis.label.set_color(TEXT_2)


def plot_adjusted(fm, out):
    fm = fm.sort_values("adjusted_precision")
    cols = [f"adj_{c}" for c in ADJ_CLASSES]
    frac = fm[cols].div(fm["both_filled"].replace(0, np.nan), axis=0).fillna(0) * 100
    y = np.arange(len(fm))
    fig, ax = plt.subplots(figsize=(10, 0.55 * len(fm) + 1.8))
    left = np.zeros(len(fm))
    for c in ADJ_CLASSES:
        v = frac[f"adj_{c}"].to_numpy()
        ax.barh(y, v, left=left, color=ADJ_COLORS[c], label=ADJ_LABELS[c], edgecolor="white", linewidth=1.5, height=0.7)
        for yy, l, w in zip(y, left, v):
            if w >= 6:
                ax.text(l + w / 2, yy, f"{w:.0f}", ha="center", va="center", fontsize=7.5,
                        color="white" if c in ("exact", "ours_supported", "conflict") else TEXT)
        left += v
    for yy, (_, row) in zip(y, fm.iterrows()):
        ax.text(101, yy, f"{row['adjusted_precision'] * 100:.1f}% correct  (n={int(row['both_filled']):,})",
                va="center", fontsize=8, color=TEXT_2)
    ax.set_yticks(y)
    ax.set_yticklabels(fm["field"], color=TEXT)
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of files where both sides are filled")
    _style(ax)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=5, frameon=False, fontsize=8)
    ax.set_title("Adjusted: why filled values differ from ReDU", loc="left", color=TEXT, fontsize=12, pad=26)
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)
    fig.savefig(out / "fig5_adjusted_breakdown.png", dpi=200)
    plt.close(fig)


def plot_metrics(fm, out):
    metrics = [("precision", "Precision, strict"), ("adjusted_precision", "Precision, adjusted"),
               ("recall", "Recall, strict"), ("adjusted_recall", "Recall, adjusted")]
    fm = fm.sort_values("adjusted_f1")
    y = np.arange(len(fm))
    h = 0.19
    fig, ax = plt.subplots(figsize=(9, 0.75 * len(fm) + 1.6))
    for i, ((col, label), color) in enumerate(zip(metrics, METRIC_COLORS)):
        vals = fm[col].fillna(0).to_numpy() * 100
        ax.barh(y + (1.5 - i) * h, vals, height=h - 0.03, color=color, label=label, edgecolor="white", linewidth=0.5)
        if col.startswith("adjusted"):
            for yy, v in zip(y + (1.5 - i) * h, vals):
                ax.text(v + 1, yy, f"{v:.0f}%", va="center", fontsize=7.5, color=TEXT)
    ax.set_yticks(y)
    ax.set_yticklabels(fm["field"], color=TEXT)
    ax.set_xlim(0, 108)
    ax.set_xlabel("% of cells  (adjusted = exact + convention + specificity + evidence-backed)")
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    _style(ax)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4, frameon=False, fontsize=8.5)
    ax.set_title("Extraction accuracy vs ReDU, analytical columns", loc="left", color=TEXT, fontsize=12, pad=28)
    fig.tight_layout()
    fig.savefig(out / "fig1_precision_recall_f1.png", dpi=200)
    plt.close(fig)


def plot_outcomes(fm, out):
    fm = fm.sort_values("precision")
    frac = fm[OUTCOMES].div(fm["n_files"], axis=0) * 100
    y = np.arange(len(fm))
    fig, ax = plt.subplots(figsize=(10, 0.55 * len(fm) + 1.8))
    left = np.zeros(len(fm))
    for oc in OUTCOMES:
        v = frac[oc].to_numpy()
        ax.barh(y, v, left=left, color=OUTCOME_COLORS[oc], label=OUTCOME_LABELS[oc],
                edgecolor="white", linewidth=1.5, height=0.7)
        for yy, l, w in zip(y, left, v):
            if w >= 6:
                ax.text(l + w / 2, yy, f"{w:.0f}", ha="center", va="center", fontsize=7.5,
                        color="white" if oc in ("exact", "conflict") else TEXT)
        left += v
    ax.set_yticks(y)
    ax.set_yticklabels(fm["field"], color=TEXT)
    ax.set_xlim(0, 100)
    ax.set_xlabel(f"% of joined files (n = {int(fm['n_files'].iloc[0]):,})")
    _style(ax)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=6, frameon=False, fontsize=8)
    ax.set_title("Where each column's cells land", loc="left", color=TEXT, fontsize=12, pad=26)
    fig.tight_layout()
    fig.savefig(out / "fig2_outcome_breakdown.png", dpi=200)
    plt.close(fig)


def plot_per_deposit(dm, fm, out):
    order = fm.sort_values("precision")["field"].tolist()
    dm = dm[dm["both_filled"] > 0]
    fig, ax = plt.subplots(figsize=(9, 0.6 * len(order) + 1.6))
    rng = np.random.default_rng(0)
    for i, f in enumerate(order):
        d = dm[dm["field"] == f]
        if d.empty:
            continue
        jitter = rng.uniform(-0.22, 0.22, len(d))
        size = 12 + 60 * np.sqrt(d["both_filled"] / dm["both_filled"].max())
        ax.scatter(d["precision"] * 100, i + jitter, s=size, color="#2a78d6", alpha=0.45,
                   edgecolor="white", linewidth=0.6)
        med = d["precision"].median() * 100
        ax.plot([med, med], [i - 0.32, i + 0.32], color=TEXT, linewidth=2)
        ax.text(102, i, f"median {med:.0f}%  (n={len(d)})", va="center", fontsize=8, color=TEXT_2)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, color=TEXT)
    ax.set_xlim(-3, 101)
    ax.set_xlabel("Per-deposit precision (exact), % of both-filled cells  -  dot size = files compared")
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    _style(ax)
    ax.set_title("Accuracy spread across deposits", loc="left", color=TEXT, fontsize=12)
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)
    fig.savefig(out / "fig3_per_deposit_precision.png", dpi=200)
    plt.close(fig)


def plot_confusion(cells, field, out, top=10):
    d = cells[(cells["field"] == field) & cells["outcome"].isin(["exact", "partial", "conflict"])]
    if d.empty:
        return
    lab = lambda s: s.str.split("|").str[0]
    ours, redu = lab(d["ours"]), lab(d["redu"])
    keep_o = ours.value_counts().index[:top]
    keep_r = redu.value_counts().index[:top]
    ours = ours.where(ours.isin(keep_o), "other")
    redu = redu.where(redu.isin(keep_r), "other")
    ct = pd.crosstab(ours, redu)
    rownorm = ct.div(ct.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(1.0 + 0.62 * ct.shape[1] + 2.5, 1.0 + 0.45 * ct.shape[0] + 1.5))
    ax.imshow(rownorm.to_numpy(), cmap=matplotlib.colors.LinearSegmentedColormap.from_list(
        "seq", ["#f4f8fd", "#86b6ef", "#2a78d6", "#0d366b"]), vmin=0, vmax=1, aspect="auto")
    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            n = ct.iat[i, j]
            if n:
                ax.text(j, i, f"{n:,}", ha="center", va="center", fontsize=7,
                        color="white" if rownorm.iat[i, j] > 0.5 else TEXT)
    ax.set_xticks(range(ct.shape[1]))
    ax.set_xticklabels(ct.columns, rotation=40, ha="right", fontsize=8, color=TEXT)
    ax.set_yticks(range(ct.shape[0]))
    ax.set_yticklabels(ct.index, fontsize=8, color=TEXT)
    ax.set_xlabel("ReDU (curated)", color=TEXT_2)
    ax.set_ylabel("Our extraction", color=TEXT_2)
    ax.set_title(f"{field}: files where both sides are filled\n(cell = file count, shade = row share)",
                 loc="left", color=TEXT, fontsize=10)
    for s in ax.spines.values():
        s.set_visible(False)
    fig.tight_layout()
    fig.savefig(out / f"fig4_confusion_{field}.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ours", type=Path, default=DEFAULT_OURS)
    ap.add_argument("--redu", type=Path, default=DEFAULT_REDU)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--fields", nargs="+", default=ANALYTICAL_FIELDS)
    ap.add_argument("--exclude-blanks", action="store_true",
                    help="drop files we typed blank_QC before scoring")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    ours = pd.read_csv(args.ours, sep="\t", dtype=str, keep_default_na=False)
    if args.exclude_blanks:
        ours = ours[ours["SampleType"] != "blank_QC"].reset_index(drop=True)
    print(f"ours: {len(ours):,} files / {ours['MassiveID'].nunique()} deposits")

    redu = load_redu(args.redu, set(ours["MassiveID"]), args.fields)
    print(f"ReDU: {len(redu):,} records for {redu['ATTRIBUTE_DatasetAccession'].nunique()} of our deposits")

    pairs = join(ours, redu)
    print(f"joined: {len(pairs):,} files / {ours.iloc[pairs['_row']]['MassiveID'].nunique()} deposits "
          f"({pairs['join_method'].value_counts().to_dict()})")

    cells = score(ours, redu, pairs, args.fields)
    fm = summarise(cells, ["field"])
    fm["field"] = pd.Categorical(fm["field"], args.fields)
    fm = fm.sort_values("field").reset_index(drop=True)
    fm["field"] = fm["field"].astype(str)
    dm = summarise(cells, ["field", "MassiveID"])
    macro = dm.groupby("field")[["precision", "adjusted_precision"]].mean().rename(columns={
        "precision": "deposit_macro_precision", "adjusted_precision": "deposit_macro_adjusted_precision"})
    fm = fm.merge(macro, on="field", how="left")

    # Recall with whole-deposit abstentions set aside: deposits where ReDU fills the column and we
    # leave every file empty. Those are one missing-protocol decision, not thousands of file errors.
    abst = dm[(dm["ours_filled"] == 0) & (dm["redu_filled"] > 0)]
    abst_cells = abst.groupby("field")["redu_filled"].sum()
    fm["abstained_deposits"] = fm["field"].map(abst.groupby("field")["MassiveID"].apply(";".join)).fillna("")
    fm["adjusted_recall_excl_abstained"] = (
        fm[[f"adj_{c}" for c in ADJ_CORRECT]].sum(axis=1)
        / (fm["redu_filled"] - fm["field"].map(abst_cells).fillna(0)).replace(0, np.nan))

    overall = summarise(cells.assign(all="ALL"), ["all"]).rename(columns={"all": "field"})
    overall["deposit_macro_precision"] = dm.groupby("MassiveID")[["exact", "partial", "conflict"]].sum().pipe(
        lambda x: (x["exact"] / x.sum(axis=1).replace(0, np.nan)).mean())
    overall["adjusted_recall_excl_abstained"] = (
        overall[[f"adj_{c}" for c in ADJ_CORRECT]].sum(axis=1) / (overall["redu_filled"] - abst_cells.sum()))
    fm_out = pd.concat([fm, overall], ignore_index=True)

    adjudicated = (cells[cells["adjusted"].isin(ADJ_CLASSES) & (cells["adjusted"] != "exact")]
                   .groupby(["field", "adjusted", "MassiveID", "ours", "redu"]).size()
                   .rename("n_files").reset_index()
                   .sort_values(["adjusted", "n_files"], ascending=[True, False]))

    disagreements = (cells[cells["outcome"].isin(["conflict", "partial"])]
                     .groupby(["field", "outcome", "ours", "redu"])
                     .agg(n_files=("filename", "size"), n_deposits=("MassiveID", "nunique"),
                          deposits=("MassiveID", lambda s: ";".join(sorted(set(s))[:10])))
                     .reset_index().sort_values(["field", "n_files"], ascending=[True, False]))

    fm_out.to_csv(args.out / "field_metrics.csv", index=False)
    dm.to_csv(args.out / "deposit_metrics.csv", index=False)
    disagreements.to_csv(args.out / "disagreements.csv", index=False)
    adjudicated.to_csv(args.out / "adjudicated_mismatches.csv", index=False)
    adjudicated[adjudicated["adjusted"] == "conflict"].to_csv(args.out / "residual_conflicts.csv", index=False)
    cells.to_csv(args.out / "cell_level_comparison.tsv.gz", sep="\t", index=False)

    plot_metrics(fm, args.out)
    plot_adjusted(fm, args.out)
    plot_outcomes(fm, args.out)
    plot_per_deposit(dm, fm, args.out)
    for f in ("SampleType", "MassSpectrometer", "IonizationSourceAndPolarity",
              "ChromatographyAndPhase", "SampleExtractionMethod"):
        if f in args.fields:
            plot_confusion(cells, f, args.out)

    cols = ["field", "both_filled", "redu_filled", "precision", "adjusted_precision", "recall",
            "adjusted_recall", "adjusted_recall_excl_abstained", "adjusted_f1", "added_verified_wrong"]
    with pd.option_context("display.width", 200, "display.float_format", "{:.3f}".format):
        print(fm_out[cols].to_string(index=False))
    print("\nadjudicated mismatches (files):")
    print(cells[cells["adjusted"].isin(ADJ_CLASSES[1:])].groupby(["field", "adjusted"]).size().unstack(fill_value=0))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
