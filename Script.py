#!/usr/bin/env python3
"""
================================================================================
SOYBEAN PHENOLOGY QA & KPI REPORT —  DC PHN GS Soy BR 26
================================================================================

Production-ready agronomic QA system for soybean phenology data.

Author : João Canal- Agronomic Data analyst 
Version: 3.1
"""

# ==================== IMPORTS ====================
import os
import re
import sys
import unicodedata
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for server / CI
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from reportlab.lib import colors as rl_colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.units import inch, cm
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    Image, PageBreak, KeepTogether, ListFlowable, ListItem,
)

warnings.filterwarnings("ignore", category=UserWarning)

# ==================== 1. CONFIGURATION ====================

# ---- File paths (edit per season) ----
INPUT_TF_PATH = (path1)
SEEDING_PATH = (path2)
BASE_OUTPUT_DIR = (path3)

# ---- Protocol filter (the ONLY protocol processed) ----
PROTOCOL_FILTER = "DC PHN GS Soy BR 26"

# ---- Farm filter (season scope) ----
FARM_FILTER = "Trials 2025/26"

# ---- Brand ----
BRAND_COLOR = "Grey"
DARK_GREY = "#4A4A4A"
LIGHT_GREY = "#CCCCCC"
ALERT_RED = "#D32F2F"
WARNING_YELLOW = "#FFC107"
SUCCESS_GREEN = "#4CAF50"

PDF_TITLE = "Soybean Phenology QA & KPI Report — Brazil — Season 25/26"

# ---- Growth-stage ordered scale ----
GS_ORDER: List[str] = [
    "VS", "VE", "VC",
    "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9",
    "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8",
]
GS_INDEX: Dict[str, int] = {gs: i for i, gs in enumerate(GS_ORDER)}

GS_PHASES: Dict[str, List[str]] = {
    "Early Vegetative":   ["VS", "VE", "VC", "V1", "V2", "V3"],
    "Late Vegetative":    ["V4", "V5", "V6", "V7", "V8", "V9"],
    "Early Reproductive": ["R1", "R2", "R3", "R4"],
    "Late Reproductive":  ["R5", "R6", "R7", "R8"],
}
GS_PHASE_ORDER = list(GS_PHASES.keys())

# ---- Development-gap threshold (business rule: > 2 stages = flag) ----
DEV_GAP_THRESHOLD = 2

# ---- Required assessment types for package completeness ----
REQUIRED_ASSESSMENTS = [
    "SEV_CANOPYLOW_PLOT",
    "SEV_CANOPYMID_PLOT",
    "SEV_CANOPY-UP_PLOT",
    "ROWCLOSING_PLANT_0-9_PLOT",
]
DISEASE_ASSESSMENTS = [
    "SEV_CANOPYLOW_PLOT",
    "SEV_CANOPYMID_PLOT",
    "SEV_CANOPY-UP_PLOT",
]
RC_ASSESSMENT = "ROWCLOSING_PLANT_0-9_PLOT"

# ---- Responsible names expected in Field strings ----
RESPONSIBLE_NAMES: List[str] = [
    "Name1", "Name2", "Name3", "Name4", "Name5", "name6",
    "Name7", "name8",
]

# ---- Brazilian state codes for extraction ----
UF_CODES = [
    "MT", "MS", "GO", "SP", "PR", "RS", "BA", "MG", "MA",
    "PI", "SC", "DF", "RO", "TO", "PA",
]

# ---- KPI weights ----
KPI_WEIGHTS = {
    "gs_trio_completeness": 0.30,
    "assessment_pkg_completeness": 0.20,
    "gs_coverage": 0.20,
    "rc_completion": 0.15,
    "inverse_issue_burden": 0.15,
}

# ---- BUA columns (Basic Unit of Analysis) ----
BUA_COLS = [
    "Trialnumber", "Replicate", "Sampleno",
    "Factorlevel", "Locationcode", "Observationdate",
]
# Trajectory key: same biological trajectory over time
TRAJECTORY_KEY = [
    "Trialnumber", "Replicate", "Sampleno",
    "Factorlevel", "Locationcode",
]


# ==================== 2. TEXT & PARSING HELPERS ====================

def _strip_accents(s: str) -> str:
    """Remove diacritics from a unicode string."""
    if pd.isna(s):
        return ""
    nfkd = unicodedata.normalize("NFKD", str(s))
    return "".join(ch for ch in nfkd if not unicodedata.combining(ch))


def _norm(s: str) -> str:
    """Uppercase, strip accents, collapse whitespace."""
    return re.sub(r"\s+", " ", _strip_accents(s).upper().strip())


def normalize_gs(val) -> Optional[str]:
    """Normalize a GS code to upper, stripped. Returns None if empty."""
    if pd.isna(val):
        return None
    t = str(val).strip().upper()
    return t if t else None


def gs_numeric(gs) -> float:
    """Growth-stage to numeric index; NaN if unknown."""
    if gs is None or pd.isna(gs):
        return np.nan
    return GS_INDEX.get(str(gs).strip().upper(), np.nan)


def get_gs_phase(gs) -> str:
    if gs is None or pd.isna(gs):
        return "Unknown"
    gs_upper = str(gs).strip().upper()
    for phase, stages in GS_PHASES.items():
        if gs_upper in stages:
            return phase
    return "Unknown"


# ---- Responsible extraction from Field string ----

# Build a list of (compiled_regex, canonical_name) for each expected name.
# We match the name with or without accents, using word-boundary logic.
_RESPONSIBLE_PATTERNS: List[Tuple[re.Pattern, str]] = []
for _name in RESPONSIBLE_NAMES:
    _ascii = _strip_accents(_name).upper()
    _original_upper = _name.upper()
    # Build alternatives: original accented + ascii
    alts = {_ascii}
    if _original_upper != _ascii:
        alts.add(_original_upper)
    alt_str = "|".join(re.escape(a) for a in alts)
    # Match with flexible boundaries (space, dash, bracket, start/end)
    pat = re.compile(
        rf"(?:^|[\s\-_\[\]\(\)])({alt_str})(?:[\s\-_\[\]\(\)]|$)",
        re.IGNORECASE,
    )
    _RESPONSIBLE_PATTERNS.append((pat, _name))


def extract_responsible(field: str) -> str:
    """
    Parse the responsible person's name directly from the Field string.
    Returns the canonical name or 'UNASSIGNED' if nothing matches.
    """
    if pd.isna(field):
        return "UNASSIGNED"
    t = _norm(field)
    for pat, canonical in _RESPONSIBLE_PATTERNS:
        if pat.search(t):
            return canonical
    return "UNASSIGNED"


# ---- State extraction ----

def extract_state(field: str) -> str:
    """
    Extract Brazilian state (UF) from Field string using robust
    word-boundary matching. Returns 'UNK' if not found.

    Priority:
      1. 'MATO GROSSO DO SUL' → MS
      2. 'MATO GROSSO'        → MT
      3. Standard UF codes with word-boundary regex
    """
    if pd.isna(field):
        return "UNK"
    t = _norm(field)

    # ---- Explicit long-form guards ----
    if "MATO GROSSO DO SUL" in t or "MATOGROSSODOSUL" in t:
        return "MS"
    if re.search(r"(?:^|[\s\-_\[\]])MS(?:[\s\-_\[\]]|$)", t):
        return "MS"
    if "MATO GROSSO" in t or "MATOGROSSO" in t:
        return "MT"
    if re.search(r"(?:^|[\s\-_\[\]])MT(?:[\s\-_\[\]]|$)", t):
        return "MT"

    # ---- Other UFs ----
    for uf in UF_CODES:
        if uf in ("MT", "MS"):
            continue
        if re.search(rf"(?:^|[\s\-_\[\]]){uf}(?:[\s\-_\[\]]|$)", t):
            return uf

    # Long-form fallbacks
    _LONG = {
        "PARANA": "PR", "GOIAS": "GO", "BAHIA": "BA",
        "RIO GRANDE DO SUL": "RS", "MINAS GERAIS": "MG",
        "SAO PAULO": "SP", "SANTA CATARINA": "SC",
        "MARANHAO": "MA", "PIAUI": "PI", "TOCANTINS": "TO",
        "RONDONIA": "RO", "DISTRITO FEDERAL": "DF",
    }
    for long_name, uf in _LONG.items():
        if long_name in t:
            return uf

    return "UNK"


# ==================== 3. DATA LOADING ====================

def _norm_col(c: str) -> str:
    """Normalize an Excel column header."""
    c = str(c).strip().replace("\n", " ").replace("\r", " ")
    return re.sub(r"\s+", " ", c)


def load_data() -> pd.DataFrame:
    """
    Load the TF data file, apply farm + protocol filters, merge seeding
    dates, extract responsible/state, normalize GS columns.
    """
    print("=" * 80)
    print(f"  Loading data from: {INPUT_TF_PATH}")
    print("=" * 80)

    df = pd.read_excel(INPUT_TF_PATH)
    df.columns = [_norm_col(c) for c in df.columns]
    print(f"  Raw rows loaded: {len(df):,}")

    # Rename common alternative column names
    rename_map = {
        "FIELD Value": "Field",
        "FARM Value": "Farm",
        "Protocol Name Value": "Protocolname",
        "Protocol Type Value": "Protocolype",
        "Trial No": "Trialnumber",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

    # ---- Farm filter ----
    if "Farm" in df.columns:
        df = df[df["Farm"] == FARM_FILTER].copy()
        print(f"  After Farm filter '{FARM_FILTER}': {len(df):,}")

    # ---- Protocol filter (THE critical filter) ----
    proto_col = None
    for candidate in ("Protocolname", "Protocolype", "Protocol"):
        if candidate in df.columns:
            proto_col = candidate
            break
    if proto_col is None:
        print("  WARNING: No protocol column found. Proceeding with all data.")
    else:
        df = df[df[proto_col] == PROTOCOL_FILTER].copy()
        print(f"  After Protocol filter '{PROTOCOL_FILTER}': {len(df):,}")

    if df.empty:
        print("  ERROR: No data remaining after filters. Check protocol name.")
        sys.exit(1)

    # ---- Parse dates ----
    df["Observationdate"] = pd.to_datetime(df["Observationdate"], errors="coerce")

    # ---- Normalize GS columns ----
    for col in ("Cropstagemincode", "Cropstagemajcode", "Cropstagemaxcode"):
        if col in df.columns:
            df[col] = df[col].apply(normalize_gs)

    # ---- Ensure BUA columns exist (fill with placeholders if missing) ----
    for col in BUA_COLS:
        if col not in df.columns:
            df[col] = "N/A"

    # ---- Numeric Originalvalue ----
    if "Originalvalue" in df.columns:
        df["OriginalvalueNum"] = pd.to_numeric(df["Originalvalue"], errors="coerce")

    # ---- Extract Responsible & State ----
    df["Responsible"] = df["Field"].apply(extract_responsible)
    df["State"] = df["Field"].apply(extract_state)

    # ---- Merge seeding dates (optional) ----
    try:
        df_seed = pd.read_excel(SEEDING_PATH)
        df_seed.columns = [_norm_col(c) for c in df_seed.columns]
        if "Trial No" in df_seed.columns:
            df_seed = df_seed.rename(columns={"Trial No": "Trialnumber"})
        if "Seeding date Value" in df_seed.columns:
            df_seed["Seeding date Value"] = pd.to_datetime(
                df_seed["Seeding date Value"], errors="coerce"
            )
            seed_map = df_seed.set_index("Trialnumber")["Seeding date Value"].to_dict()
            df["SeedingDate"] = df["Trialnumber"].map(seed_map)
            matched = df["SeedingDate"].notna().sum()
            print(f"  Seeding dates matched: {matched:,} rows")
    except Exception as exc:
        print(f"  WARNING: Could not load seeding file: {exc}")
        df["SeedingDate"] = pd.NaT

    # ---- Compute DAP (Days After Planting) ----
    if "SeedingDate" in df.columns:
        df["DAP"] = (df["Observationdate"] - df["SeedingDate"]).dt.days
        valid_dap = df["DAP"].notna().sum()
        print(f"  DAP computed: {valid_dap:,} rows with valid DAP")
    else:
        df["DAP"] = np.nan

    # ---- GS Phase ----
    if "Cropstagemajcode" in df.columns:
        df["GS_Phase"] = df["Cropstagemajcode"].apply(get_gs_phase)

    # ---- Summaries ----
    print(f"\n  Final dataset: {len(df):,} rows")
    print(f"  Unique Fields:       {df['Field'].nunique()}")
    print(f"  Unique Trials:       {df['Trialnumber'].nunique()}")
    print(f"  Responsibles:        {sorted(df['Responsible'].unique())}")
    print(f"  States:              {sorted(df['State'].unique())}")

    # ---- Validation console summary ----
    n_unassigned = (df["Responsible"] == "UNASSIGNED").sum()
    n_unk = (df["State"] == "UNK").sum()
    unassigned_fields = df.loc[df["Responsible"] == "UNASSIGNED", "Field"].unique()
    unk_fields = df.loc[df["State"] == "UNK", "Field"].unique()
    print(f"\n  VALIDATION:")
    print(f"    UNASSIGNED responsible rows: {n_unassigned}  ({len(unassigned_fields)} fields)")
    for f in unassigned_fields[:20]:
        print(f"      -> {f}")
    if len(unassigned_fields) > 20:
        print(f"      ... and {len(unassigned_fields) - 20} more")
    print(f"    UNK state rows: {n_unk}  ({len(unk_fields)} fields)")
    for f in unk_fields[:20]:
        print(f"      -> {f}")

    return df


# ==================== 4. QA CHECKS ====================

# ---------- 4.1 GS Trio Completeness ----------

def check_gs_trio_completeness(df: pd.DataFrame) -> pd.DataFrame:
    """
    For every BUA, verify presence of Cropstagemincode, Cropstagemajcode,
    Cropstagemaxcode. Flag any record missing one or more.
    """
    print("\n  [QA 4.1] GS Trio Completeness ...")
    gs_cols = ["Cropstagemincode", "Cropstagemajcode", "Cropstagemaxcode"]
    present = [c for c in gs_cols if c in df.columns]
    if not present:
        print("    WARNING: No GS columns found.")
        return pd.DataFrame()

    issues = []
    for idx, row in df.iterrows():
        missing = [c for c in present if row.get(c) is None or pd.isna(row.get(c))]
        if missing:
            issues.append({
                **{c: row.get(c, "") for c in BUA_COLS},
                "Field": row.get("Field", ""),
                "Responsible": row.get("Responsible", ""),
                "State": row.get("State", ""),
                "MissingColumns": ", ".join(missing),
                "Issue": "Incomplete GS Trio",
            })
    result = pd.DataFrame(issues)
    print(f"    Found {len(result)} records with incomplete GS trio")
    return result


# ---------- 4.2 GS Triplet Internal Consistency ----------

def check_gs_triplet_consistency(df: pd.DataFrame) -> pd.DataFrame:
    """
    Verify min <= maj <= max within each BUA.
    """
    print("\n  [QA 4.2] GS Triplet Internal Consistency ...")
    needed = ["Cropstagemincode", "Cropstagemajcode", "Cropstagemaxcode"]
    if not all(c in df.columns for c in needed):
        return pd.DataFrame()

    issues = []
    for idx, row in df.iterrows():
        mn = gs_numeric(row["Cropstagemincode"])
        mj = gs_numeric(row["Cropstagemajcode"])
        mx = gs_numeric(row["Cropstagemaxcode"])
        if np.isnan(mn) or np.isnan(mj) or np.isnan(mx):
            continue
        problems = []
        if mn > mj:
            problems.append(f"min({row['Cropstagemincode']})>maj({row['Cropstagemajcode']})")
        if mj > mx:
            problems.append(f"maj({row['Cropstagemajcode']})>max({row['Cropstagemaxcode']})")
        if mn > mx:
            problems.append(f"min({row['Cropstagemincode']})>max({row['Cropstagemaxcode']})")
        if problems:
            issues.append({
                **{c: row.get(c, "") for c in BUA_COLS},
                "Field": row.get("Field", ""),
                "Responsible": row.get("Responsible", ""),
                "State": row.get("State", ""),
                "Cropstagemincode": row["Cropstagemincode"],
                "Cropstagemajcode": row["Cropstagemajcode"],
                "Cropstagemaxcode": row["Cropstagemaxcode"],
                "ConsistencyIssue": "; ".join(problems),
            })
    result = pd.DataFrame(issues)
    print(f"    Found {len(result)} triplet consistency issues")
    return result


# ---------- 4.3 GS Regression ----------

def check_gs_regression(df: pd.DataFrame) -> pd.DataFrame:
    """
    Check GS regression separately for min, maj, max along each trajectory.
    """
    print("\n  [QA 4.3] GS Regression ...")
    gs_cols_check = {
        "Cropstagemincode": "min",
        "Cropstagemajcode": "maj",
        "Cropstagemaxcode": "max",
    }
    present_cols = {k: v for k, v in gs_cols_check.items() if k in df.columns}
    if not present_cols:
        return pd.DataFrame()

    issues = []
    for traj_key, grp in df.groupby(TRAJECTORY_KEY, dropna=False):
        grp = grp.sort_values("Observationdate").reset_index(drop=True)
        for gs_col, label in present_cols.items():
            prev_gs = None
            prev_date = None
            prev_num = np.nan
            for i, row in grp.iterrows():
                cur = row.get(gs_col)
                cur_num = gs_numeric(cur)
                if np.isnan(cur_num):
                    continue
                if not np.isnan(prev_num) and cur_num < prev_num:
                    issues.append({
                        **{c: row.get(c, "") for c in BUA_COLS},
                        "Field": row.get("Field", ""),
                        "Responsible": row.get("Responsible", ""),
                        "State": row.get("State", ""),
                        "GS_Field": label,
                        "PrevDate": prev_date,
                        "PrevGS": prev_gs,
                        "CurrDate": row["Observationdate"],
                        "CurrGS": cur,
                        "Issue": f"Regression in {label}: {prev_gs} -> {cur}",
                    })
                prev_gs = cur
                prev_date = row["Observationdate"]
                prev_num = cur_num

    result = pd.DataFrame(issues)
    print(f"    Found {len(result)} GS regression events")
    return result


# ---------- 4.4 Development Gaps (GS jump > threshold) ----------

def check_development_gaps(df: pd.DataFrame) -> pd.DataFrame:
    """
    Flag jumps > DEV_GAP_THRESHOLD stages in min, maj, or max along each trajectory.
    """
    print(f"\n  [QA 4.4] Development Gaps (jump > {DEV_GAP_THRESHOLD} stages) ...")
    gs_cols_check = {
        "Cropstagemincode": "min",
        "Cropstagemajcode": "maj",
        "Cropstagemaxcode": "max",
    }
    present_cols = {k: v for k, v in gs_cols_check.items() if k in df.columns}
    if not present_cols:
        return pd.DataFrame()

    issues = []
    for traj_key, grp in df.groupby(TRAJECTORY_KEY, dropna=False):
        grp = grp.sort_values("Observationdate").reset_index(drop=True)
        for gs_col, label in present_cols.items():
            prev_gs = None
            prev_date = None
            prev_num = np.nan
            for i, row in grp.iterrows():
                cur = row.get(gs_col)
                cur_num = gs_numeric(cur)
                if np.isnan(cur_num):
                    continue
                if not np.isnan(prev_num):
                    jump = cur_num - prev_num
                    if jump > DEV_GAP_THRESHOLD:
                        issues.append({
                            **{c: row.get(c, "") for c in BUA_COLS},
                            "Field": row.get("Field", ""),
                            "Responsible": row.get("Responsible", ""),
                            "State": row.get("State", ""),
                            "GS_Field": label,
                            "PrevDate": prev_date,
                            "PrevGS": prev_gs,
                            "CurrDate": row["Observationdate"],
                            "CurrGS": cur,
                            "StagesJumped": int(jump),
                            "Issue": f"Jump {label}: {prev_gs} -> {cur} ({int(jump)} stages)",
                        })
                prev_gs = cur
                prev_date = row["Observationdate"]
                prev_num = cur_num

    result = pd.DataFrame(issues)
    print(f"    Found {len(result)} development gap flags")
    return result


# ---------- helper: find second RC=9 date per trajectory ----------

def _compute_second_rc9(df: pd.DataFrame) -> Dict[tuple, Optional[pd.Timestamp]]:
    """
    For each trajectory, find the Observationdate of the SECOND occurrence
    of RC=9. Returns dict: trajectory_key_tuple -> date or None.
    """
    rc = df[df["Assessmenttypecode"] == RC_ASSESSMENT].copy()
    if rc.empty:
        return {}
    rc["RC_Val"] = pd.to_numeric(rc["Originalvalue"], errors="coerce")
    rc9 = rc[rc["RC_Val"] == 9].sort_values("Observationdate")

    result = {}
    for traj, grp in rc9.groupby(TRAJECTORY_KEY, dropna=False):
        key = tuple(grp.iloc[0][c] for c in TRAJECTORY_KEY)
        dates = grp["Observationdate"].tolist()
        result[key] = dates[1] if len(dates) >= 2 else None
    return result


# ---------- 4.5 Assessment Package Completeness ----------

def check_assessment_pkg_completeness(
    df: pd.DataFrame,
    second_rc9_map: Dict[tuple, Optional[pd.Timestamp]],
) -> pd.DataFrame:
    """
    Four assessments required per BUA:
      SEV_CANOPYLOW_PLOT, SEV_CANOPYMID_PLOT, SEV_CANOPY-UP_PLOT, ROWCLOSING_PLANT_0-9_PLOT

    Disease assessments are only obligatory BEFORE the second RC=9 in the trajectory.
    """
    print("\n  [QA 4.5] Assessment Package Completeness ...")
    issues = []

    # Group by BUA
    bua_groups = df.groupby(BUA_COLS, dropna=False)

    for bua_key, grp in bua_groups:
        bua_dict = dict(zip(BUA_COLS, bua_key))
        traj_tuple = tuple(bua_dict.get(c, "N/A") for c in TRAJECTORY_KEY)
        obs_date = bua_dict.get("Observationdate")
        second_rc9 = second_rc9_map.get(traj_tuple)

        present = set(grp["Assessmenttypecode"].dropna().unique())

        # Determine which assessments are required
        required = [RC_ASSESSMENT]  # RC always required
        disease_required = True
        if second_rc9 is not None and pd.notna(obs_date) and obs_date > second_rc9:
            disease_required = False
        if disease_required:
            required.extend(DISEASE_ASSESSMENTS)

        missing = [a for a in required if a not in present]
        if missing:
            field = grp["Field"].iloc[0] if "Field" in grp.columns else ""
            responsible = grp["Responsible"].iloc[0] if "Responsible" in grp.columns else ""
            state = grp["State"].iloc[0] if "State" in grp.columns else ""
            issues.append({
                **bua_dict,
                "Field": field,
                "Responsible": responsible,
                "State": state,
                "MissingAssessments": ", ".join(missing),
                "DiseaseRequired": disease_required,
                "Issue": f"Missing: {', '.join(missing)}",
            })

    result = pd.DataFrame(issues)
    print(f"    Found {len(result)} BUAs with incomplete assessment packages")
    return result


# ---------- 4.6 / 4.7 Disease Monotonicity ----------

def check_disease_monotonicity(
    df: pd.DataFrame,
    second_rc9_map: Dict[tuple, Optional[pd.Timestamp]],
) -> pd.DataFrame:
    """
    Check disease monotonicity per trajectory+Pestcode+Assessmenttypecode.
    Pre-second-RC=9: mandatory.
    Post-second-RC=9: only flag if values are actually present.
    """
    print("\n  [QA 4.6/4.7] Disease Monotonicity ...")
    disease_df = df[df["Assessmenttypecode"].isin(DISEASE_ASSESSMENTS)].copy()
    if disease_df.empty or "OriginalvalueNum" not in disease_df.columns:
        return pd.DataFrame()

    issues = []
    group_cols = TRAJECTORY_KEY + ["Pestcode", "Assessmenttypecode"]
    for gkey, grp in disease_df.groupby(group_cols, dropna=False):
        grp = grp.sort_values("Observationdate").reset_index(drop=True)
        traj_tuple = tuple(grp.iloc[0][c] for c in TRAJECTORY_KEY)
        second_rc9 = second_rc9_map.get(traj_tuple)

        prev_val = np.nan
        prev_date = None
        for _, row in grp.iterrows():
            val = row["OriginalvalueNum"]
            obs = row["Observationdate"]
            if np.isnan(val):
                continue
            # If post-second-RC=9, only check if value present (no false positive for absence)
            if not np.isnan(prev_val) and val < prev_val:
                issues.append({
                    **{c: row.get(c, "") for c in BUA_COLS},
                    "Field": row.get("Field", ""),
                    "Responsible": row.get("Responsible", ""),
                    "State": row.get("State", ""),
                    "Assessmenttypecode": row.get("Assessmenttypecode", ""),
                    "Pestcode": row.get("Pestcode", ""),
                    "PrevDate": prev_date,
                    "PrevValue": prev_val,
                    "CurrDate": obs,
                    "CurrValue": val,
                    "PostSecondRC9": (
                        "Yes" if (second_rc9 is not None and pd.notna(obs) and obs > second_rc9) else "No"
                    ),
                    "Issue": f"Disease decreased: {prev_val} -> {val}",
                })
            prev_val = val
            prev_date = obs

    result = pd.DataFrame(issues)
    print(f"    Found {len(result)} disease monotonicity violations")
    return result


# ---------- 4.8 RC Monotonicity & Progress ----------

def check_rc_monotonicity_and_progress(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns (rc_monotonicity_issues, rc_backlog).
    """
    print("\n  [QA 4.8] RC Monotonicity & Progress ...")
    rc_df = df[df["Assessmenttypecode"] == RC_ASSESSMENT].copy()
    if rc_df.empty:
        return pd.DataFrame(), pd.DataFrame()
    rc_df["RC_Val"] = pd.to_numeric(rc_df["Originalvalue"], errors="coerce")

    mono_issues = []
    backlog_rows = []
    today = pd.Timestamp.now().normalize()

    for traj_key, grp in rc_df.groupby(TRAJECTORY_KEY, dropna=False):
        grp = grp.sort_values("Observationdate").reset_index(drop=True)
        prev_val = np.nan
        prev_date = None
        first_rc9 = None
        second_rc9 = None
        rc9_count = 0

        for _, row in grp.iterrows():
            val = row["RC_Val"]
            obs = row["Observationdate"]
            if np.isnan(val):
                continue
            # Track RC=9 occurrences
            if val == 9:
                rc9_count += 1
                if rc9_count == 1:
                    first_rc9 = obs
                elif rc9_count == 2:
                    second_rc9 = obs
            # Monotonicity check
            if not np.isnan(prev_val) and val < prev_val:
                mono_issues.append({
                    **{c: row.get(c, "") for c in BUA_COLS},
                    "Field": row.get("Field", ""),
                    "Responsible": row.get("Responsible", ""),
                    "State": row.get("State", ""),
                    "PrevDate": prev_date,
                    "PrevRC": prev_val,
                    "CurrDate": obs,
                    "CurrRC": val,
                    "Issue": f"RC decreased: {prev_val} -> {val}",
                })
            prev_val = val
            prev_date = obs

        # Backlog record
        last_row = grp.iloc[-1]
        latest_rc = grp["RC_Val"].dropna().iloc[-1] if grp["RC_Val"].notna().any() else np.nan
        latest_obs = grp["Observationdate"].max()
        latest_gs = last_row.get("Cropstagemajcode", "")
        days_since = (today - latest_obs).days if pd.notna(latest_obs) else np.nan

        backlog_rows.append({
            **{c: last_row.get(c, "") for c in TRAJECTORY_KEY},
            "Field": last_row.get("Field", ""),
            "Responsible": last_row.get("Responsible", ""),
            "State": last_row.get("State", ""),
            "LatestObservationdate": latest_obs,
            "LatestGS": latest_gs,
            "LatestRC": latest_rc,
            "MaxRC": grp["RC_Val"].max(),
            "FirstRC9": first_rc9,
            "SecondRC9": second_rc9,
            "DaysSinceLatest": days_since,
            "OpenRC": 1 if (np.isnan(latest_rc) or latest_rc < 9) else 0,
        })

    df_mono = pd.DataFrame(mono_issues)
    df_backlog = pd.DataFrame(backlog_rows)
    print(f"    RC monotonicity issues: {len(df_mono)}")
    print(f"    Trajectories tracked: {len(df_backlog)}")
    print(f"    Open RC (latest < 9): {df_backlog['OpenRC'].sum() if len(df_backlog) else 0}")
    return df_mono, df_backlog


# ==================== 5. GS COVERAGE METRICS ====================

def compute_gs_coverage(df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """
    Compute GS coverage per responsible, per region, and overall.
    Both exact-stage and phase coverages.
    """
    print("\n  [METRICS] GS Coverage ...")
    results = {}
    expected_set = set(GS_ORDER)
    n_expected = len(expected_set)

    # ---- Per Responsible ----
    rows_resp = []
    for resp, grp in df.groupby("Responsible"):
        observed = set()
        for col in ("Cropstagemincode", "Cropstagemajcode", "Cropstagemaxcode"):
            if col in grp.columns:
                observed.update(grp[col].dropna().unique())
        captured = observed & expected_set
        n_captured = len(captured)
        missing = sorted(expected_set - captured, key=lambda x: GS_INDEX.get(x, 999))
        # Phase coverage
        phase_scores = {}
        for phase, stages in GS_PHASES.items():
            expected_phase = set(stages)
            observed_phase = captured & expected_phase
            phase_scores[phase] = (
                len(observed_phase) / len(expected_phase) * 100 if expected_phase else 0
            )
        phase_avg = np.mean(list(phase_scores.values()))
        rows_resp.append({
            "Responsible": resp,
            "UniqueGS": n_captured,
            "ExpectedGS": n_expected,
            "CoveragePC": round(n_captured / n_expected * 100, 1),
            "MissingGS": ", ".join(missing) if missing else "—",
            **{f"Phase_{k}": round(v, 1) for k, v in phase_scores.items()},
            "PhaseAvg": round(phase_avg, 1),
        })
    results["by_responsible"] = pd.DataFrame(rows_resp)

    # ---- Per Region ----
    rows_region = []
    for state, grp in df.groupby("State"):
        observed = set()
        for col in ("Cropstagemincode", "Cropstagemajcode", "Cropstagemaxcode"):
            if col in grp.columns:
                observed.update(grp[col].dropna().unique())
        captured = observed & expected_set
        n_captured = len(captured)
        missing = sorted(expected_set - captured, key=lambda x: GS_INDEX.get(x, 999))
        phase_scores = {}
        for phase, stages in GS_PHASES.items():
            expected_phase = set(stages)
            observed_phase = captured & expected_phase
            phase_scores[phase] = (
                len(observed_phase) / len(expected_phase) * 100 if expected_phase else 0
            )
        phase_avg = np.mean(list(phase_scores.values()))
        rows_region.append({
            "State": state,
            "UniqueGS": n_captured,
            "ExpectedGS": n_expected,
            "CoveragePC": round(n_captured / n_expected * 100, 1),
            "MissingGS": ", ".join(missing) if missing else "—",
            **{f"Phase_{k}": round(v, 1) for k, v in phase_scores.items()},
            "PhaseAvg": round(phase_avg, 1),
        })
    results["by_region"] = pd.DataFrame(rows_region)

    # ---- Heatmap matrix: region x GS ----
    mat_rows = []
    for state, grp in df.groupby("State"):
        observed = set()
        for col in ("Cropstagemincode", "Cropstagemajcode", "Cropstagemaxcode"):
            if col in grp.columns:
                observed.update(grp[col].dropna().unique())
        row = {"State": state}
        for gs in GS_ORDER:
            row[gs] = 1 if gs in observed else 0
        mat_rows.append(row)
    results["heatmap_region"] = pd.DataFrame(mat_rows).set_index("State")

    # ---- Heatmap matrix: responsible x GS ----
    mat_rows2 = []
    for resp, grp in df.groupby("Responsible"):
        observed = set()
        for col in ("Cropstagemincode", "Cropstagemajcode", "Cropstagemaxcode"):
            if col in grp.columns:
                observed.update(grp[col].dropna().unique())
        row = {"Responsible": resp}
        for gs in GS_ORDER:
            row[gs] = 1 if gs in observed else 0
        mat_rows2.append(row)
    results["heatmap_responsible"] = pd.DataFrame(mat_rows2).set_index("Responsible")

    print(f"    Coverage computed for {len(rows_resp)} responsibles, {len(rows_region)} regions")
    return results


# ==================== 6. KPI FRAMEWORK ====================

def compute_kpis(
    df: pd.DataFrame,
    qa: Dict[str, pd.DataFrame],
    gs_cov: Dict[str, pd.DataFrame],
    rc_backlog: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build KPI tables for responsible and region.

    Weights:
      30 % GS Trio Completeness
      20 % Assessment Package Completeness
      20 % GS Coverage (breadth + phase balance)
      15 % RC Completion
      15 % Inverse QA Issue Burden
    """
    print("\n  [KPI] Computing composite scores ...")

    def _kpi_for_group(df_group, group_col):
        rows = []
        for name, grp in df_group.groupby(group_col):
            n_fields = grp["Field"].nunique()
            n_trials = grp["Trialnumber"].nunique()
            n_traj = grp.groupby(TRAJECTORY_KEY, dropna=False).ngroups
            n_events = len(grp)

            # ---- GS Trio Completeness (% of BUAs with complete trio) ----
            gs_cols = ["Cropstagemincode", "Cropstagemajcode", "Cropstagemaxcode"]
            gs_present = [c for c in gs_cols if c in grp.columns]
            if gs_present:
                complete = grp[gs_present].notna().all(axis=1).sum()
                trio_pc = complete / len(grp) * 100 if len(grp) > 0 else 0
            else:
                trio_pc = 0

            # ---- Assessment Package Completeness ----
            if "assessment_pkg" in qa and len(qa["assessment_pkg"]) > 0:
                pkg_issues = qa["assessment_pkg"]
                if group_col in pkg_issues.columns:
                    n_pkg_issues = len(pkg_issues[pkg_issues[group_col] == name])
                elif group_col == "State" and "State" in pkg_issues.columns:
                    n_pkg_issues = len(pkg_issues[pkg_issues["State"] == name])
                else:
                    n_pkg_issues = 0
                total_buas = grp.groupby(BUA_COLS, dropna=False).ngroups
                pkg_pc = max(0, (1 - n_pkg_issues / max(total_buas, 1)) * 100)
            else:
                pkg_pc = 100

            # ---- GS Coverage Score ----
            if group_col == "Responsible":
                cov_df = gs_cov.get("by_responsible", pd.DataFrame())
                cov_row = cov_df[cov_df["Responsible"] == name] if len(cov_df) else pd.DataFrame()
            else:
                cov_df = gs_cov.get("by_region", pd.DataFrame())
                cov_row = cov_df[cov_df["State"] == name] if len(cov_df) else pd.DataFrame()
            if len(cov_row):
                gs_cov_score = cov_row["CoveragePC"].values[0]
                phase_score = cov_row["PhaseAvg"].values[0]
            else:
                gs_cov_score = 0
                phase_score = 0
            coverage_combined = (gs_cov_score + phase_score) / 2

            # ---- RC Completion ----
            if len(rc_backlog) > 0:
                rc_grp = rc_backlog[rc_backlog[group_col] == name] if group_col in rc_backlog.columns else pd.DataFrame()
                if len(rc_grp) > 0:
                    rc_complete = (rc_grp["OpenRC"] == 0).sum()
                    rc_total = len(rc_grp)
                    rc_pc = rc_complete / rc_total * 100 if rc_total else 0
                    open_rc_count = rc_grp["OpenRC"].sum()
                else:
                    rc_pc = 0
                    rc_total = 0
                    rc_complete = 0
                    open_rc_count = 0
            else:
                rc_pc = 0
                rc_total = 0
                rc_complete = 0
                open_rc_count = 0

            # ---- Issue Burden ----
            issue_keys = [
                "gs_regression", "dev_gaps", "gs_triplet_consistency",
                "disease_mono", "rc_mono",
            ]
            total_issues = 0
            issue_counts = {}
            for ik in issue_keys:
                if ik in qa and len(qa[ik]) > 0:
                    idf = qa[ik]
                    col_match = group_col if group_col in idf.columns else None
                    if col_match:
                        cnt = len(idf[idf[col_match] == name])
                    else:
                        cnt = 0
                else:
                    cnt = 0
                issue_counts[ik] = cnt
                total_issues += cnt
            # Inverse burden: 100 if 0 issues, decreasing.  Cap at 100.
            issue_burden_score = max(0, 100 - total_issues * 2)  # -2 pts per issue

            # ---- Composite Score ----
            overall = (
                KPI_WEIGHTS["gs_trio_completeness"] * trio_pc
                + KPI_WEIGHTS["assessment_pkg_completeness"] * pkg_pc
                + KPI_WEIGHTS["gs_coverage"] * coverage_combined
                + KPI_WEIGHTS["rc_completion"] * rc_pc
                + KPI_WEIGHTS["inverse_issue_burden"] * issue_burden_score
            )

            rows.append({
                group_col: name,
                "Fields": n_fields,
                "Trials": n_trials,
                "UniqueTrajectories": n_traj,
                "CollectionEvents": n_events,
                "GS_Trio_Completeness": round(trio_pc, 1),
                "Assessment_Pkg_Completeness": round(pkg_pc, 1),
                "GS_Coverage_Score": round(gs_cov_score, 1),
                "Phase_Coverage_Score": round(phase_score, 1),
                "GS_Regression_Count": issue_counts.get("gs_regression", 0),
                "Dev_Gap_Count": issue_counts.get("dev_gaps", 0),
                "GS_Triplet_Issues": issue_counts.get("gs_triplet_consistency", 0),
                "Disease_Mono_Issues": issue_counts.get("disease_mono", 0),
                "RC_Mono_Issues": issue_counts.get("rc_mono", 0),
                "RC_Complete": int(rc_complete),
                "RC_Total": int(rc_total),
                "RC_Completion_Rate": round(rc_pc, 1),
                "Open_RC_Backlog": int(open_rc_count),
                "Overall_Score": round(overall, 1),
            })
        return pd.DataFrame(rows).sort_values("Overall_Score", ascending=False).reset_index(drop=True)

    kpi_resp = _kpi_for_group(df, "Responsible")
    kpi_region = _kpi_for_group(df, "State")
    print(f"    KPIs computed for {len(kpi_resp)} responsibles, {len(kpi_region)} regions")
    return kpi_resp, kpi_region


# ==================== 7. FIELD & TRAJECTORY SUMMARIES ====================

def build_field_summary(
    df: pd.DataFrame,
    qa: Dict[str, pd.DataFrame],
    rc_backlog: pd.DataFrame,
) -> pd.DataFrame:
    print("\n  [SUMMARY] Building field-level summary ...")

    def _count_issues(issue_df, field_val):
        if issue_df is None or len(issue_df) == 0:
            return 0
        if "Field" in issue_df.columns:
            return len(issue_df[issue_df["Field"] == field_val])
        return 0

    rows = []
    for field, grp in df.groupby("Field"):
        resp = grp["Responsible"].iloc[0]
        state = grp["State"].iloc[0]
        trials = grp["Trialnumber"].nunique()
        n_traj = grp.groupby(TRAJECTORY_KEY, dropna=False).ngroups
        n_events = len(grp)
        first_obs = grp["Observationdate"].min()
        last_obs = grp["Observationdate"].max()
        unique_gs = set()
        for col in ("Cropstagemincode", "Cropstagemajcode", "Cropstagemaxcode"):
            if col in grp.columns:
                unique_gs.update(grp[col].dropna().unique())
        gs_cols_present = [c for c in ("Cropstagemincode", "Cropstagemajcode", "Cropstagemaxcode") if c in grp.columns]
        trio_complete = grp[gs_cols_present].notna().all(axis=1).sum() if gs_cols_present else 0
        trio_rate = trio_complete / len(grp) * 100 if len(grp) else 0

        # RC status
        if len(rc_backlog) and "Field" in rc_backlog.columns:
            rc_field = rc_backlog[rc_backlog["Field"] == field]
            open_rc = int(rc_field["OpenRC"].sum()) if len(rc_field) else 0
        else:
            open_rc = 0

        rows.append({
            "Field": field,
            "Responsible": resp,
            "Region": state,
            "Trials": trials,
            "Trajectories": n_traj,
            "CollectionEvents": n_events,
            "FirstObservation": first_obs,
            "LastObservation": last_obs,
            "UniqueGS_Captured": len(unique_gs),
            "GS_Trio_Completeness_PC": round(trio_rate, 1),
            "GS_Regression": _count_issues(qa.get("gs_regression"), field),
            "Dev_Gaps": _count_issues(qa.get("dev_gaps"), field),
            "Disease_Mono_Issues": _count_issues(qa.get("disease_mono"), field),
            "RC_Mono_Issues": _count_issues(qa.get("rc_mono"), field),
            "Open_RC": open_rc,
        })
    result = pd.DataFrame(rows).sort_values(["Region", "Responsible", "Field"])
    print(f"    Field summary: {len(result)} fields")
    return result


def build_trajectory_summary(df: pd.DataFrame) -> pd.DataFrame:
    print("\n  [SUMMARY] Building trajectory-level summary ...")
    rows = []
    for traj, grp in df.groupby(TRAJECTORY_KEY, dropna=False):
        grp_sorted = grp.sort_values("Observationdate")
        first_row = grp_sorted.iloc[0]
        last_row = grp_sorted.iloc[-1]
        unique_gs = set()
        for col in ("Cropstagemincode", "Cropstagemajcode", "Cropstagemaxcode"):
            if col in grp.columns:
                unique_gs.update(grp[col].dropna().unique())
        rows.append({
            **{c: first_row.get(c, "") for c in TRAJECTORY_KEY},
            "Field": first_row.get("Field", ""),
            "Responsible": first_row.get("Responsible", ""),
            "State": first_row.get("State", ""),
            "N_Observations": len(grp),
            "FirstDate": grp["Observationdate"].min(),
            "LastDate": grp["Observationdate"].max(),
            "FirstGS_Maj": first_row.get("Cropstagemajcode", ""),
            "LastGS_Maj": last_row.get("Cropstagemajcode", ""),
            "UniqueGS": len(unique_gs),
        })
    result = pd.DataFrame(rows)
    print(f"    Trajectory summary: {len(result)} trajectories")
    return result


# ==================== 8. UNKNOWN ASSIGNMENT SUMMARY ====================

def build_unknown_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Produce a table of fields flagged as UNASSIGNED or UNK."""
    rows = []
    for field in df["Field"].unique():
        sub = df[df["Field"] == field].iloc[0]
        resp = sub["Responsible"]
        state = sub["State"]
        if resp == "UNASSIGNED" or state == "UNK":
            rows.append({
                "Field": field,
                "Responsible": resp,
                "State": state,
                "FlaggedResponsible": resp == "UNASSIGNED",
                "FlaggedState": state == "UNK",
            })
    return pd.DataFrame(rows)


# ==================== 9. VISUALS ====================
def _priority_color(value, thresholds, metric_type="compliance"):
    """Return color based on value & thresholds."""
    if metric_type == "compliance":
        if value >= thresholds[0]:
            return SUCCESS_GREEN
        elif value >= thresholds[1]:
            return WARNING_YELLOW
        return ALERT_RED
    else:  # 'issues' – lower is better
        if value == 0:
            return SUCCESS_GREEN
        elif value <= thresholds[0]:
            return WARNING_YELLOW
        return ALERT_RED


def _priority_label(value, thresholds, metric_type="compliance"):
    """Return status label."""
    if metric_type == "compliance":
        if value >= thresholds[0]:
            return "Good"
        elif value >= thresholds[1]:
            return "Warning"
        return "Poor"
    else:
        if value == 0:
            return "Good"
        elif value <= thresholds[0]:
            return "Warning"
        return "Poor"
def _brand_cmap():
    """Custom colormap from white to brand color."""
    return LinearSegmentedColormap.from_list("xarvio", ["#FFFFFF", BRAND_COLOR])


def generate_all_plots(
    df: pd.DataFrame,
    qa: Dict[str, pd.DataFrame],
    gs_cov: Dict[str, pd.DataFrame],
    kpi_resp: pd.DataFrame,
    kpi_region: pd.DataFrame,
    rc_backlog: pd.DataFrame,
    output_dir: str,
) -> Dict[str, str]:
    """Generate all required visuals. Returns dict name -> filepath."""
    print("\n  [VISUALS] Generating plots ...")
    os.makedirs(output_dir, exist_ok=True)
    plots = {}

    def _save(fig, name):
        p = os.path.join(output_dir, f"{name}.png")
        fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        plots[name] = p
        return p

    # 1. GS coverage heatmap by region
    hm_region = gs_cov.get("heatmap_region")
    if hm_region is not None and len(hm_region):
        hm_data = hm_region.reindex(columns=GS_ORDER)
        annot_arr = hm_data.applymap(lambda v: "\u2713" if v == 1 else "\u2717")
        fig, ax = plt.subplots(figsize=(14, max(3, len(hm_data) * 0.55)))
        sns.heatmap(
            hm_data, cmap=_brand_cmap(), linewidths=0.5, linecolor="white",
            cbar_kws={"label": "Captured (1) / Missing (0)"},
            ax=ax, vmin=0, vmax=1, annot=annot_arr, fmt="",
            annot_kws={"fontsize": 9, "fontweight": "bold"},
        )
        ax.set_title("GS Coverage by Region", fontsize=14, color=BRAND_COLOR, fontweight="bold")
        ax.set_xlabel("Growth Stage")
        ax.set_ylabel("Region / State")
        _save(fig, "gs_coverage_heatmap_region")

    # 2. GS coverage heatmap by responsible
    hm_resp = gs_cov.get("heatmap_responsible")
    if hm_resp is not None and len(hm_resp):
        hm_data = hm_resp.reindex(columns=GS_ORDER)
        annot_arr = hm_data.applymap(lambda v: "\u2713" if v == 1 else "\u2717")
        fig, ax = plt.subplots(figsize=(14, max(3, len(hm_data) * 0.55)))
        sns.heatmap(
            hm_data, cmap=_brand_cmap(), linewidths=0.5, linecolor="white",
            cbar_kws={"label": "Captured (1) / Missing (0)"},
            ax=ax, vmin=0, vmax=1, annot=annot_arr, fmt="",
            annot_kws={"fontsize": 9, "fontweight": "bold"},
        )
        ax.set_title("GS Coverage by Responsible", fontsize=14, color=BRAND_COLOR, fontweight="bold")
        ax.set_xlabel("Growth Stage")
        ax.set_ylabel("Responsible")
        _save(fig, "gs_coverage_heatmap_responsible")

    # 3. Phase coverage chart
    cov_resp = gs_cov.get("by_responsible", pd.DataFrame())
    if len(cov_resp):
        phase_cols = [c for c in cov_resp.columns if c.startswith("Phase_")]
        if phase_cols:
            fig, ax = plt.subplots(figsize=(10, 5))
            x = np.arange(len(cov_resp))
            width = 0.8 / len(phase_cols)
            phase_colors = [BRAND_COLOR, "#c75ab5", "#e8a0d8", DARK_GREY]
            for i, col in enumerate(phase_cols):
                label = col.replace("Phase_", "")
                c = phase_colors[i % len(phase_colors)]
                ax.bar(x + i * width, cov_resp[col], width, label=label, color=c, edgecolor="white")
            ax.set_xticks(x + width * (len(phase_cols) - 1) / 2)
            ax.set_xticklabels(cov_resp["Responsible"], rotation=30, ha="right")
            ax.set_ylabel("Coverage %")
            ax.set_title("GS Phase Coverage by Responsible", color=BRAND_COLOR, fontweight="bold")
            ax.legend(fontsize=8)
            ax.set_ylim(0, 105)
            ax.grid(axis="y", alpha=0.3, linestyle="--")
            _save(fig, "gs_phase_coverage")

    # 4. GS completeness chart (Trio completeness by responsible)
    if len(kpi_resp):
        fig, ax = plt.subplots(figsize=(9, 5))
        sorted_kpi = kpi_resp.sort_values("GS_Trio_Completeness")
        bar_colors = [
            SUCCESS_GREEN if v >= 90 else WARNING_YELLOW if v >= 70 else ALERT_RED
            for v in sorted_kpi["GS_Trio_Completeness"]
        ]
        ax.barh(sorted_kpi["Responsible"], sorted_kpi["GS_Trio_Completeness"], color=bar_colors, edgecolor=DARK_GREY)
        ax.set_xlabel("GS Trio Completeness %")
        ax.set_title("GS Trio Completeness by Responsible", color=BRAND_COLOR, fontweight="bold")
        ax.axvline(90, color=SUCCESS_GREEN, ls="--", alpha=0.7, label="Target 90%")
        ax.set_xlim(0, 105)
        ax.legend(fontsize=8)
        for i, (v, r) in enumerate(zip(sorted_kpi["GS_Trio_Completeness"], sorted_kpi["Responsible"])):
            ax.text(v + 0.5, i, f"{v:.0f}%", va="center", fontsize=8)
        _save(fig, "gs_completeness_chart")

    # 5. GS regression distribution
    gs_reg = qa.get("gs_regression", pd.DataFrame())
    if len(gs_reg) and "Responsible" in gs_reg.columns:
        fig, ax = plt.subplots(figsize=(9, 5))
        counts = gs_reg["Responsible"].value_counts().sort_values()
        ax.barh(counts.index, counts.values, color=BRAND_COLOR, edgecolor=DARK_GREY)
        ax.set_xlabel("Regression Events")
        ax.set_title("GS Regression Distribution by Responsible", color=BRAND_COLOR, fontweight="bold")
        for i, (idx, v) in enumerate(counts.items()):
            ax.text(v + 0.3, i, str(int(v)), va="center", fontsize=9, fontweight="bold", color=BRAND_COLOR)
        _save(fig, "gs_regression_dist")

    # 6. Development gap distribution
    dg = qa.get("dev_gaps", pd.DataFrame())
    if len(dg) and "Responsible" in dg.columns:
        fig, ax = plt.subplots(figsize=(9, 6))
        counts = dg["Responsible"].value_counts().sort_values()
        ax.barh(counts.index, counts.values, color=BRAND_COLOR, edgecolor=DARK_GREY)
        ax.set_xlabel("Development Gap Events")
        for i, (idx, v) in enumerate(counts.items()):
            ax.text(v + 0.3, i, str(int(v)), va="center", fontsize=9, fontweight="bold", color=BRAND_COLOR)
        # Build subtitle with highest jump per responsible
        subtitle_parts = []
        if "Issue" in dg.columns and "Field" in dg.columns:
            for resp in counts.index:
                resp_dg = dg[dg["Responsible"] == resp]
                if len(resp_dg):
                    max_row = resp_dg.iloc[0]
                    field_short = str(max_row.get("Field", ""))[:25]
                    subtitle_parts.append(f"{resp}: {field_short}")
        subtitle = "Largest gaps — " + "; ".join(subtitle_parts[:4]) if subtitle_parts else ""
        ax.set_title("Development Gap Distribution by Responsible\n" + subtitle,
                     color=BRAND_COLOR, fontweight="bold", fontsize=11)
        _save(fig, "dev_gap_dist")

    # 7. Assessment package completeness chart
    if len(kpi_resp):
        fig, ax = plt.subplots(figsize=(9, 5))
        sorted_kpi = kpi_resp.sort_values("Assessment_Pkg_Completeness")
        bar_colors = [
            SUCCESS_GREEN if v >= 90 else WARNING_YELLOW if v >= 70 else ALERT_RED
            for v in sorted_kpi["Assessment_Pkg_Completeness"]
        ]
        ax.barh(sorted_kpi["Responsible"], sorted_kpi["Assessment_Pkg_Completeness"], color=bar_colors, edgecolor=DARK_GREY)
        ax.set_xlabel("Assessment Package Completeness %")
        ax.set_title("Assessment Package Completeness by Responsible", color=BRAND_COLOR, fontweight="bold")
        ax.set_xlim(0, 105)
        for i, (v, r) in enumerate(zip(sorted_kpi["Assessment_Pkg_Completeness"], sorted_kpi["Responsible"])):
            ax.text(v + 0.5, i, f"{v:.0f}%", va="center", fontsize=8)
        _save(fig, "assessment_completeness_chart")

    # 8. RC progress summary
    if len(kpi_resp):
        fig, ax = plt.subplots(figsize=(9, 5))
        sorted_kpi = kpi_resp.sort_values("RC_Completion_Rate")
        bar_colors = [
            SUCCESS_GREEN if v >= 80 else WARNING_YELLOW if v >= 50 else ALERT_RED
            for v in sorted_kpi["RC_Completion_Rate"]
        ]
        ax.barh(sorted_kpi["Responsible"], sorted_kpi["RC_Completion_Rate"], color=bar_colors, edgecolor=DARK_GREY)
        ax.set_xlabel("RC Completion Rate (RC=9) %")
        ax.set_title("Row Closing Completion by Responsible", color=BRAND_COLOR, fontweight="bold")
        ax.set_xlim(0, 105)
        for i, (v, r) in enumerate(zip(sorted_kpi["RC_Completion_Rate"], sorted_kpi["Responsible"])):
            ax.text(v + 0.5, i, f"{v:.0f}%", va="center", fontsize=8)
        _save(fig, "rc_progress_chart")

    # 9. RC Open Backlog summary
    if len(kpi_resp):
        fig, ax = plt.subplots(figsize=(9, 5))
        sorted_kpi = kpi_resp.sort_values("Open_RC_Backlog", ascending=False)
        ax.barh(sorted_kpi["Responsible"], sorted_kpi["Open_RC_Backlog"], color=BRAND_COLOR, edgecolor=DARK_GREY)
        ax.set_xlabel("Open RC Trajectories (latest RC < 9)")
        ax.set_title("Open RC Backlog by Responsible", color=BRAND_COLOR, fontweight="bold")
        _save(fig, "rc_open_backlog_chart")

    # 10. RC vs GS by Region (PRIMARY RC CHART)
    rc_df = df[df["Assessmenttypecode"] == RC_ASSESSMENT].copy()
    if len(rc_df) and "Cropstagemajcode" in rc_df.columns:
        rc_df["RC_Val"] = pd.to_numeric(rc_df["Originalvalue"], errors="coerce")
        rc_df = rc_df.dropna(subset=["RC_Val", "Cropstagemajcode"])
        rc_df["GS_Num"] = rc_df["Cropstagemajcode"].apply(gs_numeric)
        rc_df = rc_df.dropna(subset=["GS_Num"])
        if len(rc_df):
            rc_agg = rc_df.groupby(["State", "Cropstagemajcode"]).agg(
                MeanRC=("RC_Val", "mean"), Count=("RC_Val", "size"), GS_Num=("GS_Num", "first")
            ).reset_index().sort_values("GS_Num")
            fig, ax = plt.subplots(figsize=(12, 6))
            states = rc_agg["State"].unique()
            palette = sns.color_palette("husl", len(states))
            for i, st in enumerate(sorted(states)):
                sdf = rc_agg[rc_agg["State"] == st].sort_values("GS_Num")
                ax.plot(sdf["Cropstagemajcode"], sdf["MeanRC"], marker="o", label=st, color=palette[i], linewidth=2)
            ax.set_xlabel("Growth Stage", fontweight="bold")
            ax.set_ylabel("Mean Row Closing (0-9)", fontweight="bold")
            ax.set_title("Row Closing vs Growth Stage by Region", color=BRAND_COLOR, fontweight="bold", fontsize=14)
            ax.legend(fontsize=8, title="Region")
            ax.grid(axis="y", alpha=0.3, linestyle="--")
            ax.set_ylim(-0.5, 9.5)
            plt.xticks(rotation=45, ha="right")
            _save(fig, "rc_vs_gs_by_region")

    # 11. RC vs GS by Responsible
    if len(rc_df):
        rc_agg2 = rc_df.groupby(["Responsible", "Cropstagemajcode"]).agg(
            MeanRC=("RC_Val", "mean"), GS_Num=("GS_Num", "first")
        ).reset_index().sort_values("GS_Num")
        if len(rc_agg2):
            fig, ax = plt.subplots(figsize=(12, 6))
            resps = rc_agg2["Responsible"].unique()
            palette = sns.color_palette("husl", len(resps))
            for i, rp in enumerate(sorted(resps)):
                sdf = rc_agg2[rc_agg2["Responsible"] == rp].sort_values("GS_Num")
                ax.plot(sdf["Cropstagemajcode"], sdf["MeanRC"], marker="o", label=rp, color=palette[i], linewidth=2)
            ax.set_xlabel("Growth Stage", fontweight="bold")
            ax.set_ylabel("Mean Row Closing (0-9)", fontweight="bold")
            ax.set_title("Row Closing vs Growth Stage by Responsible", color=BRAND_COLOR, fontweight="bold", fontsize=14)
            ax.legend(fontsize=8, title="Responsible")
            ax.grid(axis="y", alpha=0.3, linestyle="--")
            ax.set_ylim(-0.5, 9.5)
            plt.xticks(rotation=45, ha="right")
            _save(fig, "rc_vs_gs_by_responsible")

    # 12. KPI overview (Overall Score)
    if len(kpi_resp):
        fig, ax = plt.subplots(figsize=(10, 5))
        sorted_kpi = kpi_resp.sort_values("Overall_Score")
        bar_colors = [
            SUCCESS_GREEN if v >= 80 else WARNING_YELLOW if v >= 60 else ALERT_RED
            for v in sorted_kpi["Overall_Score"]
        ]
        ax.barh(sorted_kpi["Responsible"], sorted_kpi["Overall_Score"], color=bar_colors, edgecolor=DARK_GREY)
        ax.set_xlabel("Overall Composite Score (0-100)")
        ax.set_title("KPI Overall Score by Responsible", color=BRAND_COLOR, fontweight="bold", fontsize=14)
        ax.axvline(80, color=SUCCESS_GREEN, ls="--", alpha=0.7, label="Good (≥80)")
        ax.axvline(60, color=WARNING_YELLOW, ls="--", alpha=0.7, label="Acceptable (≥60)")
        ax.legend(fontsize=8)
        ax.set_xlim(0, 105)
        for i, (v, r) in enumerate(zip(sorted_kpi["Overall_Score"], sorted_kpi["Responsible"])):
            ax.text(v + 0.5, i, f"{v:.0f}", va="center", fontsize=9, fontweight="bold")
        _save(fig, "kpi_overall_score")

    # 12b. Compliance Dual Panel (GS Trio vs Assessment Pkg side-by-side)
    if len(kpi_resp):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
        sorted_kpi = kpi_resp.sort_values("Overall_Score")
        # Left: GS Trio
        bar_c1 = [SUCCESS_GREEN if v >= 90 else WARNING_YELLOW if v >= 70 else ALERT_RED
                   for v in sorted_kpi["GS_Trio_Completeness"]]
        ax1.barh(sorted_kpi["Responsible"], sorted_kpi["GS_Trio_Completeness"],
                 color=bar_c1, edgecolor=DARK_GREY)
        ax1.axvline(90, color=SUCCESS_GREEN, ls="--", alpha=0.7)
        ax1.set_xlim(0, 105)
        ax1.set_xlabel("GS Trio Completeness %")
        ax1.set_title("GS Trio Completeness", color=BRAND_COLOR, fontweight="bold")
        for i, v in enumerate(sorted_kpi["GS_Trio_Completeness"]):
            ax1.text(v + 0.5, i, f"{v:.0f}%", va="center", fontsize=8)
        # Right: Assessment Pkg
        bar_c2 = [SUCCESS_GREEN if v >= 90 else WARNING_YELLOW if v >= 70 else ALERT_RED
                   for v in sorted_kpi["Assessment_Pkg_Completeness"]]
        ax2.barh(sorted_kpi["Responsible"], sorted_kpi["Assessment_Pkg_Completeness"],
                 color=bar_c2, edgecolor=DARK_GREY)
        ax2.axvline(90, color=SUCCESS_GREEN, ls="--", alpha=0.7)
        ax2.set_xlim(0, 105)
        ax2.set_xlabel("Assessment Pkg Completeness %")
        ax2.set_title("Assessment Pkg Completeness", color=BRAND_COLOR, fontweight="bold")
        for i, v in enumerate(sorted_kpi["Assessment_Pkg_Completeness"]):
            ax2.text(v + 0.5, i, f"{v:.0f}%", va="center", fontsize=8)
        fig.suptitle("Compliance Metrics Comparison", fontsize=14,
                     fontweight="bold", color=BRAND_COLOR, y=1.02)
        fig.tight_layout()
        _save(fig, "compliance_dual_panel")

    # 12c. Disease Monotonicity Chart
    dis_mono = qa.get("disease_mono", pd.DataFrame())
    if len(dis_mono) and "Responsible" in dis_mono.columns:
        fig, ax = plt.subplots(figsize=(9, 5))
        counts = dis_mono["Responsible"].value_counts().sort_values()
        ax.barh(counts.index, counts.values, color=BRAND_COLOR, edgecolor=DARK_GREY)
        ax.set_xlabel("Disease Monotonicity Violations")
        ax.set_title("Disease Monotonicity Violations by Responsible",
                     color=BRAND_COLOR, fontweight="bold")
        for i, (idx, v) in enumerate(counts.items()):
            ax.text(v + 0.3, i, str(int(v)), va="center", fontsize=9, fontweight="bold", color=BRAND_COLOR)
        _save(fig, "disease_mono_chart")

    # 12d. QA Issues by Type per Responsible (grouped bar)
    issue_types = {
        "GS Regression": qa.get("gs_regression", pd.DataFrame()),
        "Dev Gaps": qa.get("dev_gaps", pd.DataFrame()),
        "Triplet Issues": qa.get("gs_triplet_consistency", pd.DataFrame()),
        "Disease Mono": qa.get("disease_mono", pd.DataFrame()),
        "RC Mono": qa.get("rc_mono", pd.DataFrame()),
    }
    all_resps = sorted(kpi_resp["Responsible"].unique()) if len(kpi_resp) else []
    if all_resps:
        issue_matrix = {}
        for itype, idf in issue_types.items():
            if len(idf) and "Responsible" in idf.columns:
                vc = idf["Responsible"].value_counts()
                issue_matrix[itype] = {r: vc.get(r, 0) for r in all_resps}
            else:
                issue_matrix[itype] = {r: 0 for r in all_resps}
        fig, ax = plt.subplots(figsize=(14, 6))
        x = np.arange(len(all_resps))
        n_types = len(issue_types)
        width = 0.8 / n_types
        colors_issue = [BRAND_COLOR, "#c75ab5", "#CE93D8", ALERT_RED, WARNING_YELLOW]
        for i, (itype, counts_map) in enumerate(issue_matrix.items()):
            vals = [counts_map[r] for r in all_resps]
            bars = ax.bar(x + i * width, vals, width, label=itype,
                          color=colors_issue[i % len(colors_issue)], edgecolor="white")
            for bar, v in zip(bars, vals):
                if v > 0:
                    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.2,
                            str(int(v)), ha="center", va="bottom", fontsize=7, fontweight="bold")
        ax.set_xticks(x + width * (n_types - 1) / 2)
        ax.set_xticklabels(all_resps, rotation=30, ha="right")
        ax.set_ylabel("Issue Count")
        ax.set_title("QA Issues Distribution by Type per Responsible",
                     color=BRAND_COLOR, fontweight="bold", fontsize=14)
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(axis="y", alpha=0.3, linestyle="--")
        _save(fig, "issues_by_type_responsible")

    # ================================================================
    # NEW PLOTS 13-21: DAP & PHENOLOGY ANALYSIS
    # ================================================================

    # --- Prepare DAP subset (used by plots 13-19) ---
    df_dap = df[
        (df["DAP"].notna()) & (df["DAP"] > 0)
        & (df["Cropstagemajcode"].notna())
    ].copy()
    if len(df_dap):
        df_dap["GS_Num"] = df_dap["Cropstagemajcode"].apply(gs_numeric)
        df_dap["GS_Phase"] = df_dap["Cropstagemajcode"].apply(get_gs_phase)
        df_dap = df_dap[df_dap["GS_Num"].notna()]

    PHASE_COLORS = {
        "Early Vegetative":   "#7B1FA2",
        "Late Vegetative":    BRAND_COLOR,
        "Early Reproductive": "#CE93D8",
        "Late Reproductive":  "#E1BEE7",
    }

    # 13. GS Distribution by Phase (bar chart with record counts)
    if len(df_dap):
        gs_counts = df_dap.groupby(["Cropstagemajcode", "GS_Phase"]).size().reset_index(name="Count")
        gs_counts["GS_Num"] = gs_counts["Cropstagemajcode"].apply(gs_numeric)
        gs_counts = gs_counts.sort_values("GS_Num")
        fig, ax = plt.subplots(figsize=(14, 6))
        colors_bar = [PHASE_COLORS.get(p, DARK_GREY) for p in gs_counts["GS_Phase"]]
        bars = ax.bar(gs_counts["Cropstagemajcode"], gs_counts["Count"],
                       color=colors_bar, edgecolor="white", linewidth=0.5)
        for bar, cnt in zip(bars, gs_counts["Count"]):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                    str(int(cnt)), ha="center", va="bottom", fontsize=8, fontweight="bold")
        ax.set_xlabel("Growth Stage", fontweight="bold")
        ax.set_ylabel("Number of Records", fontweight="bold")
        ax.set_title("Growth Stage Distribution by Development Phase",
                      color=BRAND_COLOR, fontweight="bold", fontsize=14)
        handles = [plt.Rectangle((0, 0), 1, 1, fc=c) for c in PHASE_COLORS.values()]
        ax.legend(handles, PHASE_COLORS.keys(), fontsize=8, loc="upper right")
        ax.grid(axis="y", alpha=0.3, linestyle="--")
        plt.xticks(rotation=45, ha="right")
        _save(fig, "gs_distribution_by_phase")

    # 14. GS Records by State × Phase heatmap (count-based)
    if len(df_dap):
        ct = df_dap.groupby(["State", "GS_Phase"]).size().unstack(fill_value=0)
        ct = ct.reindex(columns=GS_PHASE_ORDER, fill_value=0)
        if len(ct):
            fig, ax = plt.subplots(figsize=(10, max(3, len(ct) * 0.6)))
            sns.heatmap(ct, cmap=_brand_cmap(), annot=True, fmt="d",
                        linewidths=0.5, linecolor="white", ax=ax)
            ax.set_title("GS Records by State and Development Phase",
                          color=BRAND_COLOR, fontweight="bold", fontsize=14)
            ax.set_xlabel("Development Phase")
            ax.set_ylabel("State / Region")
            _save(fig, "gs_records_state_phase")

    # 15. ★ GS vs DAP Scatter (THE MOST IMPORTANT PLOT)
    if len(df_dap):
        fig, ax = plt.subplots(figsize=(16, 8))
        states = sorted(df_dap["State"].unique())
        palette = sns.color_palette("husl", len(states))
        state_colors = {s: palette[i] for i, s in enumerate(states)}
        for st in states:
            sdf = df_dap[df_dap["State"] == st]
            ax.scatter(sdf["DAP"], sdf["GS_Num"], label=st,
                       color=state_colors[st], alpha=0.5, s=30, edgecolors="none")
        # Quadratic trend line
        x_all = df_dap["DAP"].values
        y_all = df_dap["GS_Num"].values
        mask = np.isfinite(x_all) & np.isfinite(y_all)
        if mask.sum() > 10:
            coeffs = np.polyfit(x_all[mask], y_all[mask], 2)
            x_fit = np.linspace(x_all[mask].min(), x_all[mask].max(), 200)
            y_fit = np.polyval(coeffs, x_fit)
            ax.plot(x_fit, y_fit, color=BRAND_COLOR, linewidth=2.5,
                    linestyle="--", label="Trend (quadratic)", zorder=5)
        # Phase background bands
        phase_boundaries = {
            "Early Vegetative":   (0, 5),     # VS=0 .. V3=5
            "Late Vegetative":    (6, 11),    # V4=6 .. V9=11
            "Early Reproductive": (12, 15),   # R1=12 .. R4=15
            "Late Reproductive":  (16, 19),   # R5=16 .. R8=19
        }
        for phase, (y0, y1) in phase_boundaries.items():
            ax.axhspan(y0 - 0.5, y1 + 0.5, alpha=0.07,
                       color=PHASE_COLORS.get(phase, LIGHT_GREY))
            ax.text(ax.get_xlim()[1] * 0.98, (y0 + y1) / 2, phase,
                    fontsize=7, ha="right", va="center", style="italic",
                    color=DARK_GREY, alpha=0.7)
        ax.set_xlabel("Days After Planting (DAP)", fontsize=12, fontweight="bold")
        ax.set_ylabel("Growth Stage (numeric index)", fontsize=12, fontweight="bold")
        ax.set_title("Growth Stage vs Days After Planting",
                      color=BRAND_COLOR, fontweight="bold", fontsize=16)
        ax.set_yticks(range(len(GS_ORDER)))
        ax.set_yticklabels(GS_ORDER, fontsize=8)
        ax.legend(fontsize=8, title="Region", loc="upper left",
                  framealpha=0.9, ncol=2)
        ax.grid(alpha=0.2, linestyle="--")
        _save(fig, "gs_vs_dap_scatter")

    # 15b. Per-Region GS vs DAP — one plot per region with individual trial lines
    if len(df_dap):
        for st in sorted(df_dap["State"].unique()):
            sdf = df_dap[df_dap["State"] == st].copy()
            if len(sdf) < 3:
                continue
            fields = sorted(sdf["Field"].unique())
            palette_f = sns.color_palette("husl", len(fields))
            fig, ax = plt.subplots(figsize=(14, 7))
            # Phase background bands
            phase_boundaries = {
                "Early Vegetative":   (0, 5),
                "Late Vegetative":    (6, 11),
                "Early Reproductive": (12, 15),
                "Late Reproductive":  (16, 19),
            }
            for phase, (y0, y1) in phase_boundaries.items():
                ax.axhspan(y0 - 0.5, y1 + 0.5, alpha=0.06,
                           color=PHASE_COLORS.get(phase, LIGHT_GREY))
            for fi, fld in enumerate(fields):
                fdf = sdf[sdf["Field"] == fld].sort_values("DAP")
                fdf_gs = fdf.dropna(subset=["GS_Num"])
                if len(fdf_gs) < 2:
                    continue
                label = str(fld)[:35]
                ax.plot(fdf_gs["DAP"], fdf_gs["GS_Num"], marker="o", markersize=4,
                        linewidth=1.2, alpha=0.7, label=label, color=palette_f[fi])
            ax.set_xlabel("Days After Planting (DAP)", fontsize=11, fontweight="bold")
            ax.set_ylabel("Growth Stage", fontsize=11, fontweight="bold")
            ax.set_yticks(range(len(GS_ORDER)))
            ax.set_yticklabels(GS_ORDER, fontsize=7)
            ax.set_title(f"GS Progression by Field — {st}",
                         color=BRAND_COLOR, fontweight="bold", fontsize=14)
            ax.legend(fontsize=6, loc="upper left", ncol=2, framealpha=0.9)
            ax.grid(alpha=0.2, linestyle="--")
            _save(fig, f"gs_vs_dap_region_{st}")

    # 16. DAP Box Plot by GS (phase-colored)
    if len(df_dap):
        fig, ax = plt.subplots(figsize=(16, 7))
        gs_present = [g for g in GS_ORDER if g in df_dap["Cropstagemajcode"].values]
        box_colors = [PHASE_COLORS.get(get_gs_phase(g), DARK_GREY) for g in gs_present]
        data_by_gs = [df_dap[df_dap["Cropstagemajcode"] == g]["DAP"].values for g in gs_present]
        bp = ax.boxplot(data_by_gs, labels=gs_present, patch_artist=True, widths=0.6)
        for patch, color in zip(bp["boxes"], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        for median in bp["medians"]:
            median.set_color("black")
            median.set_linewidth(1.5)
        ax.set_xlabel("Growth Stage", fontsize=12, fontweight="bold")
        ax.set_ylabel("Days After Planting (DAP)", fontsize=12, fontweight="bold")
        ax.set_title("DAP Distribution per Growth Stage",
                      color=BRAND_COLOR, fontweight="bold", fontsize=14)
        handles = [plt.Rectangle((0, 0), 1, 1, fc=c, alpha=0.7) for c in PHASE_COLORS.values()]
        ax.legend(handles, PHASE_COLORS.keys(), fontsize=8, loc="upper left")
        ax.grid(axis="y", alpha=0.3, linestyle="--")
        plt.xticks(rotation=45, ha="right")
        _save(fig, "dap_boxplot_by_gs")

    # 17. Regional DAP Heatmap (Mean DAP × State × GS)
    if len(df_dap):
        dap_pivot = df_dap.groupby(["State", "Cropstagemajcode"])["DAP"].mean().unstack()
        gs_cols = [g for g in GS_ORDER if g in dap_pivot.columns]
        dap_pivot = dap_pivot.reindex(columns=gs_cols)
        if len(dap_pivot):
            # Custom annotation: show value or "NO DATA" for NaN
            annot_arr = dap_pivot.copy()
            annot_str = annot_arr.applymap(lambda v: f"{v:.0f}" if pd.notna(v) else "NO DATA")
            # Fill NaN with -1 for display, use mask for special color
            dap_display = dap_pivot.fillna(-1)
            fig, ax = plt.subplots(figsize=(14, max(3, len(dap_pivot) * 0.65)))
            sns.heatmap(dap_display, cmap="YlGn", annot=annot_str, fmt="",
                        linewidths=0.5, linecolor="white", ax=ax,
                        cbar_kws={"label": "Mean DAP"}, vmin=0)
            # Color the "NO DATA" cells grey
            for i in range(dap_pivot.shape[0]):
                for j in range(dap_pivot.shape[1]):
                    if pd.isna(dap_pivot.iloc[i, j]):
                        ax.add_patch(plt.Rectangle((j, i), 1, 1,
                                     fill=True, facecolor="#E0E0E0", edgecolor="white", lw=0.5))
                        ax.text(j + 0.5, i + 0.5, "NO DATA", ha="center", va="center",
                                fontsize=7, color=ALERT_RED, fontweight="bold")
            ax.set_title("Mean DAP by Region and Growth Stage",
                          color=BRAND_COLOR, fontweight="bold", fontsize=14)
            ax.set_xlabel("Growth Stage")
            ax.set_ylabel("Region / State")
            plt.xticks(rotation=45, ha="right")
            _save(fig, "regional_dap_heatmap")

    # 18. Regional Progression Lines (Mean GS numeric vs DAP per state)
    if len(df_dap):
        fig, ax = plt.subplots(figsize=(14, 7))
        states = sorted(df_dap["State"].unique())
        palette = sns.color_palette("husl", len(states))
        for i, st in enumerate(states):
            sdf = df_dap[df_dap["State"] == st].copy()
            if len(sdf) < 3:
                continue
            agg = sdf.groupby("DAP").agg(MeanGS=("GS_Num", "mean")).reset_index().sort_values("DAP")
            ax.plot(agg["DAP"], agg["MeanGS"], marker="o", markersize=3,
                    label=st, color=palette[i], linewidth=1.5, alpha=0.8)
        ax.set_xlabel("Days After Planting (DAP)", fontsize=12, fontweight="bold")
        ax.set_ylabel("Mean Growth Stage (numeric)", fontsize=12, fontweight="bold")
        ax.set_title("Growth Stage Progression by Region",
                      color=BRAND_COLOR, fontweight="bold", fontsize=14)
        ax.set_yticks(range(len(GS_ORDER)))
        ax.set_yticklabels(GS_ORDER, fontsize=8)
        ax.legend(fontsize=8, title="Region", loc="upper left")
        ax.grid(alpha=0.3, linestyle="--")
        _save(fig, "regional_progression_lines")

    # 19. Development Phase Timing (mean DAP bar + min-max range)
    if len(df_dap):
        phase_stats = df_dap.groupby("GS_Phase")["DAP"].agg(["mean", "min", "max"]).reindex(GS_PHASE_ORDER)
        phase_stats = phase_stats.dropna(subset=["mean"])
        if len(phase_stats):
            fig, ax = plt.subplots(figsize=(10, 6))
            x = np.arange(len(phase_stats))
            colors_ph = [PHASE_COLORS.get(p, DARK_GREY) for p in phase_stats.index]
            bars = ax.bar(x, phase_stats["mean"], color=colors_ph, edgecolor="white",
                          width=0.6, alpha=0.85)
            ax.errorbar(x, phase_stats["mean"],
                        yerr=[phase_stats["mean"] - phase_stats["min"],
                              phase_stats["max"] - phase_stats["mean"]],
                        fmt="none", ecolor=DARK_GREY, capsize=5, capthick=1.5)
            for xi, (_, row) in zip(x, phase_stats.iterrows()):
                ax.text(xi, row["mean"] + 1, f'{row["mean"]:.0f}d',
                        ha="center", va="bottom", fontsize=10, fontweight="bold")
                ax.text(xi, row["min"] - 3, f'{row["min"]:.0f}',
                        ha="center", va="top", fontsize=7, color=DARK_GREY)
                ax.text(xi, row["max"] + 1.5, f'{row["max"]:.0f}',
                        ha="center", va="bottom", fontsize=7, color=DARK_GREY)
            ax.set_xticks(x)
            ax.set_xticklabels(phase_stats.index, fontsize=10)
            ax.set_xlabel("Development Phase", fontsize=12, fontweight="bold")
            ax.set_ylabel("Days After Planting (DAP)", fontsize=12, fontweight="bold")
            ax.set_title("Development Phase Timing (Mean DAP + Range)",
                          color=BRAND_COLOR, fontweight="bold", fontsize=14)
            ax.grid(axis="y", alpha=0.3, linestyle="--")
            _save(fig, "development_phase_timing")

    # 20. DAP by Region subplots (box plot grid)
    if len(df_dap):
        states = sorted(df_dap["State"].unique())
        n_states = len(states)
        if n_states > 0:
            ncols = min(3, n_states)
            nrows = int(np.ceil(n_states / ncols))
            fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4.5 * nrows),
                                     squeeze=False)
            for idx, st in enumerate(states):
                r, c = divmod(idx, ncols)
                ax = axes[r][c]
                sdf = df_dap[df_dap["State"] == st]
                gs_present = [g for g in GS_ORDER if g in sdf["Cropstagemajcode"].values]
                if not gs_present:
                    ax.set_visible(False)
                    continue
                data_by_gs = [sdf[sdf["Cropstagemajcode"] == g]["DAP"].values for g in gs_present]
                box_colors = [PHASE_COLORS.get(get_gs_phase(g), DARK_GREY) for g in gs_present]
                bp = ax.boxplot(data_by_gs, labels=gs_present, patch_artist=True, widths=0.6)
                for patch, col in zip(bp["boxes"], box_colors):
                    patch.set_facecolor(col)
                    patch.set_alpha(0.7)
                for median in bp["medians"]:
                    median.set_color("black")
                ax.set_title(f"{st}", fontsize=10, fontweight="bold", color=BRAND_COLOR)
                ax.tick_params(axis="x", rotation=45, labelsize=7)
                ax.set_ylabel("DAP", fontsize=8)
                ax.grid(axis="y", alpha=0.3, linestyle="--")
            for idx in range(n_states, nrows * ncols):
                r, c = divmod(idx, ncols)
                axes[r][c].set_visible(False)
            fig.suptitle("DAP Distribution per Growth Stage — by Region",
                         fontsize=14, fontweight="bold", color=BRAND_COLOR, y=1.02)
            fig.tight_layout()
            _save(fig, "dap_by_region_subplots")

    # 21. GS Progression Examples (3 most + 3 least complete fields)
    if len(df_dap):
        traj_counts = (
            df_dap.groupby(TRAJECTORY_KEY)
            .agg(N_Obs=("Cropstagemajcode", "count"),
                 UniqueGS=("Cropstagemajcode", "nunique"),
                 Field=("Field", "first"))
            .reset_index()
            .sort_values(["UniqueGS", "N_Obs"], ascending=[False, False])
        )
        best3 = traj_counts.head(3)
        worst3 = traj_counts.tail(3)

        for label, subset in [("most_complete", best3), ("least_complete", worst3)]:
            if len(subset) == 0:
                continue
            n = len(subset)
            fig, axes = plt.subplots(1, n, figsize=(6 * n, 5), squeeze=False)
            for i, (_, trow) in enumerate(subset.iterrows()):
                ax = axes[0][i]
                key_vals = {k: trow[k] for k in TRAJECTORY_KEY}
                mask = pd.Series(True, index=df_dap.index)
                for k, v in key_vals.items():
                    mask &= (df_dap[k] == v)
                tdata = df_dap[mask].sort_values("Observationdate").copy()
                if len(tdata) == 0:
                    ax.set_visible(False)
                    continue
                tdata["GS_Num"] = tdata["Cropstagemajcode"].apply(gs_numeric)
                tdata = tdata.dropna(subset=["GS_Num"])
                # Plot DAP on x-axis, GS_Num on y-axis
                phase_cols = [PHASE_COLORS.get(get_gs_phase(gs), DARK_GREY)
                              for gs in tdata["Cropstagemajcode"]]
                ax.scatter(tdata["DAP"], tdata["GS_Num"], c=phase_cols, s=60,
                           edgecolors=DARK_GREY, linewidths=0.5, zorder=3)
                ax.plot(tdata["DAP"], tdata["GS_Num"], color=BRAND_COLOR,
                        linewidth=1.5, alpha=0.6, zorder=2)
                ax.set_yticks(range(len(GS_ORDER)))
                ax.set_yticklabels(GS_ORDER, fontsize=7)
                ax.set_xlabel("DAP", fontsize=9)
                ax.set_ylabel("Growth Stage", fontsize=9)
                field_short = str(trow.get("Field", ""))[:30]
                ax.set_title(f"{field_short}\n({int(trow['N_Obs'])} obs, "
                             f"{int(trow['UniqueGS'])} GS)",
                             fontsize=9, fontweight="bold", color=BRAND_COLOR)
                ax.grid(alpha=0.2, linestyle="--")
            kind_title = "Most Complete" if label == "most_complete" else "Least Complete"
            fig.suptitle(f"GS Progression — {kind_title} Fields",
                         fontsize=13, fontweight="bold", color=BRAND_COLOR, y=1.02)
            fig.tight_layout()
            _save(fig, f"gs_progression_{label}")

    print(f"    Generated {len(plots)} plots")
    return plots


# ==================== 10. EXCEL EXPORTS ====================

def _write_sheet_safe(writer, df_data, sheet_name):
    """Write a DataFrame to a sheet only if non-empty."""
    if df_data is not None and len(df_data) > 0:
        df_data.to_excel(writer, sheet_name=sheet_name[:31], index=False)


def _legend_df():
    """Return a DataFrame explaining each QA check and methodology."""
    return pd.DataFrame({
        "Check": [
            "GS Trio Completeness",
            "GS Triplet Consistency",
            "GS Regression",
            "Development Gaps",
            "Assessment Pkg Completeness",
            "Disease Monotonicity",
            "RC Monotonicity",
            "RC Open Backlog",
            "GS Coverage",
        ],
        "Description": [
            "Checks presence of Cropstagemincode, Cropstagemajcode, Cropstagemaxcode per BUA.",
            "Verifies min ≤ maj ≤ max within each BUA.",
            "Detects when GS goes backward over time (per min/maj/max separately).",
            f"Flags when GS jumps > {DEV_GAP_THRESHOLD} stages between consecutive observations.",
            "Verifies all 4 required assessments present per BUA; disease optional after 2nd RC=9.",
            "Checks disease severity values are non-decreasing over time per trajectory.",
            "Checks Row Closing values are non-decreasing over time per trajectory.",
            "Lists trajectories where latest RC < 9 (not yet fully closed).",
            "Measures breadth and phase balance of GS codes captured.",
        ],
        "Methodology": [
            "% of BUAs with all three GS fields non-null.",
            "Numeric GS index comparison within same record.",
            "Time-ordered comparison within trajectory, per GS column.",
            "Numeric stage difference between consecutive observations.",
            "Disease optional after 2nd RC=9 occurrence per trajectory.",
            "Applies before & after 2nd RC=9; no false positives from optional absence.",
            "Values must never decrease along the trajectory.",
            "Descriptive operational backlog for open row closings.",
            "Unique stages / 20 expected; phase balance = avg of 4 phase coverages.",
        ],
    })


def export_excel_by_region(
    df: pd.DataFrame,
    qa: Dict[str, pd.DataFrame],
    gs_cov: Dict[str, pd.DataFrame],
    rc_backlog: pd.DataFrame,
    field_summary: pd.DataFrame,
    unknown_summary: pd.DataFrame,
    output_dir: str,
):
    """One Excel workbook per state."""
    print("\n  [EXPORT] Excel by Region ...")
    os.makedirs(output_dir, exist_ok=True)
    date_tag = datetime.now().strftime("%Y%m%d")

    for state in sorted(df["State"].unique()):
        safe_name = re.sub(r'[<>:"/\\|?*]', "_", state)
        path = os.path.join(output_dir, f"QA_Region_{safe_name}_{date_tag}.xlsx")
        def _filter(issue_df, col="State"):
            if issue_df is None or len(issue_df) == 0:
                return pd.DataFrame()
            if col in issue_df.columns:
                return issue_df[issue_df[col] == state].copy()
            return pd.DataFrame()

        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            # Region summary (KPI-like)
            region_cov = gs_cov.get("by_region", pd.DataFrame())
            region_row = region_cov[region_cov["State"] == state] if len(region_cov) else pd.DataFrame()
            if len(region_row):
                region_row.to_excel(writer, sheet_name="REGION_SUMMARY", index=False)

            _write_sheet_safe(writer, _filter(qa.get("gs_trio_completeness")), "GS_COMPLETENESS")

            cov_hm = gs_cov.get("heatmap_region")
            if cov_hm is not None and state in cov_hm.index:
                cov_hm.loc[[state]].to_excel(writer, sheet_name="GS_COVERAGE")

            _write_sheet_safe(writer, _filter(qa.get("gs_triplet_consistency")), "GS_TRIPLET_CONSISTENCY")
            _write_sheet_safe(writer, _filter(qa.get("gs_regression")), "GS_REGRESSION")
            _write_sheet_safe(writer, _filter(qa.get("dev_gaps")), "DEVELOPMENT_GAPS")
            _write_sheet_safe(writer, _filter(qa.get("assessment_pkg")), "ASSESSMENT_COMPLETENESS")
            _write_sheet_safe(writer, _filter(qa.get("disease_mono")), "DISEASE_MONOTONICITY")
            _write_sheet_safe(writer, _filter(qa.get("rc_mono")), "RC_MONOTONICITY")
            _write_sheet_safe(writer, _filter(rc_backlog), "RC_OPEN_BACKLOG")

            # Unknown assignments within this state
            unk = unknown_summary[
                (unknown_summary["State"] == state) | (unknown_summary["FlaggedState"])
            ] if len(unknown_summary) else pd.DataFrame()
            _write_sheet_safe(writer, unk, "UNKNOWN_ASSIGNMENTS")
            _write_sheet_safe(writer, _legend_df(), "LEGEND_METHOD")

        print(f"    {state}: {path}")


def export_excel_by_responsible(
    df: pd.DataFrame,
    qa: Dict[str, pd.DataFrame],
    gs_cov: Dict[str, pd.DataFrame],
    rc_backlog: pd.DataFrame,
    field_summary: pd.DataFrame,
    trajectory_summary: pd.DataFrame,
    output_dir: str,
):
    """One Excel workbook per responsible."""
    print("\n  [EXPORT] Excel by Responsible ...")
    os.makedirs(output_dir, exist_ok=True)
    date_tag = datetime.now().strftime("%Y%m%d")

    for resp in sorted(df["Responsible"].unique()):
        safe_name = re.sub(r'[<>:"/\\|?*]', "_", resp)
        path = os.path.join(output_dir, f"QA_Resp_{safe_name}_{date_tag}.xlsx")
        def _filter(issue_df, col="Responsible"):
            if issue_df is None or len(issue_df) == 0:
                return pd.DataFrame()
            if col in issue_df.columns:
                return issue_df[issue_df[col] == resp].copy()
            return pd.DataFrame()

        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            # Responsible summary
            resp_cov = gs_cov.get("by_responsible", pd.DataFrame())
            resp_row = resp_cov[resp_cov["Responsible"] == resp] if len(resp_cov) else pd.DataFrame()
            _write_sheet_safe(writer, resp_row, "RESPONSIBLE_SUMMARY")

            fs = field_summary[field_summary["Responsible"] == resp] if len(field_summary) else pd.DataFrame()
            _write_sheet_safe(writer, fs, "FIELD_SUMMARY")

            ts = trajectory_summary[trajectory_summary["Responsible"] == resp] if len(trajectory_summary) else pd.DataFrame()
            _write_sheet_safe(writer, ts, "TRAJECTORY_SUMMARY")

            _write_sheet_safe(writer, _filter(qa.get("gs_trio_completeness")), "GS_COMPLETENESS")
            _write_sheet_safe(writer, _legend_df(), "LEGEND_METHOD")
            cov_hm = gs_cov.get("heatmap_responsible")
            if cov_hm is not None and resp in cov_hm.index:
                cov_hm.loc[[resp]].to_excel(writer, sheet_name="GS_COVERAGE")

            _write_sheet_safe(writer, _filter(qa.get("gs_triplet_consistency")), "GS_TRIPLET_CONSISTENCY")
            _write_sheet_safe(writer, _filter(qa.get("gs_regression")), "GS_REGRESSION")
            _write_sheet_safe(writer, _filter(qa.get("dev_gaps")), "DEVELOPMENT_GAPS")
            _write_sheet_safe(writer, _filter(qa.get("assessment_pkg")), "ASSESSMENT_COMPLETENESS")
            _write_sheet_safe(writer, _filter(qa.get("disease_mono")), "DISEASE_MONOTONICITY")
            _write_sheet_safe(writer, _filter(qa.get("rc_mono")), "RC_MONOTONICITY")
            _write_sheet_safe(writer, _filter(rc_backlog), "RC_OPEN_BACKLOG")
            _write_sheet_safe(writer, _legend_df(), "LEGEND_METHOD")

        print(f"    {resp}: {path}")


def export_consolidated_kpi_excel(
    kpi_resp: pd.DataFrame,
    kpi_region: pd.DataFrame,
    unknown_summary: pd.DataFrame,
    output_dir: str,
):
    """Single consolidated KPI workbook."""
    print("\n  [EXPORT] Consolidated KPI Excel ...")
    os.makedirs(output_dir, exist_ok=True)
    date_tag = datetime.now().strftime("%Y%m%d")
    path = os.path.join(output_dir, f"KPI_Consolidated_{date_tag}.xlsx")

    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        kpi_resp.to_excel(writer, sheet_name="KPI_By_Responsible", index=False)
        kpi_region.to_excel(writer, sheet_name="KPI_By_Region", index=False)

        # Config sheet
        config = pd.DataFrame({
            "Parameter": [
                "Protocol", "Farm Filter", "GS Order",
                "Dev Gap Threshold", "Brand Color",
                "KPI Weight: GS Trio", "KPI Weight: Assess Pkg",
                "KPI Weight: GS Coverage", "KPI Weight: RC Completion",
                "KPI Weight: Inv Issue Burden",
                "Run Date",
            ],
            "Value": [
                PROTOCOL_FILTER, FARM_FILTER, " → ".join(GS_ORDER),
                f"> {DEV_GAP_THRESHOLD} stages", BRAND_COLOR,
                f"{KPI_WEIGHTS['gs_trio_completeness']:.0%}",
                f"{KPI_WEIGHTS['assessment_pkg_completeness']:.0%}",
                f"{KPI_WEIGHTS['gs_coverage']:.0%}",
                f"{KPI_WEIGHTS['rc_completion']:.0%}",
                f"{KPI_WEIGHTS['inverse_issue_burden']:.0%}",
                datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            ],
        })
        config.to_excel(writer, sheet_name="CONFIG", index=False)
        _write_sheet_safe(writer, unknown_summary, "VALIDATION_UNKNOWNS")

    print(f"    {path}")


# ==================== 11. PDF REPORT ====================

# ==================== 11. PDF REPORT ====================

def _pdf_styles():
    """Return a dict of custom ParagraphStyles."""
    base = getSampleStyleSheet()
    return {
        "title": ParagraphStyle(
            "PDFTitle", parent=base["Heading1"], fontSize=24,
            textColor=rl_colors.HexColor(BRAND_COLOR), alignment=TA_CENTER,
            spaceAfter=20, fontName="Helvetica-Bold",
        ),
        "h2": ParagraphStyle(
            "H2", parent=base["Heading2"], fontSize=15,
            textColor=rl_colors.HexColor(BRAND_COLOR), spaceAfter=10,
            spaceBefore=6, fontName="Helvetica-Bold",
        ),
        "h3": ParagraphStyle(
            "H3", parent=base["Heading3"], fontSize=11,
            textColor=rl_colors.HexColor(DARK_GREY), spaceAfter=6,
            fontName="Helvetica-Bold",
        ),
        "body": ParagraphStyle(
            "Body", parent=base["Normal"], fontSize=9,
            textColor=rl_colors.HexColor(DARK_GREY), spaceAfter=4,
            leading=12,
        ),
        "body_small": ParagraphStyle(
            "BodySmall", parent=base["Normal"], fontSize=8,
            textColor=rl_colors.HexColor(DARK_GREY), spaceAfter=3,
            leading=10,
        ),
        "subtitle": ParagraphStyle(
            "Subtitle", parent=base["Normal"], fontSize=13,
            alignment=TA_CENTER, textColor=rl_colors.HexColor(DARK_GREY),
        ),
        "kpi_good": ParagraphStyle(
            "KPIGood", parent=base["Normal"], fontSize=8,
            textColor=rl_colors.HexColor(SUCCESS_GREEN), fontName="Helvetica-Bold",
        ),
        "kpi_warn": ParagraphStyle(
            "KPIWarn", parent=base["Normal"], fontSize=8,
            textColor=rl_colors.HexColor(WARNING_YELLOW), fontName="Helvetica-Bold",
        ),
        "kpi_bad": ParagraphStyle(
            "KPIBad", parent=base["Normal"], fontSize=8,
            textColor=rl_colors.HexColor(ALERT_RED), fontName="Helvetica-Bold",
        ),
    }


def _make_pdf_table(data_rows, col_widths=None, header_color=BRAND_COLOR):
    """Utility to build a ReportLab table from list-of-lists."""
    if not data_rows:
        return Spacer(1, 0)
    styles = getSampleStyleSheet()
    hdr_style = ParagraphStyle(
        "TblHdr", parent=styles["Normal"], fontName="Helvetica-Bold",
        fontSize=7.5, leading=9, textColor=rl_colors.whitesmoke, alignment=1,
    )
    cell_style = ParagraphStyle(
        "TblCell", parent=styles["Normal"], fontName="Helvetica",
        fontSize=7.5, leading=9, textColor=rl_colors.HexColor(DARK_GREY),
        wordWrap="CJK",
    )
    wrapped = []
    for ri, row in enumerate(data_rows):
        wr = []
        for cell in row:
            txt = "" if cell is None else str(cell).replace("\n", "<br/>")
            wr.append(Paragraph(txt, hdr_style if ri == 0 else cell_style))
        wrapped.append(wr)
    n_cols = len(data_rows[0])
    if col_widths is None:
        avail = 6.5 * inch
        col_widths = [avail / n_cols] * n_cols
    tbl = Table(wrapped, colWidths=col_widths, repeatRows=1)
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), rl_colors.HexColor(header_color)),
        ("TEXTCOLOR", (0, 0), (-1, 0), rl_colors.whitesmoke),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 7.5),
        ("GRID", (0, 0), (-1, -1), 0.5, rl_colors.HexColor(LIGHT_GREY)),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [rl_colors.white, rl_colors.HexColor("#F5F5F5")]),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    return tbl


def generate_pdf_report(
    df: pd.DataFrame,
    qa: Dict[str, pd.DataFrame],
    gs_cov: Dict[str, pd.DataFrame],
    kpi_resp: pd.DataFrame,
    kpi_region: pd.DataFrame,
    rc_backlog: pd.DataFrame,
    plots: Dict[str, str],
    unknown_summary: pd.DataFrame,
    output_dir: str,
):
    """Generate the comprehensive technical PDF report with v2.1-quality styling."""
    print("\n  [PDF] Generating technical report ...")
    os.makedirs(output_dir, exist_ok=True)
    date_tag = datetime.now().strftime("%Y%m%d")
    pdf_path = os.path.join(output_dir, f"TM2429_QA_Report_{date_tag}.pdf")

    doc = SimpleDocTemplate(
        pdf_path, pagesize=A4,
        topMargin=0.5 * inch, bottomMargin=0.5 * inch,
        leftMargin=0.6 * inch, rightMargin=0.6 * inch,
    )
    story = []
    S = _pdf_styles()

    def _img(name, w=6.2 * inch, h=3.8 * inch):
        p = plots.get(name)
        if p and os.path.exists(p):
            story.append(Image(p, width=w, height=h))
            story.append(Spacer(1, 0.15 * inch))

    def _status_icon(value, thresholds, metric_type="compliance"):
        """Return colored status text for PDF cells."""
        lbl = _priority_label(value, thresholds, metric_type)
        if lbl == "Good":
            return f'<font color="{SUCCESS_GREEN}"><b>● {lbl}</b></font>'
        elif lbl == "Warning":
            return f'<font color="{WARNING_YELLOW}"><b>● {lbl}</b></font>'
        return f'<font color="{ALERT_RED}"><b>● {lbl}</b></font>'

    # ─────────── 1. TITLE PAGE ───────────
    story.append(Spacer(1, 1.2 * inch))
    story.append(Paragraph(PDF_TITLE, S["title"]))
    story.append(Spacer(1, 0.25 * inch))
    story.append(Paragraph("Technical QA &amp; KPI Report", S["subtitle"]))
    story.append(Spacer(1, 0.5 * inch))

    # Data Integrity Scorecard on title page
    total_issues = sum(len(v) for v in qa.values())
    n_open_rc = int(rc_backlog["OpenRC"].sum()) if len(rc_backlog) else 0
    avg_score = kpi_resp["Overall_Score"].mean() if len(kpi_resp) else 0
    scorecard = [
        ["Metric", "Value"],
        ["Protocol", PROTOCOL_FILTER],
        ["Season", FARM_FILTER],
        ["Generated", datetime.now().strftime("%Y-%m-%d %H:%M:%S")],
        ["Total Fields", str(df["Field"].nunique())],
        ["Total Trials", str(df["Trialnumber"].nunique())],
        ["Responsibles", str(df["Responsible"].nunique())],
        ["Regions", str(df["State"].nunique())],
        ["Total QA Issues", str(total_issues)],
        ["Open RC Backlog", str(n_open_rc)],
        ["Avg Overall Score", f"{avg_score:.1f}"],
    ]
    story.append(_make_pdf_table(scorecard, col_widths=[2.5 * inch, 3 * inch]))
    story.append(PageBreak())

    # ─────────── 2. SCOPE AND PROTOCOL ───────────
    story.append(Paragraph("1. Scope &amp; Protocol Filter", S["h2"]))
    story.append(Paragraph(
        f"This report covers <b>only</b> protocol <i>{PROTOCOL_FILTER}</i>. "
        "It evaluates phenological data quality, completeness of plant development capture, "
        "GS progression consistency, disease and RC assessment integrity, and regional "
        "data collection quality. The underlying goal is to improve ML model training "
        "by ensuring complete and consistent phenology capture.",
        S["body"],
    ))
    story.append(Spacer(1, 0.2 * inch))

    # ─────────── 3. METHODOLOGY ───────────
    story.append(Paragraph("2. Methodology &amp; Definitions", S["h2"]))
    story.append(Paragraph(
        "<b>Basic Unit of Analysis (BUA):</b> "
        "Trialnumber + Replicate + Sampleno + Factorlevel + Locationcode + Observationdate. "
        "Every flagged QA event preserves this grain.<br/>"
        f"<b>GS Scale:</b> {' → '.join(GS_ORDER)}",
        S["body"],
    ))
    story.append(Paragraph(
        "<b>Phases:</b> Early Vegetative (VS–V3), Late Vegetative (V4–V9), "
        "Early Reproductive (R1–R4), Late Reproductive (R5–R8).",
        S["body"],
    ))
    story.append(Paragraph(
        "<b>Composite KPI Score:</b> "
        "30% GS Trio Completeness + 20% Assessment Pkg Completeness + "
        "20% GS Coverage + 15% RC Completion + 15% Inverse Issue Burden. "
        "Each component normalized to 0–100.",
        S["body"],
    ))
    story.append(Paragraph(
        f"<b>Development Gap Threshold:</b> Flag if jump &gt; {DEV_GAP_THRESHOLD} stages.",
        S["body"],
    ))
    story.append(Paragraph(
        "Responsible names are parsed from the Field string using pattern matching. "
        "States are extracted via word-boundary regex.",
        S["body"],
    ))
    story.append(Paragraph(
        "<b>Disease Cutoff Rule:</b> Disease assessments are obligatory only BEFORE "
        "the second occurrence of RC=9 in the trajectory. After the second RC=9, "
        "absence is NOT penalized, but observed values must still be monotonic.",
        S["body"],
    ))
    story.append(PageBreak())

    # ─────────── 4. VALIDATION ───────────
    story.append(Paragraph("3. Responsible &amp; Region Assignment Validation", S["h2"]))
    story.append(Paragraph(
        "Responsible names are parsed from the Field string using pattern matching. "
        "States are extracted via word-boundary regex. Any unmatched entries are "
        "flagged below for manual review.",
        S["body"],
    ))
    if len(unknown_summary):
        rows = [["Field", "Responsible", "State", "Issue"]]
        for _, r in unknown_summary.iterrows():
            issue_parts = []
            if r["FlaggedResponsible"]:
                issue_parts.append("UNASSIGNED Responsible")
            if r["FlaggedState"]:
                issue_parts.append("UNK State")
            rows.append([r["Field"], r["Responsible"], r["State"], " / ".join(issue_parts)])
        story.append(_make_pdf_table(rows, col_widths=[2.5 * inch, 1 * inch, 0.6 * inch, 1.8 * inch]))
    else:
        story.append(Paragraph(
            '<font color="' + SUCCESS_GREEN + '"><b>✓ All fields successfully matched '
            "to a responsible and state.</b></font>", S["body"],
        ))
    story.append(PageBreak())

    # ─────────── 5. GS COMPLETENESS ───────────
    story.append(Paragraph("4. GS Trio Completeness", S["h2"]))
    story.append(Paragraph(
        "<b>What:</b> Checks if all three GS fields (min, maj, max) are filled per observation.<br/>"
        "<b>Why:</b> Incomplete GS trios degrade ML training data quality.<br/>"
        "<b>Interpretation:</b> Higher is better; target ≥ 90%. Green = ≥90%, "
        "Yellow = 70–89%, Red = &lt;70%.",
        S["body"],
    ))
    _img("gs_completeness_chart")
    story.append(PageBreak())

    # ─────────── 6. GS COVERAGE ───────────
    story.append(Paragraph("5. GS Coverage", S["h2"]))
    story.append(Paragraph(
        "<b>What:</b> Measures breadth of GS codes captured across the expected 20-stage scale.<br/>"
        "<b>Why:</b> Good phenology models require balanced representation across all stages.<br/>"
        "<b>Interpretation:</b> Dark cells = captured, white = missing. "
        "Vertical lines separate phenological phases.",
        S["body"],
    ))
    _img("gs_coverage_heatmap_region", w=6.5 * inch, h=3.5 * inch)
    _img("gs_coverage_heatmap_responsible", w=6.5 * inch, h=3.5 * inch)
    story.append(PageBreak())
    _img("gs_phase_coverage")
    story.append(PageBreak())

    # ─────────── 6b. GS DISTRIBUTION BY PHASE ───────────
    story.append(Paragraph("5b. Growth Stage Distribution by Phase", S["h2"]))
    story.append(Paragraph(
        "<b>What:</b> Shows the count of records for each growth stage, colored by development phase.<br/>"
        "<b>Why:</b> Reveals which stages are over- or under-represented in the data collection.",
        S["body"],
    ))
    _img("gs_distribution_by_phase", w=6.5 * inch, h=3.8 * inch)
    story.append(PageBreak())

    # ─────────── 6c. GS RECORDS BY STATE × PHASE ───────────
    story.append(Paragraph("5c. GS Records by State and Development Phase", S["h2"]))
    story.append(Paragraph(
        "<b>What:</b> Heatmap of observation counts per state and development phase.<br/>"
        "<b>Why:</b> Identifies regional gaps — regions with very low counts in certain "
        "phases need targeted data collection effort.",
        S["body"],
    ))
    _img("gs_records_state_phase", w=6 * inch, h=3.5 * inch)
    story.append(PageBreak())

    # ─────────── 6d. GS vs DAP ANALYSIS (★ MOST IMPORTANT) ───────────
    story.append(Paragraph("5d. Growth Stage vs Days After Planting", S["h2"]))
    story.append(Paragraph(
        "<b>What:</b> Scatter plot of GS (numeric) vs DAP, colored by region, "
        "with quadratic trend and phase background bands.<br/>"
        "<b>Why:</b> This is the <b>most critical visualization</b> — it validates "
        "whether observed growth stages align with expected phenological timing. "
        "Outliers indicate data entry errors or unusual field conditions.",
        S["body"],
    ))
    _img("gs_vs_dap_scatter", w=6.5 * inch, h=4 * inch)

    story.append(Paragraph("DAP Distribution per Growth Stage", S["h3"]))
    story.append(Paragraph(
        "Box plots showing the DAP range for each growth stage. "
        "Narrow boxes indicate consistent timing; wide boxes indicate high variability.",
        S["body"],
    ))
    _img("dap_boxplot_by_gs", w=6.5 * inch, h=3.5 * inch)
    story.append(PageBreak())

    # ─────────── 6d-ii. PER-REGION GS vs DAP ───────────
    story.append(Paragraph("5d-ii. GS Progression by Field — per Region", S["h2"]))
    story.append(Paragraph(
        "<b>What:</b> One plot per region showing each trial/field as a separate line "
        "plotting GS over DAP.<br/>"
        "<b>Why:</b> Allows assessing individual field progression within each region. "
        "Divergent lines indicate fields that may have data entry issues or unusual conditions.",
        S["body"],
    ))
    # Insert one image per region
    for st in sorted(df["State"].unique()):
        _img(f"gs_vs_dap_region_{st}", w=6.5 * inch, h=3.8 * inch)
    story.append(PageBreak())

    # ─────────── 6e. REGIONAL DAP ANALYSIS ───────────
    story.append(Paragraph("5e. Regional DAP Analysis", S["h2"]))
    story.append(Paragraph(
        "<b>Mean DAP Heatmap:</b> Shows the average DAP for each growth stage in each region. "
        "Differences reveal regional timing patterns — southern regions may advance slower "
        "due to later planting.",
        S["body"],
    ))
    _img("regional_dap_heatmap", w=6.5 * inch, h=3.5 * inch)

    story.append(Paragraph("GS Progression by Region", S["h3"]))
    story.append(Paragraph(
        "Mean GS numeric index plotted over DAP for each region. "
        "Steeper curves = faster crop development.",
        S["body"],
    ))
    _img("regional_progression_lines", w=6.5 * inch, h=3.5 * inch)

    # --- Regression fields by region table ---
    reg_df = qa.get("gs_regression", pd.DataFrame())
    if not reg_df.empty and "State" in reg_df.columns:
        reg_by_region = (
            reg_df.groupby(["State", "Responsible", "Field"])
            .size().reset_index(name="Regression Events")
            .sort_values(["State", "Regression Events"], ascending=[True, False])
        )
        if not reg_by_region.empty:
            story.append(Paragraph("Fields Contributing to GS Regressions by Region", S["h3"]))
            rows = [["Region", "Responsible", "Field", "Events"]]
            for _, r in reg_by_region.head(20).iterrows():
                rows.append([str(r["State"]), str(r["Responsible"]),
                             str(r["Field"]), str(r["Regression Events"])])
            story.append(_make_pdf_table(rows))
    story.append(PageBreak())

    # ─────────── 6f. DEVELOPMENT PHASE TIMING ───────────
    story.append(Paragraph("5f. Development Phase Timing", S["h2"]))
    story.append(Paragraph(
        "<b>What:</b> Mean DAP per development phase with min-max error bars.<br/>"
        "<b>Why:</b> Establishes expected phenological timing benchmarks. "
        "Phases should have clear temporal separation.",
        S["body"],
    ))
    _img("development_phase_timing", w=5.5 * inch, h=3.5 * inch)

    story.append(Paragraph("DAP by Region — Detailed Subplots", S["h3"]))
    story.append(Paragraph(
        "Per-region box plots of DAP per GS for granular regional analysis.",
        S["body"],
    ))
    _img("dap_by_region_subplots", w=6.5 * inch, h=5 * inch)
    story.append(PageBreak())

    # ─────────── 6g. GS PROGRESSION EXAMPLES ───────────
    story.append(Paragraph("5g. GS Progression — Best &amp; Worst Examples", S["h2"]))
    story.append(Paragraph(
        "<b>What:</b> Side-by-side comparison of the 3 most complete and 3 least complete "
        "field trajectories (by unique GS captured and observation count).<br/>"
        "<b>Why:</b> Illustrates the quality gap — best cases show smooth, complete "
        "phenological curves; worst cases show sparse or stalled data collection.",
        S["body"],
    ))
    _img("gs_progression_most_complete", w=6.5 * inch, h=3.2 * inch)
    story.append(Paragraph("Least Complete Field Trajectories", S["h3"]))
    _img("gs_progression_least_complete", w=6.5 * inch, h=3.2 * inch)
    story.append(PageBreak())

    # ─────────── 7. GS REGRESSION ───────────
    story.append(Paragraph("6. GS Regression", S["h2"]))
    n_reg = len(qa.get("gs_regression", []))
    story.append(Paragraph(
        "<b>What:</b> Flags when a GS value goes backward over time within a trajectory.<br/>"
        "<b>Why:</b> GS should only advance; regressions indicate data entry errors.<br/>"
        f"<b>Found:</b> <font color='{ALERT_RED}'><b>{n_reg}</b></font> regression events.",
        S["body"],
    ))
    _img("gs_regression_dist")

    # Regression summary table by region
    gs_reg_qa = qa.get("gs_regression", pd.DataFrame())
    if len(gs_reg_qa) and all(c in gs_reg_qa.columns for c in ["Responsible", "State", "Field"]):
        reg_summary = (
            gs_reg_qa.groupby(["State", "Responsible", "Field"])
            .size().reset_index(name="Regressions")
            .sort_values("Regressions", ascending=False)
        )
        if len(reg_summary):
            story.append(Paragraph("Fields with GS Regressions:", S["h3"]))
            tbl_rows = [["Region", "Responsible", "Field", "Count"]]
            for _, rr in reg_summary.head(15).iterrows():
                tbl_rows.append([rr["State"], rr["Responsible"],
                                str(rr["Field"])[:40], str(int(rr["Regressions"]))])
            story.append(_make_pdf_table(tbl_rows,
                col_widths=[0.6*inch, 1.0*inch, 3.2*inch, 0.7*inch],
                header_color=ALERT_RED))
    story.append(PageBreak())

    # ─────────── 8. DEVELOPMENT GAPS ───────────
    story.append(Paragraph("7. Development Gaps", S["h2"]))
    n_gaps = len(qa.get("dev_gaps", []))
    story.append(Paragraph(
        f"<b>What:</b> Flags jumps &gt; {DEV_GAP_THRESHOLD} stages between consecutive observations.<br/>"
        "<b>Why:</b> Large GS jumps mean missed intermediate phenological stages.<br/>"
        f"<b>Found:</b> <font color='{WARNING_YELLOW}'><b>{n_gaps}</b></font> development gap events.",
        S["body"],
    ))
    _img("dev_gap_dist")

    # Development gap detail table
    dg_qa = qa.get("dev_gaps", pd.DataFrame())
    if len(dg_qa) and all(c in dg_qa.columns for c in ["Responsible", "Field"]):
        dg_summary = (
            dg_qa.groupby(["Responsible", "Field"])
            .size().reset_index(name="Gaps")
            .sort_values("Gaps", ascending=False)
        )
        if len(dg_summary):
            story.append(Paragraph("Largest Gap Contributors:", S["h3"]))
            tbl_rows = [["Responsible", "Field", "Gap Events"]]
            for _, rr in dg_summary.head(10).iterrows():
                tbl_rows.append([rr["Responsible"], str(rr["Field"])[:40], str(int(rr["Gaps"]))])
            story.append(_make_pdf_table(tbl_rows,
                col_widths=[1.0*inch, 3.5*inch, 0.8*inch],
                header_color=WARNING_YELLOW))
    story.append(PageBreak())

    # ─────────── 9. ASSESSMENT PACKAGE ───────────
    story.append(Paragraph("8. Assessment Package Completeness", S["h2"]))
    story.append(Paragraph(
        "<b>What:</b> Checks if all 4 required assessments are present per BUA.<br/>"
        "<b>Rule:</b> Disease assessments are <i>obligatory only before the 2nd RC=9</i>. "
        "After that, absence is not penalized.",
        S["body"],
    ))
    _img("assessment_completeness_chart")

    # Assessment pending issues table
    pkg_qa = qa.get("assessment_pkg", pd.DataFrame())
    if len(pkg_qa):
        story.append(Paragraph("Pending Assessment Issues:", S["h3"]))
        tbl_rows = [["Responsible", "Field", "Missing Assessments"]]
        pkg_summary = pkg_qa.copy()
        if "MissingAssessments" in pkg_summary.columns and "Responsible" in pkg_summary.columns:
            pkg_agg = (pkg_summary.groupby(["Responsible", "Field"])
                       .agg(Missing=("MissingAssessments", lambda x: "; ".join(sorted(set(x)))))
                       .reset_index()
                       .sort_values("Responsible"))
            for _, rr in pkg_agg.head(20).iterrows():
                tbl_rows.append([rr["Responsible"], str(rr["Field"])[:35], str(rr["Missing"])[:50]])
            story.append(_make_pdf_table(tbl_rows,
                col_widths=[1.0*inch, 2.5*inch, 2.5*inch]))
    story.append(PageBreak())

    # ─────────── 10. COMPLIANCE DUAL PANEL ───────────
    story.append(Paragraph("8b. Compliance Metrics Comparison", S["h2"]))
    story.append(Paragraph(
        "Side-by-side comparison of GS Trio Completeness and Assessment Package Completeness. "
        "Green threshold at 90%.",
        S["body"],
    ))
    _img("compliance_dual_panel", w=6.8 * inch, h=4.2 * inch)
    story.append(PageBreak())

    # ─────────── 11. DISEASE QA ───────────
    story.append(Paragraph("9. Disease QA Logic Relative to RC", S["h2"]))
    n_dis = len(qa.get("disease_mono", []))
    story.append(Paragraph(
        "<b>What:</b> Checks monotonicity of disease severity values over time.<br/>"
        "<b>Rule:</b> Pre-2nd RC=9: disease required + monotonic. "
        "Post-2nd RC=9: disease optional, but if present must still be monotonic.<br/>"
        f"<b>Found:</b> <font color='{ALERT_RED}'><b>{n_dis}</b></font> violations.",
        S["body"],
    ))
    _img("disease_mono_chart")

    # Disease violations detail table
    dis_qa = qa.get("disease_mono", pd.DataFrame())
    if len(dis_qa):
        story.append(Paragraph("Disease Monotonicity Violations:", S["h3"]))
        tbl_rows = [["Responsible", "Field", "Assessment", "Prev Value", "Curr Value", "Post 2nd RC9"]]
        for _, rr in dis_qa.head(15).iterrows():
            tbl_rows.append([
                str(rr.get("Responsible", "")),
                str(rr.get("Field", ""))[:30],
                str(rr.get("Assessmenttypecode", ""))[:20],
                str(rr.get("PrevValue", "")),
                str(rr.get("CurrValue", "")),
                str(rr.get("PostSecondRC9", "")),
            ])
        story.append(_make_pdf_table(tbl_rows,
            col_widths=[0.8*inch, 1.5*inch, 1.2*inch, 0.7*inch, 0.7*inch, 0.8*inch],
            header_color=ALERT_RED))
    elif n_dis == 0:
        story.append(Paragraph(
            f'<font color="{SUCCESS_GREEN}"><b>\u2713 No disease monotonicity violations found.</b></font>',
            S["body"]))
    story.append(PageBreak())

    # ─────────── 12. RC PROGRESSION ───────────
    story.append(Paragraph("10. RC Progression &amp; Open Backlog", S["h2"]))
    n_rc_mono = len(qa.get("rc_mono", []))
    story.append(Paragraph(
        f"RC monotonicity violations: <b>{n_rc_mono}</b><br/>"
        f"Open RC trajectories (latest RC &lt; 9): "
        f"<font color='{ALERT_RED}'><b>{n_open_rc}</b></font>",
        S["body"],
    ))
    _img("rc_vs_gs_by_region", w=6.5 * inch, h=3.8 * inch)
    _img("rc_progress_chart")
    _img("rc_open_backlog_chart")

    # RC Open Backlog detail table
    if len(rc_backlog):
        open_rc = rc_backlog[rc_backlog["OpenRC"] == 1].copy()
        if len(open_rc):
            story.append(Paragraph("Open RC Backlog — Detail:", S["h3"]))
            open_rc_sorted = open_rc.sort_values("DaysSinceLatest", ascending=False)
            tbl_rows = [["Responsible", "Field", "Latest GS", "Latest RC", "Days Since"]]
            for _, rr in open_rc_sorted.head(20).iterrows():
                tbl_rows.append([
                    str(rr.get("Responsible", "")),
                    str(rr.get("Field", ""))[:35],
                    str(rr.get("LatestGS", "")),
                    str(int(rr.get("LatestRC", 0))) if pd.notna(rr.get("LatestRC")) else "—",
                    str(int(rr.get("DaysSinceLatest", 0))) if pd.notna(rr.get("DaysSinceLatest")) else "—",
                ])
            story.append(_make_pdf_table(tbl_rows,
                col_widths=[0.9*inch, 2.2*inch, 0.7*inch, 0.7*inch, 0.7*inch],
                header_color=ALERT_RED))
    story.append(PageBreak())

    # ─────────── 13. ISSUES DISTRIBUTION ───────────
    story.append(Paragraph("10b. QA Issues Distribution", S["h2"]))
    story.append(Paragraph(
        "Grouped bar chart showing the count of each issue type per responsible person. "
        "Lower counts are better.",
        S["body"],
    ))
    _img("issues_by_type_responsible", w=6.8 * inch, h=4 * inch)
    story.append(PageBreak())

    # ─────────── 14. KPI SUMMARY ───────────
    story.append(Paragraph("11. KPI Summary", S["h2"]))
    story.append(Paragraph(
        "Composite scores computed using documented weights. "
        "This is a <b>technical QA dashboard</b>, not a leaderboard.",
        S["body"],
    ))
    _img("kpi_overall_score")

    # KPI Ranking Table — Responsible (v2.1 style)
    if len(kpi_resp):
        story.append(Paragraph("Performance Ranking — by Responsible", S["h3"]))
        header = ["#", "Responsible", "Fields", "GS Trio %", "Assess %",
                   "GS Cov %", "RC %", "Score", "Status"]
        rows = [header]
        for rank, (_, r) in enumerate(kpi_resp.iterrows(), 1):
            score = r["Overall_Score"]
            if score >= 80:
                status = "Excellent"
            elif score >= 60:
                status = "Good"
            elif score >= 40:
                status = "Needs Improvement"
            else:
                status = "Critical"
            rows.append([
                str(rank), r["Responsible"], str(r["Fields"]),
                f"{r['GS_Trio_Completeness']:.0f}",
                f"{r['Assessment_Pkg_Completeness']:.0f}",
                f"{r['GS_Coverage_Score']:.0f}",
                f"{r['RC_Completion_Rate']:.0f}",
                f"{score:.0f}", status,
            ])
        story.append(_make_pdf_table(
            rows,
            col_widths=[0.3*inch, 1.0*inch, 0.45*inch, 0.65*inch, 0.65*inch,
                        0.65*inch, 0.5*inch, 0.5*inch, 1.0*inch],
        ))
        story.append(Spacer(1, 0.2 * inch))

    # KPI Table — Region
    if len(kpi_region):
        story.append(Paragraph("Performance — by Region", S["h3"]))
        header = ["Region", "Fields", "GS Trio %", "Assess %", "GS Cov %", "RC %", "Score"]
        rows = [header]
        for _, r in kpi_region.iterrows():
            rows.append([
                r["State"], str(r["Fields"]),
                f"{r['GS_Trio_Completeness']:.0f}",
                f"{r['Assessment_Pkg_Completeness']:.0f}",
                f"{r['GS_Coverage_Score']:.0f}",
                f"{r['RC_Completion_Rate']:.0f}",
                f"{r['Overall_Score']:.0f}",
            ])
        story.append(_make_pdf_table(
            rows,
            col_widths=[0.8*inch, 0.5*inch, 0.8*inch, 0.8*inch, 0.8*inch, 0.7*inch, 0.6*inch],
        ))
    story.append(PageBreak())

    # ─────────── 15. INDIVIDUAL RESPONSIBLE PAGES (v2.1 style) ───────────
    for _, kpi_row in kpi_resp.iterrows():
        resp = kpi_row["Responsible"]
        story.append(Paragraph(f"12. {resp} — Detailed Analysis", S["h2"]))

        # Summary metrics table
        score = kpi_row["Overall_Score"]
        summary_rows = [
            ["Metric", "Value", "Status"],
            ["Fields Managed", str(int(kpi_row["Fields"])), "—"],
            ["Unique Trajectories", str(int(kpi_row["UniqueTrajectories"])), "—"],
            ["Overall Score",
             f"{score:.1f}",
             _priority_label(score, [80, 60])],
            ["GS Trio Completeness",
             f"{kpi_row['GS_Trio_Completeness']:.1f}%",
             _priority_label(kpi_row["GS_Trio_Completeness"], [90, 70])],
            ["Assessment Pkg Completeness",
             f"{kpi_row['Assessment_Pkg_Completeness']:.1f}%",
             _priority_label(kpi_row["Assessment_Pkg_Completeness"], [90, 70])],
            ["GS Coverage Score",
             f"{kpi_row['GS_Coverage_Score']:.1f}%",
             _priority_label(kpi_row["GS_Coverage_Score"], [80, 50])],
            ["RC Completion Rate",
             f"{kpi_row['RC_Completion_Rate']:.1f}%",
             _priority_label(kpi_row["RC_Completion_Rate"], [80, 50])],
        ]
        story.append(_make_pdf_table(
            summary_rows,
            col_widths=[2.2 * inch, 1.5 * inch, 1.2 * inch],
        ))
        story.append(Spacer(1, 0.15 * inch))

        # Issue breakdown table
        issue_rows = [
            ["Issue Type", "Count", "Priority"],
            ["GS Regression", str(int(kpi_row.get("GS_Regression_Count", 0))),
             _priority_label(kpi_row.get("GS_Regression_Count", 0), [2, 5], "issues")],
            ["Development Gaps", str(int(kpi_row.get("Dev_Gap_Count", 0))),
             _priority_label(kpi_row.get("Dev_Gap_Count", 0), [3, 10], "issues")],
            ["GS Triplet Issues", str(int(kpi_row.get("GS_Triplet_Issues", 0))),
             _priority_label(kpi_row.get("GS_Triplet_Issues", 0), [2, 5], "issues")],
            ["Disease Monotonicity", str(int(kpi_row.get("Disease_Mono_Issues", 0))),
             _priority_label(kpi_row.get("Disease_Mono_Issues", 0), [2, 5], "issues")],
            ["RC Monotonicity", str(int(kpi_row.get("RC_Mono_Issues", 0))),
             _priority_label(kpi_row.get("RC_Mono_Issues", 0), [1, 3], "issues")],
            ["Open RC Backlog", str(int(kpi_row.get("Open_RC_Backlog", 0))),
             _priority_label(kpi_row.get("Open_RC_Backlog", 0), [3, 8], "issues")],
        ]
        story.append(_make_pdf_table(
            issue_rows,
            col_widths=[1.8 * inch, 0.8 * inch, 1.0 * inch],
            header_color=ALERT_RED,
        ))
        story.append(PageBreak())

    # ─────────── 16. RECOMMENDATIONS ───────────
    story.append(Paragraph("13. Technical Conclusions &amp; Recommendations", S["h2"]))
    # ─────────── 16. RECOMMENDATIONS ───────────
    # Global averages
    if len(kpi_resp):
        avg_trio = kpi_resp["GS_Trio_Completeness"].mean()
        avg_pkg = kpi_resp["Assessment_Pkg_Completeness"].mean()
        avg_cov = kpi_resp["GS_Coverage_Score"].mean()
        avg_rc = kpi_resp["RC_Completion_Rate"].mean()

        global_rows = [
            ["Global Metric", "Average", "Status"],
            ["GS Trio Completeness", f"{avg_trio:.1f}%",
             _priority_label(avg_trio, [90, 70])],
            ["Assessment Pkg Completeness", f"{avg_pkg:.1f}%",
             _priority_label(avg_pkg, [90, 70])],
            ["GS Coverage", f"{avg_cov:.1f}%",
             _priority_label(avg_cov, [80, 50])],
            ["RC Completion", f"{avg_rc:.1f}%",
             _priority_label(avg_rc, [80, 50])],
            ["Overall Score", f"{avg_score:.1f}",
             _priority_label(avg_score, [80, 60])],
        ]
        story.append(_make_pdf_table(
            global_rows, col_widths=[2.2 * inch, 1.2 * inch, 1 * inch],
        ))
        story.append(Spacer(1, 0.2 * inch))

    # Per-responsible actionable recommendations
    story.append(Paragraph("Actionable Recommendations per Responsible:", S["h3"]))
    has_recs = False
    for _, row in kpi_resp.iterrows():
        recs = []
        if row["GS_Trio_Completeness"] < 70:
            recs.append(
                f'<font color="{ALERT_RED}"><b>CRITICAL:</b></font> '
                f'Improve GS trio capture (currently {row["GS_Trio_Completeness"]:.0f}%)')
        elif row["GS_Trio_Completeness"] < 90:
            recs.append(f"Review incomplete GS trios — target ≥90%")
        if row["GS_Coverage_Score"] < 50:
            recs.append(
                f'<font color="{ALERT_RED}"><b>CRITICAL:</b></font> '
                "Broaden GS stage coverage — check heatmap for missing stages")
        if row.get("GS_Regression_Count", 0) > 3:
            recs.append(
                f'Review {int(row["GS_Regression_Count"])} GS regressions — likely data entry errors')
        if row.get("Dev_Gap_Count", 0) > 5:
            recs.append(
                f'Address {int(row["Dev_Gap_Count"])} development gaps — missed intermediate stages')
        if row.get("Disease_Mono_Issues", 0) > 3:
            recs.append(
                f'Check {int(row["Disease_Mono_Issues"])} disease monotonicity violations')
        if row.get("Open_RC_Backlog", 0) > 5:
            recs.append(
                f'<font color="{ALERT_RED}"><b>URGENT:</b></font> '
                f'{int(row["Open_RC_Backlog"])} open RC trajectories need closure')
        elif row.get("Open_RC_Backlog", 0) > 0:
            recs.append(
                f'{int(row["Open_RC_Backlog"])} open RC trajectories pending')
        if row["RC_Completion_Rate"] < 50 and row["RC_Total"] > 0:
            recs.append(
                f"Low RC completion ({row['RC_Completion_Rate']:.0f}%) — verify crop progress")

        if recs:
            has_recs = True
            story.append(Paragraph(f"<b>{row['Responsible']}:</b>", S["body"]))
            for rec in recs:
                story.append(Paragraph(f"&nbsp;&nbsp;&nbsp;• {rec}", S["body_small"]))
            story.append(Spacer(1, 0.08 * inch))

    if not has_recs:
        story.append(Paragraph(
            f'<font color="{SUCCESS_GREEN}"><b>✓ No critical actions required. '
            "All responsibles meeting targets.</b></font>",
            S["body"],
        ))

    story.append(Spacer(1, 0.3 * inch))
    story.append(Paragraph(
        "<i>This report is reproducible and can be regenerated with updated data. "
        "Methodology and weights are documented in the CONFIG sheet of the KPI Excel.</i>",
        S["body_small"],
    ))

    # Build
    doc.build(story)
    print(f"    PDF: {pdf_path}")
    return pdf_path


# ==================== 12. MAIN ====================

def main():
    """Main orchestrator."""
    print("=" * 80)
    print(f"  SOYBEAN PHENOLOGY QA — {PROTOCOL_FILTER}")
    print(f"  Run: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)

    # ---- Output dirs ----
    region_dir = os.path.join(BASE_OUTPUT_DIR, "By_Region")
    resp_dir = os.path.join(BASE_OUTPUT_DIR, "By_Responsible")
    plots_dir = os.path.join(BASE_OUTPUT_DIR, "_plots")
    for d in (BASE_OUTPUT_DIR, region_dir, resp_dir, plots_dir):
        os.makedirs(d, exist_ok=True)

    # ---- Load ----
    df = load_data()

    # ---- Second RC=9 map (used by several checks) ----
    second_rc9_map = _compute_second_rc9(df)

    # ---- QA checks ----
    print("\n" + "=" * 60)
    print("  RUNNING QA CHECKS")
    print("=" * 60)
    gs_cov = compute_gs_coverage(df)
    qa: Dict[str, pd.DataFrame] = {}
    qa["gs_trio_completeness"] = check_gs_trio_completeness(df)
    qa["gs_triplet_consistency"] = check_gs_triplet_consistency(df)
    qa["gs_regression"] = check_gs_regression(df)
    qa["dev_gaps"] = check_development_gaps(df)
    qa["assessment_pkg"] = check_assessment_pkg_completeness(df, second_rc9_map)
    qa["disease_mono"] = check_disease_monotonicity(df, second_rc9_map)
    rc_mono, rc_backlog = check_rc_monotonicity_and_progress(df)
    qa["rc_mono"] = rc_mono

    # ---- KPIs ----
    kpi_resp, kpi_region = compute_kpis(df, qa, gs_cov, rc_backlog)

    # ---- Summaries ----
    field_summary = build_field_summary(df, qa, rc_backlog)
    trajectory_summary = build_trajectory_summary(df)
    unknown_summary = build_unknown_summary(df)

    # ---- Plots ----
    plot_paths = generate_all_plots(
        df, qa, gs_cov, kpi_resp, kpi_region, rc_backlog, plots_dir,
    )

    # ---- Excel: by region ----
    export_excel_by_region(df, qa, gs_cov, rc_backlog, field_summary, unknown_summary, region_dir)

    # ---- Excel: by responsible ----
    export_excel_by_responsible(df, qa, gs_cov, rc_backlog, field_summary, trajectory_summary, resp_dir)

    # ---- Excel: consolidated KPI ----
    export_consolidated_kpi_excel(kpi_resp, kpi_region, unknown_summary, BASE_OUTPUT_DIR)

    # ---- PDF ----
    generate_pdf_report(
        df, qa, gs_cov, kpi_resp, kpi_region, rc_backlog,
        plot_paths, unknown_summary, BASE_OUTPUT_DIR,
    )

    # ---- Final console summary ----
    print("\n" + "=" * 80)
    print("  COMPLETE")
    print("=" * 80)
    print(f"\n  Outputs: {BASE_OUTPUT_DIR}")
    print(f"  Fields:       {df['Field'].nunique()}")
    print(f"  Trials:       {df['Trialnumber'].nunique()}")
    total_issues = sum(len(v) for v in qa.values())
    print(f"  Total QA issues flagged: {total_issues}")
    print(f"  UNASSIGNED responsibles: {(df['Responsible'] == 'UNASSIGNED').sum()} rows")
    print(f"  UNK state rows:          {(df['State'] == 'UNK').sum()}")


if __name__ == "__main__":
    main()
