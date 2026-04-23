"""Discover and map OASIS output files to named channels for plotting.

Usage:
    channels = discover_outputs(output_dir)
    # channels is a dict: { "display_name": (filepath, col_index) }
"""
from __future__ import annotations
import os
import re
from typing import Dict, Tuple, List, Optional

# type alias
Channel = Tuple[str, int]  # (filepath, 0-based column index)


def discover_outputs(output_dir: str) -> Dict[str, Channel]:
    """Scan an output directory and return all plottable channels."""
    channels: Dict[str, Channel] = {}
    if not os.path.isdir(output_dir):
        return channels

    files = os.listdir(output_dir)

    for fname in sorted(files):
        fpath = os.path.join(output_dir, fname)
        if not os.path.isfile(fpath):
            continue
        ext = os.path.splitext(fname)[1].lower()
        if ext not in (".txt", ".csv"):
            continue

        # ── Body DOF files: DOF_<dof>_Body_<id>.{txt,csv}
        m = re.match(r"DOF_(\d+)_Body_(\d+)\.(txt|csv)$", fname, re.IGNORECASE)
        if m:
            dof, bid = int(m.group(1)), int(m.group(2))
            dof_names = {1: "Surge", 2: "Sway", 3: "Heave",
                         4: "Roll", 5: "Pitch", 6: "Yaw"}
            label = dof_names.get(dof, f"DOF{dof}")
            prefix = f"Body {bid} / {label}"
            sep = _detect_separator(fpath)
            channels[f"{prefix} / Position"] = (fpath, 1)
            channels[f"{prefix} / Velocity"] = (fpath, 2)
            channels[f"{prefix} / Acceleration"] = (fpath, 3)
            continue

        # ── Body force files
        m = re.match(r"(WaveRadiationForce|HydroStiffnessForce|BCPForce|"
                     r"WaveExcitationForce|WindTurbineForce)_Body_(\d+)\.(txt|csv)$",
                     fname, re.IGNORECASE)
        if m:
            ftype, bid = m.group(1), int(m.group(2))
            n_cols = _count_cols(fpath)
            for ci in range(1, n_cols):
                channels[f"Body {bid} / {ftype} / col {ci}"] = (fpath, ci)
            continue

        # ── Line endpoint tensions: EndsTen_<id>.{txt,csv}
        m = re.match(r"EndsTen_(\d+)\.(txt|csv)$", fname, re.IGNORECASE)
        if m:
            lid = int(m.group(1))
            col_labels = ["TenX_Node1", "TenY_Node1", "TenZ_Node1",
                          "TenX_NodeN", "TenY_NodeN", "TenZ_NodeN"]
            for ci, lbl in enumerate(col_labels, start=1):
                channels[f"Line {lid} / Fairlead / {lbl}"] = (fpath, ci)
            channels[f"Line {lid} / Fairlead / |Tension| Node1"] = (fpath, (1, 2, 3))
            channels[f"Line {lid} / Fairlead / |Tension| NodeN"] = (fpath, (4, 5, 6))
            continue

        # ── Line node tensions: LineTen_<id>.{txt,csv}
        m = re.match(r"LineTen_(\d+)\.(txt|csv)$", fname, re.IGNORECASE)
        if m:
            lid = int(m.group(1))
            n_cols = _count_cols(fpath)
            for ci in range(1, n_cols):
                channels[f"Line {lid} / Node tension / node {ci}"] = (fpath, ci)
            continue

        # ── Line node positions: NodePos{X,Y,Z}_<id>.{txt,csv}
        m = re.match(r"NodePos([XYZ])_(\d+)\.(txt|csv)$", fname, re.IGNORECASE)
        if m:
            axis, lid = m.group(1).upper(), int(m.group(2))
            n_cols = _count_cols(fpath)
            for ci in range(1, n_cols):
                channels[f"Line {lid} / NodePos{axis} / node {ci}"] = (fpath, ci)
            continue

        # ── Wave time series
        if re.match(r"WaveTimeSeries\.(txt|csv)$", fname, re.IGNORECASE):
            channels["Waves / Elevation"] = (fpath, 1)
            continue

        # ── Wave spectrum
        if re.match(r"WaveSpectrum\.(txt|csv)$", fname, re.IGNORECASE):
            channels["Waves / Spectrum / SpectralDensity"] = (fpath, 1)
            channels["Waves / Spectrum / Amplitude"] = (fpath, 2)
            continue

        # ── OWC
        m = re.match(r"OWC_(\d+)\.(txt|csv)$", fname, re.IGNORECASE)
        if m:
            oid = int(m.group(1))
            for ci, lbl in enumerate(
                ["Displacement", "Velocity", "Pressure"], start=1
            ):
                channels[f"OWC {oid} / {lbl}"] = (fpath, ci)
            continue

        # ── Winches
        if re.match(r"WinchesTensions\.(txt|csv)$", fname, re.IGNORECASE):
            n_cols = _count_cols(fpath)
            for ci in range(1, n_cols):
                channels[f"Winches / Tension / winch {ci}"] = (fpath, ci)
            continue

        if re.match(r"WinchedLinesLengths\.(txt|csv)$", fname, re.IGNORECASE):
            n_cols = _count_cols(fpath)
            for ci in range(1, n_cols):
                channels[f"Winches / LineLength / winch {ci}"] = (fpath, ci)
            continue

        # ── Sinking
        m = re.match(r"SinkingFillingCOG_Body_(\d+)\.(txt|csv)$", fname, re.IGNORECASE)
        if m:
            bid = int(m.group(1))
            for ci, lbl in enumerate(["COG_X", "COG_Y", "COG_Z", "FillingMass"], start=1):
                channels[f"Sinking Body {bid} / {lbl}"] = (fpath, ci)
            continue

    return channels


def load_channel(filepath: str, col: int) -> Tuple[List[float], List[float]]:
    """Load (time, values) from an OASIS output text/CSV file.
    col is 1-based (column 0 is always time).
    """
    times: List[float] = []
    values: List[float] = []
    sep = _detect_separator(filepath)
    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Skip header lines (non-numeric first token)
            parts = line.split(sep) if sep else line.split()
            if len(parts) <= col:
                continue
            try:
                t = float(parts[0])
                v = float(parts[col])
                times.append(t)
                values.append(v)
            except ValueError:
                continue
    return times, values


def load_channel_magnitude(filepath: str, cols: Tuple[int, ...]) -> Tuple[List[float], List[float]]:
    """Load vector magnitude from multiple columns."""
    import math
    times: List[float] = []
    values: List[float] = []
    sep = _detect_separator(filepath)
    n_cols_needed = max(cols)
    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(sep) if sep else line.split()
            if len(parts) <= n_cols_needed:
                continue
            try:
                t = float(parts[0])
                mag = math.sqrt(sum(float(parts[c]) ** 2 for c in cols))
                times.append(t)
                values.append(mag)
            except ValueError:
                continue
    return times, values


def _detect_separator(filepath: str) -> Optional[str]:
    """Detect whether the file uses comma or space as separator."""
    try:
        with open(filepath, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    return "," if "," in line else None
    except Exception:
        pass
    return None


def _count_cols(filepath: str) -> int:
    """Count the number of columns in the first data row."""
    sep = _detect_separator(filepath)
    try:
        with open(filepath, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split(sep) if sep else line.split()
                try:
                    float(parts[0])
                    return len(parts)
                except ValueError:
                    continue
    except Exception:
        pass
    return 0
