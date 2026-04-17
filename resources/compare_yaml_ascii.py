"""
Compare YAML example outputs against their ASCII counterparts.
Handles CSV output (from YAML examples with output_format: csv) vs
TXT output (from ASCII examples with output_format: txt).

Mapping logic:
  - Body motion:  body_{id}_motion.csv  <->  DOF_{1-6}_Body_{id}.txt
  - Body forces:  body_{id}_forces.csv  <->  HydroStiffnessForce, WaveRadiation,
                  WaveExcitation (exc1+exc2 combined), BCPForce, WindTurbineForce
  - Line positions: line_{id}_positions.csv  <->  NodePosX/Y/Z_{id}.txt
  - Line tensions:  line_{id}_tensions.csv   <->  EndsTen_{id}.txt, LineTen_{id}.txt
  - Winch outputs:  winches_tensions.csv     <->  WinchesTensions.txt
                    winches_lengths.csv      <->  WinchedLinesLengths.txt
  - OWC:           owc_{id}.csv             <->  OWC_{id}.txt
  - Sinking:       sinking_{id}.csv         <->  SinkingFillingCOG_Body_{id}.txt
  - Common files (same name in both outputs) compared directly (.dat, LineIni, etc.)

Usage:
    python compare_yaml_ascii.py [--examples-dir PATH] [--rtol 1e-6] [--atol 1e-10] [--run]
"""

import argparse
import re
import subprocess
import sys
import numpy as np
from pathlib import Path


# YAML directory -> ASCII directory (relative to examples root)
YAML_ASCII_PAIRS = [
    ("body/fixed_yaml",              "body/fixed"),
    ("body/free_decay_yaml",         "body/free_decay"),
    ("body/radiation_yaml",          "body/radiation"),
    ("body/excitation_linear_yaml",  "body/excitation_linear"),
    ("body/imposed_motion_yaml",     "body/imposed_motion"),
    ("body/nonlinear_hs_flat_yaml",  "body/nonlinear_hs_flat"),
    ("body/partial_dofs_yaml",       "body/partial_dofs"),
    ("waves/regular_yaml",           "waves/regular"),
    ("waves/jonswap_yaml",           "waves/jonswap"),
    ("lines/single_line_yaml",       "lines/single_line"),
    ("lines/multi_line_yaml",        "lines/multi_line"),
    ("multibody/shared_hydb_yaml",   "multibody/shared_hydb"),
    ("springs/pile_connector_yaml",  "springs/pile_connector"),
    ("winches/constant_tension_yaml","winches/constant_tension"),
    ("solvers/bdf1_yaml",           "solvers/bdf1"),
]

# Direct CSV<->TXT mappings for files with identical column layout
CSV_TXT_DIRECT_MAP = {
    "winches_tensions.csv":    "WinchesTensions.txt",
    "winches_lengths.csv":     "WinchedLinesLengths.txt",
    "control_force.csv":       "ControlForce.txt",
    "reference_position.csv":  "ReferencePosition.txt",
}


# ---------------------------------------------------------------------------
# File loaders
# ---------------------------------------------------------------------------

def load_txt_file(filepath):
    """Load a whitespace-separated numeric file, skipping non-numeric lines."""
    rows = []
    with open(filepath, "r") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or stripped.startswith("ARMA_"):
                continue
            try:
                values = [float(x) for x in stripped.split()]
                rows.append(values)
            except ValueError:
                continue
    if not rows:
        return None
    max_cols = max(len(r) for r in rows)
    for r in rows:
        while len(r) < max_cols:
            r.append(float("nan"))
    return np.array(rows)


def load_csv_file(filepath):
    """Load a CSV file with header row. Returns (header_list, data_array or None)."""
    with open(filepath, "r") as f:
        header_line = f.readline().strip()
        headers = header_line.split(",") if header_line else []
        rows = []
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            try:
                values = [float(x) for x in stripped.split(",")]
                rows.append(values)
            except ValueError:
                continue
    if not rows:
        return headers, None
    return headers, np.array(rows)


# ---------------------------------------------------------------------------
# Comparison helpers
# ---------------------------------------------------------------------------

def compare_arrays(a, b, rtol, atol):
    """Compare two numpy arrays element-wise. Returns (match, message)."""
    if a.shape != b.shape:
        return False, f"shape mismatch: {a.shape} vs {b.shape}"
    nan_a, nan_b = np.isnan(a), np.isnan(b)
    if not np.array_equal(nan_a, nan_b):
        return False, "NaN positions differ"
    mask = ~nan_a
    if not mask.any():
        return True, "all NaN"
    y, r = a[mask], b[mask]
    max_abs = np.max(np.abs(y - r))
    nonzero = np.abs(r) > atol
    max_rel = (np.max(np.abs((y[nonzero] - r[nonzero]) / r[nonzero]))
               if nonzero.any() else 0.0)
    ok = np.allclose(y, r, rtol=rtol, atol=atol)
    prefix = "" if ok else "MISMATCH "
    return ok, f"{prefix}max_abs={max_abs:.2e}, max_rel={max_rel:.2e}"


def find_ids(output_dir, pattern):
    """Find numeric IDs from filenames matching *pattern* ({id} placeholder)."""
    escaped = re.escape(pattern).replace(r"\{id\}", r"(\d+)")
    regex = re.compile(escaped)
    ids = set()
    if output_dir.exists():
        for f in output_dir.iterdir():
            m = regex.fullmatch(f.name)
            if m:
                ids.add(int(m.group(1)))
    return sorted(ids)


# ---------------------------------------------------------------------------
# Per-module comparison functions
# Each returns (results_list, yaml_consumed_set, ascii_consumed_set)
# ---------------------------------------------------------------------------

def compare_body_motion(yaml_out, ascii_out, body_id, rtol, atol):
    """body_{id}_motion.csv  vs  DOF_{1..6}_Body_{id}.txt

    CSV cols: time(0), pos(1-6), vel(7-12), acc(13-18)
    Each DOF TXT: time(0), pos(1), vel(2), acc(3)
    """
    csv_path = yaml_out / f"body_{body_id}_motion.csv"
    if not csv_path.exists():
        return [], set(), set()
    _, csv_data = load_csv_file(csv_path)
    if csv_data is None:
        return [(f"body_{body_id} motion", "SKIP", "empty CSV")], {csv_path.name}, set()

    results, yaml_c, ascii_c = [], {csv_path.name}, set()
    dof_labels = ["x", "y", "z", "rl", "pt", "yw"]

    for di in range(6):
        txt_name = f"DOF_{di + 1}_Body_{body_id}.txt"
        txt_path = ascii_out / txt_name
        if not txt_path.exists():
            continue
        ascii_c.add(txt_name)
        txt_data = load_txt_file(txt_path)
        if txt_data is None:
            continue
        n = min(csv_data.shape[0], txt_data.shape[0])
        csv_cols = csv_data[:n, [1 + di, 7 + di, 13 + di]]
        txt_cols = txt_data[:n, 1:4]
        ok, msg = compare_arrays(csv_cols, txt_cols, rtol, atol)
        label = f"body_{body_id} DOF_{di + 1} ({dof_labels[di]})"
        results.append((label, "PASS" if ok else "FAIL", msg))

    return results, yaml_c, ascii_c


def compare_body_forces(yaml_out, ascii_out, body_id, rtol, atol):
    """body_{id}_forces.csv  vs  per-force-type TXT files.

    CSV column groups (6 each, but with variable per-BCP/per-WT columns):
      hs, rad, exc1, exc2, [bcp0..bcpN], bcp_total, wt_total
    ASCII WaveExcitationForce combines exc1+exc2 in one file (12 cols).
    ASCII BCPForce / WindTurbineForce have variable-width rows; the last
    6 columns are the totals that match the CSV bcp_total / wt_total.

    Uses header-based column lookup so per-BCP columns don't break the mapping.
    """
    csv_path = yaml_out / f"body_{body_id}_forces.csv"
    if not csv_path.exists():
        return [], set(), set()
    headers, csv_data = load_csv_file(csv_path)
    if csv_data is None:
        return [(f"body_{body_id} forces", "SKIP", "empty CSV")], {csv_path.name}, set()

    results, yaml_c, ascii_c = [], {csv_path.name}, set()
    bid = f"b{body_id}"

    def find_cols(prefix):
        return [i for i, h in enumerate(headers) if h.startswith(prefix)]

    def _cmp(txt_name, csv_cols, txt_extract, label):
        txt_path = ascii_out / txt_name
        if not txt_path.exists():
            return
        ascii_c.add(txt_name)
        txt_data = load_txt_file(txt_path)
        if txt_data is None:
            return
        if not csv_cols:
            return
        n = min(csv_data.shape[0], txt_data.shape[0])
        c = csv_data[:n][:, csv_cols]
        t = txt_extract(txt_data[:n])
        if c.shape[1] != t.shape[1]:
            results.append((label, "FAIL",
                            f"col mismatch: CSV {c.shape[1]} vs TXT {t.shape[1]}"))
            return
        ok, msg = compare_arrays(c, t, rtol, atol)
        results.append((label, "PASS" if ok else "FAIL", msg))

    _cmp(f"HydroStiffnessForce_Body_{body_id}.txt",
         find_cols(f"{bid}_hs_"), lambda d: d[:, 1:7],
         f"body_{body_id} hs forces")

    _cmp(f"WaveRadiationForce_Body_{body_id}.txt",
         find_cols(f"{bid}_rad_"), lambda d: d[:, 1:7],
         f"body_{body_id} rad forces")

    # Excitation: CSV has exc1 + exc2 separately, ASCII has combined (12 cols)
    _cmp(f"WaveExcitationForce_Body_{body_id}.txt",
         find_cols(f"{bid}_exc1_") + find_cols(f"{bid}_exc2_"),
         lambda d: d[:, 1:13],
         f"body_{body_id} exc forces")

    # BCP total forces: header-based lookup handles variable per-BCP columns
    _cmp(f"BCPForce_Body_{body_id}.txt",
         find_cols(f"{bid}_bcp_total_"), lambda d: d[:, -6:],
         f"body_{body_id} bcp forces")

    # Wind-turbine total forces
    _cmp(f"WindTurbineForce_Body_{body_id}.txt",
         find_cols(f"{bid}_wt_total_"), lambda d: d[:, -6:],
         f"body_{body_id} wt forces")

    return results, yaml_c, ascii_c


def compare_line_positions(yaml_out, ascii_out, line_id, rtol, atol):
    """line_{id}_positions.csv  vs  NodePosX/Y/Z_{id}.txt

    CSV has interleaved x,y,z per node: time, n0_x, n0_y, n0_z, n1_x, ...
    ASCII stores each axis in a separate file: time, val_node0, val_node1, ...
    """
    csv_path = yaml_out / f"line_{line_id}_positions.csv"
    if not csv_path.exists():
        return [], set(), set()
    _, csv_data = load_csv_file(csv_path)
    if csv_data is None:
        return [(f"line_{line_id} positions", "SKIP", "empty CSV")], {csv_path.name}, set()

    results, yaml_c, ascii_c = [], {csv_path.name}, set()

    for axis_idx, axis in enumerate(["X", "Y", "Z"]):
        txt_name = f"NodePos{axis}_{line_id}.txt"
        txt_path = ascii_out / txt_name
        if not txt_path.exists():
            continue
        ascii_c.add(txt_name)
        txt_data = load_txt_file(txt_path)
        if txt_data is None:
            continue
        n = min(csv_data.shape[0], txt_data.shape[0])
        num_nodes = txt_data.shape[1] - 1  # exclude time column
        # CSV axis columns: every 3rd starting from (1 + axis_idx)
        csv_axis = csv_data[:n, 1 + axis_idx::3][:, :num_nodes]
        txt_axis = txt_data[:n, 1:num_nodes + 1]
        ok, msg = compare_arrays(csv_axis, txt_axis, rtol, atol)
        results.append((f"line_{line_id} pos_{axis.lower()}", "PASS" if ok else "FAIL", msg))

    return results, yaml_c, ascii_c


def compare_line_tensions(yaml_out, ascii_out, line_id, rtol, atol):
    """line_{id}_tensions.csv  vs  EndsTen_{id}.txt + LineTen_{id}.txt

    CSV: time, start_x, start_y, start_z, end_x, end_y, end_z, n0_ten, ..., nN_ten
    EndsTen: time, start_x, start_y, start_z, end_x, end_y, end_z
    LineTen: row 0 = arc-length coords; subsequent rows = time, T_node0, ..., T_nodeN

    Uses header-based column lookup to find start/end/node tension columns.
    """
    csv_path = yaml_out / f"line_{line_id}_tensions.csv"
    if not csv_path.exists():
        return [], set(), set()
    headers, csv_data = load_csv_file(csv_path)
    if csv_data is None:
        return [(f"line_{line_id} tensions", "SKIP", "empty CSV")], {csv_path.name}, set()

    results, yaml_c, ascii_c = [], {csv_path.name}, set()
    lid = f"l{line_id}"

    # --- EndsTen: start + end 3D tensions ---
    txt_name = f"EndsTen_{line_id}.txt"
    txt_path = ascii_out / txt_name
    if txt_path.exists():
        ascii_c.add(txt_name)
        txt_data = load_txt_file(txt_path)
        if txt_data is not None:
            n = min(csv_data.shape[0], txt_data.shape[0])
            start_cols = [i for i, h in enumerate(headers)
                          if h.startswith(f"{lid}_ten_start_")]
            end_cols = [i for i, h in enumerate(headers)
                        if h.startswith(f"{lid}_ten_end_")]
            csv_ends = csv_data[:n][:, start_cols + end_cols]
            txt_ends = txt_data[:n, 1:7]
            ok, msg = compare_arrays(csv_ends, txt_ends, rtol, atol)
            results.append((f"line_{line_id} end tensions",
                            "PASS" if ok else "FAIL", msg))

    # --- LineTen: skip first arc-length row, compare per-node tensions ---
    txt_name = f"LineTen_{line_id}.txt"
    txt_path = ascii_out / txt_name
    if txt_path.exists():
        ascii_c.add(txt_name)
        txt_data = load_txt_file(txt_path)
        if txt_data is not None and txt_data.shape[0] > 1:
            txt_ten = txt_data[1:]  # skip arc-length row
            n = min(csv_data.shape[0], txt_ten.shape[0])
            num_nodes = txt_ten.shape[1] - 1
            node_cols = [i for i, h in enumerate(headers)
                         if re.match(rf"{lid}_n\d+_ten$", h)]
            if node_cols:
                csv_ten = csv_data[:n][:, node_cols][:, :num_nodes]
                txt_cols = txt_ten[:n, 1:num_nodes + 1]
                ok, msg = compare_arrays(csv_ten, txt_cols, rtol, atol)
                results.append((f"line_{line_id} node tensions",
                                "PASS" if ok else "FAIL", msg))

    return results, yaml_c, ascii_c


def compare_csv_vs_txt_direct(yaml_out, ascii_out, csv_name, txt_name,
                               label, rtol, atol):
    """Compare a CSV and a TXT file that share the same column layout."""
    csv_path, txt_path = yaml_out / csv_name, ascii_out / txt_name
    yc = {csv_name} if csv_path.exists() else set()
    ac = {txt_name} if txt_path.exists() else set()
    if not csv_path.exists() or not txt_path.exists():
        return [], yc, ac
    _, csv_data = load_csv_file(csv_path)
    txt_data = load_txt_file(txt_path)
    if csv_data is None or txt_data is None:
        return [(label, "SKIP", "no data")], yc, ac
    n = min(csv_data.shape[0], txt_data.shape[0])
    c = min(csv_data.shape[1], txt_data.shape[1])
    ok, msg = compare_arrays(csv_data[:n, :c], txt_data[:n, :c], rtol, atol)
    return [(label, "PASS" if ok else "FAIL", msg)], yc, ac


# ---------------------------------------------------------------------------
# Main pair comparison
# ---------------------------------------------------------------------------

def compare_pair(examples_dir, yaml_rel, ascii_rel, rtol, atol):
    """Compare all output files between a YAML (CSV) and ASCII (TXT) pair."""
    yaml_out = examples_dir / yaml_rel / "output"
    ascii_out = examples_dir / ascii_rel / "output"

    if not yaml_out.exists():
        return [(yaml_rel, "NO OUTPUT", "YAML output/ directory missing")]
    if not ascii_out.exists():
        return [(yaml_rel, "NO OUTPUT", "ASCII output/ directory missing")]

    yaml_files = {f.name for f in yaml_out.iterdir() if f.is_file()}
    ascii_files = {f.name for f in ascii_out.iterdir() if f.is_file()}
    yaml_acc, ascii_acc = set(), set()
    results = []

    # ---- 1. Common files (same name in both: .dat, LineIni, wave spectra) --
    for fname in sorted(yaml_files & ascii_files):
        yaml_acc.add(fname)
        ascii_acc.add(fname)
        if fname.endswith(".csv"):
            _, y = load_csv_file(yaml_out / fname)
        else:
            y = load_txt_file(yaml_out / fname)
        a = load_txt_file(ascii_out / fname)
        if y is None and a is None:
            results.append((fname, "PASS", "both non-numeric"))
        elif y is None or a is None:
            results.append((fname, "FAIL", "one has no numeric data"))
        else:
            ok, msg = compare_arrays(y, a, rtol, atol)
            results.append((fname, "PASS" if ok else "FAIL", msg))

    # ---- 2. Body motion & forces -----------------------------------------
    for bid in find_ids(yaml_out, "body_{id}_motion.csv"):
        for fn in (compare_body_motion, compare_body_forces):
            r, yc, ac = fn(yaml_out, ascii_out, bid, rtol, atol)
            results.extend(r)
            yaml_acc.update(yc)
            ascii_acc.update(ac)

    # ---- 3. Line positions & tensions -------------------------------------
    for lid in find_ids(yaml_out, "line_{id}_positions.csv"):
        for fn in (compare_line_positions, compare_line_tensions):
            r, yc, ac = fn(yaml_out, ascii_out, lid, rtol, atol)
            results.extend(r)
            yaml_acc.update(yc)
            ascii_acc.update(ac)
        # Account for debug file (created but never written)
        for debug in (f"LineDebug_{lid}.txt",):
            if debug in ascii_files:
                ascii_acc.add(debug)

    # ---- 4. Direct CSV<->TXT pairs (winches, controller) ------------------
    for csv_name, txt_name in CSV_TXT_DIRECT_MAP.items():
        r, yc, ac = compare_csv_vs_txt_direct(
            yaml_out, ascii_out, csv_name, txt_name,
            csv_name.replace(".csv", ""), rtol, atol)
        results.extend(r)
        yaml_acc.update(yc)
        ascii_acc.update(ac)

    # ---- 5. OWC -----------------------------------------------------------
    for oid in find_ids(yaml_out, "owc_{id}.csv"):
        r, yc, ac = compare_csv_vs_txt_direct(
            yaml_out, ascii_out, f"owc_{oid}.csv", f"OWC_{oid}.txt",
            f"owc_{oid}", rtol, atol)
        results.extend(r)
        yaml_acc.update(yc)
        ascii_acc.update(ac)

    # ---- 6. Sinking -------------------------------------------------------
    for sid in find_ids(yaml_out, "sinking_{id}.csv"):
        r, yc, ac = compare_csv_vs_txt_direct(
            yaml_out, ascii_out, f"sinking_{sid}.csv",
            f"SinkingFillingCOG_Body_{sid}.txt",
            f"sinking_{sid}", rtol, atol)
        results.extend(r)
        yaml_acc.update(yc)
        ascii_acc.update(ac)

    # ---- 7. Spring CSVs (YAML-only, new feature) -------------------------
    for sid in find_ids(yaml_out, "spring_{id}.csv"):
        yaml_acc.add(f"spring_{sid}.csv")

    # ---- Report unaccounted files -----------------------------------------
    ye = sorted(yaml_files - yaml_acc)
    ae = sorted(ascii_files - ascii_acc)
    if ye:
        results.append((yaml_rel, "INFO", f"YAML-only: {', '.join(ye)}"))
    if ae:
        results.append((yaml_rel, "INFO", f"ASCII-only: {', '.join(ae)}"))

    return results


# ---------------------------------------------------------------------------
# Running examples
# ---------------------------------------------------------------------------

def run_example(examples_dir, example_rel, oasis_exe):
    """Run an OASIS example. Returns exit code."""
    example_dir = examples_dir / example_rel
    result = subprocess.run(
        [str(oasis_exe), "."],
        cwd=str(example_dir),
        capture_output=True,
        timeout=300,
    )
    return result.returncode


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare YAML (CSV) vs ASCII (TXT) OASIS example outputs")
    parser.add_argument("--examples-dir", type=str, default=None,
                        help="Path to examples/ directory (default: auto-detect)")
    parser.add_argument("--rtol", type=float, default=1e-6,
                        help="Relative tolerance (default: 1e-6)")
    parser.add_argument("--atol", type=float, default=1e-10,
                        help="Absolute tolerance (default: 1e-10)")
    parser.add_argument("--run", action="store_true",
                        help="Run all examples before comparing (requires OASIS.exe)")
    parser.add_argument("--oasis-exe", type=str, default=None,
                        help="Path to OASIS.exe (default: auto-detect from bin/)")
    args = parser.parse_args()

    # Resolve paths
    if args.examples_dir:
        examples_dir = Path(args.examples_dir)
    else:
        script_dir = Path(__file__).resolve().parent
        examples_dir = script_dir.parent / "examples"

    if not examples_dir.exists():
        print(f"ERROR: Examples directory not found: {examples_dir}")
        sys.exit(1)

    # Optionally run all examples first
    if args.run:
        oasis_exe = (Path(args.oasis_exe) if args.oasis_exe
                     else examples_dir.parent / "bin" / "OASIS.exe")
        if not oasis_exe.exists():
            print(f"ERROR: OASIS.exe not found: {oasis_exe}")
            sys.exit(1)
        print(f"Running all examples with {oasis_exe}...")
        for yaml_rel, ascii_rel in YAML_ASCII_PAIRS:
            for rel in [ascii_rel, yaml_rel]:
                print(f"  Running {rel}...", end=" ", flush=True)
                try:
                    rc = run_example(examples_dir, rel, oasis_exe)
                    print(f"exit={rc}")
                except subprocess.TimeoutExpired:
                    print("TIMEOUT")
                except Exception as e:
                    print(f"ERROR: {e}")

    # Compare outputs
    print(f"\n{'=' * 80}")
    print(f"Comparing YAML (CSV) vs ASCII (TXT) outputs  "
          f"(rtol={args.rtol}, atol={args.atol})")
    print(f"{'=' * 80}\n")

    total_compared = 0
    total_passed = 0
    failed_pairs = []

    for yaml_rel, ascii_rel in YAML_ASCII_PAIRS:
        results = compare_pair(examples_dir, yaml_rel, ascii_rel,
                               args.rtol, args.atol)
        pair_ok = True

        for name, status, msg in results:
            if status == "PASS":
                total_passed += 1
                total_compared += 1
            elif status == "FAIL":
                total_compared += 1
                pair_ok = False
                print(f"  FAIL  {name}: {msg}")
            elif status == "INFO":
                print(f"  INFO  {name}: {msg}")
            else:
                # NO OUTPUT, SKIP, etc.
                pair_ok = False
                print(f"  {status:5s} {name}: {msg}")

        if pair_ok:
            n = sum(1 for _, s, _ in results if s == "PASS")
            print(f"  PASS  {yaml_rel} ({n} comparisons match)")
        else:
            failed_pairs.append(yaml_rel)

    # Summary
    print(f"\n{'=' * 80}")
    print(f"SUMMARY: {total_passed}/{total_compared} comparisons match")
    n_ok = len(YAML_ASCII_PAIRS) - len(failed_pairs)
    print(f"  Pairs: {n_ok}/{len(YAML_ASCII_PAIRS)} fully matching")
    if failed_pairs:
        print(f"  Failed pairs: {', '.join(failed_pairs)}")
    print(f"{'=' * 80}")

    sys.exit(1 if failed_pairs else 0)


if __name__ == "__main__":
    main()
"""
Compare YAML example outputs against their ASCII counterparts.
Verifies that YAML parsing produces identical numerical results.

Usage:
    python compare_yaml_ascii.py [--examples-dir PATH] [--rtol 1e-6] [--atol 1e-10] [--run]
"""

import argparse
import os
import subprocess
import sys
import numpy as np
from pathlib import Path


# YAML directory -> ASCII directory (relative to examples root)
YAML_ASCII_PAIRS = [
    ("body/fixed_yaml",              "body/fixed"),
    ("body/free_decay_yaml",         "body/free_decay"),
    ("body/radiation_yaml",          "body/radiation"),
    ("body/excitation_linear_yaml",  "body/excitation_linear"),
    ("body/imposed_motion_yaml",     "body/imposed_motion"),
    ("body/nonlinear_hs_flat_yaml",  "body/nonlinear_hs_flat"),
    ("body/partial_dofs_yaml",       "body/partial_dofs"),
    ("waves/regular_yaml",           "waves/regular"),
    ("waves/jonswap_yaml",           "waves/jonswap"),
    ("lines/single_line_yaml",       "lines/single_line"),
    ("lines/multi_line_yaml",        "lines/multi_line"),
    ("multibody/shared_hydb_yaml",   "multibody/shared_hydb"),
    ("springs/pile_connector_yaml",  "springs/pile_connector"),
    ("winches/constant_tension_yaml","winches/constant_tension"),
    ("solvers/bdf1_yaml",           "solvers/bdf1"),
]


def load_numeric_file(filepath):
    """Load a whitespace-separated numeric file, skipping comment/header lines."""
    rows = []
    with open(filepath, "r") as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            # Skip comment lines
            if stripped.startswith("#"):
                continue
            # Skip Armadillo headers
            if stripped.startswith("ARMA_"):
                continue
            # Try parsing as numbers
            try:
                values = [float(x) for x in stripped.split()]
                rows.append(values)
            except ValueError:
                continue  # Skip non-numeric lines (headers, labels)
    if not rows:
        return None
    # Handle ragged rows by padding with NaN
    max_cols = max(len(r) for r in rows)
    for r in rows:
        while len(r) < max_cols:
            r.append(float("nan"))
    return np.array(rows)


def compare_file(yaml_path, ascii_path, rtol, atol):
    """Compare two output files numerically. Returns (match, message)."""
    yaml_data = load_numeric_file(yaml_path)
    ascii_data = load_numeric_file(ascii_path)

    if yaml_data is None and ascii_data is None:
        return True, "both empty/non-numeric"
    if yaml_data is None or ascii_data is None:
        return False, "one file has no numeric data"
    if yaml_data.shape != ascii_data.shape:
        return False, f"shape mismatch: YAML {yaml_data.shape} vs ASCII {ascii_data.shape}"

    # Mask NaN positions (both should have NaN in same places)
    yaml_nan = np.isnan(yaml_data)
    ascii_nan = np.isnan(ascii_data)
    if not np.array_equal(yaml_nan, ascii_nan):
        return False, "NaN positions differ"

    mask = ~yaml_nan
    if not mask.any():
        return True, "all NaN"

    y = yaml_data[mask]
    a = ascii_data[mask]

    if np.allclose(y, a, rtol=rtol, atol=atol):
        max_rel = 0.0
        nonzero = np.abs(a) > atol
        if nonzero.any():
            max_rel = np.max(np.abs((y[nonzero] - a[nonzero]) / a[nonzero]))
        max_abs = np.max(np.abs(y - a))
        return True, f"max_abs={max_abs:.2e}, max_rel={max_rel:.2e}"
    else:
        diff = np.abs(y - a)
        max_abs = np.max(diff)
        idx_flat = np.argmax(diff)
        nonzero = np.abs(a) > atol
        max_rel = np.max(np.abs((y[nonzero] - a[nonzero]) / a[nonzero])) if nonzero.any() else 0.0
        return False, f"MISMATCH max_abs={max_abs:.2e}, max_rel={max_rel:.2e}"


def run_example(examples_dir, example_rel, oasis_exe):
    """Run an OASIS example. Returns exit code."""
    example_dir = examples_dir / example_rel
    result = subprocess.run(
        [str(oasis_exe), "."],
        cwd=str(example_dir),
        capture_output=True,
        timeout=300,
    )
    return result.returncode


def compare_pair(examples_dir, yaml_rel, ascii_rel, rtol, atol):
    """Compare all output files between a YAML and ASCII example pair."""
    yaml_out = examples_dir / yaml_rel / "output"
    ascii_out = examples_dir / ascii_rel / "output"

    if not yaml_out.exists():
        return [(yaml_rel, "NO OUTPUT", "YAML output/ directory missing")]
    if not ascii_out.exists():
        return [(yaml_rel, "NO OUTPUT", "ASCII output/ directory missing")]

    # Collect comparable files (skip non-numeric formats we can't compare)
    ascii_files = {f.name for f in ascii_out.iterdir() if f.is_file()}
    yaml_files = {f.name for f in yaml_out.iterdir() if f.is_file()}

    common = sorted(ascii_files & yaml_files)
    yaml_only = sorted(yaml_files - ascii_files)
    ascii_only = sorted(ascii_files - yaml_files)

    results = []

    if yaml_only:
        results.append((yaml_rel, "EXTRA", f"YAML-only files: {', '.join(yaml_only)}"))
    if ascii_only:
        results.append((yaml_rel, "MISSING", f"ASCII-only files: {', '.join(ascii_only)}"))

    for fname in common:
        match, msg = compare_file(yaml_out / fname, ascii_out / fname, rtol, atol)
        status = "PASS" if match else "FAIL"
        results.append((f"{yaml_rel}/{fname}", status, msg))

    return results


def main():
    parser = argparse.ArgumentParser(description="Compare YAML vs ASCII OASIS example outputs")
    parser.add_argument("--examples-dir", type=str, default=None,
                        help="Path to examples/ directory (default: auto-detect)")
    parser.add_argument("--rtol", type=float, default=1e-6, help="Relative tolerance (default: 1e-6)")
    parser.add_argument("--atol", type=float, default=1e-10, help="Absolute tolerance (default: 1e-10)")
    parser.add_argument("--run", action="store_true",
                        help="Run all examples before comparing (requires oasis.exe)")
    parser.add_argument("--oasis-exe", type=str, default=None,
                        help="Path to oasis.exe (default: auto-detect from bin/)")
    args = parser.parse_args()

    # Resolve paths
    if args.examples_dir:
        examples_dir = Path(args.examples_dir)
    else:
        script_dir = Path(__file__).resolve().parent
        examples_dir = script_dir.parent / "examples"

    if not examples_dir.exists():
        print(f"ERROR: Examples directory not found: {examples_dir}")
        sys.exit(1)

    # Optionally run all examples first
    if args.run:
        oasis_exe = Path(args.oasis_exe) if args.oasis_exe else examples_dir.parent / "bin" / "oasis.exe"
        if not oasis_exe.exists():
            print(f"ERROR: oasis.exe not found: {oasis_exe}")
            sys.exit(1)
        print(f"Running all examples with {oasis_exe}...")
        for yaml_rel, ascii_rel in YAML_ASCII_PAIRS:
            for rel in [ascii_rel, yaml_rel]:
                print(f"  Running {rel}...", end=" ", flush=True)
                try:
                    rc = run_example(examples_dir, rel, oasis_exe)
                    print(f"exit={rc}")
                except subprocess.TimeoutExpired:
                    print("TIMEOUT")
                except Exception as e:
                    print(f"ERROR: {e}")

    # Compare outputs
    print(f"\n{'='*80}")
    print(f"Comparing YAML vs ASCII outputs (rtol={args.rtol}, atol={args.atol})")
    print(f"{'='*80}\n")

    total_files = 0
    passed_files = 0
    failed_pairs = []

    for yaml_rel, ascii_rel in YAML_ASCII_PAIRS:
        results = compare_pair(examples_dir, yaml_rel, ascii_rel, args.rtol, args.atol)
        pair_ok = True

        for name, status, msg in results:
            if status == "PASS":
                passed_files += 1
                total_files += 1
            elif status == "FAIL":
                total_files += 1
                pair_ok = False
                print(f"  FAIL  {name}: {msg}")
            else:
                # EXTRA, MISSING, NO OUTPUT
                pair_ok = False
                print(f"  {status:5s} {name}: {msg}")

        if pair_ok:
            n = sum(1 for _, s, _ in results if s == "PASS")
            print(f"  PASS  {yaml_rel} ({n} files match)")
        else:
            failed_pairs.append(yaml_rel)

    # Summary
    print(f"\n{'='*80}")
    print(f"SUMMARY: {passed_files}/{total_files} files match")
    print(f"  Pairs: {len(YAML_ASCII_PAIRS) - len(failed_pairs)}/{len(YAML_ASCII_PAIRS)} fully matching")
    if failed_pairs:
        print(f"  Failed pairs: {', '.join(failed_pairs)}")
    print(f"{'='*80}")

    sys.exit(1 if failed_pairs else 0)


if __name__ == "__main__":
    main()
