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
