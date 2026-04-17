# Large Rotation Example

Comparison test for simplified vs full rigid body rotation dynamics.

## Description

This example runs **two identical free-decay simulations** that differ only in
the rotation dynamics flag (`rotSimpFlag`):

| Case | Directory | `rotSimpFlag` | Description |
|------|-----------|---------------|-------------|
| Simplified | `case_simplified/` | 1 | Small-angle approximation: $\boldsymbol{\omega} \approx \dot{\mathbf{q}}_{\text{rot}}$ |
| Full | `case_full/` | 0 | Exact mapping: $\boldsymbol{\omega} = A(\phi, \theta, \psi) \cdot \dot{\mathbf{q}}_{\text{rot}}$ |

Both cases use:
- **1 free body** with 6 DOFs, linear hydrostatics
- **No waves** (H=0) — pure free decay
- **Radiation forces enabled** — provides damping for decay
- **Initial displacement**: heave = 1 m, pitch = 0.1 rad
- **No mooring lines**, no BCPs, no wind turbines

For small rotations (like this test case), both approaches should produce
nearly identical results. Differences grow with rotation amplitude.

## Structure

```
large_rotation_example/
├── case_simplified/       # rotSimpFlag = 1
│   ├── input/
│   │   ├── dataProblem.dat
│   │   ├── dataBodies.dat
│   │   ├── dataWaves.dat
│   │   ├── dataBCPs.dat
│   │   ├── dataLines.dat
│   │   └── one_box.hydb.h5
│   └── output/
├── case_full/             # rotSimpFlag = 0
│   ├── input/             # (same files, only dataProblem.dat differs)
│   └── output/
├── plot_results.py        # Overlay comparison plot
├── run.bat                # Runs both cases (Windows)
├── run.sl                 # Runs both cases (Linux/SLURM)
└── README.md
```

## Running

### Windows
```bash
cd large_rotation_example
run.bat
python plot_results.py
```

### Linux (SLURM)
```bash
cd large_rotation_example
sbatch run.sl
python plot_results.py
```

## What to Check

1. Both cases complete without NaN/Inf errors
2. Heave and pitch show damped oscillation (free decay)
3. Simplified and full results overlap closely (small initial pitch = 0.1 rad)
4. No divergence or numerical instability in the full rotation case
