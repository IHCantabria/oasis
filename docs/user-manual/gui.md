# OASIS GUI

OASIS ships with an optional PyQt5-based graphical interface (`oasis_gui/`) that lets you build
input files, run simulations, and inspect results — without editing YAML by hand.

!!! note "Standalone executable — coming soon"
    A pre-built Windows `.exe` (no Python installation required) is planned.
    Compilation instructions will be added once available.
    In the meantime, follow the steps below to run from source.

---

## Requirements

| Software | Version |
|---|---|
| Python | ≥ 3.9 |
| OASIS executable | must be compiled — see [Installation](../installation/windows.md) |

Python package dependencies are listed in `oasis_gui/requirements.txt`:

| Package | Version | Purpose |
|---|---|---|
| PyQt5 | ≥ 5.15 | GUI framework |
| matplotlib | ≥ 3.5 | Schematic and results plots |
| ruamel.yaml | ≥ 0.17 | YAML input file I/O |
| numpy | ≥ 1.21 | Numerics |
| scipy | ≥ 1.7 | Signal processing |
| pyvista / pyvistaqt | ≥ 0.42 / 0.11 | 3-D mesh viewer (optional) |

---

## Setup

### 1 — Create a virtual environment

Open a terminal inside the repository root and run:

=== "Windows"

    ```bat
    cd oasis_gui
    python -m venv venv
    venv\Scripts\activate
    ```

=== "Linux / macOS"

    ```bash
    cd oasis_gui
    python3 -m venv venv
    source venv/bin/activate
    ```

### 2 — Install dependencies

With the virtual environment activated:

```bash
pip install -r requirements.txt
```

### 3 — Launch the interface

```bash
python main.py
```

The GUI window will open. On startup it automatically searches for the OASIS executable
at the following locations (relative to the repository root):

- `bin/OASIS.exe` (Windows)
- `bin/oasis` (Linux)
- `build/Release/OASIS.exe`
- `build/x64/Release/OASIS.exe`

If the executable is not found automatically, go to the **Run** tab and click **Browse…**
next to the *OASIS executable* field to locate it manually.

---

## Interface overview

The GUI is divided into tabs accessible from the left panel:

| Module | Description |
|---|---|
| **Simulation** | Global settings: time step, duration, output format, output directory |
| **Bodies** | Add and configure floating bodies (mass, inertia, hydrodynamic database, active DOFs, initial position) |
| **Lines** | Dynamic (FEM) and quasi-static mooring lines |
| **Waves** | Wave input (regular, irregular, time-series) |
| **Springs** | Linear and non-linear connectors between bodies |
| **Winches** | Winch-controlled mooring lines |
| **OWC** | Oscillating Water Column modules |
| **Wind Turbines** | OpenFAST coupling |
| **BCPs** | Boundary Coupling Points — anchors, joints, body-local attachments |
| **Run** | Set the OASIS executable path, launch and monitor a simulation |
| **Results** | Scan an output directory, plot channels, and compare bodies |

### Schematic view

The central panel shows a 2-D top-down schematic of the configured system.
Click **Refresh** to update it after editing any module.

### Results view

1. Set (or browse to) the output directory produced by a simulation.
2. Click **Scan** to discover all output files and channels.
3. Use **Plot all bodies** for a quick overview: one tab per body, 6-DOF position subplots.
4. Use **Plot custom** to choose arbitrary X and Y channels from the tree.

---

## Saving and loading cases

- **File → New** — start a blank case.
- **File → Open** — load an existing YAML input file (`dataProblem.yaml`).
- **File → Save / Save As** — write the current configuration to a YAML file ready to pass to OASIS.

---

## Troubleshooting

| Problem | Solution |
|---|---|
| `ModuleNotFoundError: PyQt5` | The virtual environment is not activated, or `pip install -r requirements.txt` was not run. |
| GUI opens but OASIS executable field is empty | Compile OASIS first (see [Installation](../installation/windows.md)), then point to `bin/OASIS.exe` in the **Run** tab. |
| Scan finds no channels | Check that the output directory contains `.txt` or `.csv` files produced by OASIS. Make sure `output_format` in the simulation settings matches what was used when running. |
| Plot shows "no data" for a DOF | That DOF was not enabled in the simulation (`dofs` list in Bodies). |
