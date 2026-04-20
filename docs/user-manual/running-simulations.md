# Running Simulations

## Command-Line Usage

OASIS takes a single argument: the path to the project directory containing the `input/` folder.

```bash
oasis <project_path>
```

The project path can be absolute or relative. If `.` is passed, OASIS looks for `input/` in the current directory.

**Examples:**

```bash
# From the example's input/ directory
cd examples/body/free_decay/input
../../../../bin/oasis .

# From the repository root
bin/oasis examples/body/free_decay
```

On Windows, most examples include a `run.bat` script for convenience.

## Folder Structure

Each simulation expects this layout:

```
<project>/
├── input/
│   ├── dataProblem.dat      # or dataProblem.yaml
│   ├── dataBodies.dat
│   ├── dataBCPs.dat
│   ├── dataWaves.dat
│   ├── dataLines.dat        # if mooring lines are used
│   ├── dataSprings.dat      # if springs are used
│   ├── dataWinches.dat      # if winches are used
│   └── ...
└── output/                  # created automatically
    ├── DOF_1_Body_0.txt
    ├── WaveTimeSeries.txt
    └── ...
```

- **`input/`** — all input files must be in this subdirectory
- **`output/`** — created automatically; all results are written here

## Input Format Auto-Detection

OASIS supports two input formats:

| Format | Main file | Description |
|---|---|---|
| **YAML** | `input/dataProblem.yaml` | Single file containing all configuration |
| **ASCII** | `input/dataProblem.dat` | Multiple `.dat` files, one per subsystem |

OASIS checks for `dataProblem.yaml` first. If found, it reads YAML format. Otherwise, it falls back to ASCII `.dat` files.

See [Input Format](input-format.md) for details on both systems.

## Simulation Flow

When OASIS runs, it performs these steps in order:

1. **Load** — reads all input files (problem definition, bodies, BCPs, waves, lines, springs, winches)
2. **Initialise** — sets up all subsystems, computes initial conditions (catenary or Newton's method for mooring lines)
3. **Run** — time-stepping loop with the selected ODE solver
4. **Finalise** — writes final state, closes output files

Progress is printed to the console:

```
OASIS - Offshore Advanced SImulation Software
----------------------------------------------
Reading input files...
  --> Bodies: 1
  --> BCPs: 4
  --> Lines: 3
  --> Waves: IRR (JONSWAP, Hs=2.0m, Tp=10.0s)
Initializing...
  --> Computing initial line configuration...
  --> ... done!
Running simulation (T=300.0s, dt=0.1s)...
  t = 0.0s
  t = 10.0s
  ...
  t = 300.0s
Simulation complete.
```

## Batch Execution

Run all examples at once:

=== "Windows"

    ```bat
    cd examples
    run_all_examples.bat
    ```

=== "Linux"

    ```bash
    cd examples
    chmod +x run_all_examples.sh
    ./run_all_examples.sh
    ```

Each functional group also has its own `run_group.bat` / `run_group.sh` script.

## Post-Processing

Python utilities in `resources/` can plot results:

```bash
python resources/plot_oasis_results.py <output_folder>
```

Available scripts:

| Script | Purpose |
|---|---|
| `plot_oasis_results.py` | Plot time-series output (DOFs, forces) |
| `plot_excitation_forces.py` | Visualise excitation force RAOs |
| `plot_spectrum_from_timeseries.py` | Compute and plot spectral density |
| `plot_examples.py` | Plot results from the functional test suite |
