# OASIS Examples & Functional Tests

This directory contains a systematic test suite for OASIS, organized by functionality group. Each group isolates a specific subsystem (bodies, waves, mooring lines, etc.) so that every feature can be tested independently. The examples also serve as user documentation showing how to set up inputs for each feature.

## Strategy

The examples are organized into **functional groups**, each in its own subdirectory. Within each group, individual examples test a single feature or flag combination in isolation — no unnecessary subsystems are active. This makes it easy to:

1. **Verify** that each feature works correctly after code changes
2. **Diagnose** failures by narrowing down to the exact subsystem
3. **Learn** how to configure each feature by studying minimal working examples

All examples can be run together via `run_all_examples.bat` (Windows) or `run_all_examples.sh` (Linux). Each group also has its own `run_group.bat` / `run_group.sh` scripts to run just that group.

## Functional Groups

| # | Group | Directory | Description | Status |
|---|---|---|---|---|
| 1 | **Body Behavior** | `body/` | Floating body dynamics: free/fixed/partial DOFs, imposed motion, radiation, excitation, non-linear hydrostatics, QTF | Done |
| 2 | **Wave Types** | `waves/` | Regular, irregular (JONSWAP), time-series, frequency-domain, multi-directional, piecewise decomposition | Done |
| 3 | **Mooring Lines** | `lines/` | Dynamic mooring lines: single/multi-line, BCP types (anchor, fairlead, joint, elastic anchor, actuator), material models (linear, viscoelastic, tabulated), tension models, seabed contact, friction (isotropic, anisotropic) | Done |
| 4 | **Springs** | `springs/` | Nonlinear springs: pile connector with friction, neoprene+wire body-to-body connector | Done |
| 5 | **Wind Turbines** | `turbines/` | OpenFAST coupling on fixed and moored bodies | Done |
| 6 | **OWCs** | `owcs/` | Oscillating Water Columns: chamber dynamics, pneumatic coupling, turbine models | WIP |
| 7 | **Sinking** | `sinking/` | Progressive flooding and sinking dynamics | Done |
| 8 | **Multi-body** | `multibody/` | Multiple bodies: shared multi-body HDB with cross-coupling, separate single-body HDBs | Done |
| 9 | **Integration & Solvers** | `solvers/` | ODE solver comparison (BDF1, BDFN order 2, BDFN order 4, ESDIRK46) on a cable with circular actuator motion | Done |
| 10 | **Winches** | `winches/` | Winch-controlled mooring: constant tension pretension and horizontal plane positioning (surge/sway/yaw) with state-space controller | Done |
| 11 | **Morison** | `morison/` | Body-level wind and current drag forces via Morison model | Done |

## Directory Structure

```
examples/
├── README.md                   # This file
├── run_all_examples.bat        # Run all tests (Windows)
├── run_all_examples.sh         # Run all tests (Linux)
├── data/                       # Shared hydro databases
│   ├── one_box.hydb.h5         # Single-body hydro (most examples)
│   ├── trl_plus.hydb.h5        # Turbine platform hydro
│   ├── two_boxes.hydb.h5       # Multi-body hydro (multibody examples)
│   └── isobara_owcs.hydb.h5    # OWC platform hydro (7 bodies: 1 floater + 6 OWCs)
├── body/                       # Group 1: Body behavior (DONE)
│   ├── README.md
│   ├── free_decay/
│   ├── fixed/
│   ├── partial_dofs/
│   ├── imposed_motion/
│   ├── radiation/
│   ├── excitation_linear/
│   ├── excitation_instantpos/
│   ├── nonlinear_hs_flat/
│   ├── nonlinear_hs_waves/
│   └── qtf/
├── waves/                      # Group 2: Wave types (DONE)
│   ├── README.md
│   ├── regular/
│   ├── jonswap/
│   ├── jonswap_piecewise/
│   ├── multidirectional/
│   ├── timeseries/
│   └── frequency_domain/
├── lines/                      # Group 3: Mooring lines (DONE)
│   ├── README.md
│   ├── single_line/
│   ├── multi_line/
│   ├── joint_connection/
│   ├── elastic_anchor/
│   ├── prescribed_motion/
│   ├── viscoelastic/
│   ├── tabulated_stiffness/
│   ├── tension_symmetric/
│   ├── seabed_contact/
│   ├── friction_isotropic/
│   └── friction_anisotropic/
├── springs/                    # Group 4: Springs (DONE)
│   ├── README.md
│   ├── pile_connector/
│   └── neoprene_connector/
├── turbines/                   # Group 5: Wind turbines (DONE)
│   ├── README.md
│   ├── fixed_turbine/
│   └── moored_turbine/
├── owcs/                       # Group 6: OWCs (WIP)
│   ├── README.md
│   ├── free_decay/
│   ├── excitation/
│   ├── hole_coupling/
│   └── turbine_coupling/
├── sinking/                    # Group 7: Sinking (DONE)
│   └── progressive_flooding/
├── multibody/                  # Group 8: Multi-body (DONE)
│   ├── shared_hydb/
│   └── separate_hydb/
├── solvers/                    # Group 9: Solvers (DONE)
│   ├── README.md
│   ├── generate_actuator_motion.py
│   ├── bdf1/
│   ├── bdfn_order2/
│   ├── bdfn_order4/
│   └── esdirk46/
├── winches/                    # Group 10: Winches (DONE)
│   ├── README.md
│   └── horizontal_control/
├── morison/                    # Group 11: Morison body drag (DONE)
│   ├── run_group.bat
│   ├── run_group.sh
│   ├── wind_drag/
│   ├── current_drag/
│   ├── wind_drag_yaml/
│   └── current_drag_yaml/
├── misc/                       # Pre-existing mixed examples
│   ├── freq_wave_example/
│   ├── generic_example/
│   ├── qtf_example/
│   ├── turbine_example/
│   ├── elastic_anchor_example/
│   └── large_rotation_example/
```

## Example Structure

Each individual example follows this layout:

```
<group>/<test_name>/
├── input/
│   ├── dataProblem.dat         # Simulation parameters
│   ├── dataBodies.dat          # Body definitions
│   └── dataWaves.dat           # Wave configuration
├── output/                     # Created by simulator
├── run.bat                     # Windows runner
└── run.sl                      # Linux/SLURM runner
```

Optional input files (only present when the feature is tested):
- `dataLines.dat`, `dataBCPs.dat` — mooring lines and boundary conditions
- `dataSprings.dat` — spring connections
- `dataWindTurbines.dat` — wind turbine definitions
- `dataOWCs.dat`, `dataOWCTurbines.dat` — OWC chambers and turbines
- `dataSinking.dat` — flooding/sinking setup
- `dataSeaFloor.dat` — seabed profile
- `dataLockBody_body0.dat` — imposed motion definition

## Running

### All examples
```bash
cd examples
run_all_examples.bat        # Windows
sh run_all_examples.sh      # Linux
```

### Single group
```bash
cd examples/<group>
run_group.bat               # Windows
sh run_group.sh             # Linux
```

### Single example
```bash
cd examples/<group>/<test_name>
run.bat                     # Windows
sh run.sl .                 # Linux
```

### On a SLURM cluster
```bash
cd examples/<group>/<test_name>
sbatch run.sl input
```

## Design Decisions

- **30s simulation time** for body examples (exceeds 20s IRF window to validate radiation damping)
- **ESDIRK46 integrator** (method=3, order=2) as default — stable and accurate
- **dt = 0.1s** for all time steps (output, max, hydro, FAST, controller)
- **Missing optional files are OK** — OASIS prints a warning and sets the count to 0
- **Centralized hydro databases**: Stored in `examples/data/` and referenced via relative paths (`../../../data/<file>`) to avoid duplication. `one_box.hydb.h5` for most examples, `trl_plus.hydb.h5` for turbines, `one_box.hydb.h5` + STL for non-linear hydrostatics
- **Nested folder structure**: `examples/<group>/<test_name>/input/` for clear organization

## Adding a New Group

To continue building the test suite in a future session:

1. Read this README to understand the overall strategy and which groups are planned
2. Pick the next group to implement (e.g., "Group 2: Wave Types")
3. Study the relevant source code to understand available flags and input formats
4. Create isolated examples that each test a single feature
5. Add a `README.md` inside the group directory with a summary table
6. Update `run_all_examples.bat` and `run_all_examples.sh` with the new tests
7. Mark the group as "Done" in the table above
