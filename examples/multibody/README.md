# Group 8: Multi-body Examples

Isolated tests for multi-body configurations. Each example tests a different way of setting up multiple floating bodies — shared vs. separate hydrodynamic databases — while keeping other subsystems inactive. Free-decay (no waves) is used to isolate multi-body dynamics.

## Examples

| Example | What it tests | HDB setup | Key configuration |
|---|---|---|---|
| `shared_hydb` | Two bodies sharing one multi-body HDB (radiation cross-coupling) | `two_boxes.hydb.h5` (2-body HDB), indices 1 & 2 | Body 1 displaced in heave; Body 2 at rest — observe cross-coupling |
| `separate_hydb` | Two bodies with independent single-body HDBs | `one_box.hydb.h5` (index 1) + `trl_plus.hydb.h5` (index 1) | Both displaced in heave — compare different natural frequencies |

## Shared Configuration

All examples use:
- **Simulation time:** 30s
- **Integration:** ESDIRK46 (method=3), order=2, dt=0.1s
- **Water depth:** -500m, density=1025 kg/m³, gravity=9.81 m/s²
- **Bodies:** RAD_DIFF, 6 DOFs, radiation ON, excitation OFF, linear hydrostatics
- **Waves:** REG H=0.0 T=10 (no wave excitation — free-decay only)
- **Mooring/Springs/Turbines:** None (multi-body HDB loading is the focus)

## Running

Run a single example:
```bash
cd <example_name>
run.bat          # Windows
sh run.sl .      # Linux
```

Run all multi-body examples (included in the master test suite):
```bash
cd examples
run_all_examples.bat    # Windows
sh run_all_examples.sh  # Linux
```

## Expected Behavior

- **shared_hydb:** Body 1 oscillates in heave free-decay. Body 2 (initially at rest) should show motion induced by radiation cross-coupling through the shared multi-body HDB. The added mass and radiation damping matrices have off-diagonal cross-coupling terms between the two bodies.
- **separate_hydb:** Both bodies oscillate independently in heave free-decay with different natural frequencies (box geometry vs. turbine platform geometry). No cross-coupling exists since each body has its own single-body HDB.

## Relevant Source Code

- `src/Hydro/HydroDatabase.cpp` — `LoadHydroDataH5()` for HDB loading; multi-body array indexing
- `src/Simulations/Simulation.cpp` — `ReadBodiesASCII()` for body parsing; `SetupCaseConfiguration()` for multi-body HDB validation (unique indices, no gaps)
- `src/Bodies/Body.cpp` — Body dynamics with radiation force computation
