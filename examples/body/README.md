# Group 1: Body Behavior Examples

Isolated tests for floating body dynamics. Each example tests a single aspect of body behavior with no mooring lines, springs, wind turbines, or other subsystems.

## Examples

| Example | What it tests | Hydro DB | Key configuration |
|---|---|---|---|
| `free_decay` | Free body decay in heave | isobara.ehydb | displacement z=1m, radiation ON |
| `fixed` | Fully locked body | isobara.ehydb | flag_blocked=1, no forces |
| `partial_dofs` | Heave+pitch only | isobara.ehydb | DOFs 3 5, displacement z=1m pitch=0.1 |
| `imposed_motion` | Prescribed harmonic surge | isobara.ehydb | flag_blocked=2, A=0.5m T=5s |
| `radiation` | Radiation damping (surge) | isobara.ehydb | displacement x=1m, radiation ON |
| `excitation_linear` | 1st-order excitation (locked) | isobara.ehydb | flag_blocked=1, excitation=1, REG wave H=2 |
| `excitation_instantpos` | Excitation at instant position | isobara.ehydb | free body, excitation=2, REG wave H=2 |
| `nonlinear_hs_flat` | Non-linear hydrostatics (flat water) | one_box.hydb.h5 + STL | hydrostatics=1, displacement z=1m |
| `nonlinear_hs_waves` | Non-linear hydrostatics with waves | one_box.hydb.h5 + STL | hydrostatics=2, excitation=1, REG wave H=2 |
| `qtf` | 2nd-order QTF forces | isobara.ehydb | excitation=1, QTF=1, IRR JONSWAP Hs=1 |

## Shared Configuration

All examples use:
- **Simulation time:** 30s (exceeds 20s IRF window to validate radiation damping)
- **Integration:** ESDIRK46 (method=3), order=2, dt=0.1s
- **Water depth:** -500m, density=1025 kg/m³, gravity=9.81 m/s²
- **Body type:** RAD_DIFF (radiation/diffraction)
- **Zero wave** (REG H=0.0 T=10) unless wave excitation is being tested

## Running

Run a single example:
```bash
cd <example_name>
run.bat          # Windows
sh run.sl .      # Linux
```

Run all body examples (included in the master test suite):
```bash
cd examples
run_all_examples.bat    # Windows
sh run_all_examples.sh  # Linux
```

## Expected Behavior

- **free_decay / radiation:** Oscillation with decaying amplitude (radiation damping)
- **fixed:** Output position equals initial position at all times
- **partial_dofs:** Only heave and pitch DOFs oscillate; surge/sway/roll/yaw remain zero
- **imposed_motion:** Position follows prescribed sinusoidal trajectory
- **excitation_linear:** Excitation force time series present in output (body locked)
- **excitation_instantpos:** Body moves in response to waves with position-corrected forces
- **nonlinear_hs_flat / nonlinear_hs_waves:** Restoring forces differ from linear hydrostatics
- **qtf:** Second-order slow-drift and sum-frequency forces present in output

## Relevant Source Code

- `src/Bodies/Bodies.cpp` — `ReadPropertiesASCII()` for input format; `UpdateLockBody()` for freedom=1/2 logic
- `src/Simulations/Simulation.cpp` — `ReadBodiesASCII()` for body loading; `ReadWavesASCII()` for wave loading
