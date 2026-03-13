# Group 5: Wind Turbine Examples

Isolated tests for OpenFAST wind turbine coupling. Each example tests a turbine mounted on a floating body with isobara hydro coefficients while varying the mooring configuration. No wave excitation or radiation forces are active — the focus is on verifying the FAST coupling interface.

## Examples

| Example | What it tests | Body state | Mooring | Key configuration |
|---|---|---|---|---|
| `fixed_turbine` | Turbine on locked body (no motion) | Locked (flag=1) | None | Body fixed, turbine loads only |
| `moored_turbine` | Turbine on moored floating body | Free (flag=0) | 1 line | Body free + single mooring line + turbine |

## Shared Configuration

All examples use:
- **Simulation time:** 5s
- **Integration:** ESDIRK46 (method=3), order=2, dt=0.1s
- **Water depth:** -500m, density=1025 kg/m³, gravity=9.81 m/s²
- **Body:** isobara.ehydb, RAD_DIFF, 6 DOFs, no radiation, no excitation
- **Waves:** REG H=0.0 T=10 (no wave excitation)
- **Turbine:** GE70 Beridi (FAST files in `FAST/` subdirectory)
- **FAST time steps:** forces=0.1s, controller=0.1s

## Running

Run a single example:
```bash
cd <example_name>
run.bat          # Windows
sh run.sl .      # Linux
```

Run all turbine examples (included in the master test suite):
```bash
cd examples
run_all_examples.bat    # Windows
sh run_all_examples.sh  # Linux
```

## Expected Behavior

- **fixed_turbine:** Body remains stationary; turbine aerodynamic loads are computed but body does not move. Verifies the FAST coupling initializes and runs without errors on a locked body.
- **moored_turbine:** Body responds to turbine-induced loads with mooring restoring forces. Verifies combined FAST + mooring line coupling on a floating body.

## Relevant Source Code

- `src/WindTurbine/WindTurbine.cpp` — `ReadPropertiesASCII()` for input format; FAST interface initialization
- `src/Bodies/Bodies.cpp` — Wind turbine index parsing in `ReadPropertiesASCII()`
- `src/Simulations/Simulation.cpp` — `ReadWindTurbinesASCII()` for turbine loading; FAST time stepping in `Run()`
