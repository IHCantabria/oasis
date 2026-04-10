# Group 6: OWC (Oscillating Water Column) Examples

Progressive tests for the OWC subsystem with increasing complexity. Each example builds on the previous one to isolate potential error sources. All cases use the 7-body Isobara HDB (1 floater + 6 OWC water columns).

## Examples

| Example | What it tests | OWC coupling | Status |
|---|---|---|---|
| `free_decay` | Multi-body free decay, no OWC pneumatics | None (type=-1) | **PASS** — 30s, 0 convergence failures |
| `excitation` | Wave excitation on floater, OWCs free | None (type=-1) | **PASS** — 60s, 0 convergence failures |
| `hole_coupling` | OWC pneumatic coupling via orifice | Hole (type=0) | **FAIL** — convergence failure at t≈14.4s |
| `turbine_coupling` | OWC with Buckingham-curve turbine | Turbine (type=1) | **FAIL** — NaN at t≈7.5s |

### Key Findings

- **free_decay** and **excitation** pass cleanly, confirming that the multi-body setup and wave excitation work correctly without pneumatic coupling.
- **hole_coupling** and **turbine_coupling** fail when the OWC pneumatic coupling is active, indicating numerical instability in the pressure dynamics. The OWC parameters (air volume, hole area, discharge coefficient) likely need tuning for numerical stability.
- A missing `rampTime` line in `dataWaves.dat` was originally causing NaN crashes in all wave-containing tests. This was fixed by adding `5.0 // Ramp time [s]` after the wave heading line.

## Shared Configuration

- **Hydro DB**: `isobara_owcs.hydb.h5` (7 bodies, shared in `../data/`)
- **Solver**: BDF1 (method=1, order=2), dt_max=0.01, dt_hydro=0.05
- **Bodies**: Body 1 = floater (6 DOFs), Bodies 2-7 = OWC water columns (heave only)
- **OWC params**: waterplane_area=1m², ref_air_volume=1m³, hole_area=0.1m², hole_Cd=0.1
- **Radiation**: OFF for all bodies (to isolate OWC coupling effects)
- **Viscous damping**: Floater 1e4/1e3 (lin/quad), OWC bodies 1e3/1e2 (heave only)
