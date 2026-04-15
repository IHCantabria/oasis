# Group 6: OWC (Oscillating Water Column) Examples

Progressive tests for the OWC subsystem with increasing complexity. Each example builds on the previous one to isolate potential error sources. All cases use the simplified 2-body `one_box_owc.hydb.h5` (1 floater + 1 OWC water column).

## Examples

| Example | What it tests | OWC coupling | Status |
|---|---|---|---|
| `free_decay` | 2-body free decay, no OWC pneumatics | None (type=-1) | Not tested yet |
| `excitation` | Wave excitation on floater, OWC free | None (type=-1) | Not tested yet |
| `hole_coupling` | OWC pneumatic coupling via orifice | Hole (type=0) | Not tested yet |
| `turbine_coupling` | OWC with Buckingham-curve turbine | Turbine (type=1) | Not tested yet |

### Notes

- Previous iterations used the 7-body `isobara_owcs.hydb.h5` (1 floater + 6 OWCs). That setup passed free_decay and excitation but failed hole_coupling (convergence at t≈14.4s) and turbine_coupling (NaN at t≈7.5s). Switched to the simplified 2-body HDB to tune OWC parameters for numerical stability.
- A missing `rampTime` line in `dataWaves.dat` was originally causing NaN crashes in all wave-containing tests. This was fixed by adding `5.0 // Ramp time [s]` after the wave heading line.

## Shared Configuration

- **Hydro DB**: `one_box_owc.hydb.h5` (2 bodies, shared in `../data/`)
- **Solver**: BDF1 (method=1, order=2), dt_max=0.01, dt_hydro=0.05
- **Bodies**: Body 1 = floater (6 DOFs), Body 2 = OWC water column (heave only)
- **OWC params**: waterplane_area=1m², ref_air_volume=1m³, hole_area=0.1m², hole_Cd=0.1
- **Radiation**: OFF for all bodies (to isolate OWC coupling effects)
- **Viscous damping**: Floater 1e4/1e3 (lin/quad), OWC body 1e3/1e2 (heave only)
