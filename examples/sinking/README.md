# Group 7: Sinking Examples

Isolated tests for the progressive flooding and sinking dynamics subsystem. Each example tests sinking behavior with no mooring lines, springs, wind turbines, or other subsystems.

## Examples

| Example | What it tests | Sinking HDBs | Key configuration |
|---|---|---|---|
| `progressive_flooding` | Two-compartment symmetric flooding with HDB interpolation | one_box_5p0, one_box_7p5, one_box_10p0 | 2 groups, linear fill ramp 50–400s, 3 HDBs at increasing draft |

## Shared Configuration

All examples use:
- **Simulation time:** 600s (allows full flooding ramp from t=50s to t=400s plus 200s steady-state)
- **Integration:** ESDIRK46 (method=3), order=2, dt=0.1s
- **Water depth:** -500m, density=1025 kg/m³, gravity=9.81 m/s²
- **Body type:** RAD_DIFF (radiation/diffraction), base HDB: one_box.hydb.h5
- **Zero wave** (REG H=0.0 T=10) — calm water to isolate sinking dynamics
- **Sinking HDB update interval:** 15s
- **Body freedom:** free (flag=0), radiation ON, no excitation

## Running

Run a single example:
```bash
cd progressive_flooding
run.bat          # Windows
sh run.sl .      # Linux
```

Run all sinking examples:
```bash
cd examples/sinking
run_group.bat       # Windows
sh run_group.sh     # Linux
```

## Expected Behavior

- **progressive_flooding:** Body remains stationary until t=50s, then gradually sinks as both compartments fill symmetrically (no heel/trim). At t=400s flooding is complete (total added mass = 500000 kg). Body settles at a deeper equilibrium draft. Hydrodynamic forces are interpolated between the three depth-variant HDBs as the total filling mass increases.

## Sinking Subsystem Overview

The sinking subsystem models progressive flooding of a floating body:

1. **Multiple HDBs:** Hydrodynamic databases at different drafts (load states) are interpolated based on current total filling mass
2. **Compartment groups:** Each group has a polygon footprint, a floor Z-coordinate, and a time-dependent filling schedule
3. **Mass & inertia update:** As compartments fill, the structural mass matrix is updated (added mass + shifted center of gravity)
4. **HDB interpolation:** The two adjacent HDBs (by load state) are selected and linearly interpolated for all hydrodynamic quantities

## Relevant Source Code

- `src/Sinking/Sinking.cpp` — `ReadPropertiesASCII()` for input format; `UpdateGroupsFillingState()` for mass interpolation; `UpdateSinkingHydrodynamics()` / `UpdateSinkingHydrostatics()` for force updates
- `src/Simulations/Simulation.cpp` — `ReadSinkingASCII()` for dataSinking.dat loading
