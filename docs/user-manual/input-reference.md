# Input Reference

This page provides a comprehensive field-by-field reference for all OASIS input files.

## dataProblem.dat

General simulation parameters. Each line contains a value followed by a comment.

| Line | Field | Type | Units | Description |
|------|-------|------|-------|-------------|
| 1 | `gravity` | float | m/s² | Gravitational acceleration |
| 2 | `waterDensity` | float | kg/m³ | Water density |
| 3 | `airAtmPresDensity` | float | kg/m³ | Atmospheric air density |
| 4 | `airAtmPres` | float | Pa | Atmospheric pressure |
| 5 | `airAdiabaticDilation` | float | — | Adiabatic index for air (≈1.4) |
| 6 | `waterDepth` | float | m | Water depth (positive downward from MWL) |
| 7 | `writeTimeStep` | float | s | Output writing interval |
| 8 | `maxTimeStep` | float | s | Maximum integration time step |
| 9 | `hydroTimeStep` | float | s | Hydrodynamic force update interval |
| 10 | `fastTimeStep` | float | s | OpenFAST time step |
| 11 | `fastControllerTimeStep` | float | s | OpenFAST controller time step |
| 12 | `timeIRF` | float | s | Impulse response function duration |
| 13 | `sinkingTimeStep` | float | s | Sinking HDB update interval |
| 14 | `winchesContTimeStep` | float | s | Winches controller time step |
| 15 | `owcsContTimeStep` | float | s | OWC turbine controller time step |
| 16 | `simulationTime` | float | s | Total simulation duration |
| 17 | `rotSimpFlag` | int | — | Rotation dynamics simplification (0=No, 1=Yes) |
| 18 | `timeIntMethod` | int | — | ODE solver: 1=BDF1, 2=BDFN, 3=ESDIRK46 |
| 19 | `timeIntOrder` | int | — | Solver order (for BDFN: 1–6) |
| 20 | `timeIntAdaptivity` | int | — | Adaptive time stepping (0=No, 1=Yes) |
| 21 | `timeIntJacNumStepsMax` | int | — | Max steps between Jacobian recomputation |
| 22 | `timeIntAbsTol` | float | — | Absolute tolerance |
| 23 | `timeIntRelTol` | float | — | Relative tolerance |
| 24 | `maxIterStep` | int | — | Max Newton iterations per time step |
| 25 | `readEquilibrium` | int | — | Read initial conditions from file (0=No, 1=Yes) |
| 26 | `writeEquilibrium` | int | — | Write equilibrium state to file (0=No, 1=Yes) |
| 27 | `flagStatic` | int | — | Initial mooring condition (0=Catenary, 1=Newton's) |

**Example:**

```
9.81            // gravity [m/s^2]
1025.0          // waterDensity [kg/m^3]
1.225           // airAtmPresDensity [kg/m^3]
101325.0        // airAtmPres [Pa]
1.4             // airAdiabaticDilation [-]
50.0            // waterDepth [m]
0.1             // writeTimeStep [s]
0.1             // maxTimeStep [s]
0.1             // hydroTimeStep [s]
0.1             // fastTimeStep [s]
0.1             // fastControllerTimeStep [s]
60.0            // timeIRF [s]
1.0             // sinkingTimeStep [s]
0.1             // winchesContTimeStep [s]
0.1             // owcsContTimeStep [s]
300.0           // simulationTime [s]
0               // rotSimpFlag
3               // timeIntMethod (ESDIRK46)
2               // timeIntOrder
0               // timeIntAdaptivity
0               // timeIntJacNumStepsMax
1e-5            // timeIntAbsTol
1e-3            // timeIntRelTol
20              // maxIterStep
0               // readEquilibrium
0               // writeEquilibrium
0               // flagStatic (catenary)
```

---

## dataBodies.dat

Body definitions. The number of bodies is determined automatically by counting `///` header blocks.

Each body block contains:

| # | Field | Type | Description |
|---|-------|------|-------------|
| 1 | Body type | string | Always `RAD_DIFF` |
| 2 | Take COG from HDB | int | 0=manual, 1=from database |
| 3 | DOFs list | ints | Active degrees of freedom (1–6) |
| 4 | BCP indexes | ints | Global BCP indices attached to this body |
| 5 | Wind turbine indexes | ints | Wind turbine indices (0 if none) |
| 6 | Initial position | 6 floats | \[x, y, z, rx, ry, rz\] in global frame |
| 7 | Initial displacement | 6 floats | \[dx, dy, dz, drx, dry, drz\] offsets |
| 8 | Hydro database filename | string | Path to `.hydb.h5` file |
| 9 | Database index | int | Body index within the hydro database (1-based) |
| 10 | Freedom flag | int | 0=free, 1=locked, 2=imposed motion |
| 11 | Imposed motion filename | string | File with prescribed motion (if freedom=2) |
| 12 | Hydrostatics flag | int | 0=linear, 1=nonlinear flat, 2=nonlinear with waves |
| 13 | Hydrostatics mesh filename | string | STL mesh (if nonlinear) |
| 14 | Radiation force flag | int | 0=off, 1=on |
| 15 | First-order excitation flag | int | 0=off, 1=on, 2=instantaneous position |
| 16 | Second-order excitation flag | int | 0=off, 1=slow-drift, 2=sum-frequency |
| 17 | Viscous added mass | 6 floats | Additional added mass per DOF |
| 18 | Viscous linear damping | 6 floats | Linear damping per DOF |
| 19 | Viscous quadratic damping | 6 floats | Quadratic damping per DOF |

**Example block:**

```
//////////////////
////////////////// New body [1]
//////////////////
RAD_DIFF
0
1 2 3 4 5 6
4 5 6 7
0
0.0 0.0 0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0 0.0 0.0
../../../data/one_box.hydb.h5
1
0
none
0
none
1
1
0
0.0 0.0 0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0 0.0 0.0
0.0 0.0 0.0 0.0 0.0 0.0
```

---

## dataBCPs.dat

Boundary Condition Points. The header declares the count of each BCP type.

### Header

```
1       // Number of Actuators (FairleadBCP)
2       // Number of Anchors (AnchorBCP)
0       // Number of Joints (JointBCP)
2       // Number of Fairleads (BodyBCP)
0       // Number of Elastic Anchors
```

!!! warning "Index Order"
    Global BCP indices are assigned sequentially: Actuators first (1..nAct), then Anchors, Joints, Fairleads, Elastic Anchors. All references in `dataBodies.dat` and `dataLines.dat` use these **global 1-based indices**.

### BCP Blocks

**AnchorBCP:**
```
//////////////////
////////////////// New BCP [N] [Anchor] [i]
//////////////////
200.0 0.0 -50.0     // Position [x, y, z] in global frame
0                    // Winch ID (0 = not connected)
```

**FairleadBCP (Actuator):**
```
//////////////////
////////////////// New BCP [N] [Actuator] [i]
//////////////////
5.0 0.0 -2.0        // Initial position [x, y, z]
0                    // Winch ID
dataActuator_0.dat   // Prescribed motion file
```

**JointBCP:**
```
//////////////////
////////////////// New BCP [N] [Joint] [i]
//////////////////
100.0 0.0 -30.0      // Position [x, y, z]
0                     // Winch ID
50.0                  // Mass [kg]
0.05                  // Volume [m³]
```

**BodyBCP:**
```
//////////////////
////////////////// New BCP [N] [Fairlead] [i]
//////////////////
10.0 0.0 -5.0        // Position relative to body COG
0                     // Winch ID
```

### Actuator Motion File Format

The `dataActuator_N.dat` file prescribes motion for FairleadBCP (actuator) points:

```
3001                                    // Number of time steps
0.0000 5.000 0.000 -2.000 0.000 0.000 0.000 0.000 0.000 0.000
0.0100 5.016 0.000 -2.000 0.000 0.000 0.000 0.000 0.000 0.000
...
```

Each data row has 10 columns: `t, x, y, z, vx, vy, vz, ax, ay, az`

---

## dataWaves.dat

See [Waves](waves.md) for detailed wave configuration.

### Basic Section

```
//////////////////
////////////////// Basic wave data
//////////////////
IRR             // Wave type: REG or IRR
2.0             // Height [m] (H for REG, Hs for IRR)
10.0            // Period [s] (T for REG, Tp for IRR)
0.0             // Heading [degrees]
5.0             // Ramp time [s]
```

### Irregular Wave Section (only if IRR)

```
//////////////////
////////////////// Irregular wave data
//////////////////
1               // Spectrum type: 1=JONSWAP, 2=Time-series, 3=Frequency-domain
0               // Piecewise flag: 0=single, 1=piecewise
//// JONSWAP
3.3             // gamma (peak enhancement)
0.0             // s (directional spreading, <1 = no spreading)
10.0            // dtheta [degrees]
0.001           // factor (min proportion of spectrum average)
0               // Read phases flag: 0=random, 1=from file
0.05            // Relative tolerance for phase check
0.1             // dt [s] for wave time series
none            // Phase filename
//// CUSTOM
none            // Data filename (wave.dat or wavefreq.dat)
```

---

## dataLines.dat

See [Mooring Lines](mooring-lines.md) for detailed configuration.

Each line block reads these fields in order:

| Field | Type | Description |
|-------|------|-------------|
| `lineType` | int | 1=dynamic, 2=quasi-static |
| `flag_tension` | int | 0=tension-only, 1=symmetric (tension+compression) |
| `nNodes` | int | Number of finite element nodes |
| `p` | int | Polynomial order per element |
| `L` | float | Unstretched line length [m] |
| `rho0` | float | Mass per unit length [kg/m] |
| `d` | float | Line diameter [m] |
| `flag_stiffness` | int | 0=linear, 1=viscoelastic, >1=tabulated |
| (stiffness data) | ... | Depends on `flag_stiffness` |
| `CB` | float | Seabed friction coefficient |
| `Cmn` | float | Normal added mass coefficient |
| `Cdn` | float | Normal drag coefficient |
| `Cdt` | float | Tangential drag coefficient |
| `GK` | float | Ground spring stiffness [N/m³] |
| `GC` | float | Ground damping coefficient |
| `indexSeaFloor` | int | SeaFloor index (0=flat at waterDepth) |
| `BCP_N` | int | BCP index at end node (global, 1-based) |
| `BCP_1` | int | BCP index at start node (global, 1-based) |
| `frictionModel` | int | 0=none, 1=isotropic, 2=anisotropic |
| `vth` | float | Threshold velocity for stick-slip [m/s] |
| `ust` / `usn` | float | Static friction coefficient (tangential/normal) |
| `ud` | float | Dynamic friction coefficient |
| `deltamax` | float | Max displacement for stick phase [m] |

---

## dataSprings.dat

See [Springs & Connectors](springs.md) for detailed configuration.

---

## dataWinches.dat / dataWinchesController.dat

See [Winches](winches.md) for detailed configuration.

---

## dataSeaFloor.dat

Header declares counts of each seabed type:

```
0       // Number of bathymetry (mesh) floors
1       // Number of inclined plane floors
1       // Number of flat floors
```

**Flat floor:** single depth value (overrides `waterDepth` for specific lines)

```
50.0            // Depth [m]
```

**Inclined plane:** three points defining the plane

```
0.0 0.0 -50.0       // Point 1 [x, y, z]
100.0 0.0 -45.0     // Point 2 [x, y, z]
0.0 100.0 -48.0     // Point 3 [x, y, z]
```

**Bathymetry mesh:** STL filename

```
seabed_mesh.stl      // STL mesh file
```

Floor indices are assigned sequentially: bathymetry first, then inclined, then flat.
