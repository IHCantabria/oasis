# Output Files

OASIS writes results to the `output/` directory within each simulation's project folder. Files are plain text, space-delimited, and generally have no headers unless noted.

## Bodies

**6 DOF files** per body (one per degree of freedom):

| File | Columns |
|---|---|
| `DOF_1_Body_<ID>.txt` | Time, Position_X, Velocity_X, Acceleration_X |
| `DOF_2_Body_<ID>.txt` | Time, Position_Y, Velocity_Y, Acceleration_Y |
| `DOF_3_Body_<ID>.txt` | Time, Position_Z, Velocity_Z, Acceleration_Z |
| `DOF_4_Body_<ID>.txt` | Time, Rotation_RX, Velocity_RX, Acceleration_RX |
| `DOF_5_Body_<ID>.txt` | Time, Rotation_RY, Velocity_RY, Acceleration_RY |
| `DOF_6_Body_<ID>.txt` | Time, Rotation_RZ, Velocity_RZ, Acceleration_RZ |

**5 force files** per body:

| File | Columns |
|---|---|
| `HydroStiffnessForce_Body_<ID>.txt` | Time + 6 force/moment components |
| `WaveRadiationForce_Body_<ID>.txt` | Time + 6 force/moment components |
| `WaveExcitationForce_Body_<ID>.txt` | Time + 12 components (2 sets of 6) |
| `BCPForce_Body_<ID>.txt` | Time + BCP forces + coupled body forces |
| `WindTurbineForce_Body_<ID>.txt` | Time + wind turbine forces |

## Mooring Lines

**Per line:**

| File | Columns | Precision |
|---|---|---|
| `NodePosX_<ID>.txt` | Time + N node X-coordinates | `%.14f` |
| `NodePosY_<ID>.txt` | Time + N node Y-coordinates | `%.14f` |
| `NodePosZ_<ID>.txt` | Time + N node Z-coordinates | `%.14f` |
| `EndsTen_<ID>.txt` | Time, Tx₁, Ty₁, Tz₁, TxN, TyN, TzN | `%f` |
| `LineTen_<ID>.txt` | Time + N node tensions | `%f` |
| `LineIni_<ID>.txt` | Header (`s x y z ten`) + N rows of initial state | `%f` |

`LineIni_<ID>.txt` is written once at initialisation and includes a header row.

## Waves

| File | Columns | Header |
|---|---|---|
| `WaveSpectrum.txt` | Frequency (Hz), Spectral density (m²/Hz), Amplitude (m) | `# Wave Spectrum` |
| `WaveTimeSeries.txt` | Time (s), Elevation (m) | `# Wave Time Series` |
| `WavePhases.txt` | Phase matrix (Armadillo ASCII format) | — |

## Oscillating Water Columns

| File | Columns | Header |
|---|---|---|
| `OWC_<ID>.txt` | Time, Displacement, Velocity, Pressure | `Time  Displacement  Velocity  Pressure` |

## Sinking

| File | Columns |
|---|---|
| `SinkingFillingCOG_Body_<ID>.txt` | Time, COG_X, COG_Y, COG_Z, TotalFillingMass |

## Winches

| File | Columns |
|---|---|
| `WinchesTensions.txt` | Time + N tension values (one per winch) |
| `WinchedLinesLengths.txt` | Time + N line length values |

**Horizontal controller only:**

| File | Columns |
|---|---|
| `ControlForce.txt` | Time + 3 force components (X, Y, Z) |
| `ReferencePosition.txt` | Time + 3 reference position values |

## Wind Turbines

| File | Format |
|---|---|
| `FASTurbW_<ID>_bodyInerMat.dat` | 6×6 inertia matrix (Armadillo raw ASCII) |
| `FASTurbW_<ID>_towrInerMat.dat` | Tower inertia matrix |
| `FASTurbW_<ID>_turbInerMat.dat` | Turbine inertia matrix |

## ODE Solvers (Debug)

| File | Columns | Condition |
|---|---|---|
| `time_adaptivity_debug.txt` | t, dt, EWT, LTE, rho, rejected | Only if debug flag enabled |

## Simulation State

| File | Format | Condition |
|---|---|---|
| `dataStaticIC.dat` | Armadillo ASCII (state vector) | If `writeEquilibrium=1` |

## Notes

- All files use space delimiters (typically 4 spaces)
- Floating-point values use `%f` format unless noted otherwise
- Files are overwritten on each run (opened in write mode)
- Output folder path is specified in `dataProblem.dat` / YAML
