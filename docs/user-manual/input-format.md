# Input Format

OASIS supports two input formats: the original **ASCII** (`.dat`) multi-file system and the newer **YAML** single-file format. Both provide the same functionality.

## ASCII Format (`.dat`)

The traditional format uses multiple text files in the `input/` directory, one per subsystem:

| File | Content | Always required |
|---|---|---|
| `dataProblem.dat` | General simulation parameters (gravity, time steps, tolerances, solver) | Yes |
| `dataBodies.dat` | Body definitions (hydrodynamic database, DOFs, initial position) | Yes |
| `dataBCPs.dat` | Boundary condition points (anchors, fairleads, joints) | Yes |
| `dataWaves.dat` | Wave definition (type, spectrum, direction) | Yes |
| `dataLines.dat` | Mooring line properties (elements, stiffness, friction) | If lines used |
| `dataSprings.dat` | Spring/connector definitions | If springs used |
| `dataWinches.dat` | Winch mechanical properties | If winches used |
| `dataWinchesController.dat` | Winch controller parameters | If winches used |
| `dataSeaFloor.dat` | Seabed geometry definitions | If non-flat seabed |
| `dataSinking.dat` | Progressive flooding compartments | If sinking used |

Each `.dat` file uses a simple line-based format where values are followed by comments:

```
9.81            // gravity [m/s^2]
1025.0          // waterDensity [kg/m^3]
```

Block separators use `//////////////////` comment lines.

See [Input Reference](input-reference.md) for the complete field-by-field specification of each file.

## YAML Format

The YAML format consolidates all configuration into a single file: `input/dataProblem.yaml`.

```yaml
problem:
  gravity: 9.81
  water_density: 1025.0
  water_depth: 50.0
  simulation_time: 300.0
  write_time_step: 0.1
  max_time_step: 0.1
  hydro_time_step: 0.1
  solver:
    method: 3          # ESDIRK46
    adaptivity: 0
    abs_tol: 1.0e-5
    rel_tol: 1.0e-3
    max_iter: 20

bodies:
  - type: RAD_DIFF
    database: ../../../data/one_box.hydb.h5
    database_index: 1
    dofs: [1, 2, 3, 4, 5, 6]
    position: [0, 0, 0, 0, 0, 0]
    freedom: 0          # 0=free, 1=locked, 2=imposed
    radiation: 1
    excitation: 1

bcps:
  anchors:
    - position: [200, 0, -50]
  body_bcps:
    - position: [10, 0, -5]

waves:
  type: IRR
  Hs: 2.0
  Tp: 10.0
  heading: 0.0
  ramp_time: 5.0
  spectrum:
    type: 1             # JONSWAP
    gamma: 3.3
```

!!! info "Format Detection"
    OASIS checks for `dataProblem.yaml` first. If found, YAML format is used automatically. No command-line flag is needed.

## BCP Type Mapping

The BCP (Boundary Condition Point) naming differs between formats. Here is the mapping:

| ASCII name | YAML key | Class | Description |
|---|---|---|---|
| "Actuators" | `fairleads` | `FairleadBCP` | Prescribed-motion point (reads actuator file) |
| "Anchors" | `anchors` | `AnchorBCP` | Fixed anchor point |
| "Joints" | `joints` | `JointBCP` | Multi-line joint with mass/volume |
| "Fairleads" | `body_bcps` | `BodyBCP` | Point attached to a body |
| "Elastic Anchors" | `elastic_anchors` | `ElasticAnchorBCP` | Anchor with elastic spring |

!!! warning
    BCP indices are assigned sequentially in this order: fairleads (actuators) first, then anchors, joints, body_bcps, elastic_anchors. The same ordering applies in both ASCII and YAML formats.

## Choosing a Format

| Aspect | ASCII | YAML |
|---|---|---|
| Maturity | Original format, all features | Newer, all features supported |
| Readability | Comments describe each value | Self-documenting keys |
| File count | 4–8 files per simulation | Single file |
| Editing | Simple text editor | Any YAML-aware editor |
| Existing examples | All examples have ASCII | ~15 examples have YAML equivalents |

Both formats can be used interchangeably. For new projects, YAML is recommended for its clarity.
