# Mooring Lines

OASIS models mooring lines as dynamic finite elements using a Spectral Element Method (SEM). Lines connect pairs of Boundary Condition Points (BCPs) and can interact with the seabed.

## Modelling Approach

Each mooring line is discretised into finite elements with high-order Lagrange polynomials (spectral elements). The governing equations solve for:

- **Position** of each node in 3D
- **Tension** along the line
- **Contact forces** with the seabed
- **Friction forces** (if enabled)

Initial line configuration is computed by either:

- **Catenary** solution (`flagStatic = 0`) — analytical, fast
- **Newton's method** (`flagStatic = 1`) — numerical, handles complex geometry

## Line Properties

Each line is defined by its physical and numerical parameters:

| Parameter | Symbol | Units | Description |
|---|---|---|---|
| `lineType` | — | — | 1 = dynamic, 2 = quasi-static |
| `nNodes` | $N$ | — | Number of finite element nodes |
| `p` | $p$ | — | Polynomial order per element |
| `L` | $L$ | m | Unstretched line length |
| `rho0` | $\rho_0$ | kg/m | Mass per unit length |
| `d` | $d$ | m | Line diameter |
| `Cmn` | $C_{mn}$ | — | Normal added mass coefficient |
| `Cdn` | $C_{dn}$ | — | Normal drag coefficient |
| `Cdt` | $C_{dt}$ | — | Tangential drag coefficient |
| `CB` | $C_B$ | — | Seabed friction coefficient |
| `GK` | $G_K$ | N/m³ | Ground spring stiffness |
| `GC` | $G_C$ | — | Ground damping coefficient |

## Stiffness Models

Selected via `flag_stiffness`:

| Value | Model | Description |
|---|---|---|
| 0 | Linear | Constant axial stiffness $EA$ |
| 1 | Viscoelastic | Spring-dashpot model ($EA$, $EA_d$, $\alpha$) |
| > 1 | Tabulated | Force-strain lookup table with $N$ points |

### Linear (`flag_stiffness = 0`)

```
0               // flag_stiffness
1.0e8           // EA [N]
```

### Viscoelastic (`flag_stiffness = 1`)

```
1               // flag_stiffness
1.0e8           // EA [N]
5.0e7           // EA_d [N] (dashpot stiffness)
0.01            // alpha [-] (relaxation parameter)
```

### Tabulated (`flag_stiffness = N`)

```
5               // flag_stiffness (N = number of points)
0.0   0.0       // strain [-], force [N]
0.01  1.0e5
0.02  2.5e5
0.05  8.0e5
0.10  2.0e6
```

## Tension Model

Selected via `flag_tension`:

| Value | Behaviour |
|---|---|
| 0 | Tension-only — line goes slack when compressed |
| 1 | Symmetric — supports both tension and compression |

## Seabed Contact

Lines interact with the seabed via a spring-damper ground contact model:

$$F_{\text{ground}} = G_K \cdot \delta + G_C \cdot \dot{\delta}$$

where $\delta$ is the penetration depth below the seabed surface.

The seabed surface is defined via `indexSeaFloor`:

| Value | Seabed type |
|---|---|
| 0 | Flat at `waterDepth` |
| > 0 | Custom surface from `dataSeaFloor.dat` |

Custom seabed surfaces can be flat, inclined planes, or bathymetric meshes (STL).

## Friction Models

Selected via `frictionModel`:

| Value | Model | Description |
|---|---|---|
| 0 | None | No seabed friction |
| 1 | Isotropic | Same friction in all directions |
| 2 | Anisotropic | Different tangential and normal friction |

Friction uses a stick-slip model:

| Parameter | Symbol | Description |
|---|---|---|
| `vth` | $v_{th}$ | Threshold velocity for stick-slip transition [m/s] |
| `ust` / `usn` | $\mu_s$ | Static friction coefficient (tangential/normal) |
| `ud` | $\mu_d$ | Dynamic friction coefficient |
| `deltamax` | $\delta_{\max}$ | Maximum displacement for stick phase [m] |

## BCP Connections

Each line connects two BCPs (start and end nodes). The BCP types determine how the endpoints behave:

| BCP type | Class | Behaviour |
|---|---|---|
| Anchor | `AnchorBCP` | Fixed point in space |
| Fairlead (Body) | `BodyBCP` | Moves with a floating body |
| Joint | `JointBCP` | Multi-line junction with mass and volume |
| Elastic Anchor | `ElasticAnchorBCP` | Anchor with elastic spring |
| Actuator | `FairleadBCP` | Prescribed motion from file |

!!! info "BCP Index Order"
    BCP indices are assigned globally in order: actuators, anchors, joints, body fairleads, elastic anchors. Lines reference BCPs by these 1-based global indices.

## Input Example

=== "ASCII"

    In `dataLines.dat`:

    ```
    //////////////////
    ////////////////// New line [1]
    //////////////////
    1               // lineType (dynamic)
    0               // flag_tension (tension-only)
    21              // nNodes
    4               // p (polynomial order)
    200.0           // L [m]
    100.0           // rho0 [kg/m]
    0.1             // d [m]
    0               // flag_stiffness (linear)
    1.0e8           // EA [N]
    0.5             // CB
    1.0             // Cmn
    1.0             // Cdn
    0.1             // Cdt
    1.0e5           // GK [N/m³]
    100.0           // GC
    0               // indexSeaFloor (flat)
    1               // BCP_N (end node)
    4               // BCP_1 (start node)
    0               // frictionModel (none)
    0.001           // vth
    0.3             // ust
    0.2             // ud
    0.01            // deltamax
    ```

=== "YAML"

    In `dataProblem.yaml`:

    ```yaml
    lines:
      - line_type: 1
        tension_flag: 0
        nodes: 21
        polynomial_order: 4
        length: 200.0
        mass_per_length: 100.0
        diameter: 0.1
        stiffness:
          type: 0
          EA: 1.0e8
        friction_coef: 0.5
        Cmn: 1.0
        Cdn: 1.0
        Cdt: 0.1
        ground_stiffness: 1.0e5
        ground_damping: 100.0
        seafloor_index: 0
        BCP_end: 1
        BCP_start: 4
        friction_model: 0
        vth: 0.001
        ust: 0.3
        ud: 0.2
        deltamax: 0.01
    ```

## Output Files

| File | Content |
|---|---|
| `NodePosX_<ID>.txt` | Time + N node X positions |
| `NodePosY_<ID>.txt` | Time + N node Y positions |
| `NodePosZ_<ID>.txt` | Time + N node Z positions |
| `EndsTen_<ID>.txt` | Time + tension at both ends (6 components) |
| `LineTen_<ID>.txt` | Time + N node tensions |
| `LineIni_<ID>.txt` | Initial configuration (s, x, y, z, tension) |

## Related Examples

| Example | Feature |
|---|---|
| `lines/single_line` | Single catenary mooring |
| `lines/multi_line` | Multi-line spread mooring |
| `lines/joint_connection` | Lines connected via JointBCP |
| `lines/elastic_anchor` | Elastic anchor point |
| `lines/prescribed_motion` | Actuator-driven fairlead |
| `lines/viscoelastic` | Viscoelastic material |
| `lines/tabulated_stiffness` | Tabulated force-strain curve |
| `lines/tension_symmetric` | Symmetric tension/compression |
| `lines/seabed_contact` | Ground contact model |
| `lines/friction_isotropic` | Isotropic seabed friction |
| `lines/friction_anisotropic` | Anisotropic seabed friction |
