# Examples

OASIS includes a systematic test suite organised into 10 functional groups. Each group isolates a specific subsystem so that every feature can be tested independently. The examples also serve as working configurations you can adapt for your own projects.

## Running Examples

=== "All examples"

    ```bash
    cd examples
    # Windows
    run_all_examples.bat
    # Linux
    ./run_all_examples.sh
    ```

=== "Single group"

    ```bash
    cd examples/<group>
    # Windows
    run_group.bat
    # Linux
    ./run_group.sh
    ```

=== "Individual example"

    ```bash
    cd examples/<group>/<test_name>
    # Windows
    run.bat
    # Linux — from repo root
    bin/oasis examples/<group>/<test_name>
    ```

## Shared Data

Hydrodynamic databases are stored in `examples/data/` and shared across examples:

| File | Bodies | Used by |
|---|---|---|
| `one_box.hydb.h5` | 1 box | Most body/lines/springs examples |
| `two_boxes.hydb.h5` | 2 boxes | Multi-body examples |
| `trl_plus.hydb.h5` | 1 platform | Turbine examples |
| `isobara_owcs.hydb.h5` | 7 (1 floater + 6 OWCs) | OWC examples |

These `.hydb.h5` files are produced by [**SeaMotions**](https://github.com/IHCantabria/SeaMotions), an open-source BEM solver by IHCantabria. Refer to the [SeaMotions repository](https://github.com/IHCantabria/SeaMotions) if you need to compute hydrodynamic databases for your own geometry.

## Functional Groups

### 1. Body Behavior (`body/`)

Floating body dynamics: free decay, fixed/partial DOFs, imposed motion, radiation, excitation, nonlinear hydrostatics, and QTF (Quadratic Transfer Function).

| Example | Feature tested |
|---|---|
| `free_decay` | Body released from offset, no waves |
| `fixed` | All DOFs locked, wave forces measured |
| `partial_dofs` | Subset of DOFs active |
| `imposed_motion` | Prescribed motion from file |
| `radiation` | Radiation force convolution |
| `excitation_linear` | First-order wave excitation (linear) |
| `excitation_instantpos` | Excitation at instantaneous position |
| `nonlinear_hs_flat` | Nonlinear hydrostatics (flat water) |
| `nonlinear_hs_waves` | Nonlinear hydrostatics with wave surface |
| `qtf` | Second-order slow-drift forces |

### 2. Wave Types (`waves/`)

Regular, irregular, and multi-directional seas.

| Example | Feature tested |
|---|---|
| `regular` | Monochromatic regular wave |
| `jonswap` | JONSWAP spectrum |
| `jonswap_piecewise` | Piecewise-decomposed JONSWAP |
| `multidirectional` | Multi-directional spreading |
| `timeseries` | Wave from time-series file |
| `frequency_domain` | Wave from frequency-domain file |

### 3. Mooring Lines (`lines/`)

Dynamic mooring lines with varied configurations.

| Example | Feature tested |
|---|---|
| `single_line` | Single catenary line |
| `multi_line` | Multi-line spread mooring |
| `joint_connection` | Lines connected via JointBCP |
| `elastic_anchor` | Anchor with elastic spring |
| `prescribed_motion` | Actuator-driven fairlead |
| `viscoelastic` | Viscoelastic material model |
| `tabulated_stiffness` | Tabulated force-strain curve |
| `tension_symmetric` | Symmetric tension/compression |
| `seabed_contact` | Ground contact model |
| `friction_isotropic` | Isotropic seabed friction |
| `friction_anisotropic` | Anisotropic seabed friction |

### 4. Springs & Connectors (`springs/`)

| Example | Feature tested |
|---|---|
| `pile_connector` | Pile-to-anchor spring with friction |
| `neoprene_connector` | Body-to-body neoprene+wire connection |

### 5. Wind Turbines (`turbines/`)

OpenFAST coupling.

| Example | Feature tested |
|---|---|
| `fixed_turbine` | Turbine on fixed-foundation body |
| `moored_turbine` | Turbine on moored floating platform |

### 6. Oscillating Water Columns (`owcs/`)

!!! warning "Work in Progress"
    OWC examples are under active development.

| Example | Feature tested |
|---|---|
| `free_decay` | OWC chamber free oscillation |
| `excitation` | Wave excitation on OWC |
| `hole_coupling` | Pneumatic coupling through orifice |
| `turbine_coupling` | OWC with Wells/impulse turbine |

### 7. Sinking (`sinking/`)

| Example | Feature tested |
|---|---|
| `progressive_flooding` | Compartment flooding and sinking |

### 8. Multi-Body (`multibody/`)

| Example | Feature tested |
|---|---|
| `shared_hydb` | Two bodies from one multi-body HDB |
| `separate_hydb` | Two bodies from individual HDBs |

### 9. ODE Solvers (`solvers/`)

Comparison of integration methods on a cable with circular actuator motion.

| Example | Solver |
|---|---|
| `bdf1` | BDF order 1 (backward Euler) |
| `bdfn_order2` | BDF order 2 |
| `bdfn_order4` | BDF order 4 |
| `esdirk46` | ESDIRK 4th order, 6 stages |

### 10. Winches (`winches/`)

| Example | Feature tested |
|---|---|
| `horizontal_control` | Surge/sway/yaw positioning with state-space controller |

## Example Layout

Each example follows a standard structure:

```
<group>/<test_name>/
├── input/
│   ├── dataProblem.dat      # Simulation parameters
│   ├── dataBodies.dat       # Body definitions
│   ├── dataBCPs.dat         # Boundary condition points
│   ├── dataWaves.dat        # Wave configuration
│   └── ...                  # Optional subsystem files
├── output/                  # Created by OASIS
├── run.bat                  # Windows runner
└── run.sl                   # Linux/SLURM runner
```

YAML equivalents are available for some examples (e.g., `free_decay_yaml/`, `fixed_yaml/`).
