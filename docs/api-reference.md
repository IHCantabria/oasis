# API Reference

The OASIS C++ API documentation is auto-generated from source code using [Doxygen](https://www.doxygen.nl/).

!!! info "Doxygen Output"
    The Doxygen HTML output is generated into `docs/api/` and is available alongside this site when deployed via CI. To browse locally, run `doxygen Doxyfile` and open `docs/api/index.html`.

## Generating Locally

To generate the API documentation on your machine:

1. Install [Doxygen](https://www.doxygen.nl/download.html)
2. Run from the repository root:
   ```bash
   doxygen Doxyfile
   ```
3. Open `docs/api/index.html` in your browser.

## Key Classes

| Class | Header | Description |
|---|---|---|
| `Simulation` | `src/Simulations/Simulation.hpp` | Central orchestrator — loads input, initialises subsystems, runs the time loop |
| `Body` | `src/Bodies/Bodies.hpp` | Floating rigid body with 6-DOF dynamics and hydrodynamic loading |
| `Line` | `src/Lines/Lines.hpp` | Mooring/towing line — dynamic (SEM) or quasi-static (catenary) |
| `Wave` | `src/Waves/Wave.hpp` | Wave generation and kinematics (regular, irregular, custom) |
| `Spring` | `src/Spring/Spring.hpp` | General spring/connector with linear/nonlinear stiffness and friction |
| `HydroDatabase` | `src/Hydro/HydroDatabase.hpp` | HDF5 reader for hydrodynamic databases |
| `HydroForce` | `src/Hydro/HydroForce.hpp` | Radiation, excitation and hydrostatic force computation |
| `Morison` | `src/Hydro/Morison.hpp` | Morison drag element for viscous damping |
| `SeaFloor` | `src/SeaFloor/SeaFloor.hpp` | Seabed geometry (flat, inclined, bathymetry mesh) |
| `OWC` | `src/OWC/OWC.hpp` | Oscillating water column chamber model |
| `WindTurbine` | `src/WindTurbine/WindTurbine.hpp` | OpenFAST wind turbine coupling |
| `BCP` | `src/BCPs/BCPs.hpp` | Boundary condition point (anchor, fairlead, joint, etc.) |
| `WinchieController` | `src/BCPs/WinchiesController.hpp` | Winch controller base class |
| `BDF2` / `BDFN` / `ESDIRK` | `src/ODE_solvers/ODE_solvers.hpp` | Implicit time integration solvers |
| `Sinking` | `src/Sinking/Sinking.hpp` | Progressive flooding model |
