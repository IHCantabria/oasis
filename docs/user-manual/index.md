# User Manual

This section covers everything you need to set up, run, and analyse OASIS simulations.

## Getting Started

1. **[Graphical Interface (GUI)](gui.md)** — launch the PyQt5 GUI to build cases, run simulations and plot results without editing files by hand
2. **[Running Simulations](running-simulations.md)** — command-line usage, folder conventions, input format auto-detection
2. **[Input Format](input-format.md)** — overview of ASCII (`.dat`) and YAML input systems
3. **[Input Reference](input-reference.md)** — comprehensive field-by-field reference for all input files

## Understanding Results

4. **[Output Files](output-files.md)** — every output file produced by each module, with column descriptions

## Learning from Examples

5. **[Examples](examples.md)** — 10 functional test groups with descriptions and setup instructions

## Feature Guides

Detailed guides for each simulation subsystem:

- **[Waves](waves.md)** — regular, irregular (JONSWAP), time-series and frequency-domain wave input
- **[Mooring Lines](mooring-lines.md)** — dynamic FEM and quasi-static catenary solvers
- **[Springs & Connectors](springs.md)** — linear/nonlinear connectors with friction
- **[Winches](winches.md)** — winch-controlled mooring with tension and position controllers
- **[Oscillating Water Columns](owc.md)** — OWC chamber dynamics and turbine coupling
- **[Wind Turbines](wind-turbines.md)** — OpenFAST aero-servo-elastic coupling
- **[ODE Solvers](solvers.md)** — BDF and ESDIRK time integration methods
