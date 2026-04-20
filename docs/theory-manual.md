# Theory Manual

The OASIS Theory Manual provides the mathematical formulations behind each simulation module. It is maintained as a LaTeX document and distributed as a PDF.

[:material-download: **Download Theory Manual (PDF)**](assets/OASIS_TheoryManual.pdf){ .md-button .md-button--primary }

## Contents

The theory manual covers:

| Chapter | Topic |
|---|---|
| §1 | **Introduction** — scope, notation and goals |
| §2 | **Coordinate systems and rigid body dynamics** — reference frames, Euler angles, rotation matrices, Newton–Euler equations |
| §3 | **Waves** — regular waves, irregular spectra (JONSWAP), directional spreading |
| §4 | **Hydrostatics** — linear restoring forces and non-linear Froude–Krylov on the instantaneous wetted surface |
| §5 | **Hydrodynamics** — Cummins equation, radiation/diffraction, first- and second-order excitation, QTFs, Morison wind/current |
| §6 | **Mooring and towing lines** — governing PDE, spectral element method (SEM), tension models (Rayleigh, viscoelastic), external forces, boundary conditions, quasi-static initial conditions |
| §7 | **Seafloor interaction** — flat, inclined and complex bathymetry projection algorithms, ground normal and friction forces |
| §8 | **Springs and connectors** — 6-DOF spring model, deformation measurement, elastic/damping/friction forces |
| §9 | **Winches** — mechanical model, controllers |
| §10 | **Oscillating water columns** — pneumatic chamber and turbine model |
| §11 | **Wind turbines** — OpenFAST coupling interface |
| §12 | **Sinking** — progressive flooding model |
| §13 | **Time integration** — adaptive BDF2, ESDIRK, Newton–Raphson iteration, error estimation |
| §14 | **Coupling** — global state vector, RHS evaluation order, weak coupling strategy |

## Building from Source

The LaTeX source is in `docs/theory-manual/OASIS_TheoryManual.tex`. To compile:

```bash
cd docs/theory-manual
pdflatex OASIS_TheoryManual.tex
pdflatex OASIS_TheoryManual.tex
```
