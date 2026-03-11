# OASIS — Offshore Advanced SImulation Software

OASIS is a time-domain numerical simulation tool for offshore floating structures, developed in C++ by [IHCantabria](https://ihcantabria.com). It couples potential-flow hydrodynamics with structural and mooring dynamics to predict the motion response of floating platforms under waves, wind and current.

## Features

- **Potential-flow hydrodynamics** — first-order radiation/diffraction forces from a pre-computed hydrodynamic database.
- **Second-order wave loads** — fast QTF-based evaluation of slow-drift and sum-frequency excitation forces.
- **Non-linear hydrostatics** — exact Froude–Krylov and restoring forces computed on the instantaneous wetted surface.
- **Dynamic mooring lines** — high-order Finite Element formulation with seabed contact, capable of adapting to irregular bathymetries.
- **Quasi-static mooring lines** — catenary-based solver for fast preliminary assessments.
- **Wind turbine coupling** — interface to OpenFAST for aero-servo-elastic simulation of floating offshore wind turbines.
- **Morison drag model** — current drag and viscous damping calibration via Morison coefficients.
- **Connectors & fenders** — versatile spring/damper model for multi-body contact and boundary coupling points.
- **Oscillating Water Column (OWC)** — simplified wave energy converter model.
- **OpenMP parallelisation** — multi-threaded execution on shared-memory machines.

## Project structure

```
oasis/
├── CMakeLists.txt                  Main build system
├── CMakeUserConfig.cmake.template  Template for local build paths
├── compile.bat / compile.sh        Quick build scripts (Windows / Linux)
├── src/                            C++ source code
├── examples/                       Ready-to-run example cases
├── resources/                      Python plotting utilities
├── cmake/                          CMake helper modules
└── doc/                            Theory manual (LaTeX)
```

---

## How to compile

### Prerequisites

| Dependency | Required | Notes |
|---|---|---|
| C++17 compiler | Yes | MSVC ≥ 19, GCC ≥ 10 |
| CMake ≥ 3.14 | Yes | |
| [Armadillo](http://arma.sourceforge.net/) | Yes | Linear algebra (header-only mode) |
| HDF5 with C++ bindings | Yes | I/O for hydrodynamic databases |
| BLAS / LAPACK | Yes | OpenBLAS or Intel MKL |
| [SuperLU](https://portal.nersc.gov/project/sparse/superlu/) | Recommended | Sparse solvers for dynamic mooring lines |
| [stl_reader](https://github.com/sreiter/stl_reader) | Optional | STL mesh import for body meshes |
| [OpenFAST](https://github.com/OpenFAST/openfast) + FASTurbine wrapper | Optional | Wind turbine coupling |
| Intel oneAPI Fortran runtime | Only with OpenFAST | Required to link Fortran-compiled OpenFAST libs on Windows |

### Windows (vcpkg — recommended)

1. **Install vcpkg** (if not already present):
   ```bat
   git clone https://github.com/microsoft/vcpkg C:\vcpkg
   C:\vcpkg\bootstrap-vcpkg.bat
   ```

2. **Install libraries through vcpkg**:
   ```bat
   C:\vcpkg\vcpkg install armadillo hdf5[cpp] openblas superlu metis:x64-windows
   ```

3. **Create your user configuration** by copying the template and editing paths:
   ```bat
   copy CMakeUserConfig.cmake.template CMakeUserConfig.cmake
   ```
   At a minimum, set `CMAKE_TOOLCHAIN_FILE` to your vcpkg toolchain:
   ```cmake
   set(CMAKE_TOOLCHAIN_FILE "C:/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
   ```
   See the template for all available options (Armadillo, BLAS/LAPACK, SuperLU, OpenFAST, stl_reader, etc.).

4. **Build** from a *Developer Command Prompt for VS* (or any shell with CMake and MSVC on PATH):
   ```bat
   compile.bat
   ```
   This configures, compiles, and copies the executable plus all required DLLs to `bin/`.

   Alternatively, configure and build manually:
   ```bat
   cmake -B build -S .
   cmake --build build --config Release
   ```

### Linux / HPC cluster

On EasyBuild-based clusters, OASIS auto-detects library paths from environment modules (`EBROOTGCCCORE`, `EBROOTArmadillo`, `EBROOTOPENBLAS`, etc.). Load the required modules, then:

```bash
chmod +x compile.sh
./compile.sh
```

Or manually:

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

The resulting executable is placed in `bin/oasis`.

### Build options

These CMake options can be set in `CMakeUserConfig.cmake` or passed on the command line (`-D...`):

| Option | Default | Description |
|---|---|---|
| `OASIS_USE_OPENFAST` | `OFF` | Enable OpenFAST wind turbine coupling |
| `OASIS_USE_SUPERLU` | `ON` | Enable SuperLU sparse solver via Armadillo |
| `OASIS_USE_STL_READER` | `ON` | Enable STL mesh reader for BodyMesh |

---

## How to run

### Basic usage

OASIS reads a problem definition file (typically `dataProblem.dat`) from the current working directory:

```bash
cd examples/generic_example/input
../../../bin/oasis dataProblem.dat
```

On Windows you can also use the `run.bat` script inside each example folder.

### Example cases

Five example cases are included under `examples/`. Run them all at once with:

```bat
cd examples
run_all_examples.bat        &:: Windows
./run_all_examples.sh       # Linux
```

| Example | Description |
|---|---|
| `generic_example` | Basic floating platform — JONSWAP waves, 6-DOF rigid body, linear hydrostatics |
| `freq_wave_example` | Frequency-domain wave input with spectral reconstruction |
| `qtf_example` | Second-order QTF wave loads — slow-drift and sum-frequency forces |
| `elastic_anchor_example` | Dynamic mooring lines with elastic anchor springs |
| `turbine_example` | Floating wind turbine coupled with OpenFAST |

Each example contains a `README.md` with detailed setup instructions and a `plot_results.py` script for post-processing.

### Python plotting utilities

The `resources/` directory provides general-purpose plotting scripts:

- `plot_oasis_results.py` — plot time-series output
- `plot_excitation_forces.py` — visualise excitation force RAOs
- `plot_spectrum_from_timeseries.py` — compute and plot spectral density from output

---

## License

Developed by IHCantabria. See license terms for details.
