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
├── docs/                           MkDocs documentation source
│   └── theory-manual/              Theory manual (LaTeX source + assets)
└── resources/                      Python plotting utilities
```

---

## How to compile

### Prerequisites

| Dependency | Required | Notes |
|---|---|---|
| [Visual Studio](https://visualstudio.microsoft.com/) 2019 or later | Yes (Windows) | Install with the **"Desktop development with C++"** workload. Provides MSVC compiler, linker and Windows SDK. |
| C++17 compiler | Yes | MSVC ≥ 19 (from Visual Studio), GCC ≥ 10 |
| CMake ≥ 3.14 | Yes | |
| [Armadillo](http://arma.sourceforge.net/) | Yes | Linear algebra (header-only mode) |
| HDF5 with C++ bindings | Yes | I/O for hydrodynamic databases |
| BLAS / LAPACK | Yes | OpenBLAS or Intel MKL |
| [SuperLU](https://portal.nersc.gov/project/sparse/superlu/) | Recommended | Sparse solvers for dynamic mooring lines |
| [stl_reader](https://github.com/sreiter/stl_reader) | Optional | STL mesh import for body meshes |
| [OpenFAST](https://github.com/OpenFAST/openfast) + FASTurbine wrapper | Optional | Wind turbine coupling (see [Building OpenFAST and FASTurbine wrapper](#building-openfast-and-fasturbine-wrapper-optional) below) |
| [Intel oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html) | Only with OpenFAST | Provides the Intel Fortran compiler (`ifx`) and MKL, required to compile OpenFAST and link its libraries on Windows |

> **Note — MinGW vs MSVC:** If you have both MinGW and Visual Studio installed, CMake may pick MinGW by default. To force a specific Visual Studio generator, pass `-G` when configuring:
> ```bat
> cmake -B build -S . -G "Visual Studio 17 2022"
> ```
> Replace `17 2022` with your version (`16 2019`, `17 2022`, etc.). You can list available generators with `cmake --help`.

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
   cmake -B build -S . -G "Visual Studio 17 2022"
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

## Building OpenFAST and FASTurbine wrapper (optional)

These steps are only needed if you want to enable wind turbine coupling (`OASIS_USE_OPENFAST=ON`). The components must be built in this order, since each one depends on the previous:

1. **OpenFAST** — the aero-servo-elastic engine (standalone, no OASIS dependency).
2. **FASTurbine wrapper** — a thin DLL interface around OpenFAST (requires a compiled OpenFAST installation).
3. **OASIS** — linked against the FASTurbine wrapper (and transitively against OpenFAST).

### 1. Install Intel oneAPI HPC Toolkit

Download and install the [Intel oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html). This provides the Intel Fortran compiler (`ifx`) and MKL, both required to build OpenFAST on Windows.

### 2. Build OpenFAST 3.0.0

1. **Download the source code** from the [OpenFAST v3.0.0 release](https://github.com/OpenFAST/openfast/releases/tag/v3.0.0) on GitHub (choose *Source code (zip)*).

2. **Extract** the archive (e.g. to `D:\openfast-3.0.0`).

3. **Copy the build script** `resources/build_openfast.bat` into the OpenFAST root directory:
   ```bat
   copy resources\build_openfast.bat D:\openfast-3.0.0\
   ```

4. **Build** from a regular command prompt (the script auto-detects Visual Studio and Intel oneAPI):
   ```bat
   cd D:\openfast-3.0.0
   build_openfast.bat Release
   ```
   The script configures, compiles, and installs OpenFAST into `D:\openfast-3.0.0\install\`.

### 3. Build the FASTurbine wrapper

1. **Clone the wrapper** from IHCantabria's GitHub:
   ```bat
   git clone https://github.com/IHCantabria/FASTurbine_wrapper.git
   ```

2. **Edit `compile.bat`** inside the cloned repository and set `OPENFAST_INSTALL` to your OpenFAST install path:
   ```bat
   set OPENFAST_INSTALL=D:\openfast-3.0.0\install
   ```

3. **Compile**:
   ```bat
   cd FASTurbine_wrapper
   compile.bat
   ```
   This produces `install\lib\fasturbwrapper.dll` and `install\lib\fasturbwrapper.lib`.

### 4. Configure OASIS to use OpenFAST

In your `CMakeUserConfig.cmake`, add:
```cmake
set(OASIS_USE_OPENFAST ON)
set(OPENFAST_INCLUDE_DIR "D:/openfast-3.0.0/install/include")
set(FASTURBINE_INCLUDE_DIR "D:/FASTurbine_wrapper/install/include")
set(FASTURBINE_LIBRARY     "D:/FASTurbine_wrapper/install/lib/fasturbwrapper.lib")
set(FASTURBINE_DLL         "D:/FASTurbine_wrapper/install/lib/fasturbwrapper.dll")
```

Then rebuild OASIS with `compile.bat`.

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
