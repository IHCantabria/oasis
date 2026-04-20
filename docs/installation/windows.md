# Installation — Windows

## Prerequisites

| Dependency | Required | Notes |
|---|---|---|
| [Visual Studio](https://visualstudio.microsoft.com/) 2019+ | Yes | Install with the **"Desktop development with C++"** workload |
| C++17 compiler (MSVC ≥ 19) | Yes | Included with Visual Studio |
| CMake ≥ 3.14 | Yes | Included with Visual Studio or install separately |
| [Armadillo](http://arma.sourceforge.net/) | Yes | Linear algebra (header-only mode) |
| HDF5 with C++ bindings | Yes | I/O for hydrodynamic databases |
| BLAS / LAPACK | Yes | OpenBLAS or Intel MKL |
| [SuperLU](https://portal.nersc.gov/project/sparse/superlu/) | Recommended | Sparse solvers for dynamic mooring lines |
| [stl_reader](https://github.com/sreiter/stl_reader) | Optional | STL mesh import for body meshes |
| [OpenFAST](https://github.com/OpenFAST/openfast) + FASTurbine wrapper | Optional | Wind turbine coupling |

## Using vcpkg (Recommended)

### 1. Install vcpkg

```bat
git clone https://github.com/microsoft/vcpkg C:\vcpkg
C:\vcpkg\bootstrap-vcpkg.bat
```

### 2. Install dependencies

```bat
C:\vcpkg\vcpkg install armadillo hdf5[cpp] openblas superlu metis:x64-windows
```

### 3. Configure the project

Copy the configuration template and set the vcpkg toolchain path:

```bat
copy CMakeUserConfig.cmake.template CMakeUserConfig.cmake
```

Edit `CMakeUserConfig.cmake` and set at minimum:

```cmake
set(CMAKE_TOOLCHAIN_FILE "C:/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
```

See the template file for all available options (Armadillo paths, BLAS/LAPACK, SuperLU, OpenFAST, stl_reader, etc.).

### 4. Build

From a **Developer Command Prompt for VS** (or any shell with CMake and MSVC on PATH):

```bat
compile.bat
```

This configures CMake, compiles in Release mode, and copies the executable plus all required DLLs to `bin/`.

Alternatively, configure and build manually:

```bat
cmake -B build -S . -G "Visual Studio 17 2022"
cmake --build build --config Release
```

!!! tip "MinGW vs MSVC"
    If you have both MinGW and Visual Studio installed, CMake may pick MinGW by default. Force Visual Studio with:
    ```bat
    cmake -B build -S . -G "Visual Studio 17 2022"
    ```
    Replace `17 2022` with your version (`16 2019`, etc.). List available generators with `cmake --help`.

## Build Options

These CMake options can be set in `CMakeUserConfig.cmake` or passed on the command line (`-D...`):

| Option | Default | Description |
|---|---|---|
| `OASIS_USE_OPENFAST` | `OFF` | Enable OpenFAST wind turbine coupling |
| `OASIS_USE_SUPERLU` | `ON` | Enable SuperLU sparse solver via Armadillo |
| `OASIS_USE_STL_READER` | `ON` | Enable STL mesh reader for BodyMesh |

---

## OpenFAST Setup (Optional)

These steps are only needed if you want wind turbine coupling (`OASIS_USE_OPENFAST=ON`). Build in this order — each component depends on the previous one.

### 1. Install Intel oneAPI HPC Toolkit

Download and install the [Intel oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html). This provides the Intel Fortran compiler (`ifx`) and MKL, both required for OpenFAST on Windows.

### 2. Build OpenFAST 3.0.0

1. Download [OpenFAST v3.0.0 source](https://github.com/OpenFAST/openfast/releases/tag/v3.0.0) (zip)
2. Extract (e.g. to `D:\openfast-3.0.0`)
3. Copy the build script:
   ```bat
   copy resources\build_openfast.bat D:\openfast-3.0.0\
   ```
4. Build from a regular command prompt:
   ```bat
   cd D:\openfast-3.0.0
   build_openfast.bat Release
   ```
   This installs OpenFAST into `D:\openfast-3.0.0\install\`.

### 3. Build the FASTurbine wrapper

```bat
git clone https://github.com/IHCantabria/FASTurbine_wrapper.git
```

Edit `compile.bat` inside the cloned repository and set the OpenFAST path:

```bat
set OPENFAST_INSTALL=D:\openfast-3.0.0\install
```

Then compile:

```bat
cd FASTurbine_wrapper
compile.bat
```

### 4. Configure OASIS for OpenFAST

Add to `CMakeUserConfig.cmake`:

```cmake
set(OASIS_USE_OPENFAST ON)
set(OPENFAST_INCLUDE_DIR    "D:/openfast-3.0.0/install/include")
set(FASTURBINE_INCLUDE_DIR  "D:/FASTurbine_wrapper/install/include")
set(FASTURBINE_LIBRARY      "D:/FASTurbine_wrapper/install/lib/fasturbwrapper.lib")
set(FASTURBINE_DLL          "D:/FASTurbine_wrapper/install/lib/fasturbwrapper.dll")
```

Rebuild OASIS with `compile.bat`.
