# Building OASIS

## Prerequisites

| Dependency | Version | Notes |
|---|---|---|
| CMake | ≥ 3.14 | Build system generator |
| C++ compiler | C++17 support | GCC ≥ 7, Clang ≥ 5, MSVC ≥ 2017 |
| Armadillo | — | Linear algebra library |
| HDF5 | — | With C++ bindings |
| BLAS/LAPACK | — | OpenBLAS recommended on Linux |
| yaml-cpp | — | YAML parser |

See [Installation](../installation/windows.md) for platform-specific setup.

## Quick Build

=== "Windows"

    ```bat
    compile.bat
    ```

=== "Linux"

    ```bash
    chmod +x compile.sh
    ./compile.sh
    ```

## Manual CMake Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
```

The executable is placed in `bin/oasis` (Linux) or `bin/oasis.exe` (Windows).

## Build Options

| Option | Default | Description |
|---|---|---|
| `OASIS_USE_OPENFAST` | OFF | OpenFAST wind turbine coupling |
| `OASIS_USE_SUPERLU` | ON | SuperLU sparse solver |
| `OASIS_USE_STL_READER` | ON | STL mesh reader |

Example with options:

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DOASIS_USE_OPENFAST=ON \
         -DOASIS_USE_SUPERLU=OFF
```

## User Configuration

Create `CMakeUserConfig.cmake` (git-ignored) for local library paths:

```cmake
# Example: custom Armadillo location
set(ARMADILLO_INCLUDE_DIR "/opt/armadillo/include")
set(ARMADILLO_LIBRARY "/opt/armadillo/lib/libarmadillo.so")

# Example: custom HDF5 location
set(HDF5_ROOT "/opt/hdf5")
```

This file is loaded before `project()` in CMakeLists.txt.

## Compile Definitions

The build system sets these automatically:

| Define | Purpose |
|---|---|
| `ARMA_DONT_USE_WRAPPER` | Direct BLAS/LAPACK calls |
| `ARMA_USE_HDF5` | HDF5 support in Armadillo |
| `ARMA_USE_SUPERLU` | SuperLU support (if enabled) |
| `OASIS_USE_OPENFAST` | Wind turbine coupling code |
| `OASIS_USE_STL_READER` | STL mesh import code |

## Versioning

Version information is derived from git tags via `GetGitRevisionDescription.cmake` and injected into `version.cpp` at build time.

## Release Build Optimisations

Release builds automatically enable:

- `-O3 -march=native -ffast-math` (GCC/Clang) or `/O2 /fp:fast` (MSVC)
- Link-Time Optimisation (LTO)
- OpenMP threading (if available)
