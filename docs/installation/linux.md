# Installation — Linux / HPC

## Prerequisites

Install the required development libraries using your distribution's package manager.

=== "Ubuntu / Debian"

    ```bash
    sudo apt install build-essential cmake libarmadillo-dev libhdf5-dev \
        libopenblas-dev libsuperlu-dev
    ```

=== "CentOS / RHEL"

    ```bash
    sudo yum install gcc-c++ cmake3 armadillo-devel hdf5-devel \
        openblas-devel SuperLU-devel
    ```

=== "Fedora"

    ```bash
    sudo dnf install gcc-c++ cmake armadillo-devel hdf5-devel \
        openblas-devel SuperLU-devel
    ```

Minimum compiler requirement: **GCC ≥ 10** (C++17 support).

## Building

### Quick Build

```bash
chmod +x compile.sh
./compile.sh
```

### Manual Build

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

The executable is placed in `bin/oasis`.

### Build Options

| Option | Default | Description |
|---|---|---|
| `OASIS_USE_OPENFAST` | `OFF` | Enable OpenFAST wind turbine coupling |
| `OASIS_USE_SUPERLU` | `ON` | Enable SuperLU sparse solver |
| `OASIS_USE_STL_READER` | `ON` | Enable STL mesh reader |

Pass options on the command line:

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DOASIS_USE_SUPERLU=ON
```

## HPC Clusters (EasyBuild)

On EasyBuild-based clusters, OASIS auto-detects library paths from environment modules. Load the required modules first:

```bash
module load GCCcore Armadillo OpenBLAS HDF5 SuperLU
```

Then build normally:

```bash
./compile.sh
```

The CMake configuration reads the `EBROOTGCCCORE`, `EBROOTArmadillo`, `EBROOTOPENBLAS`, etc. environment variables set by the modules.

### Deploying to a Cluster

Use the provided sync script to copy the project to a remote cluster:

```bash
# Edit sync_to_cluster.bat with your cluster hostname and path
sync_to_cluster.bat
```

## Verifying the Installation

Run a quick example to verify the build:

```bash
cd examples/body/free_decay/input
../../../../bin/oasis .
```

The simulation should complete without errors and produce output files in the `output/` directory.
