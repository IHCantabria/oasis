# OASIS Plotting and Utility Tools

This directory contains Python scripts for visualizing and analyzing OASIS simulation results, as well as utilities for generating input files.

## Core Plotting Tool

### `plot_oasis_results.py`
Main plotting library for OASIS output data. Used by all examples.

**Features:**
- 6 DOF body motion time series plots
- Mooring line tension plots
- 3D animated visualization of floating structure
- Wind turbine data visualization (if available)

**Note:** This is a library module. For plotting specific examples, use the `plot_results.py` script in each example directory.

## Specialized Plotting Scripts

### `plot_spectrum_from_timeseries.py`
Plots wave spectrum computed from time-series data with smoothing for noisy spectra.

**Usage:**
```bash
cd examples/<example_name>/output
python ../../../resources/plot_spectrum_from_timeseries.py
```

### `plot_excitation_forces.py`
Plots excitation forces on bodies.

**Usage:**
```bash
python plot_excitation_forces.py path/to/output
```

## Input Generation Tools

### `generate_wave_freq.py`
Template for generating frequency domain wave input files. 

**Note:** A customized version (`generate_wave_input.py`) is included in the `freq_wave_example` directory.

### `generate_datosPosicionFairlead.py`
Generates fairlead position data for mooring systems.

## Examples

All OASIS examples are located in the `examples/` directory with the following structure:

### Available Examples

1. **generic_example** - Basic floating platform simulation
   - Standard JONSWAP waves or regular waves
   - First-order hydrodynamics
   - Use: `cd examples/generic_example && python plot_results.py`

2. **qtf_example** - Second-order wave forces (QTF)
   - Time-series wave input with piecewise decomposition
   - First and second-order forces
   - Use: `cd examples/qtf_example && python plot_results.py`

3. **freq_wave_example** - Frequency domain wave input
   - Multi-directional wave systems
   - Custom spectral shapes
   - Use: `cd examples/freq_wave_example && python plot_results.py`

Each example contains:
- `input/` - Input data files
- `output/` - Simulation results (created after run)
- `plot_results.py` - Example-specific plotting script
- `README.md` - Detailed documentation
- `run.sl` - SLURM batch script for cluster execution

## Dependencies

Required Python packages:
```bash
pip install numpy matplotlib
```

Optional (for video export):
```bash
# Install ffmpeg for video generation
# Linux: sudo apt install ffmpeg
```

## Usage Workflow

### 1. Run simulation
```bash
cd examples/<example_name>/input
oasis dataProblem.dat
```

### 2. Plot results
```bash
cd ..
python plot_results.py
```

### 3. (Optional) Customize plots
Edit `plot_results.py` in the example directory to:
- Change plot settings
- Enable/disable video generation
- Modify body visualization parameters
- Adjust wave parameters for animation

## Creating Custom Examples

To create a new example based on existing ones:

1. Copy an existing example directory
2. Modify input files (dataProblem.dat, dataWaves.dat, etc.)
3. Update `plot_results.py` with appropriate parameters
4. Update README.md with example-specific information

## Tips

1. **Large datasets**: For long simulations, adjust animation frame rate in `plot_results.py`:
   ```python
   plotter.dt_video = 0.5  # Increase time step between frames
   ```

2. **Custom body shapes**: Modify `create_box_vertices()` in plot_oasis_results.py

3. **Batch processing**: Process multiple examples programmatically:
   ```python
   from plot_oasis_results import OASISResultsPlotter
   import glob
   
   for output_dir in glob.glob('examples/*/output'):
       plotter = OASISResultsPlotter(output_dir)
       plotter.flag_save = True
       plotter.run()
   ```

## Migration from MATLAB

Previous MATLAB scripts have been replaced with Python:

| Old MATLAB | Replacement | Location |
|-----------|-------------|----------|
| `plotBox.m` | Integrated into `plot_oasis_results.py` | resources/ |
| `plot_case.m` | `plot_results.py` for each example | examples/*/plot_results.py |
