# OASIS Plotting Tools

This directory contains Python scripts for visualizing and analyzing OASIS simulation results.

## Available Scripts

### 1. `plot_oasis_results.py`
Main plotting tool for OASIS output data. Replaces previous MATLAB scripts (`plotBox.m`, `plot_case.m`).

**Features:**
- 6 DOF body motion time series plots
- Mooring line tension plots
- 3D animated visualization of floating structure
- Wind turbine data visualization (if available)

**Usage:**
```bash
# Plot results from specific output directory
python plot_oasis_results.py path/to/output/directory

# Plot results from default test case
python plot_oasis_results.py ../examples/test/output

# With custom configuration (edit script first)
python plot_oasis_results.py
```

**Configuration:**
Edit the script to customize:
- `flag_plot`: Enable/disable time series plots
- `flag_video`: Enable/disable 3D animation
- `flag_save`: Save figures to disk
- `flag_save_video`: Save animation as MP4 (requires ffmpeg)
- `flag_lines`: Include mooring lines
- `flag_wt`: Include wind turbine data
- Video parameters: time range, speed, viewing angle

### 2. `plot_spectrum_from_timeseries.py`
Plots wave spectrum computed from time-series data.

**Usage:**
```bash
python plot_spectrum_from_timeseries.py
```
Reads `WaveSpectrum.txt` from the current directory's output folder.

### 3. `plot_excitation_forces.py`
Plots excitation forces on bodies.

**Usage:**
```bash
# Specify output directory
python plot_excitation_forces.py path/to/output

# Or edit script to set default path
python plot_excitation_forces.py
```

### 4. `generate_wave_freq.py`
Generates frequency domain wave input files for OASIS.

**Usage:**
```bash
python generate_wave_freq.py
```
Creates `wavefreq.dat` with custom multi-directional wave components.

## Dependencies

Required Python packages:
```bash
pip install numpy matplotlib
```

Optional (for video export):
```bash
# Install ffmpeg for video generation
# Windows: download from https://ffmpeg.org/
# Linux: sudo apt install ffmpeg
# Mac: brew install ffmpeg
```

## Examples

### Basic plotting workflow
```bash
# 1. Run OASIS simulation
cd examples/test/input
oasis dataProblem.dat

# 2. Plot results
cd ../../..
python resources/plot_oasis_results.py examples/test/output
```

### Creating animated visualization
```python
from resources.plot_oasis_results import OASISResultsPlotter

plotter = OASISResultsPlotter('examples/test/output')
plotter.flag_video = True
plotter.flag_save_video = True
plotter.t_ini = 0
plotter.t_fin = 60
plotter.view_angle = [45, 15]
plotter.run()
```

### Analyzing wave spectrum
```bash
cd examples/test/output
python ../../resources/plot_spectrum_from_timeseries.py
```

## Migration from MATLAB

Previous MATLAB scripts have been replaced:

| MATLAB Script | Python Equivalent | Notes |
|--------------|-------------------|-------|
| `plotBox.m` | `plot_oasis_results.py` | 3D box generation integrated into main plotter |
| `plot_case.m` | `plot_oasis_results.py` | All functionality combined in OASISResultsPlotter class |

The Python version provides:
- ✓ Cross-platform compatibility
- ✓ No MATLAB license required
- ✓ Object-oriented design for easy customization
- ✓ Better integration with scientific Python ecosystem
- ✓ Command-line interface for automation

## Tips

1. **Large datasets**: For long simulations, adjust animation frame rate:
   ```python
   plotter.dt_video = 0.5  # Increase time step between frames
   ```

2. **Custom body shapes**: Modify `create_box_vertices()` method to create custom geometries

3. **Multiple bodies**: Extend `load_dof_data()` to handle multiple bodies (Body_1, Body_2, etc.)

4. **Batch processing**: Use Python scripting to process multiple cases:
   ```python
   import glob
   for case_dir in glob.glob('examples/*/output'):
       plotter = OASISResultsPlotter(case_dir)
       plotter.flag_save = True
       plotter.run()
   ```
