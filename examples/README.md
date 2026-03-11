# OASIS Examples

This directory contains example simulations demonstrating different features and capabilities of OASIS.

## Available Examples

### 1. generic_example
**Basic floating platform simulation**

Features:
- Standard JONSWAP or regular wave input
- First-order hydrodynamic forces
- 6 DOF rigid body dynamics
- Linear hydrostatics

Use case: Learning basic OASIS workflow and setup

**Quick start:**
```bash
cd generic_example/input
oasis dataProblem.dat
cd ..
python plot_results.py
```

---

### 2. qtf_example
**Second-order wave forces (QTF - Quadratic Transfer Functions)**

Features:
- Time-series wave input
- Piecewise wave decomposition for long simulations
- First and second-order excitation forces
- Slow-drift and wave-frequency response

Use case: Accurate long-term simulations with mean drift forces

**Quick start:**
```bash
cd qtf_example/input
oasis dataProblem.dat
cd ..
python plot_results.py
```

---

### 3. freq_wave_example
**Frequency domain wave input (specType_flag=3)**

Features:
- Direct frequency domain wave definition
- Multi-directional wave systems
- Combining multiple spectra (e.g., wind sea + swell)
- Custom spectral shapes

Use case: Complex sea states with multiple wave systems from different directions

**Quick start:**
```bash
# Optional: Generate custom wave input
cd freq_wave_example
python generate_wave_input.py
mv wavefreq.dat input/

# Run simulation
cd input
oasis dataProblem.dat
cd ..
python plot_results.py
```

---

### 4. turbine_example
**OpenFAST wind turbine integration test**

Features:
- Fixed-position floating platform (no wave forces)
- OpenFAST/FASTurbine wind turbine coupling
- Wind turbine aerodynamic and structural dynamics
- Minimal hydrodynamic computation for performance testing

Use case: Testing OpenFAST integration, wind turbine only simulations, debugging turbine setup

**Quick start:**
```bash
cd turbine_example/input
oasis dataProblem.dat
cd ..
python plot_results.py
```

**Note:** Requires OASIS compiled with OpenFAST support (`OASIS_USE_OPENFAST=ON`)

---

### 5. elastic_anchor_example
**Elastic anchor mooring test**

Features:
- Elastic (compliant) mooring anchors with exponential restoring force
- 4 mooring lines with elastic anchor boundary conditions
- Linear hydrostatics, no wind/turbines
- Tests anchor dragging and restoring forces

Use case: Validating elastic anchor force models for mooring systems

**Quick start:**
```bash
cd elastic_anchor_example/input
oasis dataProblem.dat
```

---

### 6. large_rotation_example
**Simplified vs full rotation dynamics comparison**

Features:
- Two sub-cases: simplified (`rotSimpFlag=1`) and full (`rotSimpFlag=0`) rotation dynamics
- Free-decay test with initial heave + pitch displacement
- Radiation forces enabled for damped oscillation
- No waves, no mooring lines — isolates rotation mechanics

Use case: Validating full rotation dynamics against simplified approximation

**Quick start:**
```bash
cd large_rotation_example
run.bat
python plot_results.py
```

**Note:** Runs two simulations (`case_simplified/` and `case_full/`) and produces an overlay comparison plot.

---

## Example Structure

Each example follows this standard structure:

```
example_name/
├── input/                      # Input files
│   ├── dataProblem.dat        # Simulation parameters
│   ├── dataWaves.dat          # Wave configuration
│   ├── dataBodies.dat         # Body definitions
│   └── ...                    # Other input files
├── output/                     # Results (created after run)
│   ├── DOF_*.txt              # DOF time series
│   ├── WaveSpectrum.txt       # Wave spectrum data
│   └── ...                    # Other outputs
├── plot_results.py            # Plotting script for this example
├── run.sl                     # SLURM batch script (for cluster)
└── README.md                  # Example documentation
```

## Running Examples

### On a Cluster (Linux with SLURM)

```bash
cd <example_name>
sbatch run.sl input/dataProblem.dat
```

### Manual Execution

```bash
cd <example_name>/input
oasis dataProblem.dat
```

## Plotting Results

After running a simulation, visualize results using the example-specific plotting script:

```bash
cd <example_name>
python plot_results.py
```

This will:
- Load simulation results from `output/`
- Generate time series plots of body motions
- Create wave spectrum visualizations (if applicable)
- Save plots as PNG files in the example directory

## Customizing Examples

To create a new example:

1. **Copy an existing example:**
   ```bash
   cp -r generic_example my_new_example
   ```

2. **Modify input files:**
   - Edit `dataProblem.dat` for simulation parameters
   - Edit `dataWaves.dat` for wave conditions
   - Edit `dataBodies.dat` for body properties
   - Update hydrodynamic database file if needed

3. **Update plotting script:**
   - Modify `plot_results.py` to match your configuration
   - Adjust body dimensions, colors, wave parameters

4. **Update documentation:**
   - Edit `README.md` with example-specific information

## Tips and Best Practices

### Choosing the Right Wave Input

- **Regular waves (specType_flag=0)**: Simple testing, resonance studies
- **JONSWAP spectrum (specType_flag=1)**: Standard irregular seas, moderate durations
- **Time-series (specType_flag=2)**: From measurements, long simulations with piecewise
- **Frequency domain (specType_flag=3)**: Multi-directional seas, custom spectra

### Simulation Time

- Short runs (< 1 hr): Direct computation, single piece
- Long runs (> 1 hr): Enable piecewise decomposition
- Very long runs: Use time-series with piecewise (qtf_example)

### Second-Order Forces

Enable in `dataBodies.dat` for:
- Slow-drift motions
- Moored structures
- Long simulation times
- Requires QTF data in hydrodynamic database

### Performance Optimization

- Use piecewise decomposition for irregular waves (`piecewise_flag = 1`)
- Adjust time steps based on needed accuracy
- Balance hydrodynamic vs integration time steps

## Common Issues

### Output files not found
**Solution:** Verify simulation completed successfully. Check for error messages in terminal output.

### Plotting errors
**Solution:** Ensure Python dependencies are installed:
```bash
pip install numpy matplotlib
```

### Missing hydrodynamic database
**Solution:** Ensure `.ehydb` or `.hydb.h5` file is in the input directory.

## Further Information

- See `../resources/README.md` for plotting tool documentation
- See individual `README.md` in each example for detailed information
- Refer to OASIS documentation for complete simulation setup guide
