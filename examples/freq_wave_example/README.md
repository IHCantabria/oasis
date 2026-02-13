# OASIS freq_wave_example

Example demonstrating **frequency domain wave input** (specType_flag = 3) in OASIS.

## Description

This example shows how to define waves directly in the frequency domain. The wave input combines two wave systems:

- **System 1**: JONSWAP spectrum with Tp=10s, Hs=1.5m, heading=0°
- **System 2**: JONSWAP spectrum with Tp=7s, Hs=0.8m, heading=30°

This approach is ideal for:
- Multi-directional sea states
- Combining wind sea and swell from different directions
- Custom spectral shapes not available in analytical forms
- Direct output from wave measurement analysis

## Wave Input Format

The frequency domain wave file (`wavefreq.dat`) has the following format:

```
<number_of_components>  // Number of frequency components
<freq1> <heading1> <amplitude1> <phase1>
<freq2> <heading2> <amplitude2> <phase2>
...
```

Where:
- **frequency**: [Hz]
- **heading**: [degrees] (0° = +X axis, counterclockwise positive)
- **amplitude**: [m]
- **phase**: [radians]

## Files

### Input Files
- `dataProblem.dat`: General simulation parameters
- `dataWaves.dat`: Wave configuration (specType_flag=3)
- `dataBodies.dat`: Body definition
- `wavefreq.dat`: Frequency domain wave data (48 components)
- `wavefreq_spectrum.png`: Visualization of input spectrum
- `isobara.ehydb`: Hydrodynamic database

### Output Files
- `DOF_X_Body_Y.txt`: Time series of body motions
- `WaveSpectrum.txt`: Wave spectrum data
- Various force output files

### Scripts
- `generate_wave_input.py`: Generate custom frequency domain wave files
- `plot_results.py`: Visualize simulation results

## Generating Custom Wave Input

To create your own frequency domain wave file:

```bash
python generate_wave_input.py
```

This will create `wavefreq.dat` in the current directory. Move it to the `input/` folder:

```bash
mv wavefreq.dat input/
```

Edit the script to customize:
- Wave heights and periods for each system
- Number of wave systems to combine
- Heading directions
- Frequency range and resolution
- Random phases (set seed for reproducibility)

## Running the Example

### Linux (with SLURM)
```bash
sbatch run.sl input/dataProblem.dat
```

### Manual execution
```bash
cd input
oasis dataProblem.dat
cd ..
```

## Plotting Results

After running the simulation:

```bash
python plot_results.py
```

This will generate:
- 6 DOF time series plots
- Wave spectrum visualization (input vs computed)
- Amplitude spectrum
- All plots saved to the example directory

## Important Notes

- Frequency domain waves (specType_flag=3) **do not support piecewise decomposition**
- All wave components are computed for the entire simulation duration
- Ensure frequency components cover the range of interest for your structure
- Random phases are pre-generated for reproducibility
- Combined Hs for this example: 1.70m (RSS of both systems)
