# OASIS qtf_example

Example demonstrating second-order wave forces (QTF - Quadratic Transfer Functions) in OASIS.

## Description

This example simulates a floating platform under irregular waves with second-order excitation forces included. This is important for:
- Slow-drift motions
- Wave frequency response
- Long-period dynamics

The example uses:
- JONSWAP spectrum or custom time-series waves
- Piecewise wave decomposition for computational efficiency
- First and second order wave excitation forces
- ISOBARA floating platform

## Files

### Input Files
- `dataProblem.dat`: General simulation parameters (long simulation time)
- `dataWaves.dat`: Wave configuration (specType_flag=2 for time-series)
- `dataBodies.dat`: Body definition with 2nd order forces enabled
- `isobara.ehydb`: Hydrodynamic database with QTF data
- `wave_long.dat`: Long time-series wave data
- `2ndorder_long.csv`: Second-order force coefficients

### Output Files
- `DOF_X_Body_Y.txt`: Time series of body motions
- `WaveSpectrum.txt`: Wave spectrum data
- `WaveTimeSeries.txt`: Generated wave elevation time series

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

After running the simulation, visualize the results:

```bash
python plot_results.py
```

This will generate:
- 6 DOF time series plots
- Wave spectrum plots (if available)
- All plots saved to the example directory

## Notes

- Second-order forces require QTF data in the hydrodynamic database
- Piecewise decomposition is recommended for long simulations
- This example may take longer to run due to 2nd order force computations
