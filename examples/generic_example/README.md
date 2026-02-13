# OASIS generic_example

Generic example demonstrating basic OASIS simulation with a floating platform.

## Description

This example simulates a floating platform (ISOBARA-like) under regular or irregular waves. It includes:
- Potential flow hydrodynamics
- 6 DOF rigid body dynamics  
- First-order wave excitation forces
- Linear hydrostatics

## Files

### Input Files
- `dataProblem.dat`: General simulation parameters
- `dataWaves.dat`: Wave configuration
- `dataBodies.dat`: Body definition
- `dique_cadenas.hydb.h5`: Hydrodynamic database

### Output Files
- `DOF_X_Body_Y.txt`: Time series of body motions (6 DOF)
- `BCPForce_Body_Y.txt`: Boundary condition forces

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
- 6 DOF time series plots (surge, sway, heave, roll, pitch, yaw)
- Plots saved as PNG files in the example directory

## Configuration

Edit the input files to modify:
- Wave parameters (height, period, direction)
- Simulation time and time steps
- Body properties and initial conditions
- Force calculation options
