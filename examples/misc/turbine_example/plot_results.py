"""
Plot results for generic_example

This script plots the DOF time series for the generic OASIS example.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Add parent resources directory to path if needed
script_dir = os.path.dirname(os.path.abspath(__file__))
resources_dir = os.path.join(script_dir, '..', '..', 'resources')
sys.path.insert(0, resources_dir)

from plot_oasis_results import OASISResultsPlotter


def main():
    """Plot results for generic_example"""
    
    # Set output path relative to this script
    output_path = os.path.join(script_dir, 'output')
    
    if not os.path.exists(output_path):
        print(f"Error: Output directory not found: {output_path}")
        print("Run OASIS simulation first!")
        return
    
    # Create plotter instance
    plotter = OASISResultsPlotter(output_path)
    
    # Configure for this example
    plotter.flag_plot = True
    plotter.flag_video = False
    plotter.flag_save = True
    plotter.flag_lines = False
    plotter.flag_wt = False
    
    # Platform configuration
    plotter.body_dim = [20.0, 4.98, 2.25]
    plotter.body_color = [0.8, 0.8, 0.8]
    plotter.z_cog = 0.465
    plotter.z0 = 0.105
    plotter.water_depth = 10.0
    
    # Run plotting
    plotter.run()
    
    print("\nPlots saved to:", os.path.dirname(output_path))


if __name__ == "__main__":
    main()
