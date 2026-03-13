"""
Plot results for freq_wave_example

This script plots the DOF time series and wave spectrum for the 
frequency domain wave input example.
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


def plot_wave_spectrum():
    """Plot the input wave spectrum"""
    
    spectrum_file = os.path.join(script_dir, 'output', 'WaveSpectrum.txt')
    
    if not os.path.exists(spectrum_file):
        print(f"Warning: Wave spectrum file not found: {spectrum_file}")
        return
    
    print("\nPlotting wave spectrum...")
    
    # Read spectrum data
    data = np.loadtxt(spectrum_file, skiprows=1)
    frequency = data[:, 0]
    spectral_density = data[:, 1]
    amplitude = data[:, 2]
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Plot spectral density
    ax1.plot(frequency, spectral_density, 'b-', linewidth=2)
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('Spectral Density [m²/Hz]')
    ax1.set_title('Wave Spectrum - Frequency Domain Input')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([0, max(frequency)])
    
    # Mark peak
    peak_idx = np.argmax(spectral_density)
    peak_freq = frequency[peak_idx]
    ax1.axvline(peak_freq, color='r', linestyle='--', label=f'Peak: {peak_freq:.4f} Hz')
    ax1.legend()
    
    # Plot amplitude
    ax2.plot(frequency, amplitude, 'g-', linewidth=2)
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Amplitude [m]')
    ax2.set_title('Wave Amplitude Spectrum')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim([0, max(frequency)])
    
    plt.tight_layout()
    
    # Save figure
    output_dir = os.path.dirname(spectrum_file)
    plt.savefig(os.path.join(output_dir, 'WaveSpectrum_plot.png'), dpi=150)
    print(f"  Saved to: {os.path.join(output_dir, 'WaveSpectrum_plot.png')}")
    
    plt.show()


def main():
    """Plot results for freq_wave_example"""
    
    # Set output path relative to this script
    output_path = os.path.join(script_dir, 'output')
    
    if not os.path.exists(output_path):
        print(f"Error: Output directory not found: {output_path}")
        print("Run OASIS simulation first!")
        return
    
    print("="*70)
    print("Frequency Domain Wave Example - Results Plotting")
    print("="*70)
    
    # Plot wave spectrum first
    plot_wave_spectrum()
    
    # Create plotter instance for body motions
    plotter = OASISResultsPlotter(output_path)
    
    # Configure for this example
    plotter.flag_plot = True
    plotter.flag_video = False
    plotter.flag_save = True
    plotter.flag_lines = False
    plotter.flag_wt = False
    
    # Platform configuration (ISOBARA)
    plotter.body_dim = [20.0, 4.98, 2.25]
    plotter.body_color = [0.8, 0.8, 0.8]
    plotter.z_cog = 0.465
    plotter.z0 = 0.105
    plotter.water_depth = 500.0
    
    # Wave parameters (combined systems)
    plotter.wave_height = [1.5, 0.8]  # Two wave systems
    plotter.wave_period = [10.0, 7.0]
    plotter.wave_phase = [0, 0]
    plotter.wave_direction = 0
    plotter.ramp_time = 0.0
    
    # Run plotting
    plotter.run()
    
    print("\nAll plots saved to:", os.path.dirname(output_path))


if __name__ == "__main__":
    main()
