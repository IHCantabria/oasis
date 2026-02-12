import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os

def plot_spectrum():
    """
    Plot the wave spectrum from WaveSpectrum.txt output file.
    
    The WaveSpectrum.txt file format is:
    # Wave Spectrum
    # Frequency(Hz)  SpectralDensity(m^2/Hz)  Amplitude(m)
    frequency spectral_density amplitude
    """
    
    # Define file paths relative to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    
    # Try multiple possible locations
    possible_paths = [
        os.path.join(base_dir, 'examples', 'test', 'output', 'WaveSpectrum.txt'),
        os.path.join(base_dir, 'examples', 'test_qtf', 'output', 'WaveSpectrum.txt'),
        os.path.join(base_dir, 'output', 'WaveSpectrum.txt'),
        os.path.join(base_dir, 'WaveSpectrum.txt'),
    ]
    
    spectrum_file = None
    for path in possible_paths:
        if os.path.exists(path):
            spectrum_file = path
            break
    
    if spectrum_file is None:
        print(f"ERROR: WaveSpectrum.txt not found in any of these locations:")
        for path in possible_paths:
            print(f"  - {path}")
        return
    
    print(f"Reading wave spectrum from: {spectrum_file}")
    
    # Read spectrum data
    try:
        data = pd.read_csv(spectrum_file, delim_whitespace=True, comment='#', header=None)
        frequency = data.iloc[:, 0].values.astype(np.float64)
        spectral_density = data.iloc[:, 1].values.astype(np.float64)
        amplitude = data.iloc[:, 2].values.astype(np.float64)
        
        print(f"  Number of frequency points: {len(frequency)}")
        print(f"  Frequency range: {frequency[0]:.4f} to {frequency[-1]:.4f} Hz")
        print(f"  Max spectral density: {np.max(spectral_density):.6f} m^2/Hz")
        print(f"  Max amplitude: {np.max(amplitude):.6f} m")
        
        # Apply smoothing to reduce noise
        window_size = max(5, len(frequency) // 100)
        if window_size % 2 == 0:
            window_size += 1
        print(f"  Applying smoothing filter (window size: {window_size} points)...")
        spectral_density_smooth = np.convolve(spectral_density, np.ones(window_size)/window_size, mode='same')
        
        # Find spectrum peak (using smoothed spectrum)
        max_idx_smooth = np.argmax(spectral_density_smooth)
        peak_freq_smooth = frequency[max_idx_smooth]
        peak_density_smooth = spectral_density_smooth[max_idx_smooth]
        print(f"  Peak frequency (smoothed): {peak_freq_smooth:.4f} Hz (Period: {1/peak_freq_smooth:.4f} s)")
        
        # Compute FWHM from raw spectrum
        max_idx = np.argmax(spectral_density)
        peak_freq = frequency[max_idx]
        peak_density = spectral_density[max_idx]
        half_max = peak_density / 2.0
        above_half = spectral_density > half_max
        if np.any(above_half):
            indices = np.where(above_half)[0]
            f1_raw = frequency[indices[0]]
            f2_raw = frequency[indices[-1]]
            fwhm_raw = f2_raw - f1_raw
            print(f"  FWHM (raw): {fwhm_raw:.6f} Hz (f1={f1_raw:.4f} Hz, f2={f2_raw:.4f} Hz)")
        else:
            print(f"  WARNING: Could not compute FWHM from raw spectrum")
            fwhm_raw = None
        
        # Compute FWHM from smoothed spectrum
        half_max_smooth = peak_density_smooth / 2.0
        above_half_smooth = spectral_density_smooth > half_max_smooth
        if np.any(above_half_smooth):
            indices_smooth = np.where(above_half_smooth)[0]
            f1_smooth = frequency[indices_smooth[0]]
            f2_smooth = frequency[indices_smooth[-1]]
            fwhm_smooth = f2_smooth - f1_smooth
            print(f"  FWHM (smoothed): {fwhm_smooth:.6f} Hz (f1={f1_smooth:.4f} Hz, f2={f2_smooth:.4f} Hz)")
        else:
            print(f"  WARNING: Could not compute FWHM from smoothed spectrum")
            fwhm_smooth = None
        
    except Exception as e:
        print(f"ERROR reading spectrum file: {e}")
        return
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Plot 1: Spectral Density (Raw and Smoothed)
    ax1.plot(frequency, spectral_density, 'b-', linewidth=0.7, alpha=0.5, label='Raw Spectrum')
    ax1.plot(frequency, spectral_density_smooth, 'b-', linewidth=2, label='Smoothed Spectrum')
    ax1.axvline(peak_freq_smooth, color='r', linestyle='--', linewidth=1.5, label=f'Peak: {peak_freq_smooth:.4f} Hz')
    if fwhm_smooth is not None:
        ax1.axhline(half_max_smooth, color='g', linestyle=':', linewidth=1, label=f'Half Max: {half_max_smooth:.4f}')
        ax1.axvline(f1_smooth, color='orange', linestyle=':', linewidth=1.5, alpha=0.7)
        ax1.axvline(f2_smooth, color='orange', linestyle=':', linewidth=1.5, alpha=0.7)
        ax1.fill_between([f1_smooth, f2_smooth], 0, ax1.get_ylim()[1], alpha=0.2, color='yellow', 
                         label=f'FWHM: {fwhm_smooth:.4f} Hz')
    
    ax1.set_xlabel('Frequency (Hz)', fontsize=12)
    ax1.set_ylabel('Spectral Density (m²/Hz)', fontsize=12)
    ax1.set_title('Wave Spectrum (Raw vs Smoothed)', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='best')
    
    # Plot 2: Amplitude Spectrum
    ax2.plot(frequency, amplitude, 'g-', linewidth=1.5, label='Amplitude')
    ax2.axvline(peak_freq_smooth, color='r', linestyle='--', linewidth=1, label=f'Peak: {peak_freq_smooth:.4f} Hz')
    
    ax2.set_xlabel('Frequency (Hz)', fontsize=12)
    ax2.set_ylabel('Amplitude (m)', fontsize=12)
    ax2.set_title('Wave Amplitude Spectrum', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best')
    
    plt.tight_layout()
    
    # Save figure in the same directory as the spectrum file
    output_dir = os.path.dirname(spectrum_file)
    output_plot = os.path.join(output_dir, 'WaveSpectrum.png')
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_plot}")
    
    # Show plot
    plt.show()


if __name__ == "__main__":
    print("=" * 70)
    print("Wave Spectrum Plotting Tool")
    print("=" * 70)
    plot_spectrum()
    print("=" * 70)
