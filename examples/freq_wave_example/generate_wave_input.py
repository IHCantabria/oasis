"""
Generate frequency domain wave data for OASIS test_freq example
This script creates a custom wave combining two wave systems:
- System 1: JONSWAP with Tp=10s, Hs=1.5m, heading=0°
- System 2: JONSWAP with Tp=7s, Hs=0.8m, heading=30°
"""

import numpy as np
import matplotlib.pyplot as plt

def jonswap_spectrum(freq, Hs, Tp, gamma=3.3):
    """
    Compute JONSWAP spectrum according to IEC 61400-3
    
    Parameters:
    -----------
    freq : array
        Frequency array [Hz]
    Hs : float
        Significant wave height [m]
    Tp : float
        Peak period [s]
    gamma : float
        Peak enhancement factor (default=3.3)
    
    Returns:
    --------
    S : array
        Spectral density [m^2/Hz]
    """
    fp = 1.0 / Tp
    
    # Pierson-Moskowitz spectrum
    S_PM = 0.3125 * Hs**2 * fp**4 * freq**(-5) * np.exp(-1.25 * (fp / freq)**4)
    S_PM[freq == 0] = 0.0
    
    # JONSWAP enhancement
    sigma = np.where(freq <= fp, 0.07, 0.09)
    alpha = np.exp(-0.5 * ((freq - fp) / (sigma * fp))**2)
    gamma_alpha = gamma**alpha
    
    # Normalization factor
    norm_factor = 1.0 - 0.287 * np.log(gamma)
    
    S = norm_factor * S_PM * gamma_alpha
    
    return S


def generate_frequency_wave(filename):
    """
    Generate frequency domain wave data file
    
    File format:
    - Line 1: Number of components
    - Following lines: frequency[Hz] heading[deg] amplitude[m] phase[rad]
    """
    
    # Define frequency range
    f_min = 0.02  # 50s period
    f_max = 0.30  # 3.33s period
    df = 0.01     # Frequency step
    freqs = np.arange(f_min, f_max + df, df)
    
    # System 1: Main sea state
    Hs1 = 1.5
    Tp1 = 10.0
    heading1 = 0.0  # degrees
    S1 = jonswap_spectrum(freqs, Hs1, Tp1)
    
    # System 2: Secondary sea state (swell from different direction)
    Hs2 = 0.8
    Tp2 = 7.0
    heading2 = 30.0  # degrees
    S2 = jonswap_spectrum(freqs, Hs2, Tp2)
    
    # Convert spectral density to amplitudes: A = sqrt(2 * S * df)
    amp1 = np.sqrt(2.0 * S1 * df)
    amp2 = np.sqrt(2.0 * S2 * df)
    
    # Generate random phases
    np.random.seed(42)  # For reproducibility
    phases1 = np.random.uniform(-np.pi, np.pi, len(freqs))
    phases2 = np.random.uniform(-np.pi, np.pi, len(freqs))
    
    # Create list of components (frequency, heading, amplitude, phase)
    components = []
    
    # Add components from system 1
    for i, f in enumerate(freqs):
        if amp1[i] > 0.001:  # Only include significant components
            components.append([f, heading1, amp1[i], phases1[i]])
    
    # Add components from system 2
    for i, f in enumerate(freqs):
        if amp2[i] > 0.001:  # Only include significant components
            components.append([f, heading2, amp2[i], phases2[i]])
    
    # Sort by frequency then heading
    components = sorted(components, key=lambda x: (x[0], x[1]))
    
    # Write to file
    with open(filename, 'w') as f:
        f.write(f"{len(components)}  // Number of frequency components\n")
        for comp in components:
            f.write(f"{comp[0]:.6f}  {comp[1]:.3f}  {comp[2]:.6f}  {comp[3]:.6f}\n")
    
    print(f"Generated {len(components)} frequency components")
    print(f"System 1: Hs={Hs1}m, Tp={Tp1}s, heading={heading1}°")
    print(f"System 2: Hs={Hs2}m, Tp={Tp2}s, heading={heading2}°")
    print(f"Total Hs (RSS): {np.sqrt(Hs1**2 + Hs2**2):.2f}m")
    
    # Create visualization
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Plot spectra
    ax1.plot(freqs, S1, 'b-', linewidth=2, label=f'System 1 (Tp={Tp1}s, θ={heading1}°)')
    ax1.plot(freqs, S2, 'r-', linewidth=2, label=f'System 2 (Tp={Tp2}s, θ={heading2}°)')
    ax1.plot(freqs, S1 + S2, 'k--', linewidth=1.5, label='Combined')
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('Spectral Density [m²/Hz]')
    ax1.set_title('Wave Spectrum - Frequency Domain')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_xlim([0, 0.3])
    
    # Plot amplitudes
    ax2.plot(freqs, amp1, 'b-', linewidth=2, label=f'System 1 (Tp={Tp1}s, θ={heading1}°)')
    ax2.plot(freqs, amp2, 'r-', linewidth=2, label=f'System 2 (Tp={Tp2}s, θ={heading2}°)')
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Amplitude [m]')
    ax2.set_title('Wave Amplitude - Frequency Domain')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_xlim([0, 0.3])
    
    plt.tight_layout()
    plt.savefig(filename.replace('.dat', '_spectrum.png'), dpi=150)
    print(f"Saved spectrum plot to {filename.replace('.dat', '_spectrum.png')}")
    
    return components


if __name__ == "__main__":
    output_file = "wavefreq.dat"
    components = generate_frequency_wave(output_file)
    print(f"\nFile '{output_file}' created successfully!")
    print(f"Total components: {len(components)}")
