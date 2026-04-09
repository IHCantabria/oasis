"""
Generate time-domain wave data for OASIS waves/timeseries example.

This script creates a JONSWAP-based irregular wave time series and writes it
in the format expected by OASIS (specType=2):
    Line 1:  N  0          (number of time points, unused flag)
    Lines 2+: time  elevation

Parameters:
    Hs    = 1.5 m   (significant wave height)
    Tp    = 10.0 s   (peak period)
    gamma = 3.3      (JONSWAP peak enhancement)
    dt    = 0.1 s    (time step)
    T     = 1000.0 s (total duration)
    seed  = 42       (random seed for reproducibility)

Usage:
    python generate_wave_timeseries.py [output_file]
    Default output: wave.dat
"""

import sys
import numpy as np


def jonswap_spectrum(freq, Hs, Tp, gamma=3.3):
    """
    Compute JONSWAP spectrum according to IEC 61400-3.

    Parameters
    ----------
    freq : array
        Frequency array [Hz] (must not contain 0).
    Hs : float
        Significant wave height [m].
    Tp : float
        Peak period [s].
    gamma : float
        Peak enhancement factor.

    Returns
    -------
    S : array
        Spectral density [m^2/Hz].
    """
    fp = 1.0 / Tp

    # Pierson-Moskowitz base spectrum
    S_PM = 0.3125 * Hs**2 * fp**4 * freq**(-5) * np.exp(-1.25 * (fp / freq)**4)

    # JONSWAP enhancement
    sigma = np.where(freq <= fp, 0.07, 0.09)
    alpha = np.exp(-0.5 * ((freq - fp) / (sigma * fp))**2)
    gamma_alpha = gamma**alpha

    # Normalization
    norm_factor = 1.0 - 0.287 * np.log(gamma)

    return norm_factor * S_PM * gamma_alpha


def generate_wave_timeseries(filename, Hs=1.5, Tp=10.0, gamma=3.3,
                             duration=1000.0, dt=0.1, seed=42):
    """
    Generate a JONSWAP irregular wave time series and write to file.

    Parameters
    ----------
    filename : str
        Output file path.
    Hs : float
        Significant wave height [m].
    Tp : float
        Peak period [s].
    gamma : float
        JONSWAP peak enhancement factor.
    duration : float
        Total duration [s].
    dt : float
        Time step [s].
    seed : int
        Random seed for reproducibility.
    """
    np.random.seed(seed)

    # Frequency discretisation
    f_min = 0.02   # 50 s period
    f_max = 0.50   # 2 s period
    df = 0.005     # fine resolution for smooth time series
    freqs = np.arange(f_min, f_max + df / 2, df)

    # Spectrum and amplitudes
    S = jonswap_spectrum(freqs, Hs, Tp, gamma)
    amplitudes = np.sqrt(2.0 * S * df)

    # Random phases
    phases = np.random.uniform(-np.pi, np.pi, len(freqs))

    # Time vector
    times = np.arange(0.0, duration + dt / 2, dt)
    npoints = len(times)

    # Superpose sinusoidal components
    omega = 2.0 * np.pi * freqs  # angular frequencies
    elevation = np.zeros(npoints)
    for i, t in enumerate(times):
        elevation[i] = np.sum(amplitudes * np.cos(omega * t + phases))

    # Write file
    with open(filename, 'w') as f:
        f.write(f"{npoints} 0\n")
        for i in range(npoints):
            f.write(f"{times[i]:.1f} {elevation[i]:.15E}\n")

    # Print statistics
    Hs_actual = 4.0 * np.sqrt(np.var(elevation))
    print(f"Generated wave time series: {filename}")
    print(f"  Points:  {npoints}")
    print(f"  Duration: {duration} s  (dt = {dt} s)")
    print(f"  Target Hs: {Hs} m")
    print(f"  Actual Hs (4*std): {Hs_actual:.3f} m")
    print(f"  Tp: {Tp} s, gamma: {gamma}")


if __name__ == "__main__":
    output_file = sys.argv[1] if len(sys.argv) > 1 else "wave.dat"
    generate_wave_timeseries(output_file)
    print(f"\nFile '{output_file}' created successfully!")
