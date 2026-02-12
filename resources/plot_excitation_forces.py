import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os

def plot_surge_excitation():
    """
    Compare surge 2nd order excitation forces from OASIS simulation with theoretical values.
    
    The WaveExcitationForce_Body_0.txt file format is:
    Time | 1st_order_DOF1 | 1st_order_DOF2 | 1st_order_DOF3 | 1st_order_DOF4 | 1st_order_DOF5 | 1st_order_DOF6 |
    2nd_order_DOF1 | 2nd_order_DOF2 | 2nd_order_DOF3 | ... | 2nd_order_DOF6
    
    DOF1 (column index 0 for 1st order, column index 7 for 2nd order) is surge (horizontal force).
    """
    
    # Define file paths relative to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    
    output_file = os.path.join(base_dir, 'examples', 'test_qtf', 'output', 'WaveExcitationForce_Body_0.txt')
    theoretical_file = os.path.join(base_dir, 'examples', 'test_qtf', 'input', '2ndorder_long.csv')
    
    # Check if files exist
    if not os.path.exists(output_file):
        print(f"ERROR: Output file not found: {output_file}")
        return
    
    if not os.path.exists(theoretical_file):
        print(f"ERROR: Theoretical file not found: {theoretical_file}")
        return
    
    # Read OASIS simulation output
    print("Reading OASIS simulation output...")
    try:
        oasis_data = pd.read_csv(output_file, delim_whitespace=True, header=None)
        time_oasis = oasis_data.iloc[:, 0].values.astype(np.float64)
        # 1st order DOF1 (surge) is at column index 1 + 0 = 1 (0-indexed)
        surge_1st_order_oasis = oasis_data.iloc[:, 1].values.astype(np.float64)
        # 2nd order DOF1 (surge) is at column index 7 + 0 = 7 (0-indexed)
        surge_2nd_order_oasis = oasis_data.iloc[:, 7].values.astype(np.float64)
        print(f"  Read {len(time_oasis)} time steps from OASIS output")
    except Exception as e:
        print(f"ERROR reading OASIS output: {e}")
        return
    
    # Read theoretical values
    print("Reading theoretical values...")
    try:
        # Skip header lines in the CSV file and convert to float
        theoretical_data = pd.read_csv(theoretical_file, skiprows=2, delimiter=';', header=None)
        # Convert columns to numeric, coercing errors to NaN
        theoretical_data.iloc[:, 0] = pd.to_numeric(theoretical_data.iloc[:, 0], errors='coerce')
        theoretical_data.iloc[:, 1] = pd.to_numeric(theoretical_data.iloc[:, 1], errors='coerce')
        # Drop rows with NaN values
        theoretical_data = theoretical_data.dropna()
        time_theory = theoretical_data.iloc[:, 0].values.astype(np.float64)
        surge_theory = theoretical_data.iloc[:, 1].values.astype(np.float64)
        print(f"  Read {len(time_theory)} time steps from theoretical data")
    except Exception as e:
        print(f"ERROR reading theoretical data: {e}")
        return
    
    # Trim data to common time range (shorter timeseries)
    t_max_oasis = float(time_oasis[-1])
    t_max_theory = float(time_theory[-1])
    t_max_common = min(t_max_oasis, t_max_theory)
    
    # Filter OASIS data
    mask_oasis = time_oasis <= t_max_common
    time_oasis = time_oasis[mask_oasis]
    surge_1st_order_oasis = surge_1st_order_oasis[mask_oasis]
    surge_2nd_order_oasis = surge_2nd_order_oasis[mask_oasis]
    surge_oasis = surge_2nd_order_oasis #+ surge_1st_order_oasis
    
    # Filter theoretical data
    mask_theory = time_theory <= t_max_common
    time_theory = time_theory[mask_theory]
    surge_theory = surge_theory[mask_theory]
    
    print(f"\nPlotting data up to t = {t_max_common:.2f} s")
    print(f"  OASIS data: {len(time_oasis)} points")
    print(f"  Theoretical data: {len(time_theory)} points")
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Plot both datasets
    ax.plot(time_oasis, surge_oasis, 'b-', linewidth=2, label='OASIS', alpha=0.8)
    ax.plot(time_theory, surge_theory, 'r--', linewidth=2, label='AQWA', alpha=0.8)
    
    # Formatting
    ax.set_xlabel('Time (s)', fontsize=12)
    ax.set_ylabel('Surge Force (N)', fontsize=12)
    ax.set_title('Surge Excitation Force Comparison', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11, loc='best')
    
    plt.tight_layout()
    
    # Save the figure
    output_plot = os.path.join(base_dir, 'examples', 'test_qtf', 'output', 'surge_force_comparison.png')
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_plot}")
    
    # Calculate and print statistics
    print("\n=== Comparison Statistics ===")
    print(f"OASIS data range: {surge_oasis.min():.2f} to {surge_oasis.max():.2f} N")
    print(f"Theoretical data range: {surge_theory.min():.2f} to {surge_theory.max():.2f} N")
    
    # Interpolate theoretical values to OASIS time points for comparison
    theoretical_interp = np.interp(time_oasis, time_theory, surge_theory)
    error = surge_oasis - theoretical_interp
    abs_error = np.abs(error)
    
    print(f"\nAbsolute error statistics (at OASIS time points):")
    print(f"  Mean absolute error: {np.mean(abs_error):.2f} N")
    print(f"  Max absolute error: {np.max(abs_error):.2f} N")
    print(f"  RMS error: {np.sqrt(np.mean(error**2)):.2f} N")
    
    # Show the plot
    plt.show()

if __name__ == "__main__":
    plot_surge_excitation()