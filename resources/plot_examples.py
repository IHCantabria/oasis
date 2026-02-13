"""
Example usage of OASIS plotting tools

This script demonstrates common usage patterns for visualizing OASIS results.
"""

import os
import sys

# Add resources directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

from plot_oasis_results import OASISResultsPlotter


def example_basic_plots():
    """Example 1: Basic DOF time series plots"""
    print("\n" + "="*70)
    print("Example 1: Basic Time Series Plots")
    print("="*70)
    
    output_path = os.path.join(script_dir, '..', 'examples', 'generic_example', 'output')
    
    plotter = OASISResultsPlotter(output_path)
    plotter.flag_plot = True
    plotter.flag_video = False
    plotter.flag_save = False
    
    plotter.run()


def example_with_lines():
    """Example 2: Include mooring line tensions"""
    print("\n" + "="*70)
    print("Example 2: With Mooring Lines")
    print("="*70)
    
    output_path = os.path.join(script_dir, '..', 'examples', 'generic_example', 'output')
    
    plotter = OASISResultsPlotter(output_path)
    plotter.flag_plot = True
    plotter.flag_lines = True
    plotter.n_lines = 4
    plotter.lines_plot = [1, 2, 3, 4]
    plotter.flag_save = False
    
    plotter.run()


def example_animation():
    """Example 3: Create 3D animation"""
    print("\n" + "="*70)
    print("Example 3: 3D Animation")
    print("="*70)
    print("Note: This may take a while. Close the plot window to continue.")
    
    output_path = os.path.join(script_dir, '..', 'examples', 'generic_example', 'output')
    
    plotter = OASISResultsPlotter(output_path)
    plotter.flag_plot = False
    plotter.flag_video = True
    plotter.flag_save_video = False  # Set to True to save MP4 (requires ffmpeg)
    
    # Animation settings
    plotter.t_ini = 0
    plotter.t_fin = 20  # First 20 seconds
    plotter.dt_video = 0.2  # Frame every 0.2 seconds
    plotter.speed = 5  # Play at 5x speed
    plotter.view_angle = [45, 15]  # Azimuth, elevation
    
    # Wave parameters (adjust to match your simulation)
    plotter.wave_height = [2.0]
    plotter.wave_period = [10.0]
    plotter.wave_phase = [0]
    plotter.wave_direction = 0
    plotter.ramp_time = 10
    
    plotter.run()


def example_save_all():
    """Example 4: Generate and save all plots"""
    print("\n" + "="*70)
    print("Example 4: Save All Plots")
    print("="*70)
    
    output_path = os.path.join(script_dir, '..', 'examples', 'generic_example', 'output')
    
    plotter = OASISResultsPlotter(output_path)
    plotter.flag_plot = True
    plotter.flag_lines = False
    plotter.flag_save = True
    
    plotter.run()
    
    print("\nPlots saved to:", os.path.join(os.path.dirname(output_path), plotter.case_name))


def example_custom_configuration():
    """Example 5: Custom configuration for specific platform"""
    print("\n" + "="*70)
    print("Example 5: Custom Platform Configuration")
    print("="*70)
    
    output_path = os.path.join(script_dir, '..', 'examples', 'freq_wave_example', 'output')
    
    plotter = OASISResultsPlotter(output_path)
    
    # Custom body dimensions (semi-submersible)
    plotter.body_dim = [30.0, 30.0, 10.0]  # Large platform
    plotter.body_color = [0.3, 0.5, 0.8]   # Blue-gray color
    plotter.z_cog = 5.0
    plotter.z0 = 0.0
    
    # Custom wave parameters
    plotter.wave_height = [1.5, 0.8]  # Two wave systems
    plotter.wave_period = [10.0, 7.0]
    plotter.wave_phase = [0, 0]
    plotter.wave_direction = 0
    
    plotter.flag_plot = True
    plotter.flag_save = False
    
    plotter.run()


def example_batch_processing():
    """Example 6: Batch process multiple cases"""
    print("\n" + "="*70)
    print("Example 6: Batch Processing")
    print("="*70)
    
    import glob
    
    # Find all output directories in examples
    examples_dir = os.path.join(script_dir, '..', 'examples')
    output_dirs = glob.glob(os.path.join(examples_dir, '*/output'))
    
    print(f"Found {len(output_dirs)} cases to process:")
    
    for output_path in output_dirs:
        case_name = os.path.basename(os.path.dirname(output_path))
        print(f"\n  Processing case: {case_name}")
        
        # Check if DOF files exist
        dof_file = os.path.join(output_path, 'DOF_1_Body_0.txt')
        if not os.path.exists(dof_file):
            print(f"    Skipping (no DOF data found)")
            continue
        
        try:
            plotter = OASISResultsPlotter(output_path)
            plotter.flag_plot = True
            plotter.flag_save = True
            plotter.flag_video = False
            plotter.run()
            print(f"    ✓ Success")
        except Exception as e:
            print(f"    ✗ Error: {e}")


def main():
    """Main menu for examples"""
    examples = {
        '1': ('Basic time series plots', example_basic_plots),
        '2': ('With mooring line tensions', example_with_lines),
        '3': ('3D animation', example_animation),
        '4': ('Save all plots to disk', example_save_all),
        '5': ('Custom platform configuration', example_custom_configuration),
        '6': ('Batch process multiple cases', example_batch_processing),
    }
    
    print("\n" + "="*70)
    print("OASIS Plotting Examples")
    print("="*70)
    print("\nAvailable examples:")
    for key, (desc, _) in examples.items():
        print(f"  {key}. {desc}")
    print("  q. Quit")
    
    while True:
        choice = input("\nSelect example (1-6, q to quit): ").strip().lower()
        
        if choice == 'q':
            print("Goodbye!")
            break
        elif choice in examples:
            _, func = examples[choice]
            try:
                func()
            except Exception as e:
                print(f"\nError running example: {e}")
                import traceback
                traceback.print_exc()
        else:
            print("Invalid choice. Please select 1-6 or q.")


if __name__ == "__main__":
    main()
