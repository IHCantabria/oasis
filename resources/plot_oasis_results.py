"""
OASIS Results Plotting Tool

Python script to visualize OASIS simulation results including:
- Body motions (6 DOF time series)
- Mooring line tensions
- Wind turbine data
- 3D animated visualization

Replaces MATLAB scripts: plotBox.m and plot_case.m
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import os
import sys
from matplotlib.animation import FuncAnimation, FFMpegWriter


class OASISResultsPlotter:
    """Class to handle OASIS simulation results plotting"""
    
    def __init__(self, output_path):
        """
        Initialize plotter with output directory path
        
        Parameters:
        -----------
        output_path : str
            Path to OASIS output directory
        """
        self.output_path = output_path
        self.case_name = os.path.basename(os.path.dirname(output_path.rstrip('\\/')))
        
        # Configuration flags
        self.flag_plot = True
        self.flag_video = False
        self.flag_save = False
        self.flag_save_video = False
        
        # DOF configuration
        self.dofs_plot = [1, 2, 3, 4, 5, 6]
        
        # Lines configuration
        self.flag_lines = False
        self.n_lines = 4
        self.lines_plot = [1, 2, 3, 4]
        
        # Wind turbine configuration
        self.flag_wt = False
        self.n_wt = 1
        
        # Body dimensions and properties (ISOBARA-like platform)
        self.body_dim = [20.0, 4.98, 2.25]  # [length, width, height] in meters
        self.body_color = [0.8, 0.8, 0.8]
        self.z_cog = 0.465  # Center of gravity Z coordinate
        self.z0 = 0.105      # Base offset
        self.water_depth = 10.0
        self.pos0 = [0, 0, self.z_cog, 0, 0, 0]  # Initial position
        
        # Wave parameters
        self.wave_height = [6]
        self.wave_period = [10]
        self.wave_phase = [0]
        self.wave_direction = 0  # degrees
        self.ramp_time = 10
        
        # Video/animation settings
        self.t_ini = 0
        self.t_fin = 100
        self.dt_video = 0.1
        self.speed = 10
        self.view_angle = [45, 10]  # [azimuth, elevation]
        
        # Data containers
        self.dof_data = None
        self.time = None
        self.line_data = None
        self.wt_data = None
        
    def load_dof_data(self):
        """Load DOF (Degree of Freedom) time series data"""
        print("Loading DOF data...")
        self.dof_data = {}
        
        for i in range(1, 7):
            filename = os.path.join(self.output_path, f'DOF_{i}_Body_0.txt')
            if os.path.exists(filename):
                data = np.loadtxt(filename)
                self.dof_data[i] = data
            else:
                print(f"  Warning: {filename} not found")
                
        if self.dof_data:
            self.time = self.dof_data[1][:, 0]
            print(f"  Loaded {len(self.time)} time steps")
            print(f"  Time range: {self.time[0]:.2f} to {self.time[-1]:.2f} s")
    
    def load_line_data(self):
        """Load mooring line data"""
        if not self.flag_lines:
            return
            
        print("Loading line data...")
        self.line_data = {}
        self.line_positions = {}
        
        for i in range(self.n_lines):
            # Load tension data
            filename = os.path.join(self.output_path, f'EndsTen_{i}.txt')
            if os.path.exists(filename):
                self.line_data[i] = np.loadtxt(filename)
                print(f"  Loaded line {i+1} tension data")
            
            # Load position data for animation
            if self.flag_video:
                pos_files = {
                    'x': os.path.join(self.output_path, f'NodePosX_{i}.txt'),
                    'y': os.path.join(self.output_path, f'NodePosY_{i}.txt'),
                    'z': os.path.join(self.output_path, f'NodePosZ_{i}.txt')
                }
                if all(os.path.exists(f) for f in pos_files.values()):
                    self.line_positions[i] = {
                        'x': np.loadtxt(pos_files['x']),
                        'y': np.loadtxt(pos_files['y']),
                        'z': np.loadtxt(pos_files['z'])
                    }
                    print(f"  Loaded line {i+1} position data")
    
    def plot_dof_timeseries(self):
        """Plot 6 DOF time series"""
        if not self.flag_plot or not self.dof_data:
            return
            
        print("Plotting DOF time series...")
        
        fig, axes = plt.subplots(6, 1, figsize=(10, 12))
        fig.suptitle(f'Body Motions - {self.case_name}', fontsize=14, fontweight='bold')
        
        dof_labels = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']
        dof_units = ['m', 'm', 'm', '°', '°', '°']
        
        for i in range(6):
            ax = axes[i]
            dof_idx = i + 1
            
            if dof_idx in self.dof_data:
                time = self.dof_data[dof_idx][:, 0]
                value = self.dof_data[dof_idx][:, 1] - self.pos0[i]
                
                # Convert rotation DOFs to degrees
                if dof_idx > 3:
                    value = np.rad2deg(value)
                
                ax.plot(time, value, 'b-', linewidth=1.5)
                ax.set_ylabel(f'{dof_labels[i]} [{dof_units[i]}]', fontsize=10)
                ax.grid(True, alpha=0.3)
                ax.set_xlim([time[0], time[-1]])
                
                # Set reasonable y-limits
                if dof_idx <= 3:
                    ylim = ax.get_ylim()
                    ax.set_ylim([min(ylim[0], -0.1), max(ylim[1], 0.1)])
                else:
                    ylim = ax.get_ylim()
                    ax.set_ylim([min(ylim[0], -0.5), max(ylim[1], 0.5)])
        
        axes[-1].set_xlabel('Time [s]', fontsize=10)
        plt.tight_layout()
        
        if self.flag_save:
            save_dir = os.path.join(os.path.dirname(self.output_path), self.case_name)
            os.makedirs(save_dir, exist_ok=True)
            plt.savefig(os.path.join(save_dir, f'{self.case_name}_movements.png'), dpi=150)
            print(f"  Saved to {save_dir}")
        
        plt.show()
    
    def plot_line_tensions(self):
        """Plot mooring line tensions"""
        if not self.flag_lines or not self.line_data:
            return
            
        print("Plotting line tensions...")
        
        n_lines_plot = len(self.lines_plot)
        fig, axes = plt.subplots(n_lines_plot, 1, figsize=(10, 2*n_lines_plot))
        if n_lines_plot == 1:
            axes = [axes]
        
        fig.suptitle(f'Mooring Line Tensions - {self.case_name}', fontsize=14, fontweight='bold')
        
        for k, line_idx in enumerate(self.lines_plot):
            i = line_idx - 1  # Convert to 0-based index
            
            if i in self.line_data:
                ax = axes[k]
                time_line = self.line_data[i][:, 0]
                
                # Compute tension magnitude and convert to tons
                tension = np.sqrt(
                    self.line_data[i][:, 4]**2 + 
                    self.line_data[i][:, 5]**2 + 
                    self.line_data[i][:, 6]**2
                ) / 9819.0  # Convert to tons
                
                ax.plot(time_line, tension, 'k-', linewidth=1.5)
                ax.set_ylabel(f'Line {line_idx} [ton]', fontsize=10)
                ax.grid(True, alpha=0.3)
                ax.set_xlim([time_line[0], time_line[-1]])
        
        axes[-1].set_xlabel('Time [s]', fontsize=10)
        plt.tight_layout()
        
        if self.flag_save:
            save_dir = os.path.join(os.path.dirname(self.output_path), self.case_name)
            os.makedirs(save_dir, exist_ok=True)
            plt.savefig(os.path.join(save_dir, f'{self.case_name}_tensions.png'), dpi=150)
            print(f"  Saved to {save_dir}")
        
        plt.show()
    
    @staticmethod
    def rotation_matrix(roll, pitch, yaw):
        """
        Compute 3D rotation matrix from Euler angles
        
        Parameters:
        -----------
        roll, pitch, yaw : float
            Rotation angles in degrees
            
        Returns:
        --------
        M : ndarray (3x3)
            Rotation matrix
        """
        roll_rad = np.deg2rad(roll)
        pitch_rad = np.deg2rad(pitch)
        yaw_rad = np.deg2rad(yaw)
        
        cr, sr = np.cos(roll_rad), np.sin(roll_rad)
        cp, sp = np.cos(pitch_rad), np.sin(pitch_rad)
        cy, sy = np.cos(yaw_rad), np.sin(yaw_rad)
        
        M = np.array([
            [cp*cy, cy*sp*sr - cr*sy, sr*sy + cr*cy*sp],
            [cp*sy, cr*cy + sp*sr*sy, cr*sp*sy - cy*sr],
            [-sp,   cp*sr,             cp*cr]
        ])
        
        return M
    
    def create_box_vertices(self, dim, pos, z0=0):
        """
        Create vertices for a 3D box
        
        Parameters:
        -----------
        dim : list [length, width, height]
            Box dimensions
        pos : list [x, y, z, roll, pitch, yaw]
            Position and orientation
        z0 : float
            Base offset
            
        Returns:
        --------
        vertices : ndarray
            Box vertex coordinates
        """
        a, b, c = dim
        
        # Define box vertices (8 corners)
        half_vertices = np.array([
            [a/2, b/2, c/2],   # 0: top-front-right
            [-a/2, b/2, c/2],  # 1: top-front-left
            [-a/2, -b/2, c/2], # 2: top-back-left
            [a/2, -b/2, c/2],  # 3: top-back-right
            [a/2, b/2, -c/2],  # 4: bottom-front-right
            [-a/2, b/2, -c/2], # 5: bottom-front-left
            [-a/2, -b/2, -c/2],# 6: bottom-back-left
            [a/2, -b/2, -c/2]  # 7: bottom-back-right
        ])
        
        # Apply z0 offset
        half_vertices[:, 2] += z0
        
        # Apply rotation
        M_rot = self.rotation_matrix(pos[3], pos[4], pos[5])
        vertices = half_vertices @ M_rot.T
        
        # Apply translation
        vertices += pos[:3]
        
        return vertices
    
    def create_box_faces(self, vertices):
        """
        Create face definitions for a box from vertices
        
        Parameters:
        -----------
        vertices : ndarray (8x3)
            Box vertices
            
        Returns:
        --------
        faces : list of ndarrays
            List of face vertex coordinates
        """
        # Define the 6 faces of the box (each face has 4 vertices)
        face_indices = [
            [0, 1, 2, 3],  # top
            [4, 5, 6, 7],  # bottom
            [0, 3, 7, 4],  # right
            [1, 2, 6, 5],  # left
            [0, 1, 5, 4],  # front
            [2, 3, 7, 6]   # back
        ]
        
        faces = [vertices[indices] for indices in face_indices]
        
        return faces
    
    def compute_wave_surface(self, X, Y, t):
        """
        Compute wave surface elevation
        
        Parameters:
        -----------
        X, Y : ndarray
            Grid coordinates
        t : float
            Time
            
        Returns:
        --------
        Z : ndarray
            Wave elevation
        """
        Z = np.zeros_like(X)
        ramp = min(1.0, t / self.ramp_time)
        
        for H, T, P in zip(self.wave_height, self.wave_period, self.wave_phase):
            lambda_wave = 9.81 * T**2 / (2 * np.pi)
            k = 2 * np.pi / lambda_wave
            omega = 2 * np.pi / T
            
            dir_rad = np.deg2rad(self.wave_direction)
            kx = k * np.cos(dir_rad)
            ky = k * np.sin(dir_rad)
            
            Z += 0.5 * H * np.cos(omega * t - kx * X - ky * Y + P)
        
        return Z * ramp
    
    def create_animation(self):
        """Create 3D animated visualization"""
        if not self.flag_video or not self.dof_data:
            return
            
        print("Creating animation...")
        print("  This may take a while...")
        
        # Set up figure
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')
        
        # Create seabed
        X_bed, Y_bed = np.meshgrid(np.arange(-120, 121, 40), np.arange(-120, 121, 40))
        Z_bed = -self.water_depth * np.ones_like(X_bed)
        ax.plot_surface(X_bed * 10, Y_bed * 10, Z_bed, color='tan', 
                       alpha=0.8, edgecolors='none', shade=True)
        
        # Create wave grid
        X_wave, Y_wave = np.meshgrid(np.arange(-30, 31, 5), np.arange(-30, 31, 5))
        
        # Find time indices
        i1 = np.argmax(self.time >= self.t_ini)
        i2 = np.argmax(self.time > self.t_fin) if self.t_fin < self.time[-1] else len(self.time) - 1
        di = max(1, int(self.dt_video / (self.time[1] - self.time[0])))
        
        time_indices = range(i1, i2, di)
        
        def update(frame_idx):
            """Update function for animation"""
            ax.cla()  # Clear axes
            
            it = time_indices[frame_idx]
            t = self.time[it]
            
            # Seabed
            ax.plot_surface(X_bed * 10, Y_bed * 10, Z_bed, color='tan', 
                           alpha=0.8, edgecolors='none', shade=True)
            
            # Mooring lines
            if self.flag_lines and self.line_positions:
                for i in range(self.n_lines):
                    if i in self.line_positions:
                        x = self.line_positions[i]['x'][it, 1:]
                        y = self.line_positions[i]['y'][it, 1:]
                        z = self.line_positions[i]['z'][it, 1:]
                        ax.plot(x, y, z, 'kx-', linewidth=1.5, markersize=3)
            
            # Body box
            pos = [
                self.dof_data[1][it, 1],
                self.dof_data[2][it, 1],
                self.dof_data[3][it, 1],
                np.rad2deg(self.dof_data[4][it, 1]),
                np.rad2deg(self.dof_data[5][it, 1]),
                np.rad2deg(self.dof_data[6][it, 1])
            ]
            
            vertices = self.create_box_vertices(self.body_dim, pos, self.z0)
            faces = self.create_box_faces(vertices)
            
            box = Poly3DCollection(faces, alpha=1.0, facecolor=self.body_color, 
                                  edgecolor='k', linewidth=0.5)
            ax.add_collection3d(box)
            
            # Wave surface
            Z_wave = self.compute_wave_surface(X_wave, Y_wave, t)
            ax.plot_surface(X_wave, Y_wave, Z_wave, color='skyblue', 
                           alpha=0.6, edgecolors='none', shade=True)
            
            # Set view and limits
            ax.set_xlabel('X [m]')
            ax.set_ylabel('Y [m]')
            ax.set_zlabel('Z [m]')
            ax.set_xlim([-30, 30])
            ax.set_ylim([-30, 30])
            ax.set_zlim([-30, 10])
            ax.view_init(elev=self.view_angle[1], azim=self.view_angle[0])
            
            ax.set_box_aspect([1, 1, 2/3])
            ax.set_title(f'Time = {t:.1f} s', fontsize=12, fontweight='bold')
            
            return ax,
        
        # Create animation
        n_frames = len(time_indices)
        print(f"  Generating {n_frames} frames...")
        
        anim = FuncAnimation(fig, update, frames=n_frames, 
                            interval=int(1000 * self.dt_video / self.speed),
                            blit=False, repeat=True)
        
        if self.flag_save_video:
            save_dir = os.path.join(os.path.dirname(self.output_path), self.case_name)
            os.makedirs(save_dir, exist_ok=True)
            video_file = os.path.join(save_dir, f'{self.case_name}_video.mp4')
            
            print(f"  Saving video to {video_file}")
            writer = FFMpegWriter(fps=int(1.0 / (self.dt_video / self.speed)), 
                                 bitrate=1800)
            anim.save(video_file, writer=writer)
            print(f"  Video saved!")
        else:
            plt.show()
    
    def run(self):
        """Run all plotting routines"""
        print("=" * 70)
        print("OASIS Results Plotting Tool")
        print("=" * 70)
        print(f"Case: {self.case_name}")
        print(f"Output path: {self.output_path}")
        print()
        
        # Load data
        self.load_dof_data()
        self.load_line_data()
        
        # Generate plots
        if self.flag_plot:
            self.plot_dof_timeseries()
            self.plot_line_tensions()
        
        if self.flag_video:
            self.create_animation()
        
        print()
        print("=" * 70)
        print("Plotting complete!")
        print("=" * 70)


def main():
    """Main function with command line interface"""
    
    # Default configuration
    if len(sys.argv) > 1:
        output_path = sys.argv[1]
    else:
        # Default to test case
        output_path = os.path.join(os.path.dirname(__file__), 
                                   '..', 'examples', 'test', 'output')
        output_path = os.path.abspath(output_path)
    
    if not os.path.exists(output_path):
        print(f"Error: Output directory not found: {output_path}")
        print(f"\nUsage: python plot_oasis_results.py [output_directory]")
        sys.exit(1)
    
    # Create plotter instance
    plotter = OASISResultsPlotter(output_path)
    
    # Configure options (can be modified by user)
    plotter.flag_plot = True   # Plot time series
    plotter.flag_video = False  # Create 3D animation
    plotter.flag_save = False   # Save figures
    plotter.flag_save_video = False  # Save video file
    
    plotter.flag_lines = False  # Plot mooring lines
    plotter.flag_wt = False     # Plot wind turbine data
    
    # Video settings (if enabled)
    plotter.t_ini = 0
    plotter.t_fin = 100
    plotter.dt_video = 0.1
    plotter.speed = 10
    plotter.view_angle = [45, 10]
    
    # Run plotting
    plotter.run()


if __name__ == "__main__":
    main()
