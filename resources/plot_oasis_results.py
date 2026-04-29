"""
OASIS Results Plotting Tool

Python script to visualize OASIS simulation results including:
- Body motions (6 DOF time series)
- Mooring line tensions
- Wind turbine data
- 3D animated visualization

"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.tri import Triangulation
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401 — registers 3D projection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ---------------------------------------------------------------------------
# Configure these variables before running
# ---------------------------------------------------------------------------
CASE_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    '..', 'examples', 'lines', 'multi_line_yaml'
)
BODY_ID = 0  # Body index (0-based) used in CSV column names
# ---------------------------------------------------------------------------


class OASISResultsPlotter:
    """Class to handle OASIS simulation results plotting"""
    
    def __init__(self, case_path):
        """
        Initialize plotter with case directory path.

        Parameters:
        -----------
        case_path : str
            Path to OASIS case directory (parent of input/ and output/)
        """
        self.case_path   = case_path
        self.output_path = os.path.join(case_path, 'output')
        self.input_path  = os.path.join(case_path, 'input')
        self.body_id     = BODY_ID
        self.case_name   = os.path.basename(case_path.rstrip('\\/'))

        # Configuration flags
        self.flag_plot       = True
        self.flag_video      = False
        self.flag_save       = False
        self.flag_save_video = False

        # DOF configuration
        self.dofs_plot = [1, 2, 3, 4, 5, 6]

        # Lines configuration
        self.flag_lines = False
        self.n_lines    = 4
        self.lines_plot = [1, 2, 3, 4]

        # Wind turbine configuration
        self.flag_wt = False
        self.n_wt    = 1

        # Body dimensions and properties
        self.body_dim   = [10.0, 10.0, 4.0]  # [length, width, height] in metres
        self.body_color = [0.8, 0.8, 0.8]
        self.z_cog      = 0.0   # Centre-of-gravity Z offset
        self.z0         = 0.0   # Base offset for box rendering
        self.water_depth = 20.0  # Visual seabed depth for animation (m, positive)
        self.domain      = 50.0  # Domain half-extent for animation (m); updated by load_input_yaml
        self.pos0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # Initial position (YAML initial_position)

        # Wave parameters (overwritten by load_input_yaml)
        self.wave_height    = [0.0]
        self.wave_period    = [10.0]
        self.wave_phase     = [0.0]
        self.wave_direction = 0.0   # degrees
        self.ramp_time      = 0.0

        # Video/animation settings
        self.t_ini       = 0
        self.t_fin       = 100
        self.dt_video    = 0.1
        self.speed       = 10
        self.view_angle  = [45, 10]  # [azimuth, elevation]

        # Data containers
        self.dof_data       = None   # dict: dof_idx(1–6) → 1-D ndarray
        self.time           = None   # 1-D ndarray
        self.line_data      = None   # dict: line_idx(0-based) → DataFrame
        self.line_positions = None   # dict: line_idx(0-based) → DataFrame
        self.wt_data        = None
        self.bathymetry_xyz = None   # Nx3 ndarray or None
        
    def load_input_yaml(self):
        """Load case parameters from dataProblem.yaml (waves, body BCPs, seafloor)."""
        yaml_file = os.path.join(self.input_path, 'dataProblem.yaml')
        if not os.path.exists(yaml_file):
            print(f"  Warning: Input YAML not found at {yaml_file}")
            return

        print("Loading input YAML...")
        with open(yaml_file, 'r') as fh:
            data = yaml.safe_load(fh)

        # ── Wave parameters ────────────────────────────────────────────
        waves      = data.get('waves', {})
        wave_type  = waves.get('type', 'REG')
        H          = float(waves.get('height', 0.0))
        T          = float(waves.get('period', 10.0))
        self.wave_height    = [H]
        self.wave_period    = [T if T > 0 else 10.0]
        self.wave_direction = float(waves.get('heading', 0.0))
        self.ramp_time      = float(waves.get('ramp_time', 0.0))
        print(f"  Waves: type={wave_type}, H={H:.2f} m, T={self.wave_period[0]:.2f} s, "
              f"heading={self.wave_direction:.1f} deg")

        bcps      = data.get('bcps', {})
        body_bcps = bcps.get('body_bcps', [])

        # ── Body dimensions from body BCP spread ───────────────────────
        if body_bcps:
            bcp_xs = [p['position'][0] for p in body_bcps]
            bcp_ys = [p['position'][1] for p in body_bcps]
            bcp_zs = [p['position'][2] for p in body_bcps]
            x_range  = max(bcp_xs) - min(bcp_xs)
            y_range  = max(bcp_ys) - min(bcp_ys)
            z_depth  = abs(min(bcp_zs))
            self.body_dim = [
                max(x_range, 1.0),
                max(y_range, 1.0),
                max(2.0 * z_depth, 1.0),
            ]
            self.z0 = 0.0
            print(f"  Body dim (from body BCPs): {[round(v, 1) for v in self.body_dim]} m")

        # ── Initial position from body definition ──────────────────────
        bodies = data.get('bodies', [])
        if bodies:
            body  = bodies[min(self.body_id, len(bodies) - 1)]
            ip    = body.get('initial_position', [0.0] * 6)
            self.pos0 = [float(v) for v in ip]
            print(f"  Initial position: {[round(v, 3) for v in self.pos0]}")

        # ── Domain size from anchor positions ──────────────────────────
        anchors = bcps.get('anchors', [])
        if anchors:
            anchor_pos = [a['position'] for a in anchors]
            max_r      = max(np.sqrt(p[0]**2 + p[1]**2) for p in anchor_pos)
            max_depth  = max(abs(p[2]) for p in anchor_pos)
            self.domain      = max_r * 1.1
            self.water_depth = max_depth
            print(f"  Domain: radius={self.domain:.1f} m, water depth={self.water_depth:.1f} m")

        # ── Seafloor: flat or bathymetry mesh ─────────────────────────
        seafloor          = data.get('seafloor', {})
        bathymetry_entries = seafloor.get('bathymetry', [])
        if bathymetry_entries:
            all_pts = []
            for entry in bathymetry_entries:
                path = os.path.join(self.input_path, entry.get('mesh_file', ''))
                pts  = self._load_bathymetry_dat(path)
                if pts is not None:
                    all_pts.append(pts)
            if all_pts:
                self.bathymetry_xyz = np.vstack(all_pts)
                print(f"  Bathymetry: {len(self.bathymetry_xyz)} points loaded")
        else:
            flat_entries = seafloor.get('flat', [])
            if flat_entries and not anchors:
                self.water_depth = abs(float(flat_entries[0].get('depth', self.water_depth)))
                print(f"  Seafloor: flat at {-self.water_depth:.1f} m")

    @staticmethod
    def _load_bathymetry_dat(filepath):
        """
        Load an OASIS bathymetry point cloud from a .dat file.

        Expected format::

            N // Num Puntos
            ////... (separator lines, skipped)
            x1  y1  z1
            ...

        Returns Nx3 float ndarray, or None on failure.
        """
        if not os.path.exists(filepath):
            print(f"  Warning: bathymetry file not found: {filepath}")
            return None
        try:
            with open(filepath, 'r') as fh:
                raw = fh.readlines()
            n_pts = int(raw[0].split()[0])
            pts   = []
            for line in raw[1:]:
                stripped = line.strip()
                if not stripped or stripped.startswith('/'):
                    continue
                vals = stripped.split()
                if len(vals) >= 3:
                    pts.append([float(vals[0]), float(vals[1]), float(vals[2])])
                if len(pts) >= n_pts:
                    break
            return np.array(pts) if pts else None
        except Exception as exc:
            print(f"  Warning: could not load {filepath}: {exc}")
            return None

    def load_dof_data(self):
        """Load body motion time series from body_{id}_motion.csv."""
        print("Loading DOF data...")
        filename = os.path.join(self.output_path, f'body_{self.body_id}_motion.csv')
        if not os.path.exists(filename):
            print(f"  Warning: {filename} not found")
            return

        df         = pd.read_csv(filename)
        self.time  = df['time'].values
        b          = self.body_id
        col_names  = [
            f'b{b}_x', f'b{b}_y', f'b{b}_z',
            f'b{b}_rl', f'b{b}_pt', f'b{b}_yw',
        ]
        self.dof_data = {}
        for i, col in enumerate(col_names, start=1):
            if col in df.columns:
                self.dof_data[i] = df[col].values
            else:
                print(f"  Warning: column '{col}' not found in {filename}")

        print(f"  Loaded {len(self.time)} time steps")
        print(f"  Time range: {self.time[0]:.2f} to {self.time[-1]:.2f} s")
    
    def load_line_data(self):
        """Load mooring line data from line_{i}_tensions.csv and line_{i}_positions.csv."""
        if not self.flag_lines:
            return

        print("Loading line data...")
        self.line_data      = {}
        self.line_positions = {}

        for line_idx in self.lines_plot:
            i = line_idx - 1  # 0-based index

            # Load tension data
            ten_file = os.path.join(self.output_path, f'line_{i}_tensions.csv')
            if os.path.exists(ten_file):
                self.line_data[i] = pd.read_csv(ten_file)
                print(f"  Loaded line {line_idx} tension data")
            else:
                print(f"  Warning: {ten_file} not found")

            # Load position data for animation
            if self.flag_video:
                pos_file = os.path.join(self.output_path, f'line_{i}_positions.csv')
                if os.path.exists(pos_file):
                    self.line_positions[i] = pd.read_csv(pos_file)
                    print(f"  Loaded line {line_idx} position data")
    
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
                time  = self.time
                value = self.dof_data[dof_idx] - self.pos0[i]
                
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
            save_dir = os.path.join(self.case_path, 'results')
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
                ax        = axes[k]
                df        = self.line_data[i]
                time_line = df['time'].values

                # Compute tension magnitude: look for named columns, then fall back
                if 'tension_fairlead' in df.columns:
                    tension = df['tension_fairlead'].values / 9819.0
                else:
                    force_cols = [c for c in df.columns
                                  if c != 'time' and 'anchor' not in c.lower()]
                    if len(force_cols) >= 3:
                        tension = np.sqrt(
                            df[force_cols[0]].values**2 +
                            df[force_cols[1]].values**2 +
                            df[force_cols[2]].values**2
                        ) / 9819.0
                    else:
                        tension = np.zeros(len(time_line))
                
                ax.plot(time_line, tension, 'k-', linewidth=1.5)
                ax.set_ylabel(f'Line {line_idx} [ton]', fontsize=10)
                ax.grid(True, alpha=0.3)
                ax.set_xlim([time_line[0], time_line[-1]])
        
        axes[-1].set_xlabel('Time [s]', fontsize=10)
        plt.tight_layout()
        
        if self.flag_save:
            save_dir = os.path.join(self.case_path, 'results')
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
        Z = np.zeros_like(X, dtype=float)
        ramp = 1.0 if self.ramp_time <= 0 else min(1.0, t / self.ramp_time)
        
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

        fig = plt.figure(figsize=(12, 9))
        ax  = fig.add_subplot(111, projection='3d')

        d = self.domain

        # ── Pre-compute seabed geometry (reused every frame) ───────────
        if self.bathymetry_xyz is not None:
            _bx = self.bathymetry_xyz[:, 0]
            _by = self.bathymetry_xyz[:, 1]
            _bz = self.bathymetry_xyz[:, 2]
            _bathy_tri = Triangulation(_bx, _by)
        else:
            _Xb, _Yb = np.meshgrid(np.linspace(-d, d, 7), np.linspace(-d, d, 7))
            _Zb      = -self.water_depth * np.ones_like(_Xb, dtype=float)

        # ── Wave grid – resolution based on dominant wavelength ────────
        T_ref       = max(self.wave_period) if max(self.wave_period) > 0 else 10.0
        lambda_wave = 9.81 * T_ref**2 / (2 * np.pi)
        n_wave      = int(np.clip(2 * d / lambda_wave * 8, 20, 30))
        X_wave, Y_wave = np.meshgrid(
            np.linspace(-d, d, n_wave),
            np.linspace(-d, d, n_wave),
        )

        # ── Time index range ─────────────────────────────────────
        i1 = np.argmax(self.time >= self.t_ini)
        i2 = (np.argmax(self.time > self.t_fin)
              if self.t_fin < self.time[-1] else len(self.time) - 1)
        di = max(1, int(self.dt_video / (self.time[1] - self.time[0])))
        time_indices = range(i1, i2, di)

        def update(frame_idx):
            ax.cla()
            ax.computed_zorder = False  # respect explicit zorder values

            it = time_indices[frame_idx]
            t  = self.time[it]

            # Seabed (zorder=1)
            if self.bathymetry_xyz is not None:
                ax.plot_trisurf(_bathy_tri, _bz, color='tan',
                                alpha=0.80, shade=False, zorder=1)
            else:
                ax.plot_surface(_Xb, _Yb, _Zb, color='tan',
                                alpha=0.80, shade=False, zorder=1)

            # Mooring lines (zorder=3)
            if self.flag_lines and self.line_positions:
                for line_idx in self.lines_plot:
                    i = line_idx - 1
                    if i in self.line_positions:
                        df_pos  = self.line_positions[i]
                        x_cols  = [c for c in df_pos.columns if c.startswith('x_')]
                        y_cols  = [c for c in df_pos.columns if c.startswith('y_')]
                        z_cols  = [c for c in df_pos.columns if c.startswith('z_')]
                        if x_cols and y_cols and z_cols:
                            row = df_pos.iloc[(df_pos['time'] - t).abs().argmin()]
                            ax.plot(row[x_cols].values, row[y_cols].values,
                                    row[z_cols].values,
                                    'kx-', linewidth=1.5, markersize=3, zorder=3)

            # Body box (zorder=5)
            pos = [
                self.dof_data[1][it],
                self.dof_data[2][it],
                self.dof_data[3][it],
                np.rad2deg(self.dof_data[4][it]),
                np.rad2deg(self.dof_data[5][it]),
                np.rad2deg(self.dof_data[6][it]),
            ]
            vertices = self.create_box_vertices(self.body_dim, pos, self.z0)
            faces    = self.create_box_faces(vertices)
            box = Poly3DCollection(faces, alpha=1.0, facecolor=self.body_color,
                                   edgecolor='k', linewidth=0.5, zorder=5)
            ax.add_collection3d(box)

            # Wave surface (zorder=2)
            Z_wave = self.compute_wave_surface(X_wave, Y_wave, t)
            ax.plot_surface(X_wave, Y_wave, Z_wave, color='skyblue',
                            alpha=0.5, shade=False, zorder=2)

            ax.set_xlabel('X [m]')
            ax.set_ylabel('Y [m]')
            ax.set_zlabel('Z [m]')
            ax.set_xlim(-d, d)
            ax.set_ylim(-d, d)
            ax.set_zlim(-self.water_depth - 5, max(5.0, self.water_depth * 0.1))
            ax.view_init(elev=self.view_angle[1], azim=self.view_angle[0])
            ax.set_box_aspect([1, 1, 2 / 3])
            ax.set_title(f'Time = {t:.1f} s', fontsize=12, fontweight='bold')
            return ax,

        n_frames = len(time_indices)
        print(f"  Generating {n_frames} frames...")
        anim = FuncAnimation(fig, update, frames=n_frames,
                             interval=int(1000 * self.dt_video / self.speed),
                             blit=False, repeat=True)

        if self.flag_save_video:
            save_dir   = os.path.join(self.case_path, 'results')
            os.makedirs(save_dir, exist_ok=True)
            video_file = os.path.join(save_dir, f'{self.case_name}_video.mp4')
            print(f"  Saving video to {video_file}")
            writer = FFMpegWriter(
                fps=int(1.0 / (self.dt_video / self.speed)), bitrate=1800
            )
            anim.save(video_file, writer=writer)
            print("  Video saved!")
        else:
            plt.show()
    
    def run(self):
        """Run all plotting routines"""
        print("=" * 70)
        print("OASIS Results Plotting Tool")
        print("=" * 70)
        print(f"Case: {self.case_name}")
        print(f"Case path: {self.case_path}")
        print()

        # Load case parameters from YAML
        self.load_input_yaml()

        # Load output data
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
    """Main function — configure CASE_PATH and BODY_ID at the top of this file."""
    if not os.path.exists(CASE_PATH):
        print(f"Error: Case directory not found: {CASE_PATH}")
        sys.exit(1)

    plotter = OASISResultsPlotter(CASE_PATH)
    plotter.body_id = BODY_ID

    # Configure options (set as needed)
    plotter.flag_plot        = True   # Plot DOF time series
    plotter.flag_video       = True  # Create 3D animation
    plotter.flag_save        = False  # Save figures to disk
    plotter.flag_save_video  = False  # Save animation video

    plotter.flag_lines  = True  # Plot mooring line tensions
    plotter.n_lines     = 4
    plotter.lines_plot  = [1, 2, 3, 4]

    plotter.flag_wt = False  # Plot wind turbine data

    # Video settings (only used when flag_video = True)
    plotter.t_ini       = 0
    plotter.t_fin       = 100
    plotter.dt_video    = 0.1
    plotter.speed       = 10
    plotter.view_angle  = [45, 10]

    plotter.run()


if __name__ == "__main__":
    main()
