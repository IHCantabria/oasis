"""
plot_oasis_layout.py
--------------------
Static geometry check for an OASIS input folder.

Reads dataProblem.yaml and any referenced bathymetry .dat files, then
produces a 3-panel figure (3D view, XY top, XZ side) showing:

  * seafloor (triangulated bathymetry mesh or flat plane)
  * mooring lines (catenary shape solved from length + endpoints)
  * body bounding box at its initial position
  * BCP markers: anchors (black triangles) and fairleads (red circles)

No simulation output is required — only the input folder.

Usage
-----
  1. Set CASE_PATH at the bottom of this file.
  2. Run:
       python plot_oasis_layout.py

BCP indexing convention (matches OASIS YAML)
--------------------------------------------
All BCPs share a single 1-based global index list:
  - anchors occupy indices 1 … N_anchors
  - body BCPs occupy indices N_anchors+1 … N_anchors+N_body_bcps
Line fields BCP_N (fairlead) and BCP_1 (anchor) use this global numbering.
Body-BCP positions in the YAML are in the body's local frame; they are
transformed by each body's initial_position before display.
"""

import os
import sys

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.tri import Triangulation
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401 — registers 3D projection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ---------------------------------------------------------------------------
# CONFIGURE HERE before running
# ---------------------------------------------------------------------------
CASE_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    '..', 'examples', 'lines', 'multi_line_yaml'
)
# ---------------------------------------------------------------------------

_LINE_COLORS = [
    'steelblue', 'darkorange', 'seagreen', 'firebrick',
    'mediumpurple', 'goldenrod', 'deepskyblue', 'tomato',
]


# ══════════════════════════════════════════════════════════════ math helpers

def _rotation_matrix(roll, pitch, yaw):
    """ZYX (intrinsic) rotation matrix from Euler angles in radians."""
    cr, sr = np.cos(roll),  np.sin(roll)
    cp, sp = np.cos(pitch), np.sin(pitch)
    cy, sy = np.cos(yaw),   np.sin(yaw)
    return np.array([
        [cp*cy,  cy*sp*sr - cr*sy,  sr*sy + cr*cy*sp],
        [cp*sy,  cr*cy + sp*sr*sy,  cr*sp*sy - cy*sr],
        [-sp,    cp*sr,             cp*cr            ],
    ])


def _create_box_vertices(dim, center, roll=0.0, pitch=0.0, yaw=0.0):
    """
    Compute the 8 corner vertices of a 3-D box.

    Parameters
    ----------
    dim    : [length, width, height]  (metres)
    center : [x, y, z]               (global position)
    roll, pitch, yaw : rotation in *radians* (ZYX convention)

    Returns
    -------
    vertices : ndarray, shape (8, 3)
    """
    a, b, c = dim
    half = np.array([
        [ a/2,  b/2,  c/2],
        [-a/2,  b/2,  c/2],
        [-a/2, -b/2,  c/2],
        [ a/2, -b/2,  c/2],
        [ a/2,  b/2, -c/2],
        [-a/2,  b/2, -c/2],
        [-a/2, -b/2, -c/2],
        [ a/2, -b/2, -c/2],
    ])
    R = _rotation_matrix(roll, pitch, yaw)
    return half @ R.T + np.array(center)


def _create_box_faces(vertices):
    """
    Return the 6 quad faces of a box from its 8 corner vertices.

    Returns
    -------
    faces : list of ndarray, each shape (4, 3)
    """
    face_indices = [
        [0, 1, 2, 3],  # top
        [4, 5, 6, 7],  # bottom
        [0, 3, 7, 4],  # right
        [1, 2, 6, 5],  # left
        [0, 1, 5, 4],  # front
        [2, 3, 7, 6],  # back
    ]
    return [vertices[idx] for idx in face_indices]


def _catenary_nodes(anchor, fairlead, length, seabed_z, n=60):
    """
    Compute n 3D points along the catenary between *anchor* and *fairlead*.

    The catenary is solved in the vertical plane defined by the two endpoints.
    If the line is geometrically taut (length <= 3-D distance), a straight
    segment is returned.  Z values are clamped to *seabed_z*.

    Parameters
    ----------
    anchor, fairlead : array-like, shape (3,)
    length           : unstretched line length [m]
    seabed_z         : seafloor depth [m, negative]
    n                : number of sample points

    Returns
    -------
    xs, ys, zs : ndarray, shape (n,)
    """
    xa, ya, za = float(anchor[0]),   float(anchor[1]),   float(anchor[2])
    xf, yf, zf = float(fairlead[0]), float(fairlead[1]), float(fairlead[2])

    dh = np.hypot(xf - xa, yf - ya)      # horizontal distance
    dv = zf - za                           # signed vertical rise (+ = fairlead higher)
    dist3d = np.sqrt(dh**2 + dv**2)

    # Horizontal unit vector anchor → fairlead
    if dh > 1e-8:
        ehx, ehy = (xf - xa) / dh, (yf - ya) / dh
    else:
        # Nearly vertical drop — straight line
        t = np.linspace(0.0, 1.0, n)
        return xa + t*(xf - xa), ya + t*(yf - ya), za + t*(zf - za)

    # Taut check
    if length <= dist3d * 1.0005:
        t = np.linspace(0.0, 1.0, n)
        return xa + t*(xf - xa), ya + t*(yf - ya), za + t*(zf - za)

    # ── Catenary solver ──────────────────────────────────────────────────
    # Local 2-D frame: anchor at origin, fairlead at (dh, dv).
    # Shape:  z(u) = a * cosh((u - u0)/a) − a * cosh(u0/a)
    # Constraints:
    #   (i)  height: a * (cosh((dh−u0)/a) − cosh(u0/a)) = dv
    #   (ii) length: a * (sinh((dh−u0)/a) + sinh(u0/a)) = L
    solved = False
    a_sol = u0_sol = None

    try:
        from scipy.optimize import fsolve

        def _eqs(params):
            a, u0 = params
            if a < 1e-3:
                return [1e9, 1e9]
            s1 = (dh - u0) / a
            s2 = u0 / a
            return [
                a * (np.cosh(s1) - np.cosh(s2)) - dv,
                a * (np.sinh(s1) + np.sinh(s2)) - length,
            ]

        best_res = np.inf
        for a0 in (length, dh * 2.0, length * 0.3):
            for frac in (0.5, 0.3, 0.7):
                try:
                    sol, _, ier, _ = fsolve(_eqs, [a0, dh * frac],
                                            full_output=True)[:4]
                    a_c, u0_c = sol
                    if ier == 1 and a_c > 1e-3:
                        res = max(abs(r) for r in _eqs(sol))
                        if res < best_res:
                            best_res, a_sol, u0_sol = res, a_c, u0_c
                except Exception:
                    pass

        if a_sol is not None and best_res < 1.0:    # 1 m tolerance
            u = np.linspace(0.0, dh, n)
            z_loc = a_sol * np.cosh((u - u0_sol) / a_sol) \
                  - a_sol * np.cosh(u0_sol / a_sol)
            zs = np.maximum(za + z_loc, seabed_z)
            return xa + u * ehx, ya + u * ehy, zs

    except ImportError:
        pass    # scipy optional — fall through to parabolic fallback

    # ── Parabolic fallback ───────────────────────────────────────────────
    sag = np.sqrt(max(length**2 - dh**2 - dv**2, 0.0)) * 0.4
    t = np.linspace(0.0, 1.0, n)
    xs = xa + t * (xf - xa)
    ys = ya + t * (yf - ya)
    zs = za + t * (zf - za) - 4.0 * sag * t * (1.0 - t)
    return xs, ys, np.maximum(zs, seabed_z)


# ══════════════════════════════════════════════════ bathymetry file loader

def _load_bathymetry_dat(filepath):
    """
    Load an OASIS bathymetry .dat point cloud.

    Expected format::

        N // Num Puntos
        ////... (separator lines, skipped)
        x1  y1  z1
        x2  y2  z2
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
        pts = []
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


# ══════════════════════════════════════════════════════════════ main class

class OASISLayoutPlotter:
    """
    Static layout plotter for an OASIS input directory.

    No simulation results are required; only dataProblem.yaml (and
    optionally any bathymetry .dat files referenced inside it) are read.
    """

    def __init__(self, case_path):
        self.case_path  = case_path
        self.input_path = os.path.join(case_path, 'input')
        self.case_name  = os.path.basename(case_path.rstrip('/\\'))

        # Parsed state
        self.problem    = {}
        self.bodies     = []
        self.lines_cfg  = []
        self.all_bcps   = []     # 1-based list; index 0 = None (placeholder)
        self.n_anchors  = 0
        self.bathy_xyz  = None   # Nx3 ndarray or None
        self.water_depth = -100.0
        self.domain      = 100.0

    # ─────────────────────────────────────────────────────────── loading

    def load(self):
        """Parse dataProblem.yaml and load any bathymetry meshes."""
        yaml_file = os.path.join(self.input_path, 'dataProblem.yaml')
        if not os.path.exists(yaml_file):
            raise FileNotFoundError(
                f"dataProblem.yaml not found in: {self.input_path}"
            )

        with open(yaml_file, 'r') as fh:
            data = yaml.safe_load(fh)

        self.problem    = data.get('problem', {})
        self.bodies     = data.get('bodies', [])
        self.lines_cfg  = data.get('lines', [])
        self.water_depth = float(self.problem.get('water_depth', -100.0))

        # ── BCPs ────────────────────────────────────────────────────────
        bcps_sec  = data.get('bcps', {})
        anchors   = bcps_sec.get('anchors', [])
        body_bcps = bcps_sec.get('body_bcps', [])
        self.n_anchors = len(anchors)

        # Map global BCP index → body index (for body-BCP transform)
        bcp_to_body = {}
        for bi, body in enumerate(self.bodies):
            for gi in body.get('bcps_indexes', []):
                bcp_to_body[gi] = bi

        # Index 0 is unused; BCPs are 1-based
        self.all_bcps = [None]

        # Anchors — positions given in global frame
        for anc in anchors:
            self.all_bcps.append({
                'pos':  np.array(anc['position'], dtype=float),
                'type': 'anchor',
            })

        # Body BCPs — given in body-local frame; apply initial_position transform
        for i, bcp in enumerate(body_bcps):
            gi    = self.n_anchors + i + 1          # global 1-based index
            local = np.array(bcp['position'], dtype=float)
            bi    = bcp_to_body.get(gi)
            if bi is not None:
                ip = self.bodies[bi].get('initial_position', [0] * 6)
                R  = _rotation_matrix(ip[3], ip[4], ip[5])
                t  = np.array(ip[:3], dtype=float)
                gpos = R @ local + t
            else:
                gpos = local.copy()

            self.all_bcps.append({
                'pos':       gpos,
                'type':      'fairlead',
                'body_idx':  bi,
                'local_pos': local,
            })

        # ── Seafloor ────────────────────────────────────────────────────
        seafloor      = data.get('seafloor', {})
        bathy_entries = seafloor.get('bathymetry', [])
        if bathy_entries:
            all_pts = []
            for entry in bathy_entries:
                path = os.path.join(self.input_path, entry.get('mesh_file', ''))
                pts  = _load_bathymetry_dat(path)
                if pts is not None:
                    all_pts.append(pts)
            if all_pts:
                self.bathy_xyz = np.vstack(all_pts)
        else:
            flat = seafloor.get('flat', [])
            if flat:
                self.water_depth = float(flat[0].get('depth', self.water_depth))

        # ── Domain size (from anchors and bathymetry) ────────────────────
        anchor_pos = [self.all_bcps[i]['pos']
                      for i in range(1, self.n_anchors + 1)]
        if anchor_pos:
            self.domain = max(np.linalg.norm(p[:2]) for p in anchor_pos) * 1.15
        if self.bathy_xyz is not None:
            r_bathy = np.max(np.linalg.norm(self.bathy_xyz[:, :2], axis=1))
            self.domain = max(self.domain, r_bathy)

        self._print_summary()

    def _print_summary(self):
        n_fl = len(self.all_bcps) - 1 - self.n_anchors
        sf   = ('bathymetry mesh' if self.bathy_xyz is not None
                else f'flat @ {self.water_depth:.1f} m')
        print(f"Case     : {self.case_name}")
        print(f"  Bodies   : {len(self.bodies)}")
        print(f"  Anchors  : {self.n_anchors}")
        print(f"  Fairleads: {n_fl}")
        print(f"  Lines    : {len(self.lines_cfg)}")
        print(f"  Seafloor : {sf}")
        print(f"  Domain   : ±{self.domain:.1f} m")

    # ─────────────────────────────────────────────────────────── drawing

    def _draw_seafloor(self, ax3d, ax_top, ax_side):
        d = self.domain

        if self.bathy_xyz is not None:
            bx = self.bathy_xyz[:, 0]
            by = self.bathy_xyz[:, 1]
            bz = self.bathy_xyz[:, 2]
            tri = Triangulation(bx, by)

            # 3D surface
            ax3d.plot_trisurf(tri, bz, color='tan', alpha=0.70,
                              shade=False, zorder=1)
            # Top view: filled contour + colorbar
            cf = ax_top.tricontourf(tri, bz, levels=15, cmap='terrain', alpha=0.75)
            plt.colorbar(cf, ax=ax_top, label='z [m]', shrink=0.75, pad=0.02)
            ax_top.tricontour(tri, bz, levels=8, colors='k',
                              linewidths=0.4, alpha=0.4)
            # Side view: scatter  X vs Z
            ax_side.scatter(bx, bz, s=1, c='tan', alpha=0.5,
                            label=f'bathymetry')
        else:
            zbd = self.water_depth
            # 3D flat plane
            xg = np.array([[-d, d], [-d, d]])
            yg = np.array([[-d, -d], [d, d]])
            zg = np.full_like(xg, zbd)
            ax3d.plot_surface(xg, yg, zg, color='tan', alpha=0.55,
                              edgecolors='none', zorder=1)
            # Top view: filled rectangle
            ax_top.fill([-d, d, d, -d], [-d, -d, d, d],
                        color='tan', alpha=0.35,
                        label=f'seafloor ({zbd:.0f} m)')
            # Side view: horizontal line
            ax_side.axhline(zbd, color='#b8860b', linewidth=2.0,
                            label=f'seafloor ({zbd:.0f} m)')

        # Water surface — side view dashed line
        ax_side.axhline(0.0, color='steelblue', linewidth=1.5,
                        linestyle='--', label='waterline')
        # Water surface — 3D semi-transparent plane
        xg = np.array([[-d, d], [-d, d]])
        yg = np.array([[-d, -d], [d, d]])
        ax3d.plot_surface(xg, yg, np.zeros_like(xg),
                          color='skyblue', alpha=0.18,
                          edgecolors='none', zorder=2)

    def _draw_lines(self, ax3d, ax_top, ax_side):
        issues = []
        for li, line in enumerate(self.lines_cfg):
            bcp_n = line.get('BCP_N')   # fairlead (body side)
            bcp_1 = line.get('BCP_1')   # anchor
            L     = float(line.get('length', 0.0))
            color = _LINE_COLORS[li % len(_LINE_COLORS)]

            if bcp_n is None or bcp_1 is None:
                print(f"  Warning: line {li+1} missing BCP_N or BCP_1 — skipped")
                continue
            if bcp_n >= len(self.all_bcps) or self.all_bcps[bcp_n] is None:
                print(f"  Warning: line {li+1} BCP_N={bcp_n} out of range — skipped")
                continue
            if bcp_1 >= len(self.all_bcps) or self.all_bcps[bcp_1] is None:
                print(f"  Warning: line {li+1} BCP_1={bcp_1} out of range — skipped")
                continue

            fairlead = self.all_bcps[bcp_n]['pos']
            anchor   = self.all_bcps[bcp_1]['pos']
            xs, ys, zs = _catenary_nodes(anchor, fairlead, L,
                                         self.water_depth, n=80)

            label = f'Line {li+1}  (L={L:.0f} m)'
            ax3d.plot(xs, ys, zs, '-', color=color, linewidth=1.8,
                      label=label, zorder=4)
            ax_top.plot(xs, ys, '-', color=color, linewidth=1.5, label=label)
            ax_side.plot(xs, zs, '-', color=color, linewidth=1.5, label=label)

            # Seabed contact check
            if zs.min() <= self.water_depth + 0.1:
                issues.append(
                    f"Line {li+1}: touches seabed  (min z = {zs.min():.2f} m)"
                )

        if issues:
            print("\n  !! Layout warnings:")
            for w in issues:
                print(f"     {w}")

    def _draw_body(self, ax3d, ax_top, ax_side):
        for bi, body in enumerate(self.bodies):
            bcp_idxs = body.get('bcps_indexes', [])

            # Body centre = initial_position + initial_displacement (translational)
            ip = body.get('initial_position', [0.0] * 6)
            idp = body.get('initial_displacement', [0.0] * 6)
            cx = float(ip[0]) + float(idp[0])
            cy_b = float(ip[1]) + float(idp[1])
            cz = float(ip[2]) + float(idp[2])
            roll  = float(ip[3])
            pitch = float(ip[4])
            yaw   = float(ip[5])

            # Collect fairlead global positions to estimate box size
            if bcp_idxs:
                pts = np.array([
                    self.all_bcps[gi]['pos']
                    for gi in bcp_idxs
                    if gi < len(self.all_bcps) and self.all_bcps[gi] is not None
                ])
            else:
                pts = np.empty((0, 3))

            if pts.size > 0:
                # Use fairlead spread for XY dimensions; estimate height
                dx = max(pts[:, 0].max() - pts[:, 0].min(), 1.0)
                dy = max(pts[:, 1].max() - pts[:, 1].min(), 1.0)
                # Freeboard: body extends from deepest fairlead up to waterline
                dz = max(abs(cz) * 0.30, 2.0)
            else:
                # No BCP data — use a sensible default box
                dx, dy, dz = 10.0, 10.0, 4.0

            # Add a small margin around fairlead span
            mg = max(dx * 0.10, 0.5)
            body_dim = [dx + 2 * mg, dy + 2 * mg, dz]

            # ── 3D box using proper rotation ─────────────────────────────
            verts = _create_box_vertices(body_dim, [cx, cy_b, cz + dz / 2],
                                         roll, pitch, yaw)
            faces = _create_box_faces(verts)
            ax3d.add_collection3d(
                Poly3DCollection(faces, alpha=0.50, facecolor='#b8b8b8',
                                 edgecolor='k', linewidth=0.8, zorder=5)
            )

            # ── 2-D projections ─────────────────────────────────────────
            # Project 8 vertices onto XY and XZ planes
            xv, yv, zv = verts[:, 0], verts[:, 1], verts[:, 2]
            ax_top.add_patch(mpatches.Polygon(
                np.column_stack([xv[[0,1,2,3]], yv[[0,1,2,3]]]),
                closed=True, facecolor='#b8b8b8', edgecolor='k',
                alpha=0.55, linewidth=1.5, label=f'Body {bi+1}',
            ))
            ax_side.add_patch(mpatches.Polygon(
                np.column_stack([xv[[0,3,7,4]], zv[[0,3,7,4]]]),
                closed=True, facecolor='#b8b8b8', edgecolor='k',
                alpha=0.55, linewidth=1.5,
            ))

    def _draw_bcps(self, ax3d, ax_top, ax_side):
        # Anchors
        for i in range(1, self.n_anchors + 1):
            p = self.all_bcps[i]['pos']
            ax3d.scatter(*p, color='k', s=50, marker='^',
                         zorder=7, depthshade=False)
            ax_top.plot(p[0], p[1], 'k^', markersize=9, zorder=7,
                        label='Anchor' if i == 1 else None)
            ax_side.plot(p[0], p[2], 'k^', markersize=9, zorder=7,
                         label='Anchor' if i == 1 else None)

        # Fairleads
        n_fl = len(self.all_bcps) - 1 - self.n_anchors
        for j, i in enumerate(range(self.n_anchors + 1, len(self.all_bcps))):
            if self.all_bcps[i] is None:
                continue
            p = self.all_bcps[i]['pos']
            ax3d.scatter(*p, color='crimson', s=40, marker='o',
                         zorder=8, depthshade=False)
            ax_top.plot(p[0], p[1], 'o', color='crimson', markersize=7,
                        zorder=8, label='Fairlead' if j == 0 else None)
            ax_side.plot(p[0], p[2], 'o', color='crimson', markersize=7,
                         zorder=8, label='Fairlead' if j == 0 else None)

    # ───────────────────────────────────────────────────── main entry

    def plot(self, save_path=None):
        """
        Build and show (or save) the 3-panel layout figure.

        Parameters
        ----------
        save_path : str or None
            If given, the figure is saved to this path instead of displayed.
        """
        fig = plt.figure(figsize=(18, 7))
        ax3d  = fig.add_subplot(131, projection='3d')
        ax_top  = fig.add_subplot(132)
        ax_side = fig.add_subplot(133)

        d   = self.domain
        zbd = self.water_depth
        z_top = max(5.0, abs(zbd) * 0.08)

        self._draw_seafloor(ax3d, ax_top, ax_side)
        self._draw_lines(ax3d, ax_top, ax_side)
        self._draw_body(ax3d, ax_top, ax_side)
        self._draw_bcps(ax3d, ax_top, ax_side)

        # ── 3D axis ──────────────────────────────────────────────────────
        ax3d.set_xlabel('X [m]', labelpad=4)
        ax3d.set_ylabel('Y [m]', labelpad=4)
        ax3d.set_zlabel('Z [m]', labelpad=4)
        ax3d.set_xlim(-d, d)
        ax3d.set_ylim(-d, d)
        ax3d.set_zlim(zbd - 2, z_top)
        ax3d.set_title('3-D view', fontsize=11)
        ax3d.view_init(elev=20, azim=45)
        ax3d.legend(
            handles=[
                mpatches.Patch(color='#b8b8b8', label='Body'),
                plt.Line2D([0], [0], marker='^', color='k',
                           markersize=8, linestyle='none', label='Anchor'),
                plt.Line2D([0], [0], marker='o', color='crimson',
                           markersize=7, linestyle='none', label='Fairlead'),
            ],
            loc='upper right', fontsize=8,
        )

        # ── Top view ─────────────────────────────────────────────────────
        ax_top.set_aspect('equal', adjustable='box')
        ax_top.set_xlabel('X [m]')
        ax_top.set_ylabel('Y [m]')
        ax_top.set_xlim(-d, d)
        ax_top.set_ylim(-d, d)
        ax_top.set_title('Top view  (XY)', fontsize=11)
        ax_top.grid(True, alpha=0.3)
        # Deduplicate legend entries
        handles, labels = ax_top.get_legend_handles_labels()
        seen = {}
        for h, l in zip(handles, labels):
            seen.setdefault(l, h)
        ax_top.legend(seen.values(), seen.keys(), fontsize=7, loc='upper right')

        # ── Side view ────────────────────────────────────────────────────
        ax_side.set_xlabel('X [m]')
        ax_side.set_ylabel('Z [m]')
        ax_side.set_xlim(-d, d)
        ax_side.set_ylim(zbd - 5, z_top + 5)
        ax_side.set_title('Side view  (XZ)', fontsize=11)
        ax_side.grid(True, alpha=0.3)
        handles, labels = ax_side.get_legend_handles_labels()
        seen = {}
        for h, l in zip(handles, labels):
            seen.setdefault(l, h)
        ax_side.legend(seen.values(), seen.keys(), fontsize=7, loc='upper right')

        fig.suptitle(f'OASIS Layout Check — {self.case_name}',
                     fontsize=13, fontweight='bold')
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Figure saved to: {save_path}")
        else:
            plt.show()


# ══════════════════════════════════════════════════════════════ entry point

def main():
    if not os.path.exists(CASE_PATH):
        print(f"Error: case directory not found: {CASE_PATH}")
        sys.exit(1)

    plotter = OASISLayoutPlotter(CASE_PATH)
    plotter.load()
    plotter.plot()
    # To save instead of showing:
    # plotter.plot(save_path='layout.png')


if __name__ == '__main__':
    main()
