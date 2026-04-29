"""
generate_complex_bathymetry.py
------------------------------
Generates a sloped bathymetry mesh (.dat) for OASIS using a structured
nx × ny grid.  The seabed is a planar incline described by:

    z(x, y) = max_depth + slope_x * x + slope_y * y

Defaults produce a 33 × 33 grid over ±200 m with a slope of 0.1 m/m in
the x-direction (seabed rises from −120 m at x=−200 to −80 m at x=+200,
passing through −100 m at the origin).

Geometry  (nx=33 rows, ny=33 columns)
--------------------------------------
Points are arranged in a regular row-major grid:

  row 33 (y = +hy)  ●—●— … —●
  ...               |╲|╲   |╲|
  row 1  (y = -hy)  ●—●— … —●

Each rectangular cell (i, j) is split into two triangles:
    T_A: (p_bottom-left, p_bottom-right, p_top-left)
    T_B: (p_bottom-right, p_top-right,   p_top-left)

Both use counter-clockwise vertex ordering → cross product
(V1-V0) × (V2-V0) has positive z → normals face upward.

Output format (matches triangularMesh_Steep_allchain_Extrapolation.dat)
-----------------------------------------------------------------------
N // Num Puntos
////////////////////   (x3 separator block)
x1 y1 z1
...
////////////////////   (x3 separator block)
M // Num Triangulos
////////////////////   (x3 separator block)
v1 v2 v3              (1-based node indices, counter-clockwise = normals up)
...
////////////////////   (x3 separator block)
0 // [0: Normal method, 1: Barycenter]

Defaults:
  max_depth : -100.0 m   (depth at x=0, y=0)
  slope_x   :    0.1     (m/m; seabed rises toward +x)
  slope_y   :    0.0     (no y-slope)
  half-x    :  200.0 m   (domain spans x ∈ [-200, +200])
  half-y    :  200.0 m   (domain spans y ∈ [-200, +200])
  nx        :   33       (rows)
  ny        :   33       (columns; ~12.5 m cell spacing)
  output    : ../examples/lines/complex_bathymetry/input/complex_bathymetry.dat
"""

import argparse
import os

SEPARATOR = '////////////////////\n'

def depth(x, y, max_depth=-100.0, slope_x=0.1, slope_y=0.0, half_x=100.0, half_y=100.0,
          end_slope_x=None, end_slope_y=None):
    """Linear slope depth function (m, negative below MSL).

    max_depth is the deepest point of the domain (the absolute floor).
    The slope starts at the deep border and ends at end_slope_x / end_slope_y.
    Beyond those coordinates the seabed is flat (depth held constant).
    If end_slope_x / end_slope_y is None the slope extends to the opposite border.

        z(x, y) = max_depth + slope_x*(clamp(x) - x_deep) + slope_y*(clamp(y) - y_deep)

    where (x_deep, y_deep) is the domain corner with the lowest seabed
    (most negative z).  For slope_x > 0 that corner is at x = -half_x;
    for slope_x < 0 it is at x = +half_x; similarly for y.

    This guarantees z >= max_depth everywhere in the domain.
    """
    x_deep = -half_x if slope_x >= 0 else half_x
    y_deep = -half_y if slope_y >= 0 else half_y

    # Resolve defaults: slope covers the full domain when end is not specified
    if end_slope_x is None:
        end_slope_x = half_x if slope_x >= 0 else -half_x
    if end_slope_y is None:
        end_slope_y = half_y if slope_y >= 0 else -half_y

    # Clamp coordinates to the sloped region; beyond end_slope_* the bed is flat
    eff_x = min(x, end_slope_x) if slope_x >= 0 else max(x, end_slope_x)
    eff_y = min(y, end_slope_y) if slope_y >= 0 else max(y, end_slope_y)

    return max_depth + slope_x * (eff_x - x_deep) + slope_y * (eff_y - y_deep)

def build_mesh(max_depth, slope_x, slope_y, half_x, half_y, nx, ny,
               end_slope_x=None, end_slope_y=None):
    """
    Return (points, triangles) for a sloped rectangular nx × ny grid.

    Node numbering: row-major, row 1 at y=-half_y, row nx at y=+half_y.
    Node index (1-based) for row i (0-based), column j (0-based):
        idx = i * ny + j + 1

    Each rectangular cell (i, j) → two CCW triangles (normals up):
        T_A: (i*ny+j+1,  i*ny+j+2,      (i+1)*ny+j+1)  bottom-left
        T_B: (i*ny+j+2,  (i+1)*ny+j+2,  (i+1)*ny+j+1)  top-right
    """
    import numpy as np

    x = np.linspace(-half_x, half_x, ny)
    y = np.linspace(-half_y, half_y, nx)

    points = []
    for i in range(nx):
        for j in range(ny):
            points.append((x[j], y[i], depth(x[j], y[i], max_depth, slope_x, slope_y, half_x, half_y,
                                                end_slope_x, end_slope_y)))

    triangles = []
    for i in range(nx - 1):
        for j in range(ny - 1):
            p00 = i * ny + j + 1        # bottom-left
            p10 = i * ny + j + 2        # bottom-right
            p01 = (i + 1) * ny + j + 1  # top-left
            p11 = (i + 1) * ny + j + 2  # top-right
            # CCW winding → cross product z > 0 (normal faces up)
            triangles.append((p00, p10, p01))   # T_A
            triangles.append((p10, p11, p01))   # T_B

    return points, triangles


def write_dat(output_path, points, triangles, max_depth, slope_x, slope_y, half_x, half_y, nx, ny,
              end_slope_x=None, end_slope_y=None):
    """Write the .dat mesh file in OASIS bathymetry format."""
    num_points = len(points)
    num_triangles = len(triangles)

    dirname = os.path.dirname(output_path)
    if dirname:
        os.makedirs(dirname, exist_ok=True)

    with open(output_path, 'w') as f:
        # --- point cloud ---
        f.write(f'{num_points} // Num Puntos \n')
        f.write(SEPARATOR * 3)
        for px, py, pz in points:
            f.write(f'{px:22.10e} {py:22.10e} {pz:22.10e} \n')
        f.write(SEPARATOR * 3)

        # --- triangles ---
        f.write(f'{num_triangles} // Num Triangulos \n')
        f.write(SEPARATOR * 3)
        for v1, v2, v3 in triangles:
            f.write(f'{v1} {v2} {v3} \n')
        f.write(SEPARATOR * 3)

        # --- flag ---
        f.write('0 // [0: Normal method, 1: Barycenter]\n')

    dx = 2 * half_x / (ny - 1)
    dy = 2 * half_y / (nx - 1)
    depth_max_z = depth( half_x if slope_x >= 0 else -half_x,
                          half_y if slope_y >= 0 else -half_y,
                          max_depth, slope_x, slope_y, half_x, half_y,
                          end_slope_x, end_slope_y)
    print(f'Written: {output_path}')
    print(f'  Points   : {num_points}  ({nx} rows x {ny} cols)')
    print(f'  Triangles: {num_triangles}  ({nx-1} x {ny-1} cells x 2)')
    print(f'  Cell size: dx={dx:.2f} m  dy={dy:.2f} m')
    print(f'  Extent   : x=[{-half_x:.1f}, {half_x:.1f}]  y=[{-half_y:.1f}, {half_y:.1f}]')
    print(f'  Depth range     : [{max_depth:.2f} (floor), {depth_max_z:.2f}] m')
    print(f'  Slope           : slope_x={slope_x}  slope_y={slope_y} (m/m)')


if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_output = os.path.join(
        script_dir,
        '..',
        'examples', 'lines', 'complex_bathymetry', 'input',
        'complex_bathymetry.dat',
    )

    max_depth = -100.0  # seabed depth at origin (m)
    slope_x   =    0.2  # seabed rises 0.2 m per metre toward +x
    slope_y   =    0.0  # no y-slope
    half_x    =  75.0  # domain ±100 m in x  (2× anchor radius; 50 m margin)
    half_y    =  75.0  # domain ±100 m in y
    nx          =   10    # rows    (≈22 m spacing)
    ny          =   10    # columns (≈22 m spacing)
    end_slope_x = 0  # x-coordinate where slope ends (None = full domain)
    end_slope_y = None  # y-coordinate where slope ends (None = full domain)

    points, triangles = build_mesh(max_depth, slope_x, slope_y, half_x, half_y, nx, ny,
                                   end_slope_x, end_slope_y)
    write_dat(os.path.normpath(default_output), points, triangles, max_depth, slope_x, slope_y,
              half_x, half_y, nx, ny, end_slope_x, end_slope_y)

    # Anchor z note: anchors must be placed at water_depth (-100 m), NOT at the
    # seabed surface depth.  The quasi-static line initialiser uses water_depth
    # as a flat floor; placing the anchor at the sloped seabed surface depth
    # raises it above -100 m and causes the catenary to sag below the floor
    # (exception 6).  The anchor is treated as embedded in the substrate at
    # water_depth; the slope only acts on interior nodes via FEM projection.
    print('\nAnchor positions for dataProblem.yaml:')
    print('  All anchors must be placed at water_depth = -100.0 m')
    anchor_xy = [(50.0, 50.0), (-50.0, 50.0), (50.0, -50.0), (-50.0, -50.0)]
    for ax, ay in anchor_xy:
        az_seabed = depth(ax, ay, max_depth, slope_x, slope_y, half_x, half_y,
                          end_slope_x, end_slope_y)
        print(f'  - position: [{ax:.1f}, {ay:.1f}, -100.0]   (seabed surface at this xy: {az_seabed:.1f} m)')
