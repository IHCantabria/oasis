"""
generate_flat_bathymetry.py
---------------------------
Generates a flat bathymetry mesh (.dat) for OASIS using a structured
nx × ny grid.  The default produces exactly 200 points (10 rows × 20
columns) uniformly covering a ±1000 m square at a constant depth,
matching the surface of the 4-point reference rectangle:

    (-1000, -1000, -100)  …  (1000, 1000, -100)

Geometry  (default: nx=10 rows, ny=20 columns)
----------------------------------------------
Points are arranged in a regular row-major grid:

  row nx  (y = +hy)  ●—●— … —●
  ...                |╲|╲   |╲|
  row 1   (y = -hy)  ●—●— … —●

Each rectangular cell (i, j) is split into two triangles:
    T_A: (p_bottom-left, p_bottom-right, p_top-left)
    T_B: (p_bottom-right, p_top-right,   p_top-left)

Both use counter-clockwise vertex ordering so that the cross product
(V1-V0) × (V2-V0) has a positive z-component — normals face upward,
as required by the OASIS bathymetry projection algorithm.  The cells
tile the domain exactly with no gaps or overlaps.

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

Usage
-----
python generate_flat_bathymetry.py [--output PATH] [--depth D]
                                   [--half-x HX] [--half-y HY]
                                   [--nx NX] [--ny NY]

Defaults:
  depth  : -100.0 m
  half-x : 1000.0 m  (rectangle spans -1000 to +1000 in x)
  half-y : 1000.0 m  (rectangle spans -1000 to +1000 in y)
  nx     : 10        (rows;    10 × 20 = 200 points total)
  ny     : 20        (columns; dx = dy = 105 m approx)
  output : ../examples/lines/complex_bathymetry/input/flat_bathymetry.dat

Note: for a uniform flat seabed use --minimal (or nx=2, ny=2) to produce
only 4 points / 2 triangles.  Because OASIS searches all triangles with a
brute-force O(T × N_nodes) loop, extra triangles on a flat mesh add
cost without any geometric benefit.
"""

import argparse
import os

SEPARATOR = '////////////////////\n'


def build_mesh(depth, half_x, half_y, nx, ny):
    """
    Return (points, triangles) for a flat rectangular nx × ny grid.

    Node numbering: row-major, row 1 at y=-half_y, row nx at y=+half_y.
    Node index (1-based) for row i (0-based), column j (0-based):
        idx = i * ny + j + 1

    Each rectangular cell (i, j) → two CCW triangles (normals up):
        T_A: (i*ny+j+1,  i*ny+j+2,      (i+1)*ny+j+1)  bottom-left  triangle
        T_B: (i*ny+j+2,  (i+1)*ny+j+2,  (i+1)*ny+j+1)  top-right triangle

    Cells tile the domain exactly — no gaps, no overlaps.
    For nx=10, ny=20: N=200 points, T=2*(nx-1)*(ny-1)=342 triangles.
    """
    import numpy as np

    x = np.linspace(-half_x, half_x, ny)
    y = np.linspace(-half_y, half_y, nx)

    points = []
    for i in range(nx):
        for j in range(ny):
            points.append((x[j], y[i], depth))

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


def write_dat(output_path, points, triangles, depth, half_x, half_y, nx, ny):
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
    print(f'Written: {output_path}')
    print(f'  Points   : {num_points}  ({nx} rows x {ny} cols)')
    print(f'  Triangles: {num_triangles}  ({nx-1} x {ny-1} cells x 2)')
    print(f'  Cell size: dx={dx:.2f} m  dy={dy:.2f} m')
    print(f'  Extent   : x=[{-half_x:.1f}, {half_x:.1f}]  y=[{-half_y:.1f}, {half_y:.1f}]')
    print(f'  Depth    : {depth:.2f} m (uniform)')


if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_output = os.path.join(
        script_dir,
        '..',
        'examples', 'lines', 'complex_bathymetry', 'input',
        'flat_bathymetry.dat',
    )

    depth = -100
    half_x = 50
    half_y = 50
    nx = 5
    ny = 5    

    points, triangles = build_mesh(depth, half_x, half_y, nx, ny)
    write_dat(os.path.normpath(default_output), points, triangles, depth, half_x, half_y, nx, ny)
