"""
Converts dataSprings.dat files from the old format (all properties per-spring)
to the new format (spring types section + spring instances section).
"""
import os
import sys
import glob

def read_clean_lines(filepath):
    """Return list of non-empty lines with comments stripped."""
    result = []
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            clean = line.split('//')[0].strip()
            if clean:
                result.append(clean)
    return result

def parse_spring_entry(clean_lines, start_idx):
    """
    Parse one spring entry starting at start_idx.
    Skips header lines starting with '/'.
    Returns (instance_data, type_data, next_idx).
    """
    idx = start_idx

    # Skip header lines starting with '/'
    while idx < len(clean_lines) and clean_lines[idx].startswith('/'):
        idx += 1

    def next_val():
        nonlocal idx
        v = clean_lines[idx].split()[0]
        idx += 1
        return v

    def next_line():
        nonlocal idx
        v = clean_lines[idx]
        idx += 1
        return v

    def skip_comment():
        nonlocal idx
        # Skip the next line if it starts with '/' (comment)
        if idx < len(clean_lines) and clean_lines[idx].startswith('/'):
            idx += 1

    # Per-spring fields (old format order)
    stressModelFlag = next_val()
    dampingFlag     = next_val()
    frictionFlag    = next_val()
    frameFlag       = next_val()
    BCP_1           = next_val()
    BCP_2           = next_val()
    BCP_1_type      = next_val()
    BCP_2_type      = next_val()

    # 6 spring vectors (3 for BCP1, 3 for BCP2): each is "x y z"
    vectors = []
    for _ in range(6):
        vectors.append(next_line())

    # Skip stiffness matrix comment
    skip_comment()
    K_rows = []
    for _ in range(6):
        K_rows.append(next_line())

    # Skip friction coefficients comment
    skip_comment()
    mu_d = next_val()
    mu_s = next_val()
    vt   = next_val()
    Dt   = next_val()

    # Skip friction matrix comment
    skip_comment()
    M_rows = []
    for _ in range(6):
        M_rows.append(next_line())

    # Skip damping coefficients comment
    skip_comment()
    D_vals = []
    for _ in range(6):
        D_vals.append(next_val())

    # 6 DOF stress-strain blocks
    dof_blocks = []
    for _ in range(6):
        skip_comment()
        n = next_val()
        strain = next_line()
        stress_rows = []
        for _ in range(6):
            stress_rows.append(next_line())
        dof_blocks.append((n, strain, stress_rows))

    type_data = dict(
        stressModelFlag=stressModelFlag,
        dampingFlag=dampingFlag,
        frictionFlag=frictionFlag,
        frameFlag=frameFlag,
        K_rows=K_rows,
        mu_d=mu_d, mu_s=mu_s, vt=vt, Dt=Dt,
        M_rows=M_rows,
        D_vals=D_vals,
        dof_blocks=dof_blocks,
    )
    instance_data = dict(
        BCP_1=BCP_1, BCP_2=BCP_2,
        BCP_1_type=BCP_1_type, BCP_2_type=BCP_2_type,
        vectors=vectors,
    )
    return instance_data, type_data, idx

def make_type_key(td):
    """Hashable key for type deduplication."""
    dof_key = tuple(
        (b[0], b[1], tuple(b[2]))
        for b in td['dof_blocks']
    )
    return (
        td['stressModelFlag'], td['dampingFlag'],
        td['frictionFlag'], td['frameFlag'],
        tuple(td['K_rows']),
        td['mu_d'], td['mu_s'], td['vt'], td['Dt'],
        tuple(td['M_rows']),
        tuple(td['D_vals']),
        dof_key,
    )

def write_type_block(f, td, n):
    f.write(f"//////////////////\n")
    f.write(f"////////////////// New spring type [{n}]\n")
    f.write(f"//////////////////\n")
    f.write(f"{td['stressModelFlag']} // Flag for spring stress model [1: linear (matrix), 2: non-linear (stress-strain)]\n")
    f.write(f"{td['dampingFlag']} // Flag spring damping [0: no, 1: yes]\n")
    f.write(f"{td['frictionFlag']} // Flag spring friction [0: no, 1: velocity-based, 2: stick-slip]\n")
    f.write(f"{td['frameFlag']} // Flag frame for forces/deformations [0: average, 1: BCP1, 2: BCP2]\n")
    f.write(f"// STIFFNESS MATRIX\n")
    for row in td['K_rows']:
        f.write(f"{row} // stiffness row\n")
    f.write(f"// FRICTION COEFFICIENTS\n")
    f.write(f"{td['mu_d']} // mu_d - dynamic friction coefficient\n")
    f.write(f"{td['mu_s']} // mu_s - static friction coefficient\n")
    f.write(f"{td['vt']} // vt - transition velocity threshold [m/s]\n")
    f.write(f"{td['Dt']} // Dt - maximum stick-slip deformation [m]\n")
    f.write(f"// FRICTION MATRIX\n")
    for row in td['M_rows']:
        f.write(f"{row} // friction matrix row\n")
    f.write(f"// DAMPING COEFFICIENTS\n")
    for i, v in enumerate(td['D_vals']):
        f.write(f"{v} // damping DOF {i+1}\n")
    for j, (n_pts, strain, stress_rows) in enumerate(td['dof_blocks']):
        f.write(f"// Stress-strain curve DOF {j+1}\n")
        f.write(f"{n_pts} // Number of points\n")
        f.write(f"{strain} // Strain data\n")
        for k, sr in enumerate(stress_rows):
            f.write(f"{sr} // Stress in DOF {k+1}\n")

def write_instance_block(f, inst, type_idx, n):
    f.write(f"//////////////////\n")
    f.write(f"////////////////// New spring [{n}]\n")
    f.write(f"//////////////////\n")
    f.write(f"{type_idx} // Spring type index (1-based)\n")
    f.write(f"{inst['BCP_1']} // BCP_1 index\n")
    f.write(f"{inst['BCP_2']} // BCP_2 index\n")
    f.write(f"{inst['BCP_1_type']} // BCP_1 type [0: fixed, 1: on line, 2: on plane]\n")
    f.write(f"{inst['BCP_2_type']} // BCP_2 type [0: fixed, 1: on line, 2: on plane]\n")
    vec_labels = [
        "normal vector (x) BCP1", "tangent vector 1 (y) BCP1", "tangent vector 2 (z) BCP1",
        "normal vector (x) BCP2", "tangent vector 1 (y) BCP2", "tangent vector 2 (z) BCP2",
    ]
    for v, lbl in zip(inst['vectors'], vec_labels):
        f.write(f"{v} // {lbl}\n")

def convert_file(filepath):
    clean = read_clean_lines(filepath)
    if not clean:
        print(f"  WARNING: empty file {filepath}")
        return

    idx = 0
    num_springs = int(clean[idx].split()[0])
    idx += 1

    if num_springs == 0:
        with open(filepath, 'w', newline='\n') as f:
            f.write("0 // Number of spring types\n")
            f.write("0 // Number of springs\n")
        print(f"  Converted (0 springs): {filepath}")
        return

    all_instances = []
    all_type_datas = []
    for _ in range(num_springs):
        inst, tdata, idx = parse_spring_entry(clean, idx)
        all_instances.append(inst)
        all_type_datas.append(tdata)

    # Deduplicate types
    type_key_to_idx = {}
    unique_types = []
    instance_type_indices = []
    for td in all_type_datas:
        k = make_type_key(td)
        if k not in type_key_to_idx:
            type_key_to_idx[k] = len(unique_types) + 1
            unique_types.append(td)
        instance_type_indices.append(type_key_to_idx[k])

    with open(filepath, 'w', newline='\n') as f:
        f.write(f"{len(unique_types)} // Number of spring types\n")
        for i, td in enumerate(unique_types):
            write_type_block(f, td, i + 1)
        f.write(f"{num_springs} // Number of springs\n")
        for i, (inst, tidx) in enumerate(zip(all_instances, instance_type_indices)):
            write_instance_block(f, inst, tidx, i + 1)

    ntypes = len(unique_types)
    print(f"  Converted: {filepath}  ({num_springs} spring(s), {ntypes} type(s))")

if __name__ == '__main__':
    workspace = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    pattern = os.path.join(workspace, 'examples', '**', 'dataSprings.dat')
    files = glob.glob(pattern, recursive=True)
    if not files:
        print("No dataSprings.dat files found.")
        sys.exit(1)
    print(f"Converting {len(files)} dataSprings.dat files...")
    for fp in sorted(files):
        try:
            convert_file(fp)
        except Exception as e:
            import traceback
            print(f"  ERROR in {fp}: {e}")
            traceback.print_exc()
    print("Done.")
