"""
Converts dataLines.dat files from the old format (all properties per-line)
to the new format (line types section + line instances section).
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

def parse_line_entry(clean_lines, start_idx):
    """
    Parse one line entry starting at start_idx.
    Skips any lines that start with '/' (header lines).
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

    line_type     = next_val()
    flag_tension  = next_val()
    nNodos        = next_val()
    p             = next_val()
    L             = next_val()
    rho0          = next_val()
    d             = next_val()
    flag_stiffness_str = next_val()
    fs = int(flag_stiffness_str)

    stiffness_lines = []
    if fs == 0:
        stiffness_lines.append(next_val())  # EA
        stiffness_lines.append(next_val())  # beta
    elif fs == 1:
        stiffness_lines.append(next_line())  # kernel_lin_coef
        stiffness_lines.append(next_line())  # kernel_exp_coef
        stiffness_lines.append(next_line())  # elastic_coef
    else:  # fs > 1
        stiffness_lines.append(next_line())  # strain_data
        stiffness_lines.append(next_line())  # stress_data
        stiffness_lines.append(next_val())   # beta

    CB            = next_val()
    Cmn           = next_val()
    Cdn           = next_val()
    Cdt           = next_val()
    GK            = next_val()
    GC            = next_val()
    indexSeaFloor = next_val()
    BCP_N         = next_val()
    BCP_1         = next_val()
    smoothstep    = next_val()
    frictionModel = next_val()
    vth           = next_val()
    ust           = next_val()
    usn           = next_val()
    ud            = next_val()
    deltamax      = next_val()

    type_data = dict(
        flag_tension=flag_tension, rho0=rho0, d=d,
        flag_stiffness=flag_stiffness_str, fs=fs,
        stiffness_lines=stiffness_lines,
        CB=CB, Cmn=Cmn, Cdn=Cdn, Cdt=Cdt, GK=GK, GC=GC,
        smoothstep=smoothstep, frictionModel=frictionModel,
        vth=vth, ust=ust, usn=usn, ud=ud, deltamax=deltamax,
    )
    instance_data = dict(
        lineType=line_type, nNodos=nNodos, p=p, L=L,
        indexSeaFloor=indexSeaFloor, BCP_N=BCP_N, BCP_1=BCP_1,
    )
    return instance_data, type_data, idx

def make_type_key(td):
    """Hashable key for type deduplication."""
    fs = td['fs']
    if fs == 0:
        stiff = (float(td['stiffness_lines'][0]), float(td['stiffness_lines'][1]))
    elif fs == 1:
        kl = tuple(float(x) for x in td['stiffness_lines'][0].split())
        ke = tuple(float(x) for x in td['stiffness_lines'][1].split())
        ec = tuple(float(x) for x in td['stiffness_lines'][2].split())
        stiff = (kl, ke, ec)
    else:
        st  = tuple(float(x) for x in td['stiffness_lines'][0].split())
        ss  = tuple(float(x) for x in td['stiffness_lines'][1].split())
        bt  = float(td['stiffness_lines'][2])
        stiff = (st, ss, bt)
    return (
        float(td['flag_tension']), float(td['rho0']), float(td['d']), fs,
        stiff,
        float(td['CB']), float(td['Cmn']), float(td['Cdn']), float(td['Cdt']),
        float(td['GK']), float(td['GC']),
        int(td['smoothstep']), int(td['frictionModel']),
        float(td['vth']), float(td['ust']), float(td['usn']),
        float(td['ud']), float(td['deltamax']),
    )

def write_type_block(f, td, n):
    fs = td['fs']
    f.write(f"//////////////////\n")
    f.write(f"////////////////// New line type [{n}]\n")
    f.write(f"//////////////////\n")
    f.write(f"{td['flag_tension']} // Tension model [1: symmetric, 2: tension-only]\n")
    f.write(f"{td['rho0']} // Weight over water [kg/m]\n")
    f.write(f"{td['d']} // Diameter [m]\n")
    f.write(f"{td['flag_stiffness']} // Flag for axial stiffness [0: constant EA, 1: viscoelastic, >1: tabulated]\n")
    if fs == 0:
        f.write(f"{td['stiffness_lines'][0]} // EA - axial stiffness [N]\n")
        f.write(f"{td['stiffness_lines'][1]} // beta - structural damping coefficient\n")
    elif fs == 1:
        f.write(f"{td['stiffness_lines'][0]} // Kernel linear coefficients\n")
        f.write(f"{td['stiffness_lines'][1]} // Kernel exponential coefficients\n")
        f.write(f"{td['stiffness_lines'][2]} // Elastic polynomial coefficients\n")
    else:
        f.write(f"{td['stiffness_lines'][0]} // Strain data ({fs} points)\n")
        f.write(f"{td['stiffness_lines'][1]} // Stress data ({fs} points) [Pa]\n")
        f.write(f"{td['stiffness_lines'][2]} // beta - structural damping coefficient\n")
    f.write(f"{td['CB']} // CB - seabed static friction coefficient\n")
    f.write(f"{td['Cmn']} // Cmn - hydrodynamic mass coefficient\n")
    f.write(f"{td['Cdn']} // Cdn - normal drag coefficient\n")
    f.write(f"{td['Cdt']} // Cdt - tangential drag coefficient\n")
    f.write(f"{td['GK']} // GK - ground normal stiffness [N/m]\n")
    f.write(f"{td['GC']} // GC - fraction of critical damping of ground\n")
    f.write(f"{td['smoothstep']} // Smoothstep flag\n")
    f.write(f"{td['frictionModel']} // Friction model [0: None, 1: Isotropic, 2: Anisotropic]\n")
    f.write(f"{td['vth']} // vth - velocity threshold for stick-to-slip [m/s]\n")
    f.write(f"{td['ust']} // ust - static tangential friction coefficient\n")
    f.write(f"{td['usn']} // usn - static normal friction coefficient\n")
    f.write(f"{td['ud']} // ud - dynamic friction coefficient\n")
    f.write(f"{td['deltamax']} // deltamax - maximum displacement [m]\n")

def write_instance_block(f, inst, type_idx, n):
    f.write(f"//////////////////\n")
    f.write(f"////////////////// New line [{n}]\n")
    f.write(f"//////////////////\n")
    f.write(f"{inst['lineType']} // Type of line [1: Mooring, 2: Towing, 3: Tensor]\n")
    f.write(f"{type_idx} // Line type index (1-based)\n")
    f.write(f"{inst['nNodos']} // Number of nodes\n")
    f.write(f"{inst['p']} // Polynomial order\n")
    f.write(f"{inst['L']} // Length [m]\n")
    f.write(f"{inst['indexSeaFloor']} // Floor Id [from dataSeaFloor.dat]\n")
    f.write(f"{inst['BCP_N']} // BCP index of node N\n")
    f.write(f"{inst['BCP_1']} // BCP index of node 1\n")

def convert_file(filepath):
    clean = read_clean_lines(filepath)
    if not clean:
        print(f"  WARNING: empty file {filepath}")
        return

    idx = 0
    num_lines = int(clean[idx].split()[0])
    idx += 1

    if num_lines == 0:
        with open(filepath, 'w', newline='\n') as f:
            f.write("0 // Number of line types\n")
            f.write("0 // Number of lines\n")
        print(f"  Converted (0 lines): {filepath}")
        return

    # Parse all line entries
    all_instances = []
    all_type_datas = []
    for _ in range(num_lines):
        inst, tdata, idx = parse_line_entry(clean, idx)
        all_instances.append(inst)
        all_type_datas.append(tdata)

    # Deduplicate types
    type_key_to_idx = {}  # key -> 1-based index
    unique_types = []     # list of type_data in order of first appearance
    instance_type_indices = []
    for td in all_type_datas:
        k = make_type_key(td)
        if k not in type_key_to_idx:
            type_key_to_idx[k] = len(unique_types) + 1
            unique_types.append(td)
        instance_type_indices.append(type_key_to_idx[k])

    # Write new file
    with open(filepath, 'w', newline='\n') as f:
        f.write(f"{len(unique_types)} // Number of line types\n")
        for i, td in enumerate(unique_types):
            write_type_block(f, td, i + 1)

        f.write(f"{num_lines} // Number of lines\n")
        for i, (inst, tidx) in enumerate(zip(all_instances, instance_type_indices)):
            write_instance_block(f, inst, tidx, i + 1)

    ntypes = len(unique_types)
    print(f"  Converted: {filepath}  ({num_lines} lines, {ntypes} type(s))")

if __name__ == '__main__':
    workspace = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    pattern = os.path.join(workspace, 'examples', '**', 'dataLines.dat')
    files = glob.glob(pattern, recursive=True)
    if not files:
        print("No dataLines.dat files found.")
        sys.exit(1)
    print(f"Converting {len(files)} dataLines.dat files...")
    for fp in sorted(files):
        try:
            convert_file(fp)
        except Exception as e:
            print(f"  ERROR in {fp}: {e}")
    print("Done.")
