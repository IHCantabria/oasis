"""
Converts dataProblem.yaml files that have 'lines' or 'springs' sections
from the old format (all properties per-instance) to the new format
(line_types/spring_types section + instances with type_index).
"""
import os
import sys
import glob

try:
    import yaml
except ImportError:
    print("PyYAML not installed. Run: pip install pyyaml")
    sys.exit(1)

# ── YAML key lists ────────────────────────────────────────────────────────────

LINE_TYPE_KEYS_BASE = [
    'flag_tension', 'density', 'diameter', 'flag_stiffness',
]
LINE_TYPE_KEYS_FS0 = ['EA', 'beta']
LINE_TYPE_KEYS_FS1 = ['kernel_lin_coef', 'kernel_exp_coef', 'elastic_coef']
LINE_TYPE_KEYS_TAIL = [
    'CB', 'Cmn', 'Cdn', 'Cdt', 'GK', 'GC',
    'smoothstep', 'friction_model',
    'vth', 'ust', 'usn', 'ud', 'deltamax',
]
LINE_INSTANCE_KEYS = [
    'line_type', 'type_index',
    'num_nodes', 'polynomial_order', 'length',
    'seafloor_index', 'BCP_N', 'BCP_1',
]

SPRING_TYPE_KEYS = [
    'stress_model_flag', 'damping_flag', 'friction_flag', 'frame_flag',
    'stiffness_matrix', 'mu_d', 'mu_s', 'vt', 'Dt',
    'mass_matrix', 'damping_vector', 'stress_strain',
]
SPRING_INSTANCE_KEYS = [
    'type_index', 'BCP_1', 'BCP_2', 'BCP_1_type', 'BCP_2_type', 'vectors',
]


def extract_line_type(entry):
    """Extract type fields from a line entry dict."""
    td = {}
    for k in LINE_TYPE_KEYS_BASE:
        td[k] = entry[k]
    fs = entry['flag_stiffness']
    if fs == 0:
        for k in LINE_TYPE_KEYS_FS0:
            td[k] = entry[k]
    elif fs == 1:
        for k in LINE_TYPE_KEYS_FS1:
            td[k] = entry[k]
    else:
        td['strain_data'] = entry['strain_data']
        td['stress_data'] = entry['stress_data']
        td['beta'] = entry['beta']
    for k in LINE_TYPE_KEYS_TAIL:
        td[k] = entry[k]
    return td


def make_line_type_key(td):
    """Hashable key for line type deduplication."""
    fs = td['flag_stiffness']
    if fs == 0:
        stiff = (td.get('EA'), td.get('beta'))
    elif fs == 1:
        stiff = (
            tuple(td.get('kernel_lin_coef', [])),
            tuple(td.get('kernel_exp_coef', [])),
            tuple(td.get('elastic_coef', [])),
        )
    else:
        stiff = (tuple(td.get('strain_data', [])), tuple(td.get('stress_data', [])), td.get('beta'))
    return (
        td['flag_tension'], td['density'], td['diameter'], fs, stiff,
        td['CB'], td['Cmn'], td['Cdn'], td['Cdt'], td['GK'], td['GC'],
        td['smoothstep'], td['friction_model'],
        td['vth'], td['ust'], td['usn'], td['ud'], td['deltamax'],
    )


def build_line_instance(entry, type_idx):
    """Build instance dict from line entry."""
    return {
        'line_type': entry['line_type'],
        'type_index': type_idx,
        'num_nodes': entry['num_nodes'],
        'polynomial_order': entry['polynomial_order'],
        'length': entry['length'],
        'seafloor_index': entry['seafloor_index'],
        'BCP_N': entry['BCP_N'],
        'BCP_1': entry['BCP_1'],
    }


def convert_lines(lines_list):
    """
    Convert a list of line entries (old format) to
    (line_types list, new_lines list).
    """
    type_key_to_idx = {}
    unique_types = []
    new_lines = []

    for entry in lines_list:
        td = extract_line_type(entry)
        k = make_line_type_key(td)
        if k not in type_key_to_idx:
            type_key_to_idx[k] = len(unique_types) + 1
            unique_types.append(td)
        tidx = type_key_to_idx[k]
        new_lines.append(build_line_instance(entry, tidx))

    return unique_types, new_lines


def extract_spring_type(entry):
    """Extract type fields from a spring entry dict."""
    return {k: entry[k] for k in SPRING_TYPE_KEYS if k in entry}


def make_spring_type_key(td):
    """Hashable key for spring type deduplication (simple string repr)."""
    import json
    return json.dumps(td, sort_keys=True, default=str)


def build_spring_instance(entry, type_idx):
    """Build instance dict from spring entry."""
    inst = {'type_index': type_idx}
    for k in ['BCP_1', 'BCP_2', 'BCP_1_type', 'BCP_2_type', 'vectors']:
        if k in entry:
            inst[k] = entry[k]
    return inst


def convert_springs(springs_list):
    """
    Convert a list of spring entries (old format) to
    (spring_types list, new_springs list).
    """
    type_key_to_idx = {}
    unique_types = []
    new_springs = []

    for entry in springs_list:
        td = extract_spring_type(entry)
        k = make_spring_type_key(td)
        if k not in type_key_to_idx:
            type_key_to_idx[k] = len(unique_types) + 1
            unique_types.append(td)
        tidx = type_key_to_idx[k]
        new_springs.append(build_spring_instance(entry, tidx))

    return unique_types, new_springs


# ── Representer for nice float formatting ────────────────────────────────────

class FloatDumper(yaml.Dumper):
    pass

def float_representer(dumper, value):
    # Format floats cleanly (avoid scientific notation for small values if not needed)
    text = repr(value)
    return dumper.represent_scalar('tag:yaml.org,2002:float', text)

FloatDumper.add_representer(float, float_representer)


def convert_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        data = yaml.safe_load(f)

    if data is None:
        print(f"  SKIP (empty): {filepath}")
        return

    changed = False

    if 'lines' in data and data['lines']:
        unique_types, new_lines = convert_lines(data['lines'])
        data['line_types'] = unique_types
        data['lines'] = new_lines
        changed = True
        print(f"  Lines:   {len(new_lines)} instance(s), {len(unique_types)} type(s)  in {filepath}")

    if 'springs' in data and data['springs']:
        unique_types, new_springs = convert_springs(data['springs'])
        data['spring_types'] = unique_types
        data['springs'] = new_springs
        changed = True
        print(f"  Springs: {len(new_springs)} instance(s), {len(unique_types)} type(s) in {filepath}")

    if not changed:
        print(f"  SKIP (no lines/springs): {filepath}")
        return

    # Re-order top-level keys to put line_types before lines, spring_types before springs
    ordered = {}
    for k in data:
        if k not in ('line_types', 'lines', 'spring_types', 'springs'):
            ordered[k] = data[k]
    if 'line_types' in data:
        ordered['line_types'] = data['line_types']
    if 'lines' in data:
        ordered['lines'] = data['lines']
    if 'spring_types' in data:
        ordered['spring_types'] = data['spring_types']
    if 'springs' in data:
        ordered['springs'] = data['springs']

    with open(filepath, 'w', encoding='utf-8', newline='\n') as f:
        yaml.dump(ordered, f, Dumper=FloatDumper,
                  default_flow_style=False, allow_unicode=True, sort_keys=False)


if __name__ == '__main__':
    workspace = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    pattern = os.path.join(workspace, 'examples', '**', 'dataProblem.yaml')
    files = glob.glob(pattern, recursive=True)
    if not files:
        print("No dataProblem.yaml files found.")
        sys.exit(1)
    print(f"Processing {len(files)} dataProblem.yaml files...")
    for fp in sorted(files):
        try:
            convert_file(fp)
        except Exception as e:
            import traceback
            print(f"  ERROR in {fp}: {e}")
            traceback.print_exc()
    print("Done.")
