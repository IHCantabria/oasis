"""OASIS case data model — Python dataclasses mirroring the dataProblem.yaml structure.

YAML key names are preserved as field names where possible, with Python-friendly
aliases noted in comments when they differ.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional


# ──────────────────────────────────────────────────────────────────────────────
# Problem settings
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class ProblemData:
    gravity: float = 9.81
    water_density: float = 1025.0
    air_density: float = 1.225
    atmospheric_pressure: float = 101325.0
    air_adiabatic_dilation: float = 1.4
    water_depth: float = -100.0
    write_time_step: float = 0.1
    max_time_step: float = 0.1
    hydro_time_step: float = 0.1
    fast_time_step: float = 0.1
    fast_controller_time_step: float = 0.1
    irf_time: float = 20.0
    sinking_time_step: float = 15.0
    winches_controller_time_step: float = 0.1
    owcs_controller_time_step: float = 0.1
    simulation_time: float = 100.0
    rotation_simplification: int = 1
    time_integration_method: int = 3   # 1=BDF1, 2=BDFN, 3=ESDIRK46
    time_integration_order: int = 2
    time_step_adaptivity: int = 0
    jacobian_recomputation_steps: int = 0
    absolute_tolerance: float = 1e-5
    relative_tolerance: float = 1e-3
    max_iterations_per_step: int = 20
    read_equilibrium: int = 0
    write_equilibrium: int = 0
    mooring_initial_condition: int = 1  # 0=catenary, 1=Newton
    output_format: str = "csv"          # "csv" or "txt"


# ──────────────────────────────────────────────────────────────────────────────
# Sea floor
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class FlatFloorData:
    depth: float = -100.0


@dataclass
class InclinedFloorData:
    point1: List[float] = field(default_factory=lambda: [0.0, 0.0, -100.0])
    point2: List[float] = field(default_factory=lambda: [50.0, 0.0, -100.0])
    point3: List[float] = field(default_factory=lambda: [0.0, 50.0, -100.0])


@dataclass
class BathymetryFloorData:
    mesh_file: str = ""


@dataclass
class SeaFloorData:
    flat: List[FlatFloorData] = field(default_factory=list)
    inclined: List[InclinedFloorData] = field(default_factory=list)
    bathymetry: List[BathymetryFloorData] = field(default_factory=list)


# ──────────────────────────────────────────────────────────────────────────────
# Waves
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class WaveData:
    # YAML key: "type"
    wave_type: str = "REG"         # "REG" or "IRR"
    height: float = 0.0
    period: float = 10.0
    heading: float = 0.0
    ramp_time: float = 0.0
    # Irregular only — YAML keys differ from field names (noted below)
    spectrum_type: int = 1          # YAML: spectrum_type  (1=JONSWAP, 2=timeseries, 3=freq-domain)
    piecewise_flag: int = 0         # YAML: piecewise_flag
    gamma: float = 3.3
    spreading: float = 0.0          # YAML: spreading
    dtheta: float = 10.0
    factor: float = 0.001
    read_phases_flag: int = 0       # YAML: read_phases_flag
    wave_relative_tolerance: float = 0.05  # YAML: relative_tolerance (wave-specific)
    dt: float = 0.1
    phases_file: str = "WavePhases.dat"    # YAML: phases_file
    database_file: str = "unused.dat"      # YAML: database_file


# ──────────────────────────────────────────────────────────────────────────────
# BCPs
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class AnchorBCPData:
    position: List[float] = field(default_factory=lambda: [0.0, 0.0, -100.0])
    winch_id: int = 0


@dataclass
class FairleadBCPData:
    position: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    winch_id: int = 0
    actuator_file: str = ""


@dataclass
class JointBCPData:
    position: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    winch_id: int = 0
    mass: float = 0.0
    volume: float = 0.0


@dataclass
class BodyBCPData:
    position: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    winch_id: int = 0


@dataclass
class ElasticAnchorBCPData:
    position: List[float] = field(default_factory=lambda: [0.0, 0.0, -100.0])
    winch_id: int = 0
    anchor_mass: float = 0.0
    anchor_volume: float = 0.0
    c_param: float = 0.0
    k_param: float = 0.0


@dataclass
class BCPsData:
    # YAML allocation order: fairleads → anchors → joints → body_bcps → elastic_anchors
    fairleads: List[FairleadBCPData] = field(default_factory=list)
    anchors: List[AnchorBCPData] = field(default_factory=list)
    joints: List[JointBCPData] = field(default_factory=list)
    body_bcps: List[BodyBCPData] = field(default_factory=list)
    elastic_anchors: List[ElasticAnchorBCPData] = field(default_factory=list)

    def all_positions(self):
        """Returns list of (global_1based_index, type_str, [x,y,z]) for all BCPs."""
        result = []
        g = 1
        for b in self.fairleads:
            result.append((g, "fairlead", list(b.position))); g += 1
        for b in self.anchors:
            result.append((g, "anchor", list(b.position))); g += 1
        for b in self.joints:
            result.append((g, "joint", list(b.position))); g += 1
        for b in self.body_bcps:
            result.append((g, "body_bcp", list(b.position))); g += 1
        for b in self.elastic_anchors:
            result.append((g, "elastic_anchor", list(b.position))); g += 1
        return result

    def total(self):
        return (len(self.fairleads) + len(self.anchors) + len(self.joints)
                + len(self.body_bcps) + len(self.elastic_anchors))


# ──────────────────────────────────────────────────────────────────────────────
# Bodies
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class BodyData:
    body_type: str = "RAD_DIFF"             # YAML: type
    take_cog_from_hdb: int = 1              # YAML: take_cog_hydro_database
    dofs: List[int] = field(default_factory=lambda: [1, 2, 3, 4, 5, 6])
    bcps_indexes: List[int] = field(default_factory=list)
    wind_turbines_indexes: List[int] = field(default_factory=list)
    initial_position: List[float] = field(default_factory=lambda: [0.0] * 6)
    initial_displacement: List[float] = field(default_factory=lambda: [0.0] * 6)
    hdb_file: str = ""                      # YAML: hydro_database
    hdb_index: int = 1                      # YAML: hydro_database_index
    freedom_flag: int = 0                   # 0=free, 1=locked, 2=imposed
    imposed_motion_file: str = "dummy.dat"
    hydrostatics_flag: int = 0              # 0=linear, 1=nonlinear flat, 2=nonlinear+waves
    hydrostatics_mesh: str = "dummy.stl"
    radiation: int = 1                      # YAML: radiation_flag
    excitation_1st: int = 1                 # YAML: excitation_1st_flag
    excitation_2nd: int = 0                 # YAML: excitation_2nd_flag
    viscous_added_mass: List[float] = field(default_factory=lambda: [0.0] * 6)
    viscous_linear_damping: List[float] = field(default_factory=lambda: [0.0] * 6)
    viscous_quadratic_damping: List[float] = field(default_factory=lambda: [0.0] * 6)
    # GUI-only (not written to YAML)
    name: str = "Body"
    schematic_width: float = 10.0
    schematic_height: float = 5.0
    stl_3d_file: str = ""
    cog: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    mass_file: str = ""


# ──────────────────────────────────────────────────────────────────────────────
# Lines
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class LineTypeData:
    flag_tension: int = 2         # 1=symmetric, 2=tension-only
    density: float = 50.0         # kg/m
    diameter: float = 0.09        # m
    flag_stiffness: int = 0       # 0=constant EA, 1=viscoelastic, >1=tabulated
    # Mode 0
    EA: float = 1e7
    beta: float = 0.002
    # Mode 1 (viscoelastic)
    kernel_lin_coef: List[float] = field(default_factory=list)
    kernel_exp_coef: List[float] = field(default_factory=list)
    elastic_coef: List[float] = field(default_factory=list)
    # Mode >1 (tabulated strain-stress)
    strain_data: List[float] = field(default_factory=list)
    stress_data: List[float] = field(default_factory=list)
    # Hydrodynamic / seabed
    CB: float = 0.0
    Cmn: float = 3.8
    Cdn: float = 2.5
    Cdt: float = 0.5
    GK: float = 1e6
    GC: float = 0.1
    smoothstep: int = 1
    friction_model: int = 0
    vth: float = 0.01
    ust: float = 0.6
    usn: float = 0.4
    ud: float = 0.3
    deltamax: float = 0.03


@dataclass
class LineData:
    line_type: int = 1          # 1=mooring, 2=towing, 3=tensor
    type_index: int = 1         # references line_types[type_index-1]
    num_nodes: int = 11
    polynomial_order: int = 1
    length: float = 120.0
    seafloor_index: int = 1
    BCP_N: int = 1              # end-node BCP (global 1-based)
    BCP_1: int = 2              # start-node BCP (global 1-based)


# ──────────────────────────────────────────────────────────────────────────────
# Springs (type-instance pattern)
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class StressStrainCurve:
    displacements: List[float] = field(default_factory=lambda: [-1.0, 0.0, 1.0])
    # forces[i] = list of force values for DOF i+1 at each displacement point
    forces: List[List[float]] = field(
        default_factory=lambda: [[0.0, 0.0, 0.0] for _ in range(6)]
    )


@dataclass
class SpringTypeData:
    stress_model_flag: int = 1      # 1=linear matrix, 2=nonlinear stress-strain
    damping_flag: int = 0
    friction_flag: int = 0
    frame_flag: int = 0
    stiffness_matrix: List[List[float]] = field(
        default_factory=lambda: [[0.0] * 6 for _ in range(6)]
    )
    mu_d: float = 0.0
    mu_s: float = 0.0
    vt: float = 0.001
    Dt: float = 0.01
    mass_matrix: List[List[float]] = field(
        default_factory=lambda: [[0.0] * 6 for _ in range(6)]
    )
    damping_vector: List[float] = field(default_factory=lambda: [0.0] * 6)
    # 6 curves, one per DOF
    stress_strain: List[StressStrainCurve] = field(
        default_factory=lambda: [StressStrainCurve() for _ in range(6)]
    )


@dataclass
class SpringData:
    type_index: int = 1
    BCP_1: int = 1
    BCP_2: int = 2
    BCP_1_type: int = 0     # 0=fixed point, 1=on line, 2=on plane
    BCP_2_type: int = 0
    vectors: List[List[float]] = field(default_factory=lambda: [
        [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
    ])


# ──────────────────────────────────────────────────────────────────────────────
# Morison
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class MorisonData:
    body_index: int = 1
    flow_type: int = 1          # 1=constant, 2=variable HDF5
    wind_speed: float = 0.0
    wind_direction: float = 0.0
    current_speed: float = 0.0
    current_direction: float = 0.0
    symmetry_order: int = 2
    wind_fk_x: List[List[float]] = field(
        default_factory=lambda: [[0.0, 0.0] for _ in range(6)]
    )
    wind_drag_x: List[List[float]] = field(
        default_factory=lambda: [[0.0, 0.0] for _ in range(6)]
    )
    wind_fk_y: List[List[float]] = field(
        default_factory=lambda: [[0.0, 0.0] for _ in range(6)]
    )
    wind_drag_y: List[List[float]] = field(
        default_factory=lambda: [[0.0, 0.0] for _ in range(6)]
    )
    current_fk_x: List[List[float]] = field(
        default_factory=lambda: [[0.0, 0.0] for _ in range(6)]
    )
    current_drag_x: List[List[float]] = field(
        default_factory=lambda: [[0.0, 0.0] for _ in range(6)]
    )
    current_fk_y: List[List[float]] = field(
        default_factory=lambda: [[0.0, 0.0] for _ in range(6)]
    )
    current_drag_y: List[List[float]] = field(
        default_factory=lambda: [[0.0, 0.0] for _ in range(6)]
    )
    hdf5_file: str = ""     # only for flow_type=2


# ──────────────────────────────────────────────────────────────────────────────
# Winches
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class WinchItemData:
    line: int = 1
    line_bcp: int = 2       # 1=first node, 2=last node
    inertia: float = 50.0   # kg·m²
    radius: float = 0.25    # m
    drag: float = 0.5       # N·s/rad


@dataclass
class WinchControllerData:
    controller_type: int = 1        # YAML: type  (1=constant tension, 2=horizontal)
    target_tension: float = 30000.0  # for type 1 only


@dataclass
class WinchesData:
    winches: List[WinchItemData] = field(default_factory=list)
    controller: WinchControllerData = field(default_factory=WinchControllerData)


# ──────────────────────────────────────────────────────────────────────────────
# OWC
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class OWCData:
    body_index: int = 2
    floater_index: int = 1
    waterplane_area: float = 1.0
    ref_air_volume: float = 1.0
    hole_area: float = 0.1
    discharge_coefficient: float = 0.6
    turbine_type: int = -1          # -1=none, 0=hole, >0=turbine id
    initial_omega: float = 0.0
    position: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])


# ──────────────────────────────────────────────────────────────────────────────
# Sinking
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class CompartmentGroupData:
    polygon_x: List[float] = field(default_factory=lambda: [0.0, 1.0, 1.0, 0.0])
    polygon_y: List[float] = field(default_factory=lambda: [0.0, 0.0, 1.0, 1.0])
    floor_z: float = -5.0
    filling_times: List[float] = field(default_factory=lambda: [0.0, 100.0])
    filling_states: List[float] = field(default_factory=lambda: [0.0, 0.0])


@dataclass
class SinkingData:
    body_index: int = 1
    structural_mass_diagonal: List[float] = field(
        default_factory=lambda: [1.0] * 6
    )
    interp_masses: List[float] = field(default_factory=lambda: [0.0])
    hdb_files: List[str] = field(default_factory=list)
    groups: List[CompartmentGroupData] = field(default_factory=list)


# ──────────────────────────────────────────────────────────────────────────────
# Wind Turbines
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class WindTurbineData:
    body_index: int = 1
    fast_input_file: str = ""
    fast_library: str = ""


# ──────────────────────────────────────────────────────────────────────────────
# Root case
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class CaseData:
    problem: ProblemData = field(default_factory=ProblemData)
    seafloor: SeaFloorData = field(default_factory=SeaFloorData)
    waves: WaveData = field(default_factory=WaveData)
    bcps: BCPsData = field(default_factory=BCPsData)
    bodies: List[BodyData] = field(default_factory=list)
    line_types: List[LineTypeData] = field(default_factory=list)
    lines: List[LineData] = field(default_factory=list)
    spring_types: List[SpringTypeData] = field(default_factory=list)
    springs: List[SpringData] = field(default_factory=list)
    morison: List[MorisonData] = field(default_factory=list)
    winches: Optional[WinchesData] = None
    owcs: List[OWCData] = field(default_factory=list)
    sinking: List[SinkingData] = field(default_factory=list)
    wind_turbines: List[WindTurbineData] = field(default_factory=list)
    # Metadata (not written to YAML)
    project_path: str = ""
