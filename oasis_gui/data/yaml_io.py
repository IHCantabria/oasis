"""YAML I/O for OASIS case data.

load_yaml(path) -> CaseData
save_yaml(case, path) -> None
"""
from __future__ import annotations
import os
from typing import Any

from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq

from .case_data import (
    CaseData, ProblemData, WaveData, SeaFloorData, FlatFloorData,
    InclinedFloorData, BathymetryFloorData, BCPsData, AnchorBCPData,
    FairleadBCPData, JointBCPData, BodyBCPData, ElasticAnchorBCPData,
    BodyData, LineTypeData, LineData, SpringTypeData, SpringData,
    StressStrainCurve, MorisonData, WinchesData, WinchItemData,
    WinchControllerData, OWCData, SinkingData, CompartmentGroupData,
    WindTurbineData,
)

_yaml = YAML()
_yaml.default_flow_style = False
_yaml.width = 120


def load_yaml(path: str) -> CaseData:
    with open(path, "r", encoding="utf-8") as fh:
        raw = _yaml.load(fh)
    if raw is None:
        return CaseData()
    case = _parse_case(raw)
    case.project_path = os.path.dirname(os.path.abspath(path))
    return case


def save_yaml(case: CaseData, path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    data = _build_case(case)
    with open(path, "w", encoding="utf-8") as fh:
        _yaml.dump(data, fh)


# ──────────────────────────────────────────────────────────────────────────────
# Helper
# ──────────────────────────────────────────────────────────────────────────────
def _g(d: Any, key: str, default: Any) -> Any:
    """Get value from mapping with a default if key is absent or None."""
    if d is None:
        return default
    v = d.get(key, default)
    return default if v is None else v


def _flist(d: Any, key: str, default: list) -> list:
    v = _g(d, key, default)
    return list(v) if v is not None else list(default)


def _fseq(lst) -> CommentedSeq:
    """Return a ruamel CommentedSeq with inline flow style: [a, b, c]."""
    cs = CommentedSeq(list(lst))
    cs.fa.set_flow_style()
    return cs


# ──────────────────────────────────────────────────────────────────────────────
# PARSE (YAML → dataclasses)
# ──────────────────────────────────────────────────────────────────────────────
def _parse_case(d: Any) -> CaseData:
    case = CaseData()
    if "problem" in d:
        case.problem = _parse_problem(d["problem"])
    if "waves" in d:
        case.waves = _parse_waves(d["waves"])
    if "seafloor" in d:
        case.seafloor = _parse_seafloor(d["seafloor"])
    if "bcps" in d:
        case.bcps = _parse_bcps(d["bcps"])
    if "bodies" in d and d["bodies"]:
        case.bodies = [_parse_body(b) for b in d["bodies"]]
    if "line_types" in d and d["line_types"]:
        case.line_types = [_parse_line_type(lt) for lt in d["line_types"]]
    if "lines" in d and d["lines"]:
        case.lines = [_parse_line(ln) for ln in d["lines"]]
    if "spring_types" in d and d["spring_types"]:
        case.spring_types = [_parse_spring_type(st) for st in d["spring_types"]]
    if "springs" in d and d["springs"]:
        case.springs = [_parse_spring(s) for s in d["springs"]]
    if "morison" in d and d["morison"]:
        case.morison = [_parse_morison(m) for m in d["morison"]]
    if "winches" in d and d["winches"]:
        case.winches = _parse_winches(d["winches"])
    if "owcs" in d and d["owcs"]:
        case.owcs = [_parse_owc(o) for o in d["owcs"]]
    if "sinking" in d and d["sinking"]:
        case.sinking = [_parse_sinking(s) for s in d["sinking"]]
    if "wind_turbines" in d and d["wind_turbines"]:
        case.wind_turbines = [_parse_turbine(t) for t in d["wind_turbines"]]
    return case


def _parse_problem(d) -> ProblemData:
    p = ProblemData()
    for attr in vars(p):
        yaml_key = attr  # keys match field names for problem
        if yaml_key in d:
            setattr(p, attr, d[yaml_key])
    return p


def _parse_waves(d) -> WaveData:
    w = WaveData()
    w.wave_type = _g(d, "type", w.wave_type)
    w.height = _g(d, "height", w.height)
    w.period = _g(d, "period", w.period)
    w.heading = _g(d, "heading", w.heading)
    w.ramp_time = _g(d, "ramp_time", w.ramp_time)
    w.spectrum_type = _g(d, "spectrum_type", w.spectrum_type)
    w.piecewise_flag = _g(d, "piecewise_flag", w.piecewise_flag)
    w.gamma = _g(d, "gamma", w.gamma)
    w.spreading = _g(d, "spreading", w.spreading)
    w.dtheta = _g(d, "dtheta", w.dtheta)
    w.factor = _g(d, "factor", w.factor)
    w.read_phases_flag = _g(d, "read_phases_flag", w.read_phases_flag)
    w.wave_relative_tolerance = _g(d, "relative_tolerance", w.wave_relative_tolerance)
    w.dt = _g(d, "dt", w.dt)
    w.phases_file = _g(d, "phases_file", w.phases_file)
    w.database_file = _g(d, "database_file", w.database_file)
    return w


def _parse_seafloor(d) -> SeaFloorData:
    sf = SeaFloorData()
    if "flat" in d and d["flat"]:
        sf.flat = [FlatFloorData(depth=_g(item, "depth", -100.0)) for item in d["flat"]]
    if "inclined" in d and d["inclined"]:
        for item in d["inclined"]:
            sf.inclined.append(InclinedFloorData(
                point1=_flist(item, "point1", [0.0, 0.0, -100.0]),
                point2=_flist(item, "point2", [50.0, 0.0, -100.0]),
                point3=_flist(item, "point3", [0.0, 50.0, -100.0]),
            ))
    if "bathymetry" in d and d["bathymetry"]:
        sf.bathymetry = [BathymetryFloorData(mesh_file=_g(item, "mesh_file", ""))
                         for item in d["bathymetry"]]
    return sf


def _parse_bcps(d) -> BCPsData:
    bcps = BCPsData()
    if "fairleads" in d and d["fairleads"]:
        for item in d["fairleads"]:
            bcps.fairleads.append(FairleadBCPData(
                position=_flist(item, "position", [0.0, 0.0, 0.0]),
                winch_id=_g(item, "winch_id", 0),
                actuator_file=_g(item, "actuator_file", ""),
            ))
    if "anchors" in d and d["anchors"]:
        for item in d["anchors"]:
            bcps.anchors.append(AnchorBCPData(
                position=_flist(item, "position", [0.0, 0.0, -100.0]),
                winch_id=_g(item, "winch_id", 0),
            ))
    if "joints" in d and d["joints"]:
        for item in d["joints"]:
            bcps.joints.append(JointBCPData(
                position=_flist(item, "position", [0.0, 0.0, 0.0]),
                winch_id=_g(item, "winch_id", 0),
                mass=_g(item, "mass", 0.0),
                volume=_g(item, "volume", 0.0),
            ))
    if "body_bcps" in d and d["body_bcps"]:
        for item in d["body_bcps"]:
            bcps.body_bcps.append(BodyBCPData(
                position=_flist(item, "position", [0.0, 0.0, 0.0]),
                winch_id=_g(item, "winch_id", 0),
            ))
    if "elastic_anchors" in d and d["elastic_anchors"]:
        for item in d["elastic_anchors"]:
            bcps.elastic_anchors.append(ElasticAnchorBCPData(
                position=_flist(item, "position", [0.0, 0.0, -100.0]),
                winch_id=_g(item, "winch_id", 0),
                anchor_mass=_g(item, "anchor_mass", 0.0),
                anchor_volume=_g(item, "anchor_volume", 0.0),
                c_param=_g(item, "c_param", 0.0),
                k_param=_g(item, "k_param", 0.0),
            ))
    return bcps


def _parse_body(d) -> BodyData:
    b = BodyData()
    b.name = _g(d, "name", b.name)   # GUI-only label stored as YAML comment-key if present
    b.body_type = _g(d, "type", b.body_type)
    b.take_cog_from_hdb = _g(d, "take_cog_hydro_database", b.take_cog_from_hdb)
    b.dofs = _flist(d, "dofs", b.dofs)
    b.bcps_indexes = _flist(d, "bcps_indexes", b.bcps_indexes)
    b.wind_turbines_indexes = _flist(d, "wind_turbines_indexes", b.wind_turbines_indexes)
    b.initial_position = _flist(d, "initial_position", b.initial_position)
    b.initial_displacement = _flist(d, "initial_displacement", b.initial_displacement)
    b.hdb_file = _g(d, "hydro_database", b.hdb_file)
    b.hdb_index = _g(d, "hydro_database_index", b.hdb_index)
    b.freedom_flag = _g(d, "freedom_flag", b.freedom_flag)
    b.imposed_motion_file = _g(d, "imposed_motion_file", b.imposed_motion_file)
    b.hydrostatics_flag = _g(d, "hydrostatics_flag", b.hydrostatics_flag)
    b.hydrostatics_mesh = _g(d, "hydrostatics_mesh", b.hydrostatics_mesh)
    b.radiation = _g(d, "radiation_flag", b.radiation)
    b.excitation_1st = _g(d, "excitation_1st_flag", b.excitation_1st)
    b.excitation_2nd = _g(d, "excitation_2nd_flag", b.excitation_2nd)
    b.viscous_added_mass = _flist(d, "viscous_added_mass", b.viscous_added_mass)
    b.viscous_linear_damping = _flist(d, "viscous_linear_damping", b.viscous_linear_damping)
    b.viscous_quadratic_damping = _flist(d, "viscous_quadratic_damping", b.viscous_quadratic_damping)
    return b


def _parse_line_type(d) -> LineTypeData:
    lt = LineTypeData()
    for attr in ("flag_tension", "density", "diameter", "flag_stiffness",
                 "EA", "beta", "CB", "Cmn", "Cdn", "Cdt", "GK", "GC",
                 "smoothstep", "friction_model", "vth", "ust", "usn", "ud", "deltamax"):
        if attr in d:
            setattr(lt, attr, d[attr])
    lt.kernel_lin_coef = _flist(d, "kernel_lin_coef", lt.kernel_lin_coef)
    lt.kernel_exp_coef = _flist(d, "kernel_exp_coef", lt.kernel_exp_coef)
    lt.elastic_coef = _flist(d, "elastic_coef", lt.elastic_coef)
    lt.strain_data = _flist(d, "strain_data", lt.strain_data)
    lt.stress_data = _flist(d, "stress_data", lt.stress_data)
    return lt


def _parse_line(d) -> LineData:
    ln = LineData()
    for attr in ("line_type", "type_index", "num_nodes", "polynomial_order",
                 "length", "seafloor_index", "BCP_N", "BCP_1"):
        if attr in d:
            setattr(ln, attr, d[attr])
    return ln


def _parse_spring_type(d) -> SpringTypeData:
    st = SpringTypeData()
    for attr in ("stress_model_flag", "damping_flag", "friction_flag", "frame_flag",
                 "mu_d", "mu_s", "vt", "Dt"):
        if attr in d:
            setattr(st, attr, d[attr])
    if "stiffness_matrix" in d:
        st.stiffness_matrix = [list(row) for row in d["stiffness_matrix"]]
    if "mass_matrix" in d:
        st.mass_matrix = [list(row) for row in d["mass_matrix"]]
    st.damping_vector = _flist(d, "damping_vector", st.damping_vector)
    if "stress_strain" in d and d["stress_strain"]:
        st.stress_strain = []
        for curve in d["stress_strain"]:
            sc = StressStrainCurve()
            sc.displacements = _flist(curve, "displacements", sc.displacements)
            if "forces" in curve:
                sc.forces = [list(row) for row in curve["forces"]]
            st.stress_strain.append(sc)
    return st


def _parse_spring(d) -> SpringData:
    s = SpringData()
    for attr in ("type_index", "BCP_1", "BCP_2", "BCP_1_type", "BCP_2_type"):
        if attr in d:
            setattr(s, attr, d[attr])
    if "vectors" in d:
        s.vectors = [list(v) for v in d["vectors"]]
    return s


def _parse_morison(d) -> MorisonData:
    m = MorisonData()
    for attr in ("body_index", "flow_type", "wind_speed", "wind_direction",
                 "current_speed", "current_direction", "symmetry_order", "hdf5_file"):
        if attr in d:
            setattr(m, attr, d[attr])
    for mat_key in ("wind_fk_x", "wind_drag_x", "wind_fk_y", "wind_drag_y",
                    "current_fk_x", "current_drag_x", "current_fk_y", "current_drag_y"):
        if mat_key in d:
            setattr(m, mat_key, [list(row) for row in d[mat_key]])
    return m


def _parse_winches(d) -> WinchesData:
    wd = WinchesData()
    if "winches" in d and d["winches"]:
        for item in d["winches"]:
            wd.winches.append(WinchItemData(
                line=_g(item, "line", 1),
                line_bcp=_g(item, "line_bcp", 2),
                inertia=_g(item, "inertia", 50.0),
                radius=_g(item, "radius", 0.25),
                drag=_g(item, "drag", 0.5),
            ))
    if "controller" in d and d["controller"]:
        c = d["controller"]
        wd.controller = WinchControllerData(
            controller_type=_g(c, "type", 1),
            target_tension=_g(c, "target_tension", 30000.0),
        )
    return wd


def _parse_owc(d) -> OWCData:
    o = OWCData()
    for attr in ("body_index", "floater_index", "waterplane_area", "ref_air_volume",
                 "hole_area", "discharge_coefficient", "turbine_type", "initial_omega"):
        if attr in d:
            setattr(o, attr, d[attr])
    o.position = _flist(d, "position", o.position)
    return o


def _parse_sinking(d) -> SinkingData:
    s = SinkingData()
    s.body_index = _g(d, "body_index", s.body_index)
    s.structural_mass_diagonal = _flist(d, "structural_mass_diagonal",
                                         s.structural_mass_diagonal)
    s.interp_masses = _flist(d, "interp_masses", s.interp_masses)
    s.hdb_files = list(_g(d, "hdb_files", s.hdb_files))
    if "groups" in d and d["groups"]:
        for g in d["groups"]:
            cg = CompartmentGroupData()
            cg.polygon_x = _flist(g, "polygon_x", cg.polygon_x)
            cg.polygon_y = _flist(g, "polygon_y", cg.polygon_y)
            cg.floor_z = _g(g, "floor_z", cg.floor_z)
            cg.filling_times = _flist(g, "filling_times", cg.filling_times)
            cg.filling_states = _flist(g, "filling_states", cg.filling_states)
            s.groups.append(cg)
    return s


def _parse_turbine(d) -> WindTurbineData:
    t = WindTurbineData()
    for attr in ("body_index", "fast_input_file", "fast_library"):
        if attr in d:
            setattr(t, attr, d[attr])
    return t


# ──────────────────────────────────────────────────────────────────────────────
# BUILD (dataclasses → YAML-ready dict)
# ──────────────────────────────────────────────────────────────────────────────
def _build_case(case: CaseData) -> dict:
    d: dict = {}
    d["problem"] = _build_problem(case.problem)
    d["bodies"] = [_build_body(b) for b in case.bodies]
    d["waves"] = _build_waves(case.waves)
    # BCPs
    bcps_dict = _build_bcps(case.bcps)
    if bcps_dict:
        d["bcps"] = bcps_dict
    # Sea floor
    sf = _build_seafloor(case.seafloor)
    if sf:
        d["seafloor"] = sf
    # Lines
    if case.line_types:
        d["line_types"] = [_build_line_type(lt) for lt in case.line_types]
    if case.lines:
        d["lines"] = [_build_line(ln) for ln in case.lines]
    # Springs
    if case.spring_types:
        d["spring_types"] = [_build_spring_type(st) for st in case.spring_types]
    if case.springs:
        d["springs"] = [_build_spring(s) for s in case.springs]
    # Morison
    if case.morison:
        d["morison"] = [_build_morison(m) for m in case.morison]
    # Winches
    if case.winches and case.winches.winches:
        d["winches"] = _build_winches(case.winches)
    # OWC
    if case.owcs:
        d["owcs"] = [_build_owc(o) for o in case.owcs]
    # Sinking
    if case.sinking:
        d["sinking"] = [_build_sinking(s) for s in case.sinking]
    # Wind turbines
    if case.wind_turbines:
        d["wind_turbines"] = [_build_turbine(t) for t in case.wind_turbines]
    return d


def _build_problem(p: ProblemData) -> dict:
    return {
        "gravity": p.gravity,
        "water_density": p.water_density,
        "air_density": p.air_density,
        "atmospheric_pressure": p.atmospheric_pressure,
        "air_adiabatic_dilation": p.air_adiabatic_dilation,
        "water_depth": p.water_depth,
        "write_time_step": p.write_time_step,
        "max_time_step": p.max_time_step,
        "hydro_time_step": p.hydro_time_step,
        "fast_time_step": p.fast_time_step,
        "fast_controller_time_step": p.fast_controller_time_step,
        "irf_time": p.irf_time,
        "sinking_time_step": p.sinking_time_step,
        "winches_controller_time_step": p.winches_controller_time_step,
        "owcs_controller_time_step": p.owcs_controller_time_step,
        "simulation_time": p.simulation_time,
        "rotation_simplification": p.rotation_simplification,
        "time_integration_method": p.time_integration_method,
        "time_integration_order": p.time_integration_order,
        "time_step_adaptivity": p.time_step_adaptivity,
        "jacobian_recomputation_steps": p.jacobian_recomputation_steps,
        "absolute_tolerance": p.absolute_tolerance,
        "relative_tolerance": p.relative_tolerance,
        "max_iterations_per_step": p.max_iterations_per_step,
        "read_equilibrium": p.read_equilibrium,
        "write_equilibrium": p.write_equilibrium,
        "mooring_initial_condition": p.mooring_initial_condition,
        "output_format": p.output_format,
    }


def _build_waves(w: WaveData) -> dict:
    d: dict = {
        "type": w.wave_type,
        "height": w.height,
        "period": w.period,
        "heading": w.heading,
        "ramp_time": w.ramp_time,
    }
    if w.wave_type == "IRR":
        d.update({
            "spectrum_type": w.spectrum_type,
            "piecewise_flag": w.piecewise_flag,
            "gamma": w.gamma,
            "spreading": w.spreading,
            "dtheta": w.dtheta,
            "factor": w.factor,
            "read_phases_flag": w.read_phases_flag,
            "relative_tolerance": w.wave_relative_tolerance,
            "dt": w.dt,
            "phases_file": w.phases_file,
            "database_file": w.database_file,
        })
    return d


def _build_seafloor(sf: SeaFloorData) -> dict:
    d: dict = {}
    if sf.flat:
        d["flat"] = [{"depth": f.depth} for f in sf.flat]
    if sf.inclined:
        d["inclined"] = [{"point1": _fseq(i.point1), "point2": _fseq(i.point2),
                          "point3": _fseq(i.point3)} for i in sf.inclined]
    if sf.bathymetry:
        d["bathymetry"] = [{"mesh_file": b.mesh_file} for b in sf.bathymetry]
    return d


def _build_bcps(bcps: BCPsData) -> dict:
    d: dict = {}
    if bcps.fairleads:
        d["fairleads"] = []
        for b in bcps.fairleads:
            item: dict = {"position": _fseq(b.position), "winch_id": b.winch_id}
            if b.actuator_file:
                item["actuator_file"] = b.actuator_file
            d["fairleads"].append(item)
    if bcps.anchors:
        d["anchors"] = [{"position": _fseq(b.position), "winch_id": b.winch_id}
                         for b in bcps.anchors]
    if bcps.joints:
        d["joints"] = [{"position": _fseq(b.position), "winch_id": b.winch_id,
                         "mass": b.mass, "volume": b.volume} for b in bcps.joints]
    if bcps.body_bcps:
        d["body_bcps"] = [{"position": _fseq(b.position), "winch_id": b.winch_id}
                           for b in bcps.body_bcps]
    if bcps.elastic_anchors:
        d["elastic_anchors"] = [{
            "position": _fseq(b.position), "winch_id": b.winch_id,
            "anchor_mass": b.anchor_mass, "anchor_volume": b.anchor_volume,
            "c_param": b.c_param, "k_param": b.k_param,
        } for b in bcps.elastic_anchors]
    return d


def _build_body(b: BodyData) -> dict:
    return {
        "type": b.body_type,
        "take_cog_hydro_database": b.take_cog_from_hdb,
        "dofs": _fseq(b.dofs),
        "bcps_indexes": _fseq(b.bcps_indexes),
        "wind_turbines_indexes": _fseq(b.wind_turbines_indexes),
        "initial_position": _fseq(b.initial_position),
        "initial_displacement": _fseq(b.initial_displacement),
        "hydro_database": b.hdb_file,
        "hydro_database_index": b.hdb_index,
        "freedom_flag": b.freedom_flag,
        "imposed_motion_file": b.imposed_motion_file,
        "hydrostatics_flag": b.hydrostatics_flag,
        "hydrostatics_mesh": b.hydrostatics_mesh,
        "radiation_flag": b.radiation,
        "excitation_1st_flag": b.excitation_1st,
        "excitation_2nd_flag": b.excitation_2nd,
        "viscous_added_mass": _fseq(b.viscous_added_mass),
        "viscous_linear_damping": _fseq(b.viscous_linear_damping),
        "viscous_quadratic_damping": _fseq(b.viscous_quadratic_damping),
    }


def _build_line_type(lt: LineTypeData) -> dict:
    d: dict = {
        "flag_tension": lt.flag_tension,
        "density": lt.density,
        "diameter": lt.diameter,
        "flag_stiffness": lt.flag_stiffness,
    }
    if lt.flag_stiffness == 0:
        d["EA"] = lt.EA
        d["beta"] = lt.beta
    elif lt.flag_stiffness == 1:
        d["kernel_lin_coef"] = _fseq(lt.kernel_lin_coef)
        d["kernel_exp_coef"] = _fseq(lt.kernel_exp_coef)
        d["elastic_coef"] = _fseq(lt.elastic_coef)
    else:
        d["strain_data"] = _fseq(lt.strain_data)
        d["stress_data"] = _fseq(lt.stress_data)
        d["beta"] = lt.beta
    d.update({
        "CB": lt.CB, "Cmn": lt.Cmn, "Cdn": lt.Cdn, "Cdt": lt.Cdt,
        "GK": lt.GK, "GC": lt.GC, "smoothstep": lt.smoothstep,
        "friction_model": lt.friction_model,
        "vth": lt.vth, "ust": lt.ust, "usn": lt.usn, "ud": lt.ud,
        "deltamax": lt.deltamax,
    })
    return d


def _build_line(ln: LineData) -> dict:
    return {
        "line_type": ln.line_type,
        "type_index": ln.type_index,
        "num_nodes": ln.num_nodes,
        "polynomial_order": ln.polynomial_order,
        "length": ln.length,
        "seafloor_index": ln.seafloor_index,
        "BCP_N": ln.BCP_N,
        "BCP_1": ln.BCP_1,
    }


def _build_spring_type(st: SpringTypeData) -> dict:
    d: dict = {
        "stress_model_flag": st.stress_model_flag,
        "damping_flag": st.damping_flag,
        "friction_flag": st.friction_flag,
        "frame_flag": st.frame_flag,
        "stiffness_matrix": [_fseq(row) for row in st.stiffness_matrix],
        "mu_d": st.mu_d, "mu_s": st.mu_s, "vt": st.vt, "Dt": st.Dt,
        "mass_matrix": [_fseq(row) for row in st.mass_matrix],
        "damping_vector": _fseq(st.damping_vector),
        "stress_strain": [],
    }
    for curve in st.stress_strain:
        d["stress_strain"].append({
            "displacements": _fseq(curve.displacements),
            "forces": [_fseq(row) for row in curve.forces],
        })
    return d


def _build_spring(s: SpringData) -> dict:
    return {
        "type_index": s.type_index,
        "BCP_1": s.BCP_1,
        "BCP_2": s.BCP_2,
        "BCP_1_type": s.BCP_1_type,
        "BCP_2_type": s.BCP_2_type,
        "vectors": [_fseq(v) for v in s.vectors],
    }


def _build_morison(m: MorisonData) -> dict:
    d: dict = {
        "body_index": m.body_index,
        "flow_type": m.flow_type,
        "wind_speed": m.wind_speed,
        "wind_direction": m.wind_direction,
        "current_speed": m.current_speed,
        "current_direction": m.current_direction,
        "symmetry_order": m.symmetry_order,
    }
    if m.flow_type == 2:
        d["hdf5_file"] = m.hdf5_file
    for mat_key in ("wind_fk_x", "wind_drag_x", "wind_fk_y", "wind_drag_y",
                    "current_fk_x", "current_drag_x", "current_fk_y", "current_drag_y"):
        d[mat_key] = _fseq([_fseq(row) for row in getattr(m, mat_key)])
    return d


def _build_winches(wd: WinchesData) -> dict:
    return {
        "winches": [{"line": w.line, "line_bcp": w.line_bcp,
                     "inertia": w.inertia, "radius": w.radius, "drag": w.drag}
                    for w in wd.winches],
        "controller": {
            "type": wd.controller.controller_type,
            "target_tension": wd.controller.target_tension,
        },
    }


def _build_owc(o: OWCData) -> dict:
    return {
        "body_index": o.body_index,
        "floater_index": o.floater_index,
        "waterplane_area": o.waterplane_area,
        "ref_air_volume": o.ref_air_volume,
        "hole_area": o.hole_area,
        "discharge_coefficient": o.discharge_coefficient,
        "turbine_type": o.turbine_type,
        "initial_omega": o.initial_omega,
        "position": _fseq(o.position),
    }


def _build_sinking(s: SinkingData) -> dict:
    d: dict = {
        "body_index": s.body_index,
        "structural_mass_diagonal": _fseq(s.structural_mass_diagonal),
        "interp_masses": _fseq(s.interp_masses),
        "hdb_files": _fseq(s.hdb_files),
        "groups": [],
    }
    for g in s.groups:
        d["groups"].append({
            "polygon_x": _fseq(g.polygon_x),
            "polygon_y": _fseq(g.polygon_y),
            "floor_z": g.floor_z,
            "filling_times": _fseq(g.filling_times),
            "filling_states": _fseq(g.filling_states),
        })
    return d


def _build_turbine(t: WindTurbineData) -> dict:
    return {
        "body_index": t.body_index,
        "fast_input_file": t.fast_input_file,
        "fast_library": t.fast_library,
    }
