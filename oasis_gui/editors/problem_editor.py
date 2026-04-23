"""Problem settings editor."""
from __future__ import annotations
from PyQt5.QtWidgets import (
    QFormLayout, QGroupBox, QVBoxLayout, QDoubleSpinBox,
    QSpinBox, QComboBox, QWidget,
)
from PyQt5.QtCore import Qt

from .base_editor import BaseEditor, make_scrollable, make_double_spin, make_spin, make_combo
from data.case_data import ProblemData


class ProblemEditor(BaseEditor):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self):
        inner = QWidget()
        vbox = QVBoxLayout(inner)
        vbox.setAlignment(Qt.AlignTop)

        # ── Environment ──────────────────────────────────────────────────────
        env_box = QGroupBox("Environment")
        env_form = QFormLayout(env_box)
        self._gravity = make_double_spin(9.81, 0, 100, decimals=4)
        self._water_density = make_double_spin(1025.0, 0, 5000, decimals=2)
        self._air_density = make_double_spin(1.225, 0, 100, decimals=4)
        self._atm_pressure = make_double_spin(101325.0, 0, 1e7, decimals=1)
        self._air_adiab = make_double_spin(1.4, 0, 10, decimals=4)
        self._water_depth = make_double_spin(-100.0, -20000, 0, decimals=2)
        env_form.addRow("Gravity [m/s²]", self._gravity)
        env_form.addRow("Water density [kg/m³]", self._water_density)
        env_form.addRow("Air density [kg/m³]", self._air_density)
        env_form.addRow("Atmospheric pressure [Pa]", self._atm_pressure)
        env_form.addRow("Air adiabatic dilation [-]", self._air_adiab)
        env_form.addRow("Water depth [m] (negative)", self._water_depth)
        vbox.addWidget(env_box)

        # ── Time steps ───────────────────────────────────────────────────────
        dt_box = QGroupBox("Time Steps")
        dt_form = QFormLayout(dt_box)
        self._sim_time = make_double_spin(100.0, 0, 1e6, decimals=2)
        self._write_dt = make_double_spin(0.1, 1e-9, 1e6, decimals=6)
        self._max_dt = make_double_spin(0.1, 1e-9, 1e6, decimals=6)
        self._hydro_dt = make_double_spin(0.1, 1e-9, 1e6, decimals=6)
        self._fast_dt = make_double_spin(0.1, 1e-9, 1e6, decimals=6)
        self._fast_ctrl_dt = make_double_spin(0.1, 1e-9, 1e6, decimals=6)
        self._irf_time = make_double_spin(20.0, 0, 1e5, decimals=2)
        self._sinking_dt = make_double_spin(15.0, 1e-9, 1e6, decimals=4)
        self._winches_dt = make_double_spin(0.1, 1e-9, 1e6, decimals=6)
        self._owcs_dt = make_double_spin(0.1, 1e-9, 1e6, decimals=6)
        dt_form.addRow("Simulation time [s]", self._sim_time)
        dt_form.addRow("Output time step [s]", self._write_dt)
        dt_form.addRow("Max integration time step [s]", self._max_dt)
        dt_form.addRow("Hydro time step [s]", self._hydro_dt)
        dt_form.addRow("FAST turbine time step [s]", self._fast_dt)
        dt_form.addRow("FAST controller time step [s]", self._fast_ctrl_dt)
        dt_form.addRow("IRF time [s]", self._irf_time)
        dt_form.addRow("Sinking update time step [s]", self._sinking_dt)
        dt_form.addRow("Winches controller time step [s]", self._winches_dt)
        dt_form.addRow("OWC controller time step [s]", self._owcs_dt)
        vbox.addWidget(dt_box)

        # ── Solver ───────────────────────────────────────────────────────────
        sol_box = QGroupBox("ODE Solver")
        sol_form = QFormLayout(sol_box)
        self._method = make_combo(
            ["1 – BDF1", "2 – BDFN", "3 – ESDIRK46"], current="3 – ESDIRK46"
        )
        self._order = make_spin(2, 1, 10)
        self._adaptivity = make_combo(["0 – Fixed", "1 – Adaptive"])
        self._jac_steps = make_spin(0, 0, 9999)
        self._abs_tol = make_double_spin(1e-5, 1e-15, 1.0, decimals=10, step=1e-6)
        self._rel_tol = make_double_spin(1e-3, 1e-15, 1.0, decimals=10, step=1e-4)
        self._max_iter = make_spin(20, 1, 9999)
        self._rot_simp = make_combo(["0 – No simplification", "1 – Simplified"])
        sol_form.addRow("Integration method", self._method)
        sol_form.addRow("Integration order (BDFN)", self._order)
        sol_form.addRow("Time step adaptivity", self._adaptivity)
        sol_form.addRow("Jacobian recomputation steps", self._jac_steps)
        sol_form.addRow("Absolute tolerance", self._abs_tol)
        sol_form.addRow("Relative tolerance", self._rel_tol)
        sol_form.addRow("Max iterations per step", self._max_iter)
        sol_form.addRow("Rotation simplification", self._rot_simp)
        vbox.addWidget(sol_box)

        # ── Initial conditions / Output ───────────────────────────────────────
        out_box = QGroupBox("Output & Initial Conditions")
        out_form = QFormLayout(out_box)
        self._read_eq = make_combo(["0 – No", "1 – Yes"])
        self._write_eq = make_combo(["0 – No", "1 – Yes"])
        self._mooring_ic = make_combo(["0 – Catenary", "1 – Newton's method"])
        self._output_fmt = make_combo(["csv", "txt"])
        out_form.addRow("Read equilibrium IC", self._read_eq)
        out_form.addRow("Write equilibrium IC", self._write_eq)
        out_form.addRow("Mooring initial condition", self._mooring_ic)
        out_form.addRow("Output format", self._output_fmt)
        vbox.addWidget(out_box)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(make_scrollable(inner))

    # ── Data binding ─────────────────────────────────────────────────────────
    def load_from_case(self, case) -> None:
        p = case.problem
        self._gravity.setValue(p.gravity)
        self._water_density.setValue(p.water_density)
        self._air_density.setValue(p.air_density)
        self._atm_pressure.setValue(p.atmospheric_pressure)
        self._air_adiab.setValue(p.air_adiabatic_dilation)
        self._water_depth.setValue(p.water_depth)
        self._sim_time.setValue(p.simulation_time)
        self._write_dt.setValue(p.write_time_step)
        self._max_dt.setValue(p.max_time_step)
        self._hydro_dt.setValue(p.hydro_time_step)
        self._fast_dt.setValue(p.fast_time_step)
        self._fast_ctrl_dt.setValue(p.fast_controller_time_step)
        self._irf_time.setValue(p.irf_time)
        self._sinking_dt.setValue(p.sinking_time_step)
        self._winches_dt.setValue(p.winches_controller_time_step)
        self._owcs_dt.setValue(p.owcs_controller_time_step)
        self._method.setCurrentIndex(p.time_integration_method - 1)
        self._order.setValue(p.time_integration_order)
        self._adaptivity.setCurrentIndex(p.time_step_adaptivity)
        self._jac_steps.setValue(p.jacobian_recomputation_steps)
        self._abs_tol.setValue(p.absolute_tolerance)
        self._rel_tol.setValue(p.relative_tolerance)
        self._max_iter.setValue(p.max_iterations_per_step)
        self._rot_simp.setCurrentIndex(p.rotation_simplification)
        self._read_eq.setCurrentIndex(p.read_equilibrium)
        self._write_eq.setCurrentIndex(p.write_equilibrium)
        self._mooring_ic.setCurrentIndex(p.mooring_initial_condition)
        self._output_fmt.setCurrentText(p.output_format)

    def save_to_case(self, case) -> None:
        p = case.problem
        p.gravity = self._gravity.value()
        p.water_density = self._water_density.value()
        p.air_density = self._air_density.value()
        p.atmospheric_pressure = self._atm_pressure.value()
        p.air_adiabatic_dilation = self._air_adiab.value()
        p.water_depth = self._water_depth.value()
        p.simulation_time = self._sim_time.value()
        p.write_time_step = self._write_dt.value()
        p.max_time_step = self._max_dt.value()
        p.hydro_time_step = self._hydro_dt.value()
        p.fast_time_step = self._fast_dt.value()
        p.fast_controller_time_step = self._fast_ctrl_dt.value()
        p.irf_time = self._irf_time.value()
        p.sinking_time_step = self._sinking_dt.value()
        p.winches_controller_time_step = self._winches_dt.value()
        p.owcs_controller_time_step = self._owcs_dt.value()
        p.time_integration_method = self._method.currentIndex() + 1
        p.time_integration_order = self._order.value()
        p.time_step_adaptivity = self._adaptivity.currentIndex()
        p.jacobian_recomputation_steps = self._jac_steps.value()
        p.absolute_tolerance = self._abs_tol.value()
        p.relative_tolerance = self._rel_tol.value()
        p.max_iterations_per_step = self._max_iter.value()
        p.rotation_simplification = self._rot_simp.currentIndex()
        p.read_equilibrium = self._read_eq.currentIndex()
        p.write_equilibrium = self._write_eq.currentIndex()
        p.mooring_initial_condition = self._mooring_ic.currentIndex()
        p.output_format = self._output_fmt.currentText()
