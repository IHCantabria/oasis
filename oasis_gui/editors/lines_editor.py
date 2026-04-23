"""Lines editor — line types and instances."""
from __future__ import annotations
from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox, QWidget,
    QSplitter, QLabel,
)
from PyQt5.QtCore import Qt

from .base_editor import (
    ListFormEditor, make_double_spin, make_spin, make_combo,
    make_line_edit, make_scrollable,
)
from data.case_data import LineTypeData, LineData


# ──────────────────────────────────────────────────────────────────────────────
# Line types
# ──────────────────────────────────────────────────────────────────────────────
class LineTypesEditor(ListFormEditor):
    def __init__(self, parent=None):
        super().__init__("Line Types", parent)

    def _new_item_data(self):
        return LineTypeData()

    def _item_label(self, item) -> str:
        return f"Type {self._items.index(item)+1}  EA={item.EA:.2e}"

    def _build_form(self) -> QWidget:
        inner = QWidget()
        vbox = QVBoxLayout(inner)
        vbox.setAlignment(Qt.AlignTop)

        # Mechanical
        mech_box = QGroupBox("Mechanical properties")
        mf = QFormLayout(mech_box)
        self._t_density = make_double_spin(50.0, 0, 1e9, decimals=4)
        self._t_diameter = make_double_spin(0.09, 0, 100, decimals=6)
        self._t_flag_tension = make_combo(["1 – Symmetric", "2 – Tension-only"])
        self._t_flag_stiffness = make_combo(["0 – Constant EA/beta",
                                              "1 – Viscoelastic kernel",
                                              "2+ – Tabulated strain-stress"])
        self._t_flag_stiffness.currentIndexChanged.connect(self._on_stiffness_changed)
        mf.addRow("Density [kg/m]:", self._t_density)
        mf.addRow("Diameter [m]:", self._t_diameter)
        mf.addRow("Tension flag:", self._t_flag_tension)
        mf.addRow("Stiffness mode:", self._t_flag_stiffness)
        vbox.addWidget(mech_box)

        # Mode 0 – constant EA
        self._ea_box = QGroupBox("Mode 0 – Constant EA")
        ef = QFormLayout(self._ea_box)
        self._t_EA = make_double_spin(1e7, 0, 1e15, decimals=0)
        self._t_beta = make_double_spin(0.002, 0, 1.0, decimals=6)
        ef.addRow("EA [N]:", self._t_EA)
        ef.addRow("beta (structural damping):", self._t_beta)
        vbox.addWidget(self._ea_box)

        # Mode 1 – viscoelastic kernel (comma-separated coefficients)
        self._ve_box = QGroupBox("Mode 1 – Viscoelastic kernel")
        vef = QFormLayout(self._ve_box)
        self._t_kernel_lin = make_line_edit("1e7")
        self._t_kernel_exp = make_line_edit("-1.0")
        self._t_elastic = make_line_edit("1e7")
        vef.addRow("Linear coefs (comma-sep):", self._t_kernel_lin)
        vef.addRow("Exp coefs (comma-sep):", self._t_kernel_exp)
        vef.addRow("Elastic coefs (comma-sep):", self._t_elastic)
        self._ve_box.setVisible(False)
        vbox.addWidget(self._ve_box)

        # Mode >1 – tabulated
        self._tab_box = QGroupBox("Mode 2+ – Tabulated strain-stress")
        tabf = QFormLayout(self._tab_box)
        self._t_strain = make_line_edit("0.0,0.01,0.1")
        self._t_stress = make_line_edit("0.0,1e5,1e6")
        tabf.addRow("Strain values (comma-sep):", self._t_strain)
        tabf.addRow("Stress values [Pa] (comma-sep):", self._t_stress)
        self._tab_box.setVisible(False)
        vbox.addWidget(self._tab_box)

        # Hydrodynamics / seabed
        hyd_box = QGroupBox("Hydrodynamics & Seabed")
        hf = QFormLayout(hyd_box)
        self._t_CB = make_double_spin(0.0, 0, 10, decimals=4)
        self._t_Cmn = make_double_spin(3.8, 0, 100, decimals=4)
        self._t_Cdn = make_double_spin(2.5, 0, 100, decimals=4)
        self._t_Cdt = make_double_spin(0.5, 0, 100, decimals=4)
        self._t_GK = make_double_spin(1e6, 0, 1e15, decimals=0)
        self._t_GC = make_double_spin(0.1, 0, 1e9, decimals=4)
        self._t_smoothstep = make_combo(["0 – No", "1 – Yes"])
        self._t_friction_model = make_combo(["0 – None", "1 – Coulomb"])
        self._t_vth = make_double_spin(0.01, 0, 100, decimals=6)
        self._t_ust = make_double_spin(0.6, 0, 10, decimals=4)
        self._t_usn = make_double_spin(0.4, 0, 10, decimals=4)
        self._t_ud = make_double_spin(0.3, 0, 10, decimals=4)
        self._t_deltamax = make_double_spin(0.03, 0, 100, decimals=6)
        hf.addRow("Seabed contact coefficient CB:", self._t_CB)
        hf.addRow("Normal added mass Cmn:", self._t_Cmn)
        hf.addRow("Normal drag Cdn:", self._t_Cdn)
        hf.addRow("Tangential drag Cdt:", self._t_Cdt)
        hf.addRow("Seabed stiffness GK [N/m²]:", self._t_GK)
        hf.addRow("Seabed damping GC:", self._t_GC)
        hf.addRow("Smoothstep seabed:", self._t_smoothstep)
        hf.addRow("Friction model:", self._t_friction_model)
        hf.addRow("Threshold velocity vth [m/s]:", self._t_vth)
        hf.addRow("Static friction μ_st:", self._t_ust)
        hf.addRow("Normal static friction μ_sn:", self._t_usn)
        hf.addRow("Dynamic friction μ_d:", self._t_ud)
        hf.addRow("Max deformation delta_max [m]:", self._t_deltamax)
        vbox.addWidget(hyd_box)

        return make_scrollable(inner)

    def _on_stiffness_changed(self, idx):
        self._ea_box.setVisible(idx == 0)
        self._ve_box.setVisible(idx == 1)
        self._tab_box.setVisible(idx >= 2)

    def _populate_form(self, item: LineTypeData) -> None:
        self._t_density.setValue(item.density)
        self._t_diameter.setValue(item.diameter)
        self._t_flag_tension.setCurrentIndex(item.flag_tension - 1)
        idx = min(item.flag_stiffness, 2)
        self._t_flag_stiffness.setCurrentIndex(idx)
        self._t_EA.setValue(item.EA)
        self._t_beta.setValue(item.beta)
        self._t_kernel_lin.setText(",".join(str(x) for x in item.kernel_lin_coef))
        self._t_kernel_exp.setText(",".join(str(x) for x in item.kernel_exp_coef))
        self._t_elastic.setText(",".join(str(x) for x in item.elastic_coef))
        self._t_strain.setText(",".join(str(x) for x in item.strain_data))
        self._t_stress.setText(",".join(str(x) for x in item.stress_data))
        self._t_CB.setValue(item.CB)
        self._t_Cmn.setValue(item.Cmn)
        self._t_Cdn.setValue(item.Cdn)
        self._t_Cdt.setValue(item.Cdt)
        self._t_GK.setValue(item.GK)
        self._t_GC.setValue(item.GC)
        self._t_smoothstep.setCurrentIndex(item.smoothstep)
        self._t_friction_model.setCurrentIndex(item.friction_model)
        self._t_vth.setValue(item.vth)
        self._t_ust.setValue(item.ust)
        self._t_usn.setValue(item.usn)
        self._t_ud.setValue(item.ud)
        self._t_deltamax.setValue(item.deltamax)
        self._on_stiffness_changed(idx)

    def _read_form(self):
        item = self._items[self._current_idx]
        item.density = self._t_density.value()
        item.diameter = self._t_diameter.value()
        item.flag_tension = self._t_flag_tension.currentIndex() + 1
        item.flag_stiffness = self._t_flag_stiffness.currentIndex()
        item.EA = self._t_EA.value()
        item.beta = self._t_beta.value()
        item.kernel_lin_coef = _parse_floats(self._t_kernel_lin.text())
        item.kernel_exp_coef = _parse_floats(self._t_kernel_exp.text())
        item.elastic_coef = _parse_floats(self._t_elastic.text())
        item.strain_data = _parse_floats(self._t_strain.text())
        item.stress_data = _parse_floats(self._t_stress.text())
        item.CB = self._t_CB.value()
        item.Cmn = self._t_Cmn.value()
        item.Cdn = self._t_Cdn.value()
        item.Cdt = self._t_Cdt.value()
        item.GK = self._t_GK.value()
        item.GC = self._t_GC.value()
        item.smoothstep = self._t_smoothstep.currentIndex()
        item.friction_model = self._t_friction_model.currentIndex()
        item.vth = self._t_vth.value()
        item.ust = self._t_ust.value()
        item.usn = self._t_usn.value()
        item.ud = self._t_ud.value()
        item.deltamax = self._t_deltamax.value()
        return item


# ──────────────────────────────────────────────────────────────────────────────
# Line instances
# ──────────────────────────────────────────────────────────────────────────────
class LineInstancesEditor(ListFormEditor):
    def __init__(self, parent=None):
        super().__init__("Line Instances", parent)

    def _new_item_data(self):
        return LineData()

    def _item_label(self, item) -> str:
        return f"Line  BCP{item.BCP_1}→{item.BCP_N}"

    def _build_form(self) -> QWidget:
        inner = QWidget()
        vbox = QVBoxLayout(inner)
        vbox.setAlignment(Qt.AlignTop)
        form_box = QGroupBox("Line instance")
        ff = QFormLayout(form_box)
        self._i_line_type = make_combo(["1 – Mooring", "2 – Towing", "3 – Tensor"])
        self._i_type_index = make_spin(1, 1, 999)
        self._i_num_nodes = make_spin(11, 2, 9999)
        self._i_poly_order = make_spin(1, 1, 10)
        self._i_length = make_double_spin(120.0, 0, 1e6, decimals=4)
        self._i_seafloor_index = make_spin(1, 1, 999)
        self._i_bcp_1 = make_spin(1, 1, 9999)
        self._i_bcp_n = make_spin(2, 1, 9999)
        ff.addRow("Line type:", self._i_line_type)
        ff.addRow("Line type index (1-based):", self._i_type_index)
        ff.addRow("Number of nodes:", self._i_num_nodes)
        ff.addRow("Polynomial order:", self._i_poly_order)
        ff.addRow("Unstretched length [m]:", self._i_length)
        ff.addRow("Sea floor section index:", self._i_seafloor_index)
        ff.addRow("BCP start (global 1-based):", self._i_bcp_1)
        ff.addRow("BCP end (global 1-based):", self._i_bcp_n)
        vbox.addWidget(form_box)
        return make_scrollable(inner)

    def _populate_form(self, item: LineData) -> None:
        self._i_line_type.setCurrentIndex(item.line_type - 1)
        self._i_type_index.setValue(item.type_index)
        self._i_num_nodes.setValue(item.num_nodes)
        self._i_poly_order.setValue(item.polynomial_order)
        self._i_length.setValue(item.length)
        self._i_seafloor_index.setValue(item.seafloor_index)
        self._i_bcp_1.setValue(item.BCP_1)
        self._i_bcp_n.setValue(item.BCP_N)

    def _read_form(self):
        item = self._items[self._current_idx]
        item.line_type = self._i_line_type.currentIndex() + 1
        item.type_index = self._i_type_index.value()
        item.num_nodes = self._i_num_nodes.value()
        item.polynomial_order = self._i_poly_order.value()
        item.length = self._i_length.value()
        item.seafloor_index = self._i_seafloor_index.value()
        item.BCP_1 = self._i_bcp_1.value()
        item.BCP_N = self._i_bcp_n.value()
        return item


# ──────────────────────────────────────────────────────────────────────────────
# Top-level Lines editor
# ──────────────────────────────────────────────────────────────────────────────
def _parse_floats(text: str) -> list:
    try:
        return [float(x.strip()) for x in text.split(",") if x.strip()]
    except ValueError:
        return []


class LinesEditor(QWidget):
    """Wraps LineTypesEditor and LineInstancesEditor side by side."""

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(2, 2, 2, 2)

        splitter = QSplitter(Qt.Horizontal)
        self._types_ed = LineTypesEditor()
        self._instances_ed = LineInstancesEditor()
        splitter.addWidget(self._types_ed)
        splitter.addWidget(self._instances_ed)
        splitter.setSizes([500, 300])
        layout.addWidget(splitter)

    def load_from_case(self, case) -> None:
        self._types_ed.set_items(case.line_types)
        self._instances_ed.set_items(case.lines)

    def save_to_case(self, case) -> None:
        self._types_ed._commit_current()
        self._instances_ed._commit_current()
        case.line_types = list(self._types_ed._items)
        case.lines = list(self._instances_ed._items)
