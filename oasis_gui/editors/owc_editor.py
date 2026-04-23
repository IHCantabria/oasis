"""OWC editor."""
from __future__ import annotations
from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox, QWidget,
)
from PyQt5.QtCore import Qt

from .base_editor import (
    ListFormEditor, make_double_spin, make_spin, make_combo, make_scrollable,
)
from data.case_data import OWCData


class OWCEditor(ListFormEditor):
    def __init__(self, parent=None):
        super().__init__("OWC devices", parent)

    def _new_item_data(self):
        return OWCData()

    def _item_label(self, item) -> str:
        return f"OWC  body={item.body_index}"

    def _build_form(self) -> QWidget:
        inner = QWidget()
        vbox = QVBoxLayout(inner)
        vbox.setAlignment(Qt.AlignTop)
        form_box = QGroupBox("OWC parameters")
        ff = QFormLayout(form_box)
        self._o_body = make_spin(2, 1, 999)
        self._o_floater = make_spin(1, 1, 999)
        self._o_wp_area = make_double_spin(1.0, 0, 1e6, decimals=4)
        self._o_air_vol = make_double_spin(1.0, 0, 1e6, decimals=4)
        self._o_hole_area = make_double_spin(0.1, 0, 1e6, decimals=6)
        self._o_Cd = make_double_spin(0.6, 0, 10, decimals=4)
        self._o_turbine_type = make_spin(-1, -1, 999)
        self._o_omega = make_double_spin(0.0, 0, 1e6, decimals=4)
        self._o_pos = [make_double_spin(0.0, -1e6, 1e6, decimals=4) for _ in range(3)]
        pos_row = QHBoxLayout()
        for s in self._o_pos:
            pos_row.addWidget(s)
        ff.addRow("Body index:", self._o_body)
        ff.addRow("Floater index:", self._o_floater)
        ff.addRow("Waterplane area [m²]:", self._o_wp_area)
        ff.addRow("Reference air volume [m³]:", self._o_air_vol)
        ff.addRow("Orifice hole area [m²]:", self._o_hole_area)
        ff.addRow("Discharge coefficient Cd:", self._o_Cd)
        ff.addRow("Turbine type (-1=none, 0=hole):", self._o_turbine_type)
        ff.addRow("Initial omega [rad/s]:", self._o_omega)
        ff.addRow("Position [x,y,z]:", pos_row)
        vbox.addWidget(form_box)
        return make_scrollable(inner)

    def _populate_form(self, item: OWCData) -> None:
        self._o_body.setValue(item.body_index)
        self._o_floater.setValue(item.floater_index)
        self._o_wp_area.setValue(item.waterplane_area)
        self._o_air_vol.setValue(item.ref_air_volume)
        self._o_hole_area.setValue(item.hole_area)
        self._o_Cd.setValue(item.discharge_coefficient)
        self._o_turbine_type.setValue(item.turbine_type)
        self._o_omega.setValue(item.initial_omega)
        for i, s in enumerate(self._o_pos):
            s.setValue(item.position[i] if i < len(item.position) else 0.0)

    def _read_form(self):
        item = self._items[self._current_idx]
        item.body_index = self._o_body.value()
        item.floater_index = self._o_floater.value()
        item.waterplane_area = self._o_wp_area.value()
        item.ref_air_volume = self._o_air_vol.value()
        item.hole_area = self._o_hole_area.value()
        item.discharge_coefficient = self._o_Cd.value()
        item.turbine_type = self._o_turbine_type.value()
        item.initial_omega = self._o_omega.value()
        item.position = [s.value() for s in self._o_pos]
        return item

    def load_from_case(self, case) -> None:
        self._load_list(case.owcs)

    def save_to_case(self, case) -> None:
        self._commit_current()
        case.owcs = list(self._items)
