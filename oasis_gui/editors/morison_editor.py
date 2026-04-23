"""Morison editor — per-body wind/current Morison coefficients."""
from __future__ import annotations
from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox, QWidget,
    QTabWidget, QTableWidget, QTableWidgetItem, QHeaderView,
)
from PyQt5.QtCore import Qt

from .base_editor import (
    ListFormEditor, make_double_spin, make_spin, make_combo,
    make_scrollable,
)
from data.case_data import MorisonData


_COEF_TABLES = [
    ("wind_fk_x",    "Wind FK X"),
    ("wind_drag_x",  "Wind Drag X"),
    ("wind_fk_y",    "Wind FK Y"),
    ("wind_drag_y",  "Wind Drag Y"),
    ("current_fk_x", "Current FK X"),
    ("current_drag_x","Current Drag X"),
    ("current_fk_y", "Current FK Y"),
    ("current_drag_y","Current Drag Y"),
]


def _make_coef_table() -> QTableWidget:
    t = QTableWidget(6, 2)
    t.setHorizontalHeaderLabels(["A (amplitude)", "Phase [deg]"])
    t.setVerticalHeaderLabels([f"DOF {i+1}" for i in range(6)])
    t.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    for r in range(6):
        for c in range(2):
            t.setItem(r, c, QTableWidgetItem("0.0"))
    return t


def _fill_coef(table: QTableWidget, mat: list) -> None:
    for r in range(6):
        for c in range(2):
            val = mat[r][c] if r < len(mat) and c < len(mat[r]) else 0.0
            table.item(r, c).setText(f"{val:.6g}")


def _read_coef(table: QTableWidget) -> list:
    mat = []
    for r in range(6):
        row = []
        for c in range(2):
            item = table.item(r, c)
            try:
                row.append(float(item.text()) if item else 0.0)
            except ValueError:
                row.append(0.0)
        mat.append(row)
    return mat


class MorisonEditor(ListFormEditor):
    def __init__(self, parent=None):
        super().__init__("Morison bodies", parent)

    def _new_item_data(self):
        return MorisonData()

    def _item_label(self, item) -> str:
        return f"Body {item.body_index}"

    def _build_form(self) -> QWidget:
        inner = QWidget()
        vbox = QVBoxLayout(inner)
        vbox.setAlignment(Qt.AlignTop)

        # Basic params
        basic_box = QGroupBox("Basic parameters")
        bf = QFormLayout(basic_box)
        self._m_body_index = make_spin(1, 1, 999)
        self._m_flow_type = make_combo(["1 – Constant", "2 – Variable (HDF5)"])
        self._m_wind_speed = make_double_spin(0.0, 0, 200, decimals=3)
        self._m_wind_dir = make_double_spin(0.0, -360, 360, decimals=2)
        self._m_cur_speed = make_double_spin(0.0, 0, 20, decimals=4)
        self._m_cur_dir = make_double_spin(0.0, -360, 360, decimals=2)
        self._m_sym_order = make_spin(2, 0, 8)
        bf.addRow("Body index (1-based):", self._m_body_index)
        bf.addRow("Flow type:", self._m_flow_type)
        bf.addRow("Wind speed [m/s]:", self._m_wind_speed)
        bf.addRow("Wind direction [deg]:", self._m_wind_dir)
        bf.addRow("Current speed [m/s]:", self._m_cur_speed)
        bf.addRow("Current direction [deg]:", self._m_cur_dir)
        bf.addRow("Symmetry order:", self._m_sym_order)
        vbox.addWidget(basic_box)

        # Coefficient tables (8 tabs)
        coef_box = QGroupBox("Coefficient matrices (6 DOF × [amplitude, phase])")
        cv = QVBoxLayout(coef_box)
        self._coef_tabs = QTabWidget()
        self._coef_tbl: dict[str, QTableWidget] = {}
        for attr, label in _COEF_TABLES:
            table = _make_coef_table()
            self._coef_tbl[attr] = table
            self._coef_tabs.addTab(table, label)
        cv.addWidget(self._coef_tabs)
        vbox.addWidget(coef_box)

        return make_scrollable(inner)

    def _populate_form(self, item: MorisonData) -> None:
        self._m_body_index.setValue(item.body_index)
        self._m_flow_type.setCurrentIndex(item.flow_type - 1)
        self._m_wind_speed.setValue(item.wind_speed)
        self._m_wind_dir.setValue(item.wind_direction)
        self._m_cur_speed.setValue(item.current_speed)
        self._m_cur_dir.setValue(item.current_direction)
        self._m_sym_order.setValue(item.symmetry_order)
        for attr, _ in _COEF_TABLES:
            mat = getattr(item, attr, [[0.0, 0.0]] * 6)
            _fill_coef(self._coef_tbl[attr], mat)

    def _read_form(self):
        item = self._items[self._current_idx]
        item.body_index = self._m_body_index.value()
        item.flow_type = self._m_flow_type.currentIndex() + 1
        item.wind_speed = self._m_wind_speed.value()
        item.wind_direction = self._m_wind_dir.value()
        item.current_speed = self._m_cur_speed.value()
        item.current_direction = self._m_cur_dir.value()
        item.symmetry_order = self._m_sym_order.value()
        for attr, _ in _COEF_TABLES:
            setattr(item, attr, _read_coef(self._coef_tbl[attr]))
        return item

    def load_from_case(self, case) -> None:
        self._load_list(case.morison)

    def save_to_case(self, case) -> None:
        self._commit_current()
        case.morison = list(self._items)
