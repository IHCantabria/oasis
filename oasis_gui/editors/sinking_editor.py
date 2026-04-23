"""Sinking editor."""
from __future__ import annotations
from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox, QWidget,
    QListWidget, QPushButton, QLabel, QSplitter, QTableWidget,
    QTableWidgetItem, QHeaderView, QTabWidget,
)
from PyQt5.QtCore import Qt

from .base_editor import (
    ListFormEditor, make_double_spin, make_spin, make_line_edit,
    make_file_browse_row, make_scrollable,
)
from data.case_data import SinkingData, CompartmentGroupData


class SinkingEditor(ListFormEditor):
    def __init__(self, parent=None):
        super().__init__("Sinking scenarios", parent)

    def _new_item_data(self):
        s = SinkingData()
        s.structural_mass_diagonal = [1.0] * 6
        s.interp_masses = [0.0]
        return s

    def _item_label(self, item) -> str:
        return f"Body {item.body_index} sinking"

    def _build_form(self) -> QWidget:
        inner = QWidget()
        vbox = QVBoxLayout(inner)
        vbox.setAlignment(Qt.AlignTop)

        basic_box = QGroupBox("Basic")
        bf = QFormLayout(basic_box)
        self._sk_body = make_spin(1, 1, 999)
        self._sk_mass_diag = make_line_edit("1,1,1,1,1,1")
        self._sk_interp_masses = make_line_edit("0.0")
        self._sk_hdb_files = make_line_edit("")
        bf.addRow("Body index:", self._sk_body)
        bf.addRow("Structural mass diag (6 vals):", self._sk_mass_diag)
        bf.addRow("Interpolation masses (comma-sep):", self._sk_interp_masses)
        bf.addRow("HDB files (comma-sep paths):", self._sk_hdb_files)
        vbox.addWidget(basic_box)

        # Compartment groups sub-editor
        grp_box = QGroupBox("Compartment groups")
        gv = QVBoxLayout(grp_box)
        self._grp_list = QListWidget()
        self._grp_list.setMaximumHeight(80)
        btn_row = QHBoxLayout()
        btn_add = QPushButton("+")
        btn_add.setMaximumWidth(36)
        btn_rem = QPushButton("−")
        btn_rem.setMaximumWidth(36)
        btn_row.addWidget(btn_add)
        btn_row.addWidget(btn_rem)
        btn_row.addStretch()
        gv.addWidget(self._grp_list)
        gv.addLayout(btn_row)

        grp_form = QFormLayout()
        self._grp_poly_x = make_line_edit("0,1,1,0")
        self._grp_poly_y = make_line_edit("0,0,1,1")
        self._grp_floor_z = make_double_spin(-5.0, -1e5, 0, decimals=3)
        self._grp_fill_times = make_line_edit("0.0,100.0")
        self._grp_fill_states = make_line_edit("0.0,0.0")
        grp_form.addRow("Polygon X (comma-sep):", self._grp_poly_x)
        grp_form.addRow("Polygon Y (comma-sep):", self._grp_poly_y)
        grp_form.addRow("Floor Z [m]:", self._grp_floor_z)
        grp_form.addRow("Filling times [s] (comma-sep):", self._grp_fill_times)
        grp_form.addRow("Filling states (0–1, comma-sep):", self._grp_fill_states)
        gv.addLayout(grp_form)
        vbox.addWidget(grp_box)

        btn_add.clicked.connect(self._add_group)
        btn_rem.clicked.connect(self._rem_group)
        self._grp_list.currentRowChanged.connect(self._on_grp_select)
        self._groups: list[CompartmentGroupData] = []
        self._current_grp = -1

        return make_scrollable(inner)

    # ── Compartment group list helpers ────────────────────────────────────────
    def _add_group(self):
        self._commit_group()
        self._groups.append(CompartmentGroupData())
        self._grp_list.addItem(f"Group #{len(self._groups)}")
        self._grp_list.setCurrentRow(len(self._groups) - 1)

    def _rem_group(self):
        row = self._grp_list.currentRow()
        if 0 <= row < len(self._groups):
            self._groups.pop(row)
            self._current_grp = -1
            self._grp_list.takeItem(row)

    def _on_grp_select(self, row):
        self._commit_group()
        self._current_grp = row
        if 0 <= row < len(self._groups):
            g = self._groups[row]
            self._grp_poly_x.setText(",".join(str(v) for v in g.polygon_x))
            self._grp_poly_y.setText(",".join(str(v) for v in g.polygon_y))
            self._grp_floor_z.setValue(g.floor_z)
            self._grp_fill_times.setText(",".join(str(v) for v in g.filling_times))
            self._grp_fill_states.setText(",".join(str(v) for v in g.filling_states))

    def _commit_group(self):
        row = self._current_grp
        if 0 <= row < len(self._groups):
            g = self._groups[row]
            g.polygon_x = _pf(self._grp_poly_x.text())
            g.polygon_y = _pf(self._grp_poly_y.text())
            g.floor_z = self._grp_floor_z.value()
            g.filling_times = _pf(self._grp_fill_times.text())
            g.filling_states = _pf(self._grp_fill_states.text())

    # ── ListFormEditor contract ───────────────────────────────────────────────
    def _populate_form(self, item: SinkingData) -> None:
        self._sk_body.setValue(item.body_index)
        self._sk_mass_diag.setText(",".join(str(v) for v in item.structural_mass_diagonal))
        self._sk_interp_masses.setText(",".join(str(v) for v in item.interp_masses))
        self._sk_hdb_files.setText(",".join(item.hdb_files))
        self._groups = list(item.groups)
        self._grp_list.clear()
        for i, _ in enumerate(self._groups):
            self._grp_list.addItem(f"Group #{i+1}")
        if self._groups:
            self._grp_list.setCurrentRow(0)
        self._current_grp = -1

    def _read_form(self):
        self._commit_group()
        item = self._items[self._current_idx]
        item.body_index = self._sk_body.value()
        item.structural_mass_diagonal = _pf(self._sk_mass_diag.text())
        item.interp_masses = _pf(self._sk_interp_masses.text())
        raw_hdb = self._sk_hdb_files.text()
        item.hdb_files = [s.strip() for s in raw_hdb.split(",") if s.strip()]
        item.groups = list(self._groups)
        return item

    def load_from_case(self, case) -> None:
        self._load_list(case.sinking)

    def save_to_case(self, case) -> None:
        self._commit_current()
        case.sinking = list(self._items)


def _pf(text: str) -> list:
    try:
        return [float(x.strip()) for x in text.split(",") if x.strip()]
    except ValueError:
        return []
