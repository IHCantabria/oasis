"""Sea floor editor."""
from __future__ import annotations
from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QGroupBox, QFormLayout,
    QComboBox, QStackedWidget, QWidget, QPushButton, QListWidget, QLabel,
)
from PyQt5.QtCore import Qt

from .base_editor import BaseEditor, make_double_spin, make_scrollable, make_line_edit, make_file_browse_row
from data.case_data import (
    SeaFloorData, FlatFloorData, InclinedFloorData, BathymetryFloorData,
)


class SeaFloorEditor(BaseEditor):
    """Editor for seafloor: flat, inclined or bathymetry sections."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self):
        inner = QWidget()
        vbox = QVBoxLayout(inner)
        vbox.setAlignment(Qt.AlignTop)

        # ── Flat floors ──────────────────────────────────────────────────────
        flat_box = QGroupBox("Flat sea floor sections")
        flat_v = QVBoxLayout(flat_box)
        flat_label = QLabel("Each flat section defines a constant depth.")
        flat_v.addWidget(flat_label)
        self._flat_list = QListWidget()
        self._flat_list.setMaximumHeight(100)
        flat_v.addWidget(self._flat_list)
        flat_btns = QHBoxLayout()
        btn_add_flat = QPushButton("+")
        btn_add_flat.setMaximumWidth(36)
        btn_rem_flat = QPushButton("−")
        btn_rem_flat.setMaximumWidth(36)
        flat_btns.addWidget(btn_add_flat)
        flat_btns.addWidget(btn_rem_flat)
        flat_btns.addStretch()
        flat_v.addLayout(flat_btns)
        flat_form = QFormLayout()
        self._flat_depth = make_double_spin(-100.0, -20000, 0, decimals=3)
        flat_form.addRow("Depth [m] (negative):", self._flat_depth)
        flat_v.addLayout(flat_form)
        vbox.addWidget(flat_box)

        # ── Inclined floors ──────────────────────────────────────────────────
        inc_box = QGroupBox("Inclined sea floor sections")
        inc_v = QVBoxLayout(inc_box)
        inc_label = QLabel("Each inclined section is defined by three points.")
        inc_v.addWidget(inc_label)
        self._inc_list = QListWidget()
        self._inc_list.setMaximumHeight(80)
        inc_v.addWidget(self._inc_list)
        inc_btns = QHBoxLayout()
        btn_add_inc = QPushButton("+")
        btn_add_inc.setMaximumWidth(36)
        btn_rem_inc = QPushButton("−")
        btn_rem_inc.setMaximumWidth(36)
        inc_btns.addWidget(btn_add_inc)
        inc_btns.addWidget(btn_rem_inc)
        inc_btns.addStretch()
        inc_v.addLayout(inc_btns)
        inc_form = QFormLayout()
        self._inc_p1 = [make_double_spin(0.0, -1e6, 1e6, decimals=3) for _ in range(3)]
        self._inc_p2 = [make_double_spin(50.0, -1e6, 1e6, decimals=3) for _ in range(3)]
        self._inc_p3 = [make_double_spin(0.0, -1e6, 1e6, decimals=3) for _ in range(3)]
        self._inc_p1[2].setValue(-100.0)
        self._inc_p2[2].setValue(-100.0)
        self._inc_p3[2].setValue(-100.0)
        for label, spins in [("Point 1 [x,y,z]", self._inc_p1),
                              ("Point 2 [x,y,z]", self._inc_p2),
                              ("Point 3 [x,y,z]", self._inc_p3)]:
            row = QHBoxLayout()
            for s in spins:
                row.addWidget(s)
            inc_form.addRow(label + ":", row)
        inc_v.addLayout(inc_form)
        vbox.addWidget(inc_box)

        # ── Bathymetry floors ─────────────────────────────────────────────────
        bat_box = QGroupBox("Bathymetry sections (mesh file)")
        bat_v = QVBoxLayout(bat_box)
        self._bat_list = QListWidget()
        self._bat_list.setMaximumHeight(80)
        bat_v.addWidget(self._bat_list)
        bat_btns = QHBoxLayout()
        btn_add_bat = QPushButton("+")
        btn_add_bat.setMaximumWidth(36)
        btn_rem_bat = QPushButton("−")
        btn_rem_bat.setMaximumWidth(36)
        bat_btns.addWidget(btn_add_bat)
        bat_btns.addWidget(btn_rem_bat)
        bat_btns.addStretch()
        bat_v.addLayout(bat_btns)
        bat_form = QFormLayout()
        self._bat_file = make_line_edit("")
        browse_row = make_file_browse_row("Bathymetry mesh", self._bat_file, self,
                                          file_filter="Mesh files (*.dat *.txt *.msh);;All Files (*)")
        bat_form.addRow("Mesh file:", browse_row)
        bat_v.addLayout(bat_form)
        vbox.addWidget(bat_box)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(make_scrollable(inner))

        # Connect list signals to form updates
        self._flat_list.currentRowChanged.connect(self._on_flat_select)
        self._inc_list.currentRowChanged.connect(self._on_inc_select)
        self._bat_list.currentRowChanged.connect(self._on_bat_select)

        btn_add_flat.clicked.connect(self._add_flat)
        btn_rem_flat.clicked.connect(self._rem_flat)
        btn_add_inc.clicked.connect(self._add_inc)
        btn_rem_inc.clicked.connect(self._rem_inc)
        btn_add_bat.clicked.connect(self._add_bat)
        btn_rem_bat.clicked.connect(self._rem_bat)

        self._flat_items: list[FlatFloorData] = []
        self._inc_items: list[InclinedFloorData] = []
        self._bat_items: list[BathymetryFloorData] = []

    # ── Flat helpers ──────────────────────────────────────────────────────────
    def _add_flat(self):
        self._flat_items.append(FlatFloorData(depth=-100.0))
        self._flat_list.addItem(f"Flat #{len(self._flat_items)}")
        self._flat_list.setCurrentRow(len(self._flat_items) - 1)

    def _rem_flat(self):
        row = self._flat_list.currentRow()
        if 0 <= row < len(self._flat_items):
            self._flat_items.pop(row)
            self._flat_list.takeItem(row)

    def _on_flat_select(self, row):
        if 0 <= row < len(self._flat_items):
            self._flat_depth.setValue(self._flat_items[row].depth)

    # ── Inclined helpers ──────────────────────────────────────────────────────
    def _add_inc(self):
        self._inc_items.append(InclinedFloorData())
        self._inc_list.addItem(f"Inclined #{len(self._inc_items)}")
        self._inc_list.setCurrentRow(len(self._inc_items) - 1)

    def _rem_inc(self):
        row = self._inc_list.currentRow()
        if 0 <= row < len(self._inc_items):
            self._inc_items.pop(row)
            self._inc_list.takeItem(row)

    def _on_inc_select(self, row):
        if 0 <= row < len(self._inc_items):
            item = self._inc_items[row]
            for i, s in enumerate(self._inc_p1):
                s.setValue(item.point1[i])
            for i, s in enumerate(self._inc_p2):
                s.setValue(item.point2[i])
            for i, s in enumerate(self._inc_p3):
                s.setValue(item.point3[i])

    # ── Bathymetry helpers ────────────────────────────────────────────────────
    def _add_bat(self):
        self._bat_items.append(BathymetryFloorData(mesh_file=""))
        self._bat_list.addItem(f"Bathymetry #{len(self._bat_items)}")
        self._bat_list.setCurrentRow(len(self._bat_items) - 1)

    def _rem_bat(self):
        row = self._bat_list.currentRow()
        if 0 <= row < len(self._bat_items):
            self._bat_items.pop(row)
            self._bat_list.takeItem(row)

    def _on_bat_select(self, row):
        if 0 <= row < len(self._bat_items):
            self._bat_file.setText(self._bat_items[row].mesh_file)

    # ── Data binding ──────────────────────────────────────────────────────────
    def load_from_case(self, case) -> None:
        sf = case.seafloor
        self._flat_items = list(sf.flat)
        self._inc_items = list(sf.inclined)
        self._bat_items = list(sf.bathymetry)
        self._flat_list.clear()
        for i, _ in enumerate(self._flat_items):
            self._flat_list.addItem(f"Flat #{i+1}  depth={self._flat_items[i].depth}")
        self._inc_list.clear()
        for i, _ in enumerate(self._inc_items):
            self._inc_list.addItem(f"Inclined #{i+1}")
        self._bat_list.clear()
        for i, item in enumerate(self._bat_items):
            self._bat_list.addItem(f"Bathymetry #{i+1}  {item.mesh_file}")
        if self._flat_items:
            self._flat_list.setCurrentRow(0)
        if self._inc_items:
            self._inc_list.setCurrentRow(0)
        if self._bat_items:
            self._bat_list.setCurrentRow(0)

    def save_to_case(self, case) -> None:
        # Commit currently displayed flat item
        row = self._flat_list.currentRow()
        if 0 <= row < len(self._flat_items):
            self._flat_items[row].depth = self._flat_depth.value()
        row = self._inc_list.currentRow()
        if 0 <= row < len(self._inc_items):
            item = self._inc_items[row]
            item.point1 = [s.value() for s in self._inc_p1]
            item.point2 = [s.value() for s in self._inc_p2]
            item.point3 = [s.value() for s in self._inc_p3]
        row = self._bat_list.currentRow()
        if 0 <= row < len(self._bat_items):
            self._bat_items[row].mesh_file = self._bat_file.text()

        case.seafloor = SeaFloorData(
            flat=list(self._flat_items),
            inclined=list(self._inc_items),
            bathymetry=list(self._bat_items),
        )
