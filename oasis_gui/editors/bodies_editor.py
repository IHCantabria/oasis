"""Bodies editor."""
from __future__ import annotations
from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox, QWidget,
    QCheckBox,
)
from PyQt5.QtCore import Qt

from .base_editor import (
    ListFormEditor, make_double_spin, make_spin, make_combo,
    make_line_edit, make_file_browse_row, make_scrollable,
)
from data.case_data import BodyData


def _spin6() -> list:
    return [make_double_spin(0.0, -1e9, 1e9, decimals=4) for _ in range(6)]


class BodiesEditor(ListFormEditor):
    def __init__(self, parent=None):
        super().__init__("Bodies", parent)

    def _new_item_data(self):
        b = BodyData(name="Body")
        b.initial_position = [0.0] * 6
        b.initial_displacement = [0.0] * 6
        b.dofs = [1] * 6
        b.cog = [0.0] * 3
        return b

    def _item_label(self, item) -> str:
        return item.name or "Body"

    def _build_form(self) -> QWidget:
        inner = QWidget()
        vbox = QVBoxLayout(inner)
        vbox.setAlignment(Qt.AlignTop)

        # â”€â”€ Identity â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        basic_box = QGroupBox("Identity")
        bf = QFormLayout(basic_box)
        self._f_name = make_line_edit("Body")
        self._f_type = make_combo(["0 â€“ Rigid floating body", "1 â€“ Rigid fixed body"])
        bf.addRow("Name:", self._f_name)
        bf.addRow("Body type:", self._f_type)
        vbox.addWidget(basic_box)

        # â”€â”€ Freedom â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        dof_box = QGroupBox("DOFs & Freedom")
        dof_v = QVBoxLayout(dof_box)
        self._f_dofs = []
        dof_row = QHBoxLayout()
        for i in range(6):
            cb = QCheckBox(f"DOF {i+1}")
            cb.setChecked(True)
            self._f_dofs.append(cb)
            dof_row.addWidget(cb)
        dof_v.addLayout(dof_row)
        ff = QFormLayout()
        self._f_freedom = make_combo(["0 â€“ Free", "1 â€“ Fixed", "2 â€“ Imposed motion"])
        self._f_freedom.currentIndexChanged.connect(self._on_freedom_changed)
        ff.addRow("Freedom flag:", self._f_freedom)
        dof_v.addLayout(ff)
        # Imposed motion container (conditionally visible)
        self._imposed_group = QGroupBox("Imposed motion file")
        imp_f = QFormLayout(self._imposed_group)
        self._f_imposed_file = make_line_edit("")
        imp_f.addRow("File:", make_file_browse_row(
            "Imposed motion file", self._f_imposed_file, inner,
            file_filter="Data files (*.dat *.txt);;All (*)"))
        self._imposed_group.setVisible(False)
        dof_v.addWidget(self._imposed_group)
        vbox.addWidget(dof_box)

        # â”€â”€ BCP associations â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        bcp_box = QGroupBox("Associated BCPs")
        bcpf = QFormLayout(bcp_box)
        self._f_bcp_indices = make_line_edit("1,2,3")
        bcpf.addRow("BCP global indices (comma-sep):", self._f_bcp_indices)
        vbox.addWidget(bcp_box)

        # â”€â”€ Initial conditions â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        ic_box = QGroupBox("Initial conditions")
        icf = QFormLayout(ic_box)
        self._f_init_pos = _spin6()
        self._f_init_disp = _spin6()
        pos_row = QHBoxLayout()
        disp_row = QHBoxLayout()
        for row, spins in [(pos_row, self._f_init_pos), (disp_row, self._f_init_disp)]:
            for s in spins:
                row.addWidget(s)
        icf.addRow("Initial position [m/deg]:", pos_row)
        icf.addRow("Initial displacement [m/deg]:", disp_row)
        vbox.addWidget(ic_box)

        # â”€â”€ Hydrodynamics â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        hydro_box = QGroupBox("Hydrodynamics (HDB)")
        hydrof = QFormLayout(hydro_box)
        self._f_hdb_file = make_line_edit("")
        self._f_hdb_index = make_spin(0, 0, 99)
        self._f_take_cog = QCheckBox("Use CoG from HDB")
        self._f_radiation = QCheckBox("Radiation")
        self._f_excitation_1st = QCheckBox("1st-order excitation")
        self._f_excitation_2nd = QCheckBox("2nd-order excitation")
        hydrof.addRow("HDB file:", make_file_browse_row(
            "HDB file", self._f_hdb_file, inner,
            file_filter="HDB files (*.h5 *.hdb);;All (*)"))
        hydrof.addRow("HDB body index:", self._f_hdb_index)
        hydrof.addRow("", self._f_take_cog)
        flags_row = QHBoxLayout()
        for cb in [self._f_radiation, self._f_excitation_1st, self._f_excitation_2nd]:
            flags_row.addWidget(cb)
        hydrof.addRow("Hydro flags:", flags_row)
        self._f_hydrostatics = make_combo(["0 â€“ None", "1 â€“ Mesh STL", "2 â€“ From file"])
        self._f_hydrostatics.currentIndexChanged.connect(self._on_hydrostatics_changed)
        hydrof.addRow("Hydrostatics:", self._f_hydrostatics)
        # STL container (conditionally visible)
        self._stl_group = QGroupBox("Hydrostatics STL file")
        stlf = QFormLayout(self._stl_group)
        self._f_stl_file = make_line_edit("")
        stlf.addRow("STL file:", make_file_browse_row(
            "STL file", self._f_stl_file, inner,
            file_filter="STL files (*.stl);;All (*)"))
        self._stl_group.setVisible(False)
        hydrof.addRow(self._stl_group)
        vbox.addWidget(hydro_box)

        # â”€â”€ Mass & CoG â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        misc_box = QGroupBox("Mass & CoG (if not from HDB)")
        miscf = QFormLayout(misc_box)
        self._f_mass_file = make_line_edit("")
        self._f_cog = [make_double_spin(0.0, -1e6, 1e6, decimals=4) for _ in range(3)]
        cog_row = QHBoxLayout()
        for s in self._f_cog:
            cog_row.addWidget(s)
        miscf.addRow("Mass/inertia file:", make_file_browse_row(
            "Mass/inertia file", self._f_mass_file, inner,
            file_filter="Data files (*.dat);;All (*)"))
        miscf.addRow("CoG override [x,y,z]:", cog_row)
        vbox.addWidget(misc_box)

        # â”€â”€ Schematic appearance â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
        sch_box = QGroupBox("Schematic appearance (GUI only)")
        schf = QFormLayout(sch_box)
        self._f_sch_width = make_double_spin(10.0, 0.1, 1000, decimals=2)
        self._f_sch_height = make_double_spin(10.0, 0.1, 1000, decimals=2)
        self._f_sch_stl = make_line_edit("")
        schf.addRow("Width [m]:", self._f_sch_width)
        schf.addRow("Height [m]:", self._f_sch_height)
        schf.addRow("STL for 3D view:", make_file_browse_row(
            "Body STL", self._f_sch_stl, inner,
            file_filter="STL files (*.stl);;All (*)"))
        vbox.addWidget(sch_box)

        return make_scrollable(inner)

    def _populate_form(self, item) -> None:
        self._f_name.setText(item.name)
        self._f_type.setCurrentIndex(item.body_type)
        for i, cb in enumerate(self._f_dofs):
            cb.setChecked(bool(item.dofs[i]) if i < len(item.dofs) else True)
        self._f_freedom.setCurrentIndex(item.freedom_flag)
        self._f_imposed_file.setText(item.imposed_motion_file)
        self._f_bcp_indices.setText(",".join(str(b) for b in item.bcps))
        for i, s in enumerate(self._f_init_pos):
            s.setValue(item.initial_position[i] if i < len(item.initial_position) else 0.0)
        for i, s in enumerate(self._f_init_disp):
            s.setValue(item.initial_displacement[i] if i < len(item.initial_displacement) else 0.0)
        self._f_hdb_file.setText(item.hdb_file)
        self._f_hdb_index.setValue(item.hdb_index)
        self._f_take_cog.setChecked(item.take_cog_from_hdb)
        self._f_hydrostatics.setCurrentIndex(item.hydrostatics_flag)
        self._f_stl_file.setText(item.stl_file)
        self._f_radiation.setChecked(item.radiation)
        self._f_excitation_1st.setChecked(item.excitation_1st)
        self._f_excitation_2nd.setChecked(item.excitation_2nd)
        self._f_mass_file.setText(item.mass_file)
        for i, s in enumerate(self._f_cog):
            s.setValue(item.cog[i] if i < len(item.cog) else 0.0)
        self._f_sch_width.setValue(item.schematic_width)
        self._f_sch_height.setValue(item.schematic_height)
        self._f_sch_stl.setText(item.stl_3d_file)
        self._on_freedom_changed(item.freedom_flag)
        self._on_hydrostatics_changed(item.hydrostatics_flag)

    def _read_form(self):
        item = self._items[self._current_idx]
        item.name = self._f_name.text()
        item.body_type = self._f_type.currentIndex()
        item.dofs = [1 if cb.isChecked() else 0 for cb in self._f_dofs]
        item.freedom_flag = self._f_freedom.currentIndex()
        item.imposed_motion_file = self._f_imposed_file.text()
        raw = self._f_bcp_indices.text()
        try:
            item.bcps = [int(x.strip()) for x in raw.split(",") if x.strip()]
        except ValueError:
            item.bcps = []
        item.initial_position = [s.value() for s in self._f_init_pos]
        item.initial_displacement = [s.value() for s in self._f_init_disp]
        item.hdb_file = self._f_hdb_file.text()
        item.hdb_index = self._f_hdb_index.value()
        item.take_cog_from_hdb = self._f_take_cog.isChecked()
        item.hydrostatics_flag = self._f_hydrostatics.currentIndex()
        item.stl_file = self._f_stl_file.text()
        item.radiation = self._f_radiation.isChecked()
        item.excitation_1st = self._f_excitation_1st.isChecked()
        item.excitation_2nd = self._f_excitation_2nd.isChecked()
        item.mass_file = self._f_mass_file.text()
        item.cog = [s.value() for s in self._f_cog]
        item.schematic_width = self._f_sch_width.value()
        item.schematic_height = self._f_sch_height.value()
        item.stl_3d_file = self._f_sch_stl.text()
        return item

    def _on_freedom_changed(self, idx):
        if hasattr(self, '_imposed_group'):
            self._imposed_group.setVisible(idx == 2)

    def _on_hydrostatics_changed(self, idx):
        if hasattr(self, '_stl_group'):
            self._stl_group.setVisible(idx > 0)

    def load_from_case(self, case) -> None:
        self._load_list(case.bodies)

    def save_to_case(self, case) -> None:
        self._commit_current()
        case.bodies = list(self._items)

