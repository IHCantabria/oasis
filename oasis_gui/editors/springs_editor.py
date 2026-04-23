"""Springs editor — types and instances."""
from __future__ import annotations
from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox, QWidget,
    QTabWidget, QTableWidget, QTableWidgetItem, QHeaderView, QSplitter,
    QLabel,
)
from PyQt5.QtCore import Qt

from .base_editor import (
    ListFormEditor, make_double_spin, make_spin, make_combo,
    make_scrollable,
)
from data.case_data import SpringTypeData, SpringData, StressStrainCurve


def _make_6x6_table() -> QTableWidget:
    t = QTableWidget(6, 6)
    t.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    t.verticalHeader().setSectionResizeMode(QHeaderView.Stretch)
    t.setMaximumHeight(200)
    for r in range(6):
        for c in range(6):
            t.setItem(r, c, QTableWidgetItem("0.0"))
    return t


def _read_6x6(table: QTableWidget) -> list:
    mat = []
    for r in range(6):
        row = []
        for c in range(6):
            item = table.item(r, c)
            try:
                row.append(float(item.text()) if item else 0.0)
            except ValueError:
                row.append(0.0)
        mat.append(row)
    return mat


def _fill_6x6(table: QTableWidget, mat: list) -> None:
    for r in range(6):
        for c in range(6):
            val = mat[r][c] if r < len(mat) and c < len(mat[r]) else 0.0
            table.item(r, c).setText(f"{val:.6g}")


class SpringTypesEditor(ListFormEditor):
    def __init__(self, parent=None):
        super().__init__("Spring Types", parent)

    def _new_item_data(self):
        return SpringTypeData()

    def _item_label(self, item) -> str:
        return f"Type {self._items.index(item)+1}"

    def _build_form(self) -> QWidget:
        inner = QWidget()
        vbox = QVBoxLayout(inner)
        vbox.setAlignment(Qt.AlignTop)

        # Flags
        flag_box = QGroupBox("Flags")
        ff = QFormLayout(flag_box)
        self._st_stress_model = make_combo(["1 – Linear stiffness matrix",
                                             "2 – Nonlinear stress-strain"])
        self._st_stress_model.currentIndexChanged.connect(self._on_stress_model)
        self._st_damping_flag = make_combo(["0 – None", "1 – Linear", "2 – Rayleigh"])
        self._st_friction_flag = make_combo(["0 – No", "1 – Yes"])
        self._st_frame_flag = make_combo(["0 – Global", "1 – Local"])
        ff.addRow("Stress model:", self._st_stress_model)
        ff.addRow("Damping flag:", self._st_damping_flag)
        ff.addRow("Friction flag:", self._st_friction_flag)
        ff.addRow("Frame flag:", self._st_frame_flag)
        vbox.addWidget(flag_box)

        # Linear stiffness 6×6
        self._lin_box = QGroupBox("Stiffness matrix K (6×6)")
        lv = QVBoxLayout(self._lin_box)
        self._st_K = _make_6x6_table()
        lv.addWidget(self._st_K)
        vbox.addWidget(self._lin_box)

        # Nonlinear stress-strain (6 tabs, one per DOF)
        self._nl_box = QGroupBox("Nonlinear stress-strain curves (one per DOF)")
        nv = QVBoxLayout(self._nl_box)
        self._ss_tabs = QTabWidget()
        self._ss_tables = []
        for i in range(6):
            table = QTableWidget(3, 2)
            table.setHorizontalHeaderLabels(["Displacement [m/rad]", "Force [N/Nm]"])
            table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            for r in range(3):
                for c in range(2):
                    table.setItem(r, c, QTableWidgetItem("0.0"))
            self._ss_tables.append(table)
            self._ss_tabs.addTab(table, f"DOF {i+1}")
        nv.addWidget(self._ss_tabs)
        self._nl_box.setVisible(False)
        vbox.addWidget(self._nl_box)

        # Friction & damping params
        misc_box = QGroupBox("Friction & damping parameters")
        mf = QFormLayout(misc_box)
        self._st_mu_d = make_double_spin(0.0, 0, 100, decimals=4)
        self._st_mu_s = make_double_spin(0.0, 0, 100, decimals=4)
        self._st_vt = make_double_spin(0.001, 0, 100, decimals=6)
        self._st_Dt = make_double_spin(0.01, 0, 100, decimals=6)
        mf.addRow("Dynamic friction μ_d:", self._st_mu_d)
        mf.addRow("Static friction μ_s:", self._st_mu_s)
        mf.addRow("Velocity threshold vt [m/s]:", self._st_vt)
        mf.addRow("Dissipation Dt:", self._st_Dt)
        vbox.addWidget(misc_box)

        return make_scrollable(inner)

    def _on_stress_model(self, idx):
        self._lin_box.setVisible(idx == 0)
        self._nl_box.setVisible(idx == 1)

    def _populate_form(self, item: SpringTypeData) -> None:
        self._st_stress_model.setCurrentIndex(item.stress_model_flag - 1)
        self._st_damping_flag.setCurrentIndex(item.damping_flag)
        self._st_friction_flag.setCurrentIndex(item.friction_flag)
        self._st_frame_flag.setCurrentIndex(item.frame_flag)
        _fill_6x6(self._st_K, item.stiffness_matrix)
        for i, curve in enumerate(item.stress_strain):
            table = self._ss_tables[i]
            n_pts = len(curve.displacements)
            table.setRowCount(n_pts)
            for r in range(n_pts):
                for c in range(2):
                    if c == 0:
                        val = curve.displacements[r] if r < len(curve.displacements) else 0.0
                    else:
                        forces = curve.forces[i] if i < len(curve.forces) else []
                        val = forces[r] if r < len(forces) else 0.0
                    item_w = table.item(r, c)
                    if item_w is None:
                        item_w = QTableWidgetItem(f"{val:.6g}")
                        table.setItem(r, c, item_w)
                    else:
                        item_w.setText(f"{val:.6g}")
        self._st_mu_d.setValue(item.mu_d)
        self._st_mu_s.setValue(item.mu_s)
        self._st_vt.setValue(item.vt)
        self._st_Dt.setValue(item.Dt)
        self._on_stress_model(item.stress_model_flag - 1)

    def _read_form(self):
        item = self._items[self._current_idx]
        item.stress_model_flag = self._st_stress_model.currentIndex() + 1
        item.damping_flag = self._st_damping_flag.currentIndex()
        item.friction_flag = self._st_friction_flag.currentIndex()
        item.frame_flag = self._st_frame_flag.currentIndex()
        item.stiffness_matrix = _read_6x6(self._st_K)
        for i, table in enumerate(self._ss_tables):
            disps = []
            forces = []
            for r in range(table.rowCount()):
                d_item = table.item(r, 0)
                f_item = table.item(r, 1)
                try:
                    disps.append(float(d_item.text()) if d_item else 0.0)
                    forces.append(float(f_item.text()) if f_item else 0.0)
                except ValueError:
                    disps.append(0.0)
                    forces.append(0.0)
            while i >= len(item.stress_strain):
                item.stress_strain.append(StressStrainCurve())
            item.stress_strain[i] = StressStrainCurve(
                displacements=disps,
                forces=[forces if j == i else item.stress_strain[j].forces[0]
                        if i < len(item.stress_strain) else [0.0]
                        for j in range(6)],
            )
        item.mu_d = self._st_mu_d.value()
        item.mu_s = self._st_mu_s.value()
        item.vt = self._st_vt.value()
        item.Dt = self._st_Dt.value()
        return item


class SpringInstancesEditor(ListFormEditor):
    def __init__(self, parent=None):
        super().__init__("Spring Instances", parent)

    def _new_item_data(self):
        return SpringData()

    def _item_label(self, item) -> str:
        return f"Spring  BCP{item.BCP_1}→{item.BCP_2}"

    def _build_form(self) -> QWidget:
        inner = QWidget()
        vbox = QVBoxLayout(inner)
        vbox.setAlignment(Qt.AlignTop)
        form_box = QGroupBox("Spring instance")
        ff = QFormLayout(form_box)
        self._si_type_index = make_spin(1, 1, 999)
        self._si_bcp1 = make_spin(1, 1, 9999)
        self._si_bcp2 = make_spin(2, 1, 9999)
        self._si_bcp1_type = make_combo(["0 – Fixed point", "1 – On line", "2 – On plane"])
        self._si_bcp2_type = make_combo(["0 – Fixed point", "1 – On line", "2 – On plane"])
        ff.addRow("Spring type index (1-based):", self._si_type_index)
        ff.addRow("BCP 1 (global):", self._si_bcp1)
        ff.addRow("BCP 2 (global):", self._si_bcp2)
        ff.addRow("BCP 1 type:", self._si_bcp1_type)
        ff.addRow("BCP 2 type:", self._si_bcp2_type)
        vbox.addWidget(form_box)

        # Vectors 6×3
        vec_box = QGroupBox("Direction vectors (6 rows × 3 cols)")
        vec_v = QVBoxLayout(vec_box)
        self._si_vectors = QTableWidget(6, 3)
        self._si_vectors.setHorizontalHeaderLabels(["X", "Y", "Z"])
        self._si_vectors.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        default_vecs = [[1,0,0],[0,1,0],[0,0,1],[1,0,0],[0,1,0],[0,0,1]]
        for r in range(6):
            for c in range(3):
                self._si_vectors.setItem(r, c, QTableWidgetItem(str(default_vecs[r][c])))
        vec_v.addWidget(self._si_vectors)
        vbox.addWidget(vec_box)
        return make_scrollable(inner)

    def _populate_form(self, item: SpringData) -> None:
        self._si_type_index.setValue(item.type_index)
        self._si_bcp1.setValue(item.BCP_1)
        self._si_bcp2.setValue(item.BCP_2)
        self._si_bcp1_type.setCurrentIndex(item.BCP_1_type)
        self._si_bcp2_type.setCurrentIndex(item.BCP_2_type)
        for r in range(6):
            for c in range(3):
                val = item.vectors[r][c] if r < len(item.vectors) and c < len(item.vectors[r]) else 0.0
                self._si_vectors.item(r, c).setText(f"{val:.4g}")

    def _read_form(self):
        item = self._items[self._current_idx]
        item.type_index = self._si_type_index.value()
        item.BCP_1 = self._si_bcp1.value()
        item.BCP_2 = self._si_bcp2.value()
        item.BCP_1_type = self._si_bcp1_type.currentIndex()
        item.BCP_2_type = self._si_bcp2_type.currentIndex()
        item.vectors = []
        for r in range(6):
            row = []
            for c in range(3):
                wi = self._si_vectors.item(r, c)
                try:
                    row.append(float(wi.text()) if wi else 0.0)
                except ValueError:
                    row.append(0.0)
            item.vectors.append(row)
        return item


class SpringsEditor(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(2, 2, 2, 2)
        splitter = QSplitter(Qt.Horizontal)
        self._types_ed = SpringTypesEditor()
        self._instances_ed = SpringInstancesEditor()
        splitter.addWidget(self._types_ed)
        splitter.addWidget(self._instances_ed)
        splitter.setSizes([600, 300])
        layout.addWidget(splitter)

    def load_from_case(self, case) -> None:
        self._types_ed.set_items(case.spring_types)
        self._instances_ed.set_items(case.springs)

    def save_to_case(self, case) -> None:
        self._types_ed._commit_current()
        self._instances_ed._commit_current()
        case.spring_types = list(self._types_ed._items)
        case.springs = list(self._instances_ed._items)
