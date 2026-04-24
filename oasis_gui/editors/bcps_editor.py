"""BCPs editor — manages all five BCP types in tabbed lists."""
from __future__ import annotations
from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QTabWidget, QWidget, QPushButton, QListWidget, QLabel,
    QDoubleSpinBox, QSplitter, QComboBox, QSizePolicy,
)
from PyQt5.QtCore import Qt

from .base_editor import BaseEditor, make_double_spin, make_spin, make_line_edit, make_file_browse_row
from data.case_data import (
    BCPsData, AnchorBCPData, FairleadBCPData, JointBCPData,
    BodyBCPData, ElasticAnchorBCPData,
)


def _xyz_row(parent):
    """Returns (HBoxLayout, [spin_x, spin_y, spin_z])."""
    row = QHBoxLayout()
    spins = []
    for label, default in [("X", 0.0), ("Y", 0.0), ("Z", 0.0)]:
        l = QLabel(label)
        l.setFixedWidth(14)
        s = make_double_spin(default, -1e6, 1e6, decimals=4)
        row.addWidget(l)
        row.addWidget(s)
        spins.append(s)
    return row, spins


class _BCPTypeTab(QWidget):
    """Generic BCP sub-editor for one BCP type."""

    def __init__(self, title: str, extra_fields: list = None, parent=None):
        """extra_fields: list of (label, widget) tuples appended after position+winch_id."""
        super().__init__(parent)
        self._items: list = []
        self._current = -1
        self._extra_fields = extra_fields or []

        splitter = QSplitter(Qt.Horizontal)

        # Left: list + buttons
        left = QWidget()
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 0, 0)
        self._list = QListWidget()
        self._list.setMaximumWidth(160)
        lv.addWidget(self._list)
        btn_row = QHBoxLayout()
        self._btn_add = QPushButton("+")
        self._btn_add.setMaximumWidth(36)
        self._btn_rem = QPushButton("−")
        self._btn_rem.setMaximumWidth(36)
        btn_row.addWidget(self._btn_add)
        btn_row.addWidget(self._btn_rem)
        btn_row.addStretch()
        lv.addLayout(btn_row)
        splitter.addWidget(left)

        # Right: form
        right = QWidget()
        form = QFormLayout(right)
        pos_row, self._pos_spins = _xyz_row(self)
        form.addRow("Position [m]:", pos_row)
        self._winch_id = make_spin(0, 0, 999)
        form.addRow("Winch ID (0=none):", self._winch_id)
        # Extra fields
        self._extra_widgets = []
        for lbl, widget in self._extra_fields:
            form.addRow(lbl, widget)
            self._extra_widgets.append(widget)
        splitter.addWidget(right)
        splitter.setSizes([140, 400])

        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(4, 4, 4, 4)
        main_layout.setSpacing(4)
        lbl = QLabel(f"<b>{title}</b>")
        lbl.setMaximumHeight(24)
        lbl.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        main_layout.addWidget(lbl)
        main_layout.addWidget(splitter, 1)

        self._btn_add.clicked.connect(self._on_add)
        self._btn_rem.clicked.connect(self._on_rem)
        self._list.currentRowChanged.connect(self._on_select)

    def set_items(self, items):
        self._items = list(items)
        self._refresh_list()
        if self._items:
            self._list.setCurrentRow(0)

    def get_items(self):
        self._commit()
        return list(self._items)

    def _refresh_list(self):
        self._list.blockSignals(True)
        self._list.clear()
        for i, item in enumerate(self._items):
            p = item.position
            self._list.addItem(f"#{i+1}  ({p[0]:.1f}, {p[1]:.1f}, {p[2]:.1f})")
        self._list.blockSignals(False)

    def _on_select(self, row):
        if self._current >= 0 and self._current < len(self._items):
            self._commit()
        self._current = row
        if 0 <= row < len(self._items):
            self._populate(self._items[row])

    def _commit(self):
        if 0 <= self._current < len(self._items):
            item = self._items[self._current]
            item.position = [s.value() for s in self._pos_spins]
            item.winch_id = self._winch_id.value()
            self._commit_extra(item)
            p = item.position
            list_item = self._list.item(self._current)
            if list_item:
                list_item.setText(f"#{self._current+1}  ({p[0]:.1f}, {p[1]:.1f}, {p[2]:.1f})")

    def _populate(self, item):
        for i, s in enumerate(self._pos_spins):
            s.setValue(item.position[i] if i < len(item.position) else 0.0)
        self._winch_id.setValue(item.winch_id)
        self._populate_extra(item)

    def _on_add(self):
        new_item = self._make_new()
        self._items.append(new_item)
        p = new_item.position
        self._list.addItem(f"#{len(self._items)}  ({p[0]:.1f}, {p[1]:.1f}, {p[2]:.1f})")
        self._list.setCurrentRow(len(self._items) - 1)

    def _on_rem(self):
        row = self._list.currentRow()
        if 0 <= row < len(self._items):
            self._items.pop(row)
            self._current = -1
            self._refresh_list()
            if self._items:
                self._list.setCurrentRow(min(row, len(self._items) - 1))

    def _make_new(self):
        raise NotImplementedError

    def _populate_extra(self, item):
        pass

    def _commit_extra(self, item):
        pass


class _AnchorTab(_BCPTypeTab):
    def __init__(self, parent=None):
        super().__init__("Anchors (fixed seabed points)", parent=parent)

    def _make_new(self):
        return AnchorBCPData()


class _ActuatorTab(_BCPTypeTab):
    def __init__(self, parent=None):
        self._actuator_le = make_line_edit("")
        super().__init__("Actuators (prescribed-motion BCPs, global frame)",
                         extra_fields=[("Actuator file (if driven):", self._actuator_le)],
                         parent=parent)

    def _make_new(self):
        return FairleadBCPData()

    def _populate_extra(self, item):
        self._actuator_le.setText(item.actuator_file)

    def _commit_extra(self, item):
        item.actuator_file = self._actuator_le.text()


class _JointTab(_BCPTypeTab):
    def __init__(self, parent=None):
        self._mass_spin = make_double_spin(0.0, 0, 1e12, decimals=3)
        self._vol_spin = make_double_spin(0.0, 0, 1e9, decimals=6)
        super().__init__("Joints (free underwater connection points)",
                         extra_fields=[("Mass [kg]:", self._mass_spin),
                                        ("Volume [m³]:", self._vol_spin)],
                         parent=parent)

    def _make_new(self):
        return JointBCPData()

    def _populate_extra(self, item):
        self._mass_spin.setValue(item.mass)
        self._vol_spin.setValue(item.volume)

    def _commit_extra(self, item):
        item.mass = self._mass_spin.value()
        item.volume = self._vol_spin.value()


class _BodyBCPTab(_BCPTypeTab):
    def __init__(self, parent=None):
        super().__init__("Body BCPs (attachment points on body, body-local frame)", parent=parent)

    def _make_new(self):
        return BodyBCPData()


class _ElasticAnchorTab(_BCPTypeTab):
    def __init__(self, parent=None):
        self._anc_mass = make_double_spin(0.0, 0, 1e12, decimals=3)
        self._anc_vol = make_double_spin(0.0, 0, 1e9, decimals=6)
        self._c_param = make_double_spin(0.0, 0, 1e12, decimals=4)
        self._k_param = make_double_spin(0.0, 0, 1e12, decimals=4)
        super().__init__("Elastic Anchors (spring-mounted anchors)",
                         extra_fields=[("Anchor mass [kg]:", self._anc_mass),
                                        ("Anchor volume [m³]:", self._anc_vol),
                                        ("C damping param:", self._c_param),
                                        ("K stiffness param:", self._k_param)],
                         parent=parent)

    def _make_new(self):
        return ElasticAnchorBCPData()

    def _populate_extra(self, item):
        self._anc_mass.setValue(item.anchor_mass)
        self._anc_vol.setValue(item.anchor_volume)
        self._c_param.setValue(item.c_param)
        self._k_param.setValue(item.k_param)

    def _commit_extra(self, item):
        item.anchor_mass = self._anc_mass.value()
        item.anchor_volume = self._anc_vol.value()
        item.c_param = self._c_param.value()
        item.k_param = self._k_param.value()


class BCPsEditor(BaseEditor):
    """Top-level BCP editor: tabs for each BCP type."""

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        info = QLabel(
            "BCP global index allocation order (as seen by OASIS): "
            "<b>Actuators → Anchors → Joints → Body BCPs → Elastic Anchors</b>"
        )
        info.setWordWrap(True)
        info.setMaximumHeight(36)
        info.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        layout.addWidget(info)

        tabs = QTabWidget()
        self._bodybcp_tab = _BodyBCPTab()
        self._anchor_tab = _AnchorTab()
        self._elastic_tab = _ElasticAnchorTab()
        self._joint_tab = _JointTab()
        self._actuator_tab = _ActuatorTab()
        tabs.addTab(self._bodybcp_tab, "Body BCPs")
        tabs.addTab(self._anchor_tab, "Anchors")
        tabs.addTab(self._elastic_tab, "Elastic Anchors")
        tabs.addTab(self._joint_tab, "Joints")
        tabs.addTab(self._actuator_tab, "Actuators")
        layout.addWidget(tabs, 1)

    def load_from_case(self, case) -> None:
        b = case.bcps
        self._actuator_tab.set_items(b.fairleads)
        self._anchor_tab.set_items(b.anchors)
        self._joint_tab.set_items(b.joints)
        self._bodybcp_tab.set_items(b.body_bcps)
        self._elastic_tab.set_items(b.elastic_anchors)

    def save_to_case(self, case) -> None:
        case.bcps = BCPsData(
            fairleads=self._actuator_tab.get_items(),
            anchors=self._anchor_tab.get_items(),
            joints=self._joint_tab.get_items(),
            body_bcps=self._bodybcp_tab.get_items(),
            elastic_anchors=self._elastic_tab.get_items(),
        )
