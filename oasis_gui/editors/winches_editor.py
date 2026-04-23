"""Winches editor."""
from __future__ import annotations
from PyQt5.QtWidgets import (
    QVBoxLayout, QFormLayout, QGroupBox, QWidget, QSplitter,
)
from PyQt5.QtCore import Qt

from .base_editor import (
    ListFormEditor, make_double_spin, make_spin, make_combo,
    make_scrollable, BaseEditor,
)
from data.case_data import WinchesData, WinchItemData, WinchControllerData


class _WinchListEditor(ListFormEditor):
    def __init__(self, parent=None):
        super().__init__("Winches", parent)

    def _new_item_data(self):
        return WinchItemData()

    def _item_label(self, item) -> str:
        return f"Winch  line={item.line}"

    def _build_form(self) -> QWidget:
        inner = QWidget()
        vbox = QVBoxLayout(inner)
        vbox.setAlignment(Qt.AlignTop)
        form_box = QGroupBox("Winch")
        ff = QFormLayout(form_box)
        self._w_line = make_spin(1, 1, 999)
        self._w_bcp = make_combo(["1 – First node", "2 – Last node"])
        self._w_inertia = make_double_spin(50.0, 0, 1e9, decimals=4)
        self._w_radius = make_double_spin(0.25, 0, 100, decimals=6)
        self._w_drag = make_double_spin(0.5, 0, 1e9, decimals=4)
        ff.addRow("Line index (1-based):", self._w_line)
        ff.addRow("Attachment BCP:", self._w_bcp)
        ff.addRow("Drum inertia [kg·m²]:", self._w_inertia)
        ff.addRow("Drum radius [m]:", self._w_radius)
        ff.addRow("Drag torque [N·s/rad]:", self._w_drag)
        vbox.addWidget(form_box)
        return make_scrollable(inner)

    def _populate_form(self, item: WinchItemData) -> None:
        self._w_line.setValue(item.line)
        self._w_bcp.setCurrentIndex(item.line_bcp - 1)
        self._w_inertia.setValue(item.inertia)
        self._w_radius.setValue(item.radius)
        self._w_drag.setValue(item.drag)

    def _read_form(self):
        item = self._items[self._current_idx]
        item.line = self._w_line.value()
        item.line_bcp = self._w_bcp.currentIndex() + 1
        item.inertia = self._w_inertia.value()
        item.radius = self._w_radius.value()
        item.drag = self._w_drag.value()
        return item


class WinchesEditor(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        # Controller (top)
        ctrl_box = QGroupBox("Winch controller")
        cf = QFormLayout(ctrl_box)
        self._ctrl_type = make_combo(["1 – Constant tension", "2 – Horizontal component"])
        self._ctrl_tension = make_double_spin(30000.0, 0, 1e12, decimals=2)
        cf.addRow("Controller type:", self._ctrl_type)
        cf.addRow("Target tension [N]:", self._ctrl_tension)
        layout.addWidget(ctrl_box)

        # Winch list (below)
        self._winch_list_ed = _WinchListEditor()
        layout.addWidget(self._winch_list_ed)

    def load_from_case(self, case) -> None:
        wd = case.winches
        self._ctrl_type.setCurrentIndex(wd.controller.controller_type - 1)
        self._ctrl_tension.setValue(wd.controller.target_tension)
        self._winch_list_ed.set_items(wd.winches)

    def save_to_case(self, case) -> None:
        self._winch_list_ed._commit_current()
        case.winches = WinchesData(
            winches=list(self._winch_list_ed._items),
            controller=WinchControllerData(
                controller_type=self._ctrl_type.currentIndex() + 1,
                target_tension=self._ctrl_tension.value(),
            ),
        )
