"""Wind turbines editor (OpenFAST)."""
from __future__ import annotations
from PyQt5.QtWidgets import (
    QVBoxLayout, QFormLayout, QGroupBox, QWidget,
)
from PyQt5.QtCore import Qt

from .base_editor import (
    ListFormEditor, make_spin, make_line_edit,
    make_file_browse_row, make_scrollable,
)
from data.case_data import WindTurbineData


class TurbinesEditor(ListFormEditor):
    def __init__(self, parent=None):
        super().__init__("Wind Turbines (OpenFAST)", parent)

    def _new_item_data(self):
        return WindTurbineData()

    def _item_label(self, item) -> str:
        return f"Turbine  body={item.body_index}"

    def _build_form(self) -> QWidget:
        inner = QWidget()
        vbox = QVBoxLayout(inner)
        vbox.setAlignment(Qt.AlignTop)
        form_box = QGroupBox("OpenFAST turbine")
        ff = QFormLayout(form_box)
        self._tr_body = make_spin(1, 1, 999)
        self._tr_fst_file = make_line_edit("")
        self._tr_fst_lib = make_line_edit("")
        ff.addRow("Body index:", self._tr_body)
        ff.addRow("FAST input file (.fst):",
                  make_file_browse_row("FAST input file", self._tr_fst_file, inner,
                                       file_filter="FAST files (*.fst);;All (*)"))
        ff.addRow("FAST library (.dll/.so):",
                  make_file_browse_row("FAST library", self._tr_fst_lib, inner,
                                       file_filter="Libraries (*.dll *.so);;All (*)"))
        vbox.addWidget(form_box)
        return make_scrollable(inner)

    def _populate_form(self, item: WindTurbineData) -> None:
        self._tr_body.setValue(item.body_index)
        self._tr_fst_file.setText(item.fast_input_file)
        self._tr_fst_lib.setText(item.fast_library)

    def _read_form(self):
        item = self._items[self._current_idx]
        item.body_index = self._tr_body.value()
        item.fast_input_file = self._tr_fst_file.text()
        item.fast_library = self._tr_fst_lib.text()
        return item

    def load_from_case(self, case) -> None:
        self._load_list(case.wind_turbines)

    def save_to_case(self, case) -> None:
        self._commit_current()
        case.wind_turbines = list(self._items)
