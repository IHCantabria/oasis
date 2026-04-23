"""Base editor class and shared utilities for all OASIS module editors."""
from __future__ import annotations
from typing import Any

from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox, QCheckBox,
    QPushButton, QListWidget, QListWidgetItem, QSplitter,
    QFileDialog, QLabel, QSizePolicy, QScrollArea, QFrame,
)
from PyQt5.QtCore import Qt


class BaseEditor(QWidget):
    """Abstract base for all module editors.

    Subclasses implement:
        load_from_case(case: CaseData)  – populate widgets from data model
        save_to_case(case: CaseData)    – write widget values back to data model
    """
    data_changed = pyqtSignal()

    def load_from_case(self, case) -> None:  # noqa: ANN001
        raise NotImplementedError

    def save_to_case(self, case) -> None:  # noqa: ANN001
        raise NotImplementedError


# ──────────────────────────────────────────────────────────────────────────────
# Widget factory helpers
# ──────────────────────────────────────────────────────────────────────────────
def make_double_spin(value: float = 0.0, lo: float = -1e18, hi: float = 1e18,
                     decimals: int = 6, step: float = 0.1) -> QDoubleSpinBox:
    w = QDoubleSpinBox()
    w.setDecimals(decimals)
    w.setRange(lo, hi)
    w.setSingleStep(step)
    w.setValue(value)
    w.setMinimumWidth(120)
    return w


def make_spin(value: int = 0, lo: int = -999999, hi: int = 999999) -> QSpinBox:
    w = QSpinBox()
    w.setRange(lo, hi)
    w.setValue(value)
    return w


def make_line_edit(text: str = "") -> QLineEdit:
    w = QLineEdit(str(text))
    return w


def make_combo(choices: list, current: Any = None) -> QComboBox:
    w = QComboBox()
    for c in choices:
        w.addItem(str(c))
    if current is not None:
        idx = w.findText(str(current))
        if idx >= 0:
            w.setCurrentIndex(idx)
    return w


def make_file_browse_row(label: str, line_edit: QLineEdit,
                          parent_widget: QWidget,
                          file_filter: str = "All Files (*)",
                          save: bool = False) -> QHBoxLayout:
    """Returns an HBoxLayout with a label, line_edit and Browse button."""
    row = QHBoxLayout()
    btn = QPushButton("Browse…")

    def _browse():
        if save:
            path, _ = QFileDialog.getSaveFileName(parent_widget, f"Select {label}", "", file_filter)
        else:
            path, _ = QFileDialog.getOpenFileName(parent_widget, f"Select {label}", "", file_filter)
        if path:
            line_edit.setText(path)

    btn.clicked.connect(_browse)
    row.addWidget(line_edit)
    row.addWidget(btn)
    return row


def make_scrollable(widget: QWidget) -> QScrollArea:
    """Wrap a widget in a QScrollArea."""
    sa = QScrollArea()
    sa.setWidgetResizable(True)
    sa.setWidget(widget)
    sa.setFrameShape(QFrame.NoFrame)
    return sa


def section_box(title: str) -> QGroupBox:
    gb = QGroupBox(title)
    gb.setLayout(QFormLayout())
    gb.layout().setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
    return gb


# ──────────────────────────────────────────────────────────────────────────────
# Generic list+form editor (used by Bodies, BCPs, Lines, Springs, etc.)
# ──────────────────────────────────────────────────────────────────────────────
class ListFormEditor(BaseEditor):
    """A generic editor with a list on the left and a form on the right.

    Subclasses must implement:
        _new_item_data() -> object   – create a new default data object
        _build_form()                – build self._form_widget (QWidget)
        _populate_form(item_data)    – fill form widgets from item_data
        _read_form() -> object       – read form widgets and return item_data
        _item_label(item_data) -> str – label for list widget item
    """

    def __init__(self, title: str = "Items", parent=None):
        super().__init__(parent)
        self._items: list = []
        self._current_idx: int = -1
        self._building = False

        main = QVBoxLayout(self)
        main.setContentsMargins(4, 4, 4, 4)

        label = QLabel(f"<b>{title}</b>")
        main.addWidget(label)

        splitter = QSplitter(Qt.Horizontal)
        main.addWidget(splitter)

        # Left: list + buttons
        left = QWidget()
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 0, 0)
        self._list = QListWidget()
        self._list.setMaximumWidth(180)
        self._list.setMinimumWidth(120)
        lv.addWidget(self._list)

        btn_row = QHBoxLayout()
        self._btn_add = QPushButton("+")
        self._btn_add.setToolTip("Add item")
        self._btn_add.setMaximumWidth(36)
        self._btn_remove = QPushButton("−")
        self._btn_remove.setToolTip("Remove selected")
        self._btn_remove.setMaximumWidth(36)
        self._btn_up = QPushButton("↑")
        self._btn_up.setMaximumWidth(36)
        self._btn_down = QPushButton("↓")
        self._btn_down.setMaximumWidth(36)
        btn_row.addWidget(self._btn_add)
        btn_row.addWidget(self._btn_remove)
        btn_row.addWidget(self._btn_up)
        btn_row.addWidget(self._btn_down)
        btn_row.addStretch()
        lv.addLayout(btn_row)
        splitter.addWidget(left)

        # Right: form (subclass builds this)
        self._form_container = QWidget()
        self._form_container.setLayout(QVBoxLayout())
        self._form_container.layout().setContentsMargins(4, 0, 0, 0)
        form_scroll = make_scrollable(self._form_container)
        splitter.addWidget(form_scroll)

        splitter.setSizes([160, 600])

        # Connections
        self._list.currentRowChanged.connect(self._on_selection_changed)
        self._btn_add.clicked.connect(self._on_add)
        self._btn_remove.clicked.connect(self._on_remove)
        self._btn_up.clicked.connect(self._on_move_up)
        self._btn_down.clicked.connect(self._on_move_down)

        # Allow subclass to build form (called after __init__ to use polymorphism)
        form_widget = self._build_form()
        if form_widget is not None:
            self._form_container.layout().addWidget(form_widget)

    # ── Public API ────────────────────────────────────────────────────────────
    def set_items(self, items: list) -> None:
        self._items = list(items)
        self._refresh_list()
        if self._items:
            self._list.setCurrentRow(0)
        else:
            self._clear_form()

    # Alias
    def _load_list(self, items: list) -> None:
        self.set_items(items)

    def get_items(self) -> list:
        self._commit_current()
        return list(self._items)

    # ── Internal ──────────────────────────────────────────────────────────────
    def _refresh_list(self):
        self._list.blockSignals(True)
        self._list.clear()
        for item in self._items:
            self._list.addItem(self._item_label(item))
        self._list.blockSignals(False)

    def _on_selection_changed(self, row: int):
        if self._building:
            return
        # Save current form before switching
        if self._current_idx >= 0 and self._current_idx < len(self._items):
            self._commit_current()
        self._current_idx = row
        if 0 <= row < len(self._items):
            self._populate_form(self._items[row])
        else:
            self._clear_form()

    def _commit_current(self):
        if 0 <= self._current_idx < len(self._items):
            self._items[self._current_idx] = self._read_form()
            label = self._item_label(self._items[self._current_idx])
            item = self._list.item(self._current_idx)
            if item:
                item.setText(label)

    def _on_add(self):
        new_item = self._new_item_data()
        self._items.append(new_item)
        self._list.addItem(self._item_label(new_item))
        self._list.setCurrentRow(len(self._items) - 1)
        self.data_changed.emit()

    def _on_remove(self):
        row = self._list.currentRow()
        if 0 <= row < len(self._items):
            self._items.pop(row)
            self._current_idx = -1
            self._refresh_list()
            new_row = min(row, len(self._items) - 1)
            if new_row >= 0:
                self._list.setCurrentRow(new_row)
            else:
                self._clear_form()
            self.data_changed.emit()

    def _on_move_up(self):
        row = self._list.currentRow()
        if row > 0:
            self._commit_current()
            self._items[row], self._items[row - 1] = self._items[row - 1], self._items[row]
            self._refresh_list()
            self._list.setCurrentRow(row - 1)
            self.data_changed.emit()

    def _on_move_down(self):
        row = self._list.currentRow()
        if 0 <= row < len(self._items) - 1:
            self._commit_current()
            self._items[row], self._items[row + 1] = self._items[row + 1], self._items[row]
            self._refresh_list()
            self._list.setCurrentRow(row + 1)
            self.data_changed.emit()

    def _clear_form(self):
        pass  # subclasses can override

    # ── Subclass contract ─────────────────────────────────────────────────────
    def _new_item_data(self):
        raise NotImplementedError

    def _build_form(self) -> QWidget:
        raise NotImplementedError

    def _populate_form(self, item_data) -> None:
        raise NotImplementedError

    def _read_form(self):
        raise NotImplementedError

    def _item_label(self, item_data) -> str:
        return str(item_data)
