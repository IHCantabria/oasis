"""Results view — browse output channels and plot time series."""
from __future__ import annotations
import os
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QPushButton, QLabel, QComboBox, QSplitter, QTreeWidget,
    QTreeWidgetItem, QLineEdit,
)
from PyQt5.QtCore import Qt

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from data.output_reader import discover_outputs, load_channel


_PRESETS = {
    "Custom": None,
    "Body 1 – 6 DOF positions": ["DOF 1 – Body 1", "DOF 2 – Body 1", "DOF 3 – Body 1",
                                   "DOF 4 – Body 1", "DOF 5 – Body 1", "DOF 6 – Body 1"],
    "Line 1 end tensions": ["End tension – Line 1"],
    "Wave elevation": ["Wave elevation"],
    "OWC 1 pressure": ["OWC pressure – 1"],
}


class ResultsView(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._channels: dict = {}
        self._case = None
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        # Toolbar: open output dir + preset
        toolbar = QHBoxLayout()
        self._output_dir_le = QLineEdit()
        self._output_dir_le.setPlaceholderText("Output directory…")
        btn_browse = QPushButton("Browse…")
        btn_browse.clicked.connect(self._browse_output)
        btn_scan = QPushButton("Scan")
        btn_scan.clicked.connect(self._scan_output)
        toolbar.addWidget(QLabel("Output dir:"))
        toolbar.addWidget(self._output_dir_le)
        toolbar.addWidget(btn_browse)
        toolbar.addWidget(btn_scan)
        layout.addLayout(toolbar)

        preset_row = QHBoxLayout()
        self._preset = QComboBox()
        for name in _PRESETS:
            self._preset.addItem(name)
        self._preset.currentTextChanged.connect(self._apply_preset)
        preset_row.addWidget(QLabel("Preset:"))
        preset_row.addWidget(self._preset)
        preset_row.addStretch()
        btn_plot = QPushButton("Plot selected")
        btn_plot.clicked.connect(self._plot_selected)
        preset_row.addWidget(btn_plot)
        layout.addLayout(preset_row)

        # Splitter: channel tree | plot
        splitter = QSplitter(Qt.Horizontal)

        # Left: channel tree
        left = QWidget()
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 0, 0)
        lv.addWidget(QLabel("Available channels:"))
        self._tree = QTreeWidget()
        self._tree.setHeaderLabel("Channel")
        self._tree.setSelectionMode(QTreeWidget.MultiSelection)
        lv.addWidget(self._tree)
        splitter.addWidget(left)

        # Right: Matplotlib
        right = QWidget()
        rv = QVBoxLayout(right)
        rv.setContentsMargins(0, 0, 0, 0)
        self._fig = Figure(figsize=(8, 5), tight_layout=True)
        self._ax = self._fig.add_subplot(111)
        self._canvas = FigureCanvas(self._fig)
        rv.addWidget(self._canvas)
        splitter.addWidget(right)
        splitter.setSizes([250, 700])

        layout.addWidget(splitter)

    def set_case(self, case) -> None:
        self._case = case
        if case and case.project_path:
            output_dir = os.path.join(case.project_path, "output")
            if os.path.isdir(output_dir):
                self._output_dir_le.setText(output_dir)

    def _browse_output(self):
        from PyQt5.QtWidgets import QFileDialog
        path = QFileDialog.getExistingDirectory(self, "Select output directory")
        if path:
            self._output_dir_le.setText(path)

    def _scan_output(self):
        path = self._output_dir_le.text().strip()
        if not os.path.isdir(path):
            return
        self._channels = discover_outputs(path)
        self._tree.clear()
        # Group by prefix
        groups: dict[str, list[str]] = {}
        for name in sorted(self._channels):
            prefix = name.split("–")[0].strip() if "–" in name else name
            groups.setdefault(prefix, []).append(name)
        for grp, names in groups.items():
            parent = QTreeWidgetItem(self._tree, [grp])
            for name in names:
                child = QTreeWidgetItem(parent, [name])
            parent.setExpanded(True)
        self._tree.repaint()

    def _apply_preset(self, preset_name: str):
        channels = _PRESETS.get(preset_name)
        if not channels:
            return
        self._tree.clearSelection()
        for name in channels:
            items = self._tree.findItems(name, Qt.MatchRecursive | Qt.MatchExactly)
            for item in items:
                item.setSelected(True)

    def _plot_selected(self):
        selected_names = []
        it = self._tree.selectedItems()
        for item in it:
            if item.parent():  # leaf
                selected_names.append(item.text(0))

        if not selected_names:
            # also try top-level with no children
            it2 = self._tree.selectedItems()
            for item in it2:
                selected_names.append(item.text(0))

        self._ax.clear()
        for name in selected_names:
            ch = self._channels.get(name)
            if ch is None:
                continue
            try:
                time, values = load_channel(ch.filepath, ch.column_index)
                self._ax.plot(time, values, label=name)
            except Exception as e:
                pass
        self._ax.set_xlabel("Time [s]")
        self._ax.legend(fontsize=7)
        self._ax.grid(True, linestyle="--", linewidth=0.5)
        self._canvas.draw()
