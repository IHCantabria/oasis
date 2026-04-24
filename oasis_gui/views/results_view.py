"""Results view — browse output channels and plot time series."""
from __future__ import annotations
import os
import re
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QSplitter, QTreeWidget,
    QTreeWidgetItem, QLineEdit, QTabWidget, QComboBox,
)
from PyQt5.QtCore import Qt

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from data.output_reader import discover_outputs, load_channel, load_channel_magnitude

_TIME_LABEL = "Time [s]"
_DOF_NAMES = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]


class ResultsView(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._channels: dict = {}
        self._case = None
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        # ── Toolbar: output dir
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

        # ── Splitter: left panel | tab widget
        splitter = QSplitter(Qt.Horizontal)

        # ── Left panel ────────────────────────────────────────────────────────
        left = QWidget()
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 2, 4, 0)

        # X axis selector
        lv.addWidget(QLabel("X axis:"))
        self._x_combo = QComboBox()
        self._x_combo.addItem(_TIME_LABEL)   # time always first
        lv.addWidget(self._x_combo)

        # Y channel tree (multi-select)
        lv.addWidget(QLabel("Y channels (multi-select):"))
        self._tree = QTreeWidget()
        self._tree.setHeaderLabel("Channel")
        self._tree.setSelectionMode(QTreeWidget.MultiSelection)
        lv.addWidget(self._tree, 1)

        # Buttons
        btn_row = QHBoxLayout()
        btn_custom = QPushButton("Plot custom")
        btn_custom.clicked.connect(self._plot_custom)
        btn_all = QPushButton("Plot all bodies")
        btn_all.clicked.connect(self._plot_all_bodies)
        btn_row.addWidget(btn_custom)
        btn_row.addWidget(btn_all)
        lv.addLayout(btn_row)

        splitter.addWidget(left)

        # ── Right: tab widget
        self._tabs = QTabWidget()
        splitter.addWidget(self._tabs)
        splitter.setSizes([240, 730])

        layout.addWidget(splitter)

    # ── Public API ───────────────────────────────────────────────────────────

    def set_case(self, case) -> None:
        self._case = case
        if case and case.project_path:
            output_dir = os.path.join(case.project_path, "output")
            if os.path.isdir(output_dir):
                self._output_dir_le.setText(output_dir)

    # ── Slots ────────────────────────────────────────────────────────────────

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

        # Populate X combo: Time first, then all channel names sorted
        self._x_combo.blockSignals(True)
        prev_x = self._x_combo.currentText()
        self._x_combo.clear()
        self._x_combo.addItem(_TIME_LABEL)
        for name in sorted(self._channels):
            self._x_combo.addItem(name)
        idx = self._x_combo.findText(prev_x)
        self._x_combo.setCurrentIndex(idx if idx >= 0 else 0)
        self._x_combo.blockSignals(False)

        # Populate Y tree
        self._tree.clear()
        groups: dict[str, list[str]] = {}
        for name in sorted(self._channels):
            parts = name.split(" / ")
            prefix = parts[0] if len(parts) > 1 else name
            groups.setdefault(prefix, []).append(name)
        for grp, names in groups.items():
            parent = QTreeWidgetItem(self._tree, [grp])
            for name in names:
                QTreeWidgetItem(parent, [name])
            parent.setExpanded(True)
        self._tree.repaint()

    def _plot_custom(self):
        """Plot user-selected Y channels against the chosen X channel."""
        # Collect selected leaf names
        y_names = [item.text(0) for item in self._tree.selectedItems()
                   if item.parent() is not None]
        if not y_names:
            return

        x_label = self._x_combo.currentText()
        x_is_time = (x_label == _TIME_LABEL)
        x_ch = None if x_is_time else self._channels.get(x_label)

        fig = Figure(figsize=(9, 5), tight_layout=True)
        ax = fig.add_subplot(111)

        for y_name in y_names:
            ch = self._channels.get(y_name)
            if ch is None:
                continue
            try:
                fy, cy = ch[0], ch[1]
                if isinstance(cy, tuple):
                    t_y, y_vals = load_channel_magnitude(fy, cy)
                else:
                    t_y, y_vals = load_channel(fy, cy)

                if x_is_time:
                    ax.plot(t_y, y_vals, label=y_name, linewidth=0.9)
                else:
                    if x_ch is None:
                        continue
                    fx, cx = x_ch[0], x_ch[1]
                    if isinstance(cx, tuple):
                        _, x_vals = load_channel_magnitude(fx, cx)
                    else:
                        _, x_vals = load_channel(fx, cx)
                    # Align lengths
                    n = min(len(x_vals), len(y_vals))
                    ax.plot(x_vals[:n], y_vals[:n], label=y_name, linewidth=0.9)
            except Exception:
                pass

        ax.set_xlabel(x_label, fontsize=8)
        ax.legend(fontsize=7)
        ax.grid(True, linestyle="--", linewidth=0.4)
        ax.tick_params(labelsize=7)

        canvas = FigureCanvas(fig)

        # Replace existing Custom tab if present
        for i in range(self._tabs.count()):
            if self._tabs.tabText(i) == "Custom":
                old = self._tabs.widget(i)
                self._tabs.removeTab(i)
                old.deleteLater()
                break
        self._tabs.insertTab(0, canvas, "Custom")
        self._tabs.setCurrentIndex(0)

    def _plot_all_bodies(self):
        """Scan (if needed), then create one tab per body with 6 DOF subplots."""
        path = self._output_dir_le.text().strip()
        if not os.path.isdir(path):
            return
        if not self._channels:
            self._scan_output()

        # Find all body IDs that have position channels
        body_ids: list[int] = []
        for name in self._channels:
            m = re.match(r"Body (\d+) / \w+ / Position", name)
            if m:
                bid = int(m.group(1))
                if bid not in body_ids:
                    body_ids.append(bid)
        body_ids.sort()

        if not body_ids:
            return

        # Remove old body tabs (keep any non-body tabs untouched)
        for i in reversed(range(self._tabs.count())):
            if self._tabs.tabText(i).startswith("Body "):
                widget = self._tabs.widget(i)
                self._tabs.removeTab(i)
                widget.deleteLater()

        for bid in body_ids:
            canvas = self._make_body_tab(bid)
            self._tabs.addTab(canvas, f"Body {bid}")

        if self._tabs.count():
            self._tabs.setCurrentIndex(0)

    # ── Helpers ──────────────────────────────────────────────────────────────

    def _make_body_tab(self, bid: int) -> QWidget:
        """Return a widget containing a 2×3 subplot figure for body *bid*."""
        fig = Figure(figsize=(10, 7), tight_layout=True)
        axes = fig.subplots(nrows=2, ncols=3)  # shape (2, 3)

        dof_labels = _DOF_NAMES  # Surge Sway Heave Roll Pitch Yaw
        flat_axes = [axes[r][c] for r in range(2) for c in range(3)]

        for ax, dof_name in zip(flat_axes, dof_labels):
            ch_name = f"Body {bid} / {dof_name} / Position"
            ch = self._channels.get(ch_name)
            ax.set_title(dof_name, fontsize=9)
            ax.set_xlabel("Time [s]", fontsize=7)
            ax.grid(True, linestyle="--", linewidth=0.4)
            ax.tick_params(labelsize=7)
            if ch is not None:
                filepath, col = ch[0], ch[1]
                try:
                    if isinstance(col, tuple):
                        time, values = load_channel_magnitude(filepath, col)
                    else:
                        time, values = load_channel(filepath, col)
                    ax.plot(time, values, linewidth=0.8)
                    ax.set_ylabel("Position [m / rad]", fontsize=7)
                except Exception:
                    ax.text(0.5, 0.5, "load error", transform=ax.transAxes,
                            ha="center", va="center", fontsize=8, color="red")
            else:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=8, color="gray")

        canvas = FigureCanvas(fig)
        return canvas
