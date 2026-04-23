"""2D schematic view using Matplotlib."""
from __future__ import annotations
import math
import numpy as np
from typing import Optional

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QComboBox, QLabel, QPushButton,
)

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

# BCP type colors
_BCP_COLORS = {
    "fairlead": "royalblue",
    "anchor": "saddlebrown",
    "joint": "gray",
    "body_bcp": "limegreen",
    "elastic_anchor": "orange",
}


class SchematicView(QWidget):
    """2D schematic of bodies, BCPs, lines and springs."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._case = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        # Toolbar
        toolbar = QHBoxLayout()
        toolbar.addWidget(QLabel("Projection:"))
        self._proj = QComboBox()
        self._proj.addItems(["XZ (side view)", "XY (top view)"])
        self._proj.currentIndexChanged.connect(self.refresh)
        toolbar.addWidget(self._proj)
        toolbar.addStretch()
        btn_refresh = QPushButton("↺ Refresh")
        btn_refresh.clicked.connect(self.refresh)
        toolbar.addWidget(btn_refresh)
        layout.addLayout(toolbar)

        # Matplotlib canvas
        self._fig = Figure(figsize=(8, 6), tight_layout=True)
        self._ax = self._fig.add_subplot(111)
        self._canvas = FigureCanvas(self._fig)
        layout.addWidget(self._canvas)

        # Tooltip annotation
        self._annot = self._ax.annotate(
            "", xy=(0, 0), xytext=(10, 10), textcoords="offset points",
            bbox=dict(boxstyle="round", fc="w"),
            arrowprops=dict(arrowstyle="->"),
            fontsize=8,
        )
        self._annot.set_visible(False)

        self._bcp_positions: list = []  # (label, x, y) for hover
        self._canvas.mpl_connect("motion_notify_event", self._on_hover)

    def update_case(self, case) -> None:
        self._case = case
        self.refresh()

    def refresh(self) -> None:
        if self._case is None:
            return
        self._ax.clear()
        self._bcp_positions = []
        xz = self._proj.currentIndex() == 0

        case = self._case
        water_depth = case.problem.water_depth

        # Horizontal axis: X; vertical: Z or Y
        ix, iy = (0, 2) if xz else (0, 1)
        horiz_label = "X [m]"
        vert_label = "Z [m]" if xz else "Y [m]"

        # Water surface (z=0 in XZ)
        if xz:
            xlim = [-200, 200]
            self._ax.axhline(0, color="deepskyblue", linewidth=1.5, label="Water surface")
            self._ax.axhline(water_depth, color="sienna", linewidth=1.5, linestyle="--",
                              label=f"Seabed (z={water_depth}m)")

        # Build BCP global map: global index → (type, position)
        bcp_map: dict[int, tuple[str, list]] = {}
        idx = 1
        for fl in case.bcps.fairleads:
            bcp_map[idx] = ("fairlead", fl.position)
            idx += 1
        for an in case.bcps.anchors:
            bcp_map[idx] = ("anchor", an.position)
            idx += 1
        for jt in case.bcps.joints:
            bcp_map[idx] = ("joint", jt.position)
            idx += 1
        for bb in case.bcps.body_bcps:
            bcp_map[idx] = ("body_bcp", bb.position)
            idx += 1
        for ea in case.bcps.elastic_anchors:
            bcp_map[idx] = ("elastic_anchor", ea.position)
            idx += 1

        # Draw bodies as rectangles
        for body in case.bodies:
            pos = body.initial_position
            if len(pos) < 3:
                continue
            cx, cy = pos[ix], pos[iy]
            w = body.schematic_width
            h = body.schematic_height
            rect = plt_rectangle(cx - w / 2, cy - h / 2, w, h,
                                  edgecolor="navy", facecolor="lightsteelblue", linewidth=2)
            self._ax.add_patch(rect)
            self._ax.text(cx, cy, body.name or "Body", ha="center", va="center",
                          fontsize=8, color="navy")

        # Draw lines as segments between BCPs
        for line in case.lines:
            p1 = bcp_map.get(line.BCP_1)
            p2 = bcp_map.get(line.BCP_N)
            if p1 and p2:
                x1, y1 = p1[1][ix], p1[1][iy]
                x2, y2 = p2[1][ix], p2[1][iy]
                self._ax.plot([x1, x2], [y1, y2], color="dimgray", linewidth=1.5)

        # Draw springs as zigzag
        for spring in case.springs:
            p1 = bcp_map.get(spring.BCP_1)
            p2 = bcp_map.get(spring.BCP_2)
            if p1 and p2:
                x1, y1 = p1[1][ix], p1[1][iy]
                x2, y2 = p2[1][ix], p2[1][iy]
                _draw_zigzag(self._ax, x1, y1, x2, y2, n=6, amplitude=0.03, color="darkorange")

        # Draw BCPs as colored markers
        for gidx, (bcp_type, pos) in bcp_map.items():
            px, py = pos[ix], pos[iy]
            color = _BCP_COLORS.get(bcp_type, "black")
            self._ax.plot(px, py, "o", color=color, markersize=7, zorder=5)
            self._bcp_positions.append((f"#{gidx} {bcp_type} ({px:.1f},{py:.1f})", px, py))

        # Legend for BCP types (once)
        from matplotlib.lines import Line2D
        legend_handles = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor=c, label=t, markersize=8)
            for t, c in _BCP_COLORS.items()
        ]
        self._ax.legend(handles=legend_handles, loc="best", fontsize=7)

        self._ax.set_xlabel(horiz_label)
        self._ax.set_ylabel(vert_label)
        self._ax.set_aspect("equal", adjustable="datalim")
        self._ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)
        self._canvas.draw()

    def _on_hover(self, event):
        if event.inaxes != self._ax or not self._bcp_positions:
            self._annot.set_visible(False)
            self._canvas.draw_idle()
            return
        for label, px, py in self._bcp_positions:
            dist = math.hypot(event.xdata - px, event.ydata - py)
            if dist < 2.0:
                self._annot.set_text(label)
                self._annot.xy = (px, py)
                self._annot.set_visible(True)
                self._canvas.draw_idle()
                return
        self._annot.set_visible(False)
        self._canvas.draw_idle()


def plt_rectangle(x, y, w, h, **kwargs):
    from matplotlib.patches import FancyBboxPatch
    return FancyBboxPatch((x, y), w, h, boxstyle="square,pad=0", **kwargs)


def _draw_zigzag(ax, x1, y1, x2, y2, n=6, amplitude=0.05, color="darkorange"):
    """Draw a spring-like zigzag between two points."""
    dx = x2 - x1
    dy = y2 - y1
    length = math.hypot(dx, dy)
    if length < 1e-9:
        return
    ux, uy = dx / length, dy / length
    nx, ny = -uy, ux  # normal

    n_pts = n * 2 + 2
    xs, ys = [], []
    for i in range(n_pts):
        t = i / (n_pts - 1)
        base_x = x1 + t * dx
        base_y = y1 + t * dy
        if 0 < i < n_pts - 1:
            sign = 1 if i % 2 == 1 else -1
            off = amplitude * length * sign
            base_x += off * nx
            base_y += off * ny
        xs.append(base_x)
        ys.append(base_y)
    ax.plot(xs, ys, color=color, linewidth=1.2)
