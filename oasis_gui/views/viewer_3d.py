"""3D viewer using PyVista (optional dependency)."""
from __future__ import annotations
import os
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QSlider, QFileDialog, QGroupBox, QFormLayout, QSpinBox,
)
from PyQt5.QtCore import Qt, QTimer

try:
    import pyvista as pv
    from pyvistaqt import QtInteractor
    _PYVISTA_OK = True
except ImportError:
    _PYVISTA_OK = False


class Viewer3D(QWidget):
    """3D viewer for STL bodies with 6-DOF animation and mooring lines."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._case = None
        self._case_path: str = ""
        self._actors: dict = {}   # body_name → (mesh_actor, initial_mesh)
        self._line_actors: list = []
        self._time_data: dict = {}  # body_index → [[6 DOF time series], time_array]
        self._current_frame = 0
        self._n_frames = 0
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._next_frame)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        if not _PYVISTA_OK:
            layout.addWidget(QLabel(
                "PyVista / pyvistaqt not installed.\n"
                "Run: pip install pyvista pyvistaqt"
            ))
            return

        # ── Controls ──────────────────────────────────────────────────────────
        ctrl_row = QHBoxLayout()
        self._btn_load = QPushButton("Load STLs from case")
        self._btn_load.clicked.connect(self._load_stls)
        self._btn_play = QPushButton("▶ Play")
        self._btn_play.setEnabled(False)
        self._btn_play.clicked.connect(self._toggle_play)
        self._btn_stop = QPushButton("■ Stop")
        self._btn_stop.setEnabled(False)
        self._btn_stop.clicked.connect(self._stop_anim)
        self._fps_spin = QSpinBox()
        self._fps_spin.setRange(1, 120)
        self._fps_spin.setValue(25)
        ctrl_row.addWidget(self._btn_load)
        ctrl_row.addWidget(self._btn_play)
        ctrl_row.addWidget(self._btn_stop)
        ctrl_row.addWidget(QLabel("FPS:"))
        ctrl_row.addWidget(self._fps_spin)
        ctrl_row.addStretch()
        layout.addLayout(ctrl_row)

        # Slider
        slider_row = QHBoxLayout()
        slider_row.addWidget(QLabel("Frame:"))
        self._slider = QSlider(Qt.Horizontal)
        self._slider.setMinimum(0)
        self._slider.setMaximum(0)
        self._slider.valueChanged.connect(self._seek)
        self._frame_label = QLabel("0 / 0")
        slider_row.addWidget(self._slider)
        slider_row.addWidget(self._frame_label)
        layout.addLayout(slider_row)

        # PyVista interactor
        self._plotter = QtInteractor(self)
        self._plotter.show_axes()
        self._plotter.add_axes_at_origin()
        layout.addWidget(self._plotter.interactor)

    def set_case(self, case, case_path: str = "") -> None:
        self._case = case
        self._case_path = case_path

    def _load_stls(self):
        if not _PYVISTA_OK or self._case is None:
            return
        self._plotter.clear()
        self._actors = {}
        self._line_actors = []

        for body in self._case.bodies:
            stl_path = body.stl_3d_file
            if not stl_path:
                continue
            if not os.path.isabs(stl_path) and self._case_path:
                stl_path = os.path.join(self._case_path, stl_path)
            if os.path.isfile(stl_path):
                try:
                    mesh = pv.read(stl_path)
                    # Position at initial_position (x,y,z translation)
                    pos = body.initial_position[:3] if len(body.initial_position) >= 3 else [0, 0, 0]
                    mesh = mesh.translate(pos, inplace=False)
                    actor = self._plotter.add_mesh(mesh, color="lightsteelblue",
                                                    opacity=0.8, label=body.name)
                    self._actors[body.name] = actor
                except Exception:
                    pass

        # Add mooring lines (straight for now)
        bcp_map = self._build_bcp_map()
        for line in self._case.lines:
            p1 = bcp_map.get(line.BCP_1)
            p2 = bcp_map.get(line.BCP_N)
            if p1 and p2:
                pts = pv.Line(p1[1][:3], p2[1][:3])
                actor = self._plotter.add_mesh(pts, color="dimgray", line_width=2)
                self._line_actors.append(actor)

        # Add water plane
        wp = pv.Plane(center=(0, 0, 0), direction=(0, 0, 1), i_size=400, j_size=400)
        self._plotter.add_mesh(wp, color="deepskyblue", opacity=0.3)
        self._plotter.reset_camera()

        # Load DOF time series
        self._load_time_data()

    def _build_bcp_map(self) -> dict:
        bcp_map = {}
        idx = 1
        for fl in self._case.bcps.fairleads:
            bcp_map[idx] = ("fairlead", fl.position)
            idx += 1
        for an in self._case.bcps.anchors:
            bcp_map[idx] = ("anchor", an.position)
            idx += 1
        for jt in self._case.bcps.joints:
            bcp_map[idx] = ("joint", jt.position)
            idx += 1
        for bb in self._case.bcps.body_bcps:
            bcp_map[idx] = ("body_bcp", bb.position)
            idx += 1
        for ea in self._case.bcps.elastic_anchors:
            bcp_map[idx] = ("elastic_anchor", ea.position)
            idx += 1
        return bcp_map

    def _load_time_data(self):
        """Load DOF time series from output folder if available."""
        from data.output_reader import discover_outputs, load_channel
        if not self._case_path:
            return
        output_dir = os.path.join(self._case_path, "output")
        if not os.path.isdir(output_dir):
            return
        channels = discover_outputs(output_dir)
        self._time_data = {}
        for body_idx, body in enumerate(self._case.bodies, start=1):
            body_series = []
            time_arr = None
            all_found = True
            for dof in range(1, 7):
                key = f"DOF {dof} – Body {body_idx}"
                ch = channels.get(key)
                if ch is None:
                    all_found = False
                    break
                t, vals = load_channel(ch.filepath, ch.column_index)
                if time_arr is None:
                    time_arr = t
                body_series.append(vals)
            if all_found and body_series:
                self._time_data[body.name] = (time_arr, body_series)

        if self._time_data:
            first_name = next(iter(self._time_data))
            n = len(self._time_data[first_name][0])
            self._n_frames = n
            self._slider.setMaximum(n - 1)
            self._slider.setValue(0)
            self._btn_play.setEnabled(True)
        else:
            self._n_frames = 0
            self._btn_play.setEnabled(False)

    def _toggle_play(self):
        if self._timer.isActive():
            self._timer.stop()
            self._btn_play.setText("▶ Play")
        else:
            interval = max(1, int(1000 / self._fps_spin.value()))
            self._timer.start(interval)
            self._btn_play.setText("⏸ Pause")
            self._btn_stop.setEnabled(True)

    def _stop_anim(self):
        self._timer.stop()
        self._btn_play.setText("▶ Play")
        self._btn_stop.setEnabled(False)
        self._current_frame = 0
        self._slider.setValue(0)

    def _next_frame(self):
        if self._current_frame >= self._n_frames - 1:
            self._stop_anim()
            return
        self._current_frame += 1
        self._slider.setValue(self._current_frame)

    def _seek(self, frame: int):
        self._current_frame = frame
        self._frame_label.setText(f"{frame} / {self._n_frames}")
        # Update actor positions
        for body_name, (time_arr, body_series) in self._time_data.items():
            actor = self._actors.get(body_name)
            if actor is None or frame >= len(time_arr):
                continue
            # body_series: list of 6 arrays (x,y,z,rx,ry,rz)
            dx = body_series[0][frame] if len(body_series) > 0 else 0.0
            dy = body_series[1][frame] if len(body_series) > 1 else 0.0
            dz = body_series[2][frame] if len(body_series) > 2 else 0.0
            # Simple translation only (full 6-DOF rotation requires transform matrix)
            # Find initial position
            initial_pos = [0, 0, 0]
            for b in self._case.bodies:
                if b.name == body_name:
                    initial_pos = b.initial_position[:3]
                    break
            new_pos = [initial_pos[0] + dx, initial_pos[1] + dy, initial_pos[2] + dz]
            actor.SetPosition(*new_pos)
        self._plotter.render()
