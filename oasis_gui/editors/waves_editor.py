"""Waves editor."""
from __future__ import annotations
from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QGroupBox, QFormLayout, QWidget,
    QComboBox, QLabel, QStackedWidget,
)
from PyQt5.QtCore import Qt

from .base_editor import BaseEditor, make_double_spin, make_spin, make_scrollable, make_line_edit, make_file_browse_row
from data.case_data import WaveData


class WavesEditor(BaseEditor):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self):
        inner = QWidget()
        vbox = QVBoxLayout(inner)
        vbox.setAlignment(Qt.AlignTop)

        # ── Wave type selector ────────────────────────────────────────────────
        type_box = QGroupBox("Wave type")
        type_form = QFormLayout(type_box)
        self._wave_type = QComboBox()
        self._wave_type.addItems(["REG – Regular", "IRR – Irregular"])
        self._wave_type.currentIndexChanged.connect(self._on_type_changed)
        type_form.addRow("Type:", self._wave_type)
        vbox.addWidget(type_box)

        # ── Basic parameters (always shown) ──────────────────────────────────
        basic_box = QGroupBox("Basic parameters")
        basic_form = QFormLayout(basic_box)
        self._height = make_double_spin(0.0, 0, 100, decimals=4)
        self._period = make_double_spin(10.0, 0.01, 1000, decimals=4)
        self._heading = make_double_spin(0.0, -360, 360, decimals=2)
        self._ramp_time = make_double_spin(0.0, 0, 10000, decimals=2)
        basic_form.addRow("Height H [m]:", self._height)
        basic_form.addRow("Period T [s]:", self._period)
        basic_form.addRow("Heading [deg]:", self._heading)
        basic_form.addRow("Ramp time [s]:", self._ramp_time)
        vbox.addWidget(basic_box)

        # ── Irregular-only parameters ─────────────────────────────────────────
        self._irr_box = QGroupBox("Irregular wave parameters")
        irr_vbox = QVBoxLayout(self._irr_box)

        spec_form = QFormLayout()
        self._spec_type = QComboBox()
        self._spec_type.addItems([
            "1 – JONSWAP analytical",
            "2 – Custom time-series (wave.dat)",
            "3 – Custom frequency-domain (wavefreq.dat)",
        ])
        self._spec_type.currentIndexChanged.connect(self._on_spec_type_changed)
        self._piecewise = QComboBox()
        self._piecewise.addItems(["0 – Single piece", "1 – Piecewise"])
        spec_form.addRow("Spectrum type:", self._spec_type)
        spec_form.addRow("Piecewise flag:", self._piecewise)
        irr_vbox.addLayout(spec_form)

        # JONSWAP sub-params
        self._jonswap_box = QGroupBox("JONSWAP parameters")
        jf = QFormLayout(self._jonswap_box)
        self._gamma = make_double_spin(3.3, 1.0, 10.0, decimals=4)
        self._spreading = make_double_spin(0.0, 0, 100, decimals=4)
        self._dtheta = make_double_spin(10.0, 1, 360, decimals=2)
        self._factor = make_double_spin(0.001, 1e-10, 1.0, decimals=8)
        self._read_phases = QComboBox()
        self._read_phases.addItems(["0 – Generate random", "1 – Read from file"])
        self._rel_tol = make_double_spin(0.05, 1e-10, 1.0, decimals=6)
        self._dt = make_double_spin(0.1, 1e-6, 1000, decimals=6)
        self._phases_file = make_line_edit("WavePhases.dat")
        jf.addRow("Peak enhancement γ:", self._gamma)
        jf.addRow("Directional spreading s:", self._spreading)
        jf.addRow("Heading step dθ [deg]:", self._dtheta)
        jf.addRow("Min spectrum fraction:", self._factor)
        jf.addRow("Phase source:", self._read_phases)
        jf.addRow("Phase relative tolerance:", self._rel_tol)
        jf.addRow("Wave series dt [s]:", self._dt)
        jf.addRow("Phases filename:", self._phases_file)
        irr_vbox.addWidget(self._jonswap_box)

        # Custom data file sub-param
        self._custom_box = QGroupBox("Custom data file")
        cf = QFormLayout(self._custom_box)
        self._data_file = make_line_edit("wave.dat")
        browse_row = make_file_browse_row("wave data file", self._data_file, self,
                                          file_filter="Data files (*.dat *.txt);;All Files (*)")
        cf.addRow("Data file:", browse_row)
        irr_vbox.addWidget(self._custom_box)
        self._custom_box.setVisible(False)

        vbox.addWidget(self._irr_box)
        self._irr_box.setVisible(False)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(make_scrollable(inner))

    def _on_type_changed(self, idx):
        self._irr_box.setVisible(idx == 1)

    def _on_spec_type_changed(self, idx):
        self._jonswap_box.setVisible(idx == 0)
        self._custom_box.setVisible(idx > 0)

    # ── Data binding ──────────────────────────────────────────────────────────
    def load_from_case(self, case) -> None:
        w = case.waves
        self._wave_type.setCurrentIndex(0 if w.wave_type == "REG" else 1)
        self._height.setValue(w.height)
        self._period.setValue(w.period)
        self._heading.setValue(w.heading)
        self._ramp_time.setValue(w.ramp_time)
        self._spec_type.setCurrentIndex(w.spectrum_type - 1)
        self._piecewise.setCurrentIndex(w.piecewise_flag)
        self._gamma.setValue(w.gamma)
        self._spreading.setValue(w.spreading)
        self._dtheta.setValue(w.dtheta)
        self._factor.setValue(w.factor)
        self._read_phases.setCurrentIndex(w.read_phases_flag)
        self._rel_tol.setValue(w.wave_relative_tolerance)
        self._dt.setValue(w.dt)
        self._phases_file.setText(w.phases_file)
        self._data_file.setText(w.database_file)
        self._on_type_changed(0 if w.wave_type == "REG" else 1)
        self._on_spec_type_changed(w.spectrum_type - 1)

    def save_to_case(self, case) -> None:
        w = case.waves
        w.wave_type = "REG" if self._wave_type.currentIndex() == 0 else "IRR"
        w.height = self._height.value()
        w.period = self._period.value()
        w.heading = self._heading.value()
        w.ramp_time = self._ramp_time.value()
        w.spectrum_type = self._spec_type.currentIndex() + 1
        w.piecewise_flag = self._piecewise.currentIndex()
        w.gamma = self._gamma.value()
        w.spreading = self._spreading.value()
        w.dtheta = self._dtheta.value()
        w.factor = self._factor.value()
        w.read_phases_flag = self._read_phases.currentIndex()
        w.wave_relative_tolerance = self._rel_tol.value()
        w.dt = self._dt.value()
        w.phases_file = self._phases_file.text()
        w.database_file = self._data_file.text()
